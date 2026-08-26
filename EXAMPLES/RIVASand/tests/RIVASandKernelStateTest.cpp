/* Standalone state-contract tests for the OpenSees RIVA-Sand kernel copy. */
#include "../../../SRC/material/nD/RIVASand/RIVASandKernel.h"

#include <array>
#include <cstdlib>
#include <cstring>
#include <iostream>

static bool sameState(const riva_state_t &a, const riva_state_t &b)
{
    double av[RIVA_STATE_VALUE_COUNT] = {};
    double bv[RIVA_STATE_VALUE_COUNT] = {};
    if (riva_state_values(&a, av) != RIVA_STATE_VALUE_COUNT ||
        riva_state_values(&b, bv) != RIVA_STATE_VALUE_COUNT)
        return false;
    for (int i = 0; i < RIVA_STATE_VALUE_COUNT; ++i) {
        if (av[i] != bv[i]) {
            std::cerr << "state field " << i << " differs: "
                      << av[i] << " versus " << bv[i] << '\n';
            return false;
        }
    }
    return a.initialized == b.initialized;
}

static riva_state_t initialState(
    const riva_parameters_t &parameters,
    const riva_material_parameters_t &material, double voidRatio)
{
    const riva_tensor_t stress = {-19.4, -19.4, -40.0, 0, 0, 0};
    riva_state_t state = {};
    if (!riva_initialize_material(&parameters, &material, stress, voidRatio,
                                 &state) ||
        !riva_begin_dynamic_phase(&parameters, stress, &state)) {
        std::cerr << "RIVASand initialization failed\n";
        std::exit(1);
    }
    return state;
}

static bool advance(const riva_parameters_t &parameters,
                    const riva_material_parameters_t &material,
                    riva_tensor_t increment, int substeps,
                    riva_state_t &state)
{
    riva_tensor_t stress = state.stress;
    riva_update_info_t info = {};
    const int expected = riva_ddot(increment, increment) == 0.0 ? 0 : substeps;
    return riva_update_material(&parameters, &material, increment, substeps,
                               &state, &stress, &info) &&
           info.accepted_substeps == expected && riva_finite_tensor(stress);
}

int main()
{
    const riva_parameters_t parameters = riva_reference_parameters(1.0);
    const riva_material_parameters_t material =
        riva_reference_material_parameters(&parameters);
    if (!riva_material_parameters_valid(&parameters, &material)) {
        std::cerr << "reference material parameters are invalid\n";
        return 1;
    }

    const double Dr = (material.e_max-0.601)/(material.e_max-material.e_min);
    if (riva_void_ratio_from_material_relative_density(&material, Dr) != 0.601) {
        std::cerr << "relative-density conversion failed\n";
        return 1;
    }

    /* Components xy, yz, xz below are physical tensor shear increments.
     * The OpenSees adapter supplies half of each engineering shear increment. */
    const std::array<riva_tensor_t, 8> path = {{
        {0, 0, 0,  2.5e-4, 0, 0},
        {0, 0, 0,  2.5e-4, 0, 0},
        {0, 0, 0, -2.5e-4, 0, 0},
        {0, 0, 0, -2.5e-4, 0, 0},
        {-1e-5, -1e-5, -1e-5, -2.5e-4, 0, 0},
        {0, 0, 0,  2.5e-4, 0, 0},
        {0, 0, 0,  2.5e-4, 0, 0},
        {1e-5, 1e-5, 1e-5, -2.5e-4, 0, 0},
    }};

    riva_state_t baseline = initialState(parameters, material, 0.601);
    for (const riva_tensor_t increment : path)
        if (!advance(parameters, material, increment, 4, baseline)) return 1;

    /* The frozen wrapper and configurable entry point must agree exactly
     * when the custom row contains the frozen reference values. */
    riva_state_t frozen = initialState(parameters, material, 0.601);
    for (const riva_tensor_t increment : path) {
        riva_tensor_t stress = frozen.stress;
        riva_update_info_t info = {};
        if (!riva_update(&parameters, increment, 4, &frozen, &stress, &info))
            return 1;
    }
    if (!sameState(baseline, frozen)) {
        std::cerr << "custom reference row changed frozen V8 response\n";
        return 1;
    }

    riva_material_parameters_t changed = material;
    changed.zeta *= 2.0;
    if (!riva_material_parameters_valid(&parameters, &changed) ||
        initialState(parameters, changed, 0.601).D_ir ==
            initialState(parameters, material, 0.601).D_ir) {
        std::cerr << "zeta is not active in the custom material block\n";
        return 1;
    }

    /* Copying a committed state is the trial-state operation used by the
     * OpenSees adapter. A discarded trial must not alter the committed copy. */
    riva_state_t committed = initialState(parameters, material, 0.601);
    riva_state_t trial = committed;
    if (!advance(parameters, material, path[0], 4, trial) ||
        riva_ddot(riva_sub(committed.stress, trial.stress),
                 riva_sub(committed.stress, trial.stress)) == 0.0 ||
        !sameState(committed, initialState(parameters, material, 0.601))) {
        std::cerr << "trial/committed state isolation failed\n";
        return 1;
    }

    riva_state_t restarted = initialState(parameters, material, 0.601);
    for (std::size_t i = 0; i < path.size()/2; ++i)
        if (!advance(parameters, material, path[i], 4, restarted)) return 1;
    std::array<unsigned char, sizeof(riva_state_t)> bytes = {};
    std::memcpy(bytes.data(), &restarted, sizeof(restarted));
    restarted = {};
    std::memcpy(&restarted, bytes.data(), sizeof(restarted));
    for (std::size_t i = path.size()/2; i < path.size(); ++i)
        if (!advance(parameters, material, path[i], 4, restarted)) return 1;
    if (!sameState(baseline, restarted)) {
        std::cerr << "checkpoint/restart changed the response\n";
        return 1;
    }

    const riva_state_t beforeZero = restarted;
    if (!advance(parameters, material, riva_zero(), 4, restarted) ||
        !sameState(beforeZero, restarted)) {
        std::cerr << "zero increment changed history\n";
        return 1;
    }

    std::cout << "PASS: RIVA-Sand kernel reference/custom equivalence, "
                 "zeta exposure, trial isolation, restart, and zero increment\n";
    return 0;
}
