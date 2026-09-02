/* Standalone state and reversal-contract checks for the research kernel. */
#include "../../../SRC/material/nD/RIVASandIntermediateBiasResearch/RIVASandIntermediateBiasResearchKernel.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>

using namespace riva_ib_native;

static bool sameState(const riva_ib_state_t &a, const riva_ib_state_t &b)
{
    double av[RIVA_IB_STATE_VALUE_COUNT] = {};
    double bv[RIVA_IB_STATE_VALUE_COUNT] = {};
    if (riva_ib_state_values(&a, av) != RIVA_IB_STATE_VALUE_COUNT ||
        riva_ib_state_values(&b, bv) != RIVA_IB_STATE_VALUE_COUNT)
        return false;
    for (int i = 0; i < RIVA_IB_STATE_VALUE_COUNT; ++i) {
        if (av[i] != bv[i]) {
            std::cerr << "state value " << i << " differs: " << av[i]
                      << " versus " << bv[i] << '\n';
            return false;
        }
    }
    return a.base.initialized == b.base.initialized &&
           a.base.geostatic_admitted == b.base.geostatic_admitted;
}

static riva_ib_state_t biasedState(const riva_ib_parameters_t &parameters,
                                   const riva_material_parameters_t &material,
                                   double voidRatio, double shear)
{
    const tensor_t reference = {-19.4, -19.4, -40.0, 0, 0, 0};
    tensor_t stress = reference;
    stress.xz = shear;
    riva_ib_state_t state = {};
    if (!riva_ib_initialize_material(&parameters, &material, stress, voidRatio,
                                     &state) ||
        !riva_ib_begin_dynamic_phase(&parameters, &material, &reference,
                                     &state)) {
        std::cerr << "research-kernel initialization failed\n";
        std::exit(1);
    }
    return state;
}

static bool advance(const riva_ib_parameters_t &parameters,
                    const riva_material_parameters_t &material,
                    double engineeringShearIncrement, riva_ib_state_t &state)
{
    const tensor_t deps = {0, 0, 0, 0, 0, 0.5*engineeringShearIncrement};
    riva_update_info_t info = {};
    return riva_ib_update_material(&parameters, &material, deps, 4, &state,
                                   nullptr, &info) &&
           info.accepted_substeps == 4 && riva_finite_tensor(state.base.stress);
}

static bool driveCycles(const riva_ib_parameters_t &parameters,
                        const riva_material_parameters_t &material,
                        double amplitude, int cycles, riva_ib_state_t &state)
{
    double oldGamma = 0.0;
    for (int step = 1; step <= 32*cycles; ++step) {
        const double gamma = amplitude*std::sin(
            2.0*3.14159265358979323846*step/32.0);
        if (!advance(parameters, material, gamma-oldGamma, state)) return false;
        oldGamma = gamma;
    }
    return true;
}

int main()
{
    if (RIVA_IB_PARAMETER_COUNT != 249 || RIVA_IB_LOGICAL_STATE_COUNT != 63 ||
        RIVA_IB_STATE_VALUE_COUNT != 138 || RIVA_IB_KERNEL_REVISION != 3u) {
        std::cerr << "research kernel contract constants are inconsistent\n";
        return 1;
    }
    const riva_ib_parameters_t parameters = riva_ib_reference_parameters(1.0);
    const riva_material_parameters_t material =
        riva_reference_material_parameters(&parameters.base);
    if (!(parameters.bias_reversible_volume_ep_ref > 0.0)) return 1;

    /* Isolate the inherited bias wave and verify that a weakly plastic
     * microcycle cannot create a finite pressure-wave target. */
    riva_ib_parameters_t microParameters = parameters;
    microParameters.phase_transformation_enabled = 0;
    microParameters.mapping_backstress_enabled = 0;
    microParameters.intermediate_bias_flow_enabled = 0;
    riva_ib_state_t micro = biasedState(microParameters, material, 0.601, 10.0);
    if (!driveCycles(microParameters, material, 1.0e-5, 2, micro) ||
        !(micro.ep_half_last < microParameters.bias_reversible_volume_ep_ref) ||
        std::fabs(micro.base.bias_reversible_volume) >= 1.0e-10) {
        std::cerr << "plastic-activity gate did not suppress a microcycle\n";
        return 1;
    }

    /* Activation may initialize geometric anchors, but it must not invent
     * completed cyclic reversals or plastic activity. */
    riva_ib_state_t intermediate = biasedState(parameters, material, 0.601, 10.0);
    if (intermediate.base.amplitude_reversals != 0 ||
        intermediate.ep_half_last != 0.0 ||
        riva_ib_mapping_gate(&parameters, &intermediate) > 1.0e-14) {
        std::cerr << "dynamic activation fabricated intermediate-bias history\n";
        return 1;
    }
    if (!driveCycles(parameters, material, 0.01, 2, intermediate) ||
        intermediate.base.amplitude_reversals < 2 ||
        intermediate.ep_half_last <= parameters.bias_reversible_volume_ep_ref) {
        std::cerr << "base backbone did not commit completed-half-cycle activity\n";
        return 1;
    }

    riva_ib_state_t mapping = biasedState(parameters, material, 0.662, 15.0);
    if (riva_ib_mapping_gate(&parameters, &mapping) <= 1.0e-14 ||
        !driveCycles(parameters, material, 0.01, 2, mapping) ||
        mapping.ep_half_last <= parameters.bias_reversible_volume_ep_ref) {
        std::cerr << "mapping backbone did not commit completed-half-cycle activity\n";
        return 1;
    }

    double values[RIVA_IB_STATE_VALUE_COUNT] = {};
    if (riva_ib_state_values(&mapping, values) != RIVA_IB_STATE_VALUE_COUNT ||
        values[49] != mapping.base.pressure_anchor ||
        values[RIVA_IB_STATE_VALUE_COUNT-1] != mapping.ep_half_last) {
        std::cerr << "state-vector append contract is inconsistent\n";
        return 1;
    }

    /* A raw checkpoint copy is sufficient here to prove that the new memory
     * participates in the complete native point-state contract. */
    riva_ib_state_t baseline = mapping;
    riva_ib_state_t restarted = {};
    std::array<unsigned char, sizeof(riva_ib_state_t)> bytes = {};
    std::memcpy(bytes.data(), &mapping, sizeof(mapping));
    std::memcpy(&restarted, bytes.data(), sizeof(restarted));
    const std::array<double, 6> tail = {0.001, -0.002, 0.0015,
                                       -0.001, 0.0005, 0.0002};
    for (double increment : tail) {
        if (!advance(parameters, material, increment, baseline) ||
            !advance(parameters, material, increment, restarted))
            return 1;
    }
    if (!sameState(baseline, restarted)) {
        std::cerr << "checkpoint/restart changed the research response\n";
        return 1;
    }

    std::cout << "PASS: plastic activity captured by both backbones; no "
                 "fictitious activation history; state append and restart "
                 "contracts are stable\n";
    return 0;
}
