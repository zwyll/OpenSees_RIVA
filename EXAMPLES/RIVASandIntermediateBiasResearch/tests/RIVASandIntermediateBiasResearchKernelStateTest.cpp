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
    if (RIVA_IB_PARAMETER_COUNT != 257 || RIVA_IB_LOGICAL_STATE_COUNT != 64 ||
        RIVA_IB_STATE_VALUE_COUNT != 139 || RIVA_IB_KERNEL_REVISION != 4u) {
        std::cerr << "research kernel contract constants are inconsistent\n";
        return 1;
    }
    const riva_ib_parameters_t parameters = riva_ib_reference_parameters(1.0);
    const riva_material_parameters_t material =
        riva_reference_material_parameters(&parameters.base);
    if (!(parameters.bias_reversible_volume_ep_ref > 0.0)) return 1;

    /* The field correction is bounded and acts on only the inherited mean
     * term in its low-bias/high-confinement window. It is disabled by
     * default and therefore cannot alter the accepted research histories. */
    riva_ib_state_t field = biasedState(parameters, material, 0.601, 4.0);
    field.base.cyclic_direction = {0, 0, 0, 0, 0, 1};
    field.base.static_bias_tensor = {
        0, 0, 0, 0, 0, 0.19/std::sqrt(2.0)};
    field.base.static_bias_index = 0.19;
    field.base.pressure_anchor =
        parameters.base.bias_reversible_mean_transition_pressure*1.25;
    field.base.amplitude_reversals = 6;
    field.base.cyclic_amplitude = 0.2;
    field.ep_half_last = 2.0*parameters.bias_reversible_volume_ep_ref;
    field.phase_accumulation_lambda_anchor = field.base.lambda_total;
    riva_ib_parameters_t fieldOff = parameters;
    fieldOff.phase_transformation_enabled = 0;
    riva_ib_parameters_t fieldOn = fieldOff;
    fieldOn.field_bias_mean_correction_enabled = 1;
    const double offTarget = riva_ib_bias_volume_target(&fieldOff, &field, 0);
    const double onTarget = riva_ib_bias_volume_target(&fieldOn, &field, 0);
    if (onTarget != offTarget) {
        std::cerr << "field correction did not start continuously from zero\n";
        return 1;
    }
    field.field_bias_mean_activity = 1.0;
    const double activeTarget = riva_ib_bias_volume_target(&fieldOn, &field, 0);
    const double fieldBias = riva_ib_projected_bias(&field);
    const double meanBuildup = 1.0-std::exp(-6.0/
        parameters.base.bias_reversible_mean_buildup_reversals);
    const double expectedMean = parameters.base.bias_reversible_mean_scale*
        (parameters.base.bias_reversible_mean_transition_pressure-
         field.base.pressure_anchor)/parameters.base.bias_reference_pressure*
        (fieldBias/parameters.base.bias_reversible_volume_reference_bias)*
        meanBuildup;
    if (std::fabs((offTarget-activeTarget)-expectedMean) > 1.0e-16 ||
        std::fabs(offTarget-activeTarget) > std::fabs(expectedMean)+1.0e-16) {
        std::cerr << "field correction exceeded its bounded mean adjustment\n";
        return 1;
    }
    riva_ib_state_t dilativeField = field;
    dilativeField.base.pressure_anchor =
        parameters.base.bias_reversible_mean_transition_pressure*0.65;
    if (riva_ib_bias_volume_target(&fieldOff, &dilativeField, 0) !=
        riva_ib_bias_volume_target(&fieldOn, &dilativeField, 0)) {
        std::cerr << "field correction suppressed a helpful dilative mean shift\n";
        return 1;
    }
    riva_ib_state_t calibrated = field;
    calibrated.base.static_bias_tensor = {
        0, 0, 0, 0, 0, 0.29/std::sqrt(2.0)};
    calibrated.base.static_bias_index = 0.29;
    if (riva_ib_bias_volume_target(&fieldOff, &calibrated, 0) !=
        riva_ib_bias_volume_target(&fieldOn, &calibrated, 0)) {
        std::cerr << "field correction leaked outside its low-bias window\n";
        return 1;
    }

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
        values[RIVA_IB_STATE_VALUE_COUNT-2] != mapping.ep_half_last ||
        values[RIVA_IB_STATE_VALUE_COUNT-1] !=
            mapping.field_bias_mean_activity) {
        std::cerr << "state-vector append contract is inconsistent\n";
        return 1;
    }

    riva_ib_parameters_t activityParameters = parameters;
    activityParameters.field_bias_mean_correction_enabled = 1;
    riva_ib_state_t activity = biasedState(
        activityParameters, material, 0.601, 4.0);
    if (!advance(activityParameters, material, 0.002, activity) ||
        !(activity.field_bias_mean_activity > 0.0) ||
        !(activity.field_bias_mean_activity <= 1.0)) {
        std::cerr << "field correction activity did not evolve monotonically\n";
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

    /* The extended entry point must preserve the original detector when the
     * adapter supplies no override, and a forced reversal must be registered
     * exactly once rather than once per constitutive substep. */
    const tensor_t latchIncrement = {0, 0, 0, 0, 0, -0.0005};
    riva_ib_state_t automatic = mapping;
    riva_ib_state_t original = mapping;
    riva_update_info_t automaticInfo = {};
    riva_update_info_t originalInfo = {};
    if (!riva_ib_update_material_ex(&parameters, &material, latchIncrement, 4,
                                    &automatic, nullptr, &automaticInfo, -1) ||
        !riva_ib_update_material(&parameters, &material, latchIncrement, 4,
                                 &original, nullptr, &originalInfo) ||
        !sameState(automatic, original) ||
        automaticInfo.reversal_registered !=
            originalInfo.reversal_registered) {
        std::cerr << "automatic reversal override changed the original path\n";
        return 1;
    }

    riva_ib_state_t suppressed = mapping;
    riva_update_info_t suppressedInfo = {};
    if (!riva_ib_update_material_ex(&parameters, &material, latchIncrement, 4,
                                    &suppressed, nullptr, &suppressedInfo, 0) ||
        suppressedInfo.reversal_registered != 0 ||
        suppressed.base.reversals != mapping.base.reversals) {
        std::cerr << "forced no-reversal decision was not respected\n";
        return 1;
    }

    riva_ib_state_t forced = mapping;
    riva_ib_state_t repeatedTrial = mapping;
    riva_update_info_t forcedInfo = {};
    riva_update_info_t repeatedInfo = {};
    if (!riva_ib_update_material_ex(&parameters, &material, latchIncrement, 4,
                                    &forced, nullptr, &forcedInfo, 1) ||
        !riva_ib_update_material_ex(&parameters, &material, latchIncrement, 4,
                                    &repeatedTrial, nullptr, &repeatedInfo, 1) ||
        forcedInfo.reversal_registered != 1 ||
        forced.base.reversals != mapping.base.reversals + 1 ||
        !sameState(forced, repeatedTrial)) {
        std::cerr << "latched reversal was not deterministic and single-event\n";
        return 1;
    }

    std::cout << "PASS: bounded field-bias mean correction; plastic activity "
                 "captured by both backbones; no "
                 "fictitious activation history; state append, restart, and "
                 "reversal-override contracts are stable\n";
    return 0;
}
