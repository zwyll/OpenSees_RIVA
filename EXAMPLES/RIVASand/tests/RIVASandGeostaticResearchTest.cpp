#include <cmath>
#include <iostream>

#include "../../../SRC/material/nD/RIVASand/RIVASandKernel.h"

namespace {

bool close(double a, double b, double tolerance = 1.0e-10)
{
    return std::fabs(a-b) <= tolerance*(1.0+std::fabs(a)+std::fabs(b));
}

bool sameTensor(const riva_tensor_t &a, const riva_tensor_t &b)
{
    return close(a.xx,b.xx) && close(a.yy,b.yy) && close(a.zz,b.zz) &&
           close(a.xy,b.xy) && close(a.yz,b.yz) && close(a.xz,b.xz);
}

} // namespace

int main()
{
    riva_parameters_t frozen = riva_reference_parameters(1.0);
    const riva_material_parameters_t material =
        riva_reference_material_parameters(&frozen);
    const double voidRatio =
        material.e_max-0.641*(material.e_max-material.e_min);

    const riva_tensor_t shallow = {-0.20,-0.20,-0.20,0.0,0.0,0.0};
    riva_state_t state = {};
    if (!riva_initialize_material(&frozen,&material,shallow,voidRatio,&state)) {
        std::cerr << "FAIL: frozen positive-pressure initialization\n";
        return 1;
    }

    const riva_tensor_t tensile = {0.50,0.50,0.50,0.0,0.0,0.0};
    state = riva_state_t{};
    if (riva_initialize_material(&frozen,&material,tensile,voidRatio,&state)) {
        std::cerr << "FAIL: frozen kernel admitted tensile pressure\n";
        return 1;
    }

    riva_parameters_t research = frozen;
    research.p_residual = 1.013;
    research.geostatic_admission_enabled = 1;
    state = riva_state_t{};
    if (!riva_initialize_material(
            &research,&material,tensile,voidRatio,&state)) {
        std::cerr << "FAIL: translated cone did not admit p' > -pResidual\n";
        return 1;
    }
    if (!close(state.pressure_anchor,0.513) ||
        !close(riva_pressure(state.stress),-0.50)) {
        std::cerr << "FAIL: cone and physical pressures were not separated\n";
        return 1;
    }

    riva_tensor_t stressNew = {};
    riva_update_info_t info = {};
    if (!riva_update_material(&research,&material,riva_zero(),1,
                              &state,&stressNew,&info) ||
        !sameTensor(stressNew,tensile)) {
        std::cerr << "FAIL: zero increment changed translated physical stress\n";
        return 1;
    }

    const riva_tensor_t biased = {-2.0,-2.0,-2.0,0.0,0.0,3.0};
    riva_state_t baseState = {}, shiftedState = {};
    if (!riva_initialize_material(&frozen,&material,biased,voidRatio,&baseState) ||
        !riva_initialize_material(&research,&material,biased,voidRatio,
                                  &shiftedState)) {
        std::cerr << "FAIL: biased initialization\n";
        return 1;
    }
    if (!(riva_norm(shiftedState.alpha) < riva_norm(baseState.alpha))) {
        std::cerr << "FAIL: pResidual did not reduce low-pressure stress ratio\n";
        return 1;
    }

    shiftedState.bias_reversible_volume = 1.0e-4;
    const riva_state_t mechanical = riva_without_reversible_bias(
        &research,&material,&shiftedState);
    int reconstructionFloor = 0;
    if (!close(riva_cone_pressure(&research,mechanical.stress),
               riva_pressure_from_confining(
                   &research,&material,mechanical.eps_v_confining,
                   mechanical.pressure_anchor,&reconstructionFloor))) {
        std::cerr << "FAIL: reversible-volume reconstruction lost cone translation\n";
        return 1;
    }
    shiftedState.bias_reversible_volume = 0.0;

    /* Mirror the adapter's stress-preserving activation recentering and
     * verify that the research cap does not snap an inherited over-bound
     * stress during a zero increment. */
    shiftedState.alpha0 = shiftedState.alpha;
    shiftedState.alpha01 = shiftedState.alpha;
    stressNew = riva_zero();
    info = riva_update_info_t{};
    if (!riva_update_material(&research,&material,riva_zero(),1,
                              &shiftedState,&stressNew,&info) ||
        !sameTensor(stressNew,biased)) {
        std::cerr << "FAIL: geostatic admission changed equilibrated stress\n";
        return 1;
    }

    std::cout << "PASS: residual-pressure translation and stress-preserving "
                 "geostatic admission\n";
    return 0;
}
