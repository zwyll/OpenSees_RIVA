/* Dependency-free verification of the opt-in OpenSees geostatic-admission
 * state. The frozen constitutive parameters and ordinary V8 paths are not
 * changed by this numerical stage-transition safeguard. */
#include "../../../SRC/material/nD/RIVASand/RIVASandKernel.h"

#include <array>
#include <cmath>
#include <cstring>
#include <iostream>

namespace {

bool close(double a, double b, double tolerance = 1.0e-12)
{
    return std::fabs(a-b) <= tolerance*(1.0+std::fabs(a)+std::fabs(b));
}

bool sameTensor(const riva_tensor_t &a, const riva_tensor_t &b)
{
    return close(a.xx,b.xx) && close(a.yy,b.yy) && close(a.zz,b.zz) &&
           close(a.xy,b.xy) && close(a.yz,b.yz) && close(a.xz,b.xz);
}

bool sameState(const riva_state_t &a, const riva_state_t &b)
{
    double av[RIVA_STATE_VALUE_COUNT] = {};
    double bv[RIVA_STATE_VALUE_COUNT] = {};
    if (riva_state_values(&a,av) != RIVA_STATE_VALUE_COUNT ||
        riva_state_values(&b,bv) != RIVA_STATE_VALUE_COUNT)
        return false;
    for (int i=0;i<RIVA_STATE_VALUE_COUNT;++i)
        if (av[i] != bv[i]) return false;
    return a.initialized == b.initialized &&
           a.geostatic_admission_radius == b.geostatic_admission_radius &&
           a.geostatic_admitted == b.geostatic_admitted;
}

double boundingRadius(const riva_parameters_t &p,
                      const riva_material_parameters_t &material,
                      const riva_state_t &state)
{
    double mb=0.0,md=0.0,xi=0.0;
    riva_surfaces(&p,&material,riva_pressure(state.stress),state.void_ratio,
                  &mb,&md,&xi);
    return std::sqrt(2.0/3.0)*mb;
}

bool advance(const riva_parameters_t &p,
             const riva_material_parameters_t &material,
             riva_tensor_t increment,riva_state_t &state)
{
    riva_tensor_t stress=state.stress;
    riva_update_info_t info={};
    return riva_update_material(&p,&material,increment,1,&state,&stress,&info) &&
           riva_finite_tensor(stress);
}

} // namespace

int main()
{
    const riva_parameters_t p=riva_reference_parameters(1.0);
    const riva_material_parameters_t material=
        riva_reference_material_parameters(&p);
    const double voidRatio=material.e_max-
        0.641*(material.e_max-material.e_min);

    /* Enabling admission on an already admissible state must be an exact
     * no-op, including all subsequent constitutive history. */
    const riva_tensor_t ordinaryStress={-19.4,-19.4,-40.0,0.0,0.0,0.0};
    riva_state_t ordinary={},admissionEnabled={};
    if (!riva_initialize_material(&p,&material,ordinaryStress,voidRatio,
                                  &ordinary) ||
        !riva_initialize_material(&p,&material,ordinaryStress,voidRatio,
                                  &admissionEnabled) ||
        !riva_admit_geostatic_state(&p,&material,&admissionEnabled) ||
        admissionEnabled.geostatic_admitted ||
        !sameState(ordinary,admissionEnabled)) {
        std::cerr << "FAIL: admissible-state equivalence at activation\n";
        return 1;
    }
    if (!riva_begin_dynamic_phase(&p,ordinaryStress,&ordinary) ||
        !riva_begin_dynamic_phase(&p,ordinaryStress,&admissionEnabled))
        return 1;
    const std::array<riva_tensor_t,4> ordinaryPath={{
        {0,0,0,2.0e-4,0,0},
        {0,0,0,-4.0e-4,0,0},
        {-1.0e-5,-1.0e-5,-1.0e-5,2.0e-4,0,0},
        {1.0e-5,1.0e-5,1.0e-5,0,0,0}
    }};
    for (const riva_tensor_t increment : ordinaryPath) {
        if (!advance(p,material,increment,ordinary) ||
            !advance(p,material,increment,admissionEnabled) ||
            !sameState(ordinary,admissionEnabled)) {
            std::cerr << "FAIL: admission changed an admissible V8 path\n";
            return 1;
        }
    }

    /* Admission is forbidden at or below the positive constitutive floor;
     * no residual-pressure translation or tensile cone is available. */
    const riva_tensor_t belowFloor={-5.0e-4,-5.0e-4,-5.0e-4,0,0,0};
    riva_state_t shallow={};
    if (!riva_initialize_material(&p,&material,belowFloor,voidRatio,&shallow) ||
        riva_admit_geostatic_state(&p,&material,&shallow)) {
        std::cerr << "FAIL: noncompressive/floor admission guard\n";
        return 1;
    }

    /* This stress is compressive in every normal component but lies well
     * outside the current cone. Admission must preserve it exactly. */
    const riva_tensor_t overboundStress={-2.0,-2.0,-56.0,0.0,0.0,0.0};
    riva_state_t admitted={};
    if (!riva_initialize_material(&p,&material,overboundStress,voidRatio,
                                  &admitted)) return 1;
    const riva_state_t beforeAdmission=admitted;
    if (!riva_admit_geostatic_state(&p,&material,&admitted) ||
        !admitted.geostatic_admitted ||
        !sameTensor(admitted.stress,beforeAdmission.stress) ||
        !sameTensor(admitted.alpha0,admitted.alpha) ||
        !sameTensor(admitted.alpha01,admitted.alpha) ||
        !close(admitted.geostatic_admission_radius,
               riva_norm(admitted.alpha))) {
        std::cerr << "FAIL: stress-preserving over-bound admission\n";
        return 1;
    }

    /* A finite outward candidate and 100 repeated outward candidates may not
     * grow the admitted normalized radius by even one smoothing band. */
    riva_state_t radial=admitted;
    const double bound=boundingRadius(p,material,radial);
    const double initialRadius=riva_norm(radial.alpha);
    for (int step=0;step<100;++step) {
        const double previousRadius=riva_norm(radial.alpha);
        const riva_tensor_t candidate=riva_scale(radial.alpha,1.25);
        radial.alpha=riva_geostatic_admission_cap(
            bound,&radial,candidate,&radial.geostatic_admitted,
            &radial.geostatic_admission_radius);
        if (!radial.geostatic_admitted ||
            riva_norm(radial.alpha)>previousRadius*(1.0+1.0e-13)) {
            std::cerr << "FAIL: repeated outward admission growth at step "
                      << step << '\n';
            return 1;
        }
    }
    if (!close(riva_norm(radial.alpha),initialRadius)) {
        std::cerr << "FAIL: finite outward cap changed admitted radius\n";
        return 1;
    }

    /* Exercise the cap through the complete stress update as well. */
    if (!riva_begin_dynamic_phase(&p,overboundStress,&admitted)) return 1;
    const riva_tensor_t outwardIncrement={2.0e-4,2.0e-4,-4.0e-4,0,0,0};
    for (int step=0;step<20 && admitted.geostatic_admitted;++step) {
        const double previousRadius=riva_norm(admitted.alpha);
        if (!advance(p,material,outwardIncrement,admitted) ||
            riva_norm(admitted.alpha)>previousRadius*(1.0+1.0e-12)) {
            std::cerr << "FAIL: full update expanded admitted radius\n";
            return 1;
        }
    }

    /* Restart while admission is active, then require identical continuation. */
    riva_state_t restarted=admitted;
    std::array<unsigned char,sizeof(riva_state_t)> bytes={};
    std::memcpy(bytes.data(),&restarted,sizeof(restarted));
    restarted={};
    std::memcpy(&restarted,bytes.data(),sizeof(restarted));
    riva_state_t uninterrupted=admitted;
    if (!advance(p,material,outwardIncrement,restarted) ||
        !advance(p,material,outwardIncrement,uninterrupted) ||
        !sameState(restarted,uninterrupted)) {
        std::cerr << "FAIL: admitted-state restart\n";
        return 1;
    }

    /* An inward candidate that reaches the cone disables admission. A later
     * outward candidate then receives the ordinary hard V8 cone projection. */
    riva_state_t reentry=radial;
    riva_tensor_t direction=riva_scale(reentry.alpha,1.0/riva_norm(reentry.alpha));
    riva_tensor_t inside=riva_scale(direction,0.5*bound);
    reentry.alpha=riva_geostatic_admission_cap(
        bound,&reentry,inside,&reentry.geostatic_admitted,
        &reentry.geostatic_admission_radius);
    if (reentry.geostatic_admitted ||
        reentry.geostatic_admission_radius != 0.0 ||
        !(riva_norm(reentry.alpha)<bound)) {
        std::cerr << "FAIL: inward re-entry did not disable admission\n";
        return 1;
    }
    const riva_tensor_t outside=riva_scale(direction,1.25*bound);
    reentry.alpha=riva_geostatic_admission_cap(
        bound,&reentry,outside,&reentry.geostatic_admitted,
        &reentry.geostatic_admission_radius);
    if (reentry.geostatic_admitted ||
        !close(riva_norm(reentry.alpha),bound)) {
        std::cerr << "FAIL: post-reentry ordinary cone projection\n";
        return 1;
    }

    std::cout << "PASS: compressive stress-preserving geostatic admission, "
                 "non-expansive radius, re-entry, restart, and frozen-path "
                 "equivalence\n";
    return 0;
}
