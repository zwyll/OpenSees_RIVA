// Raw-kernel probe of the reverse-push response from a static-bias state.
// Usage: RIVASandBiasProbe [bias=0.4] [sigv=100] [dir=-1] [zeta=0.05]
//        [h=380] [kd=1.125] [M=1.25] [Dr=0.5] [nsub=20] [steps=40] [dg=2.5e-5]
// Initializes exactly like the adapter (riva_initialize_material +
// riva_begin_dynamic_phase with shear-free reference => committed shear is the
// static bias), then feeds pure shear strain increments and prints the path.
#include <cstdio>
#include <cstdlib>
#include "../../../SRC/material/nD/RIVASand/RIVASandKernel.h"

int main(int argc, char **argv)
{
    double bias = argc > 1 ? atof(argv[1]) : 0.4;
    double sigv = argc > 2 ? atof(argv[2]) : 100.0;
    double dir  = argc > 3 ? atof(argv[3]) : -1.0;
    double zeta = argc > 4 ? atof(argv[4]) : 0.05;
    double h    = argc > 5 ? atof(argv[5]) : 380.0;
    double kd   = argc > 6 ? atof(argv[6]) : 1.125;
    double M    = argc > 7 ? atof(argv[7]) : 1.25;
    double Dr   = argc > 8 ? atof(argv[8]) : 0.5;
    int nsub    = argc > 9 ? atoi(argv[9]) : 20;
    int steps   = argc > 10 ? atoi(argv[10]) : 40;
    double dg   = argc > 11 ? atof(argv[11]) : 2.5e-5;
    int admit   = argc > 12 ? atoi(argv[12]) : 0;

    riva_parameters_t p = riva_reference_parameters(1.0);
    p.zeta = zeta; p.h = h; p.kd = kd; p.M = M;
    p.admit_inherited_overbound = admit;
    riva_material_parameters_t material = riva_reference_material_parameters(&p);
    const double e0 = p.e_max - Dr*(p.e_max - p.e_min);

    const double k0 = 3.0/7.0;   // elastic K0 of the equilibrium test
    riva_tensor_t stress = {-k0*sigv, -k0*sigv, -sigv, 0.0, 0.0, bias*sigv};
    riva_state_t state = {};
    if (!riva_initialize_material(&p, &material, stress, e0, &state)) {
        printf("INIT FAILED\n"); return 1;
    }
    riva_tensor_t reference = stress;
    reference.xy = reference.yz = reference.xz = 0.0;
    if (!riva_begin_dynamic_phase(&p, reference, &state)) {
        printf("BEGIN FAILED\n"); return 1;
    }
    printf("# bias=%g sigv=%g dir=%g zeta=%g h=%g kd=%g M=%g Dr=%g nsub=%d\n",
           bias, sigv, dir, zeta, h, kd, M, Dr, nsub);
    {
        const riva_tensor_t s0 = riva_dev(state.stress);
        const double p0 = -(state.stress.xx+state.stress.yy+state.stress.zz)/3.0;
        const double a0 = riva_norm(riva_scale(s0, 1.0/p0));
        double mb, md, xi;
        riva_surfaces(&p, &material, p0, state.void_ratio, &mb, &md, &xi);
        printf("# init: p=%g |alpha|=%g bound=sqrt(2/3)*mb=%g (mb=%g) %s\n",
               p0, a0, sqrt(2.0/3.0)*mb, mb,
               a0 > sqrt(2.0/3.0)*mb ? "OVER-BOUND -> first call will snap-project"
                                     : "inside bound");
    }
    printf("%10s %12s %12s %12s %6s %4s %4s\n",
           "gamma", "tau_xz", "p", "D", "rev", "ok", "amp");
    double gamma = 0.0;
    for (int i = 0; i < steps; ++i) {
        riva_tensor_t deps = {0,0,0,0,0, dir*dg*0.5};  // tensor shear = gamma/2
        riva_tensor_t snew = {};
        riva_update_info_t info = {};
        int ok = riva_update_material(&p, &material, deps, nsub, &state,
                                      &snew, &info);
        gamma += dir*dg;
        const double pr = -(snew.xx+snew.yy+snew.zz)/3.0;
        printf("%10.3e %12.5g %12.5g %12.5g %6d %4d %4g\n",
               gamma, snew.xz, pr, state.D,
               (int)info.reversal_registered, ok, state.amplitude_factor);
        if (!ok) { printf("KERNEL REJECT at step %d\n", i+1); return 2; }
    }
    return 0;
}
