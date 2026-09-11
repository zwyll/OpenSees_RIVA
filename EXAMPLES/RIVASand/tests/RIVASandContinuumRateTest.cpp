// Standalone check of the smooth backbone rate, not a claim that this is
// the algorithmic derivative of the complete substepped/overlay update.
#include "../../../SRC/material/nD/RIVASand/RIVASandKernel.h"
#include "../../../SRC/material/nD/RIVASand/RIVASandContinuumTangent.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>

static void require(bool yes, const char *message) {
    if (!yes) throw std::runtime_error(message);
}
static riva_tensor_t engineeringBasis(int i) {
    riva_tensor_t t={};
    if(i==0)t.xx=1; if(i==1)t.yy=1; if(i==2)t.zz=1;
    if(i==3)t.xy=.5; if(i==4)t.yz=.5; if(i==5)t.xz=.5;
    return t;
}
int main() {
    for(double scale: {1.,1000.}) for(double d: {-.12,0.,.2}) {
        auto p=riva_reference_parameters(scale);
        auto m=riva_reference_material_parameters(&p);
        riva_tensor_t n={.4,-.3,-.1,.2,.3,-.15};
        n=riva_scale(n,1/riva_norm(n));
        const double pressure=100*scale;
        const auto stress=riva_sub(riva_scale(n,.3*pressure),riva_iso(pressure));
        riva_state_t s={};
        require(riva_initialize_material(&p,&m,stress,.601,&s),"initialize");
        s.n=n; s.beta=1.7; s.D=d; s.D_ir=std::max(0.,d);
        s.D_re=std::min(0.,d);
        double g,k; riva_moduli_for_state(&p,&m,pressure,&s,&g,&k);
        double h=pressure*riva_hardening_for_state(&p,&m,pressure,&s)*std::pow(s.beta,m.m);
        double c[6][6];
        require(riva_continuum::build(g,k,n,n,riva_ddot(s.alpha,n),d,
            (2./3.)*h,p.denominator_floor_ratio,c),"build");
        double last=0;
        for(double eps: {1.e-7,1.e-8,1.e-9}) {
            double maxError=0;
            // Differentiate about a positive loading direction so all six
            // central perturbations stay on the same active branch.
            for(int j=0;j<6;++j) {
                const auto center=riva_scale(n,eps);
                const auto delta=riva_scale(engineeringBasis(j),eps*.01);
                const auto plus=riva_backbone_forward_euler(&p,&m,&s,riva_add(center,delta),0,0);
                const auto minus=riva_backbone_forward_euler(&p,&m,&s,riva_sub(center,delta),0,0);
                require(plus.lambda_total>s.lambda_total && minus.lambda_total>s.lambda_total,"active branch");
                require(plus.pressure_floor_hits==s.pressure_floor_hits &&
                    minus.denominator_floor_hits==s.denominator_floor_hits,"smooth neighborhood");
                double ds[6];
                riva_continuum::components(riva_scale(riva_sub(plus.stress,minus.stress),1/(.02*eps)),ds);
                for(int i=0;i<6;++i)maxError=std::max(maxError,std::abs(ds[i]-c[i][j])/k);
            }
            if(eps==1.e-9) require(maxError<2.e-5,"continuum rate/FD mismatch");
            last=maxError;
        }
        // Nonassociated flow must retain the nonsymmetry, including the
        // volumetric-to-shear coupling and all engineering shear columns.
        if(d!=-.3)require(std::abs(c[0][3]-c[3][0])>1.e-3*scale,"lost nonsymmetry");
        double untouched[6][6]={{123}};
        require(!riva_continuum::build(g,k,riva_zero(),n,.3,d,h,.05,untouched),"zero normal");
        require(untouched[0][0]==123,"failed build modified output");
        require(!riva_continuum::build(g,k,n,n,.3,d,-1.e30,.05,untouched),"invalid denominator");
        std::cout<<"PASS scale="<<scale<<" D="<<d<<" relative FD error="<<last<<'\n';
    }
}
