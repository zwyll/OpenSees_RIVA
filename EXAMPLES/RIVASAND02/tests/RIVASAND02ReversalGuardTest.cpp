// Synthetic eligibility and transaction contracts for the type-3 research guard.
#include "../../../SRC/material/nD/RIVASAND02/RIVASAND02Kernel.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <cstdlib>
using namespace riva_ib_native;

static void require(bool value, const char *message) {
    if (!value) { std::cerr << message << '\n'; std::exit(1); }
}
static bool same(const riva_ib_state_t &a, const riva_ib_state_t &b) {
    double av[RIVA_IB_STATE_VALUE_COUNT], bv[RIVA_IB_STATE_VALUE_COUNT];
    riva_ib_state_values(&a,av); riva_ib_state_values(&b,bv);
    for (int i=0;i<RIVA_IB_STATE_VALUE_COUNT;++i) if (av[i]!=bv[i]) return false;
    return a.base.initialized==b.base.initialized &&
        a.base.geostatic_admitted==b.base.geostatic_admitted;
}
int main() {
    const auto p=riva_ib_reference_parameters(1.0);
    const auto m=riva_reference_material_parameters(&p.base);
    riva_ib_state_t initial={};
    const tensor_t stress={-100,-100,-100,0,0,0};
    require(riva_ib_initialize_material(&p,&m,stress,.601,&initial),"initialize");
    require(riva_ib_begin_dynamic_phase(&p,&m,&stress,&initial),"activate");

    // Arrange an unambiguous reversal so the eligibility boundaries can be
    // tested independently of a solver's choice of strain history.
    const tensor_t deps={0,0,0,0,0,1e-6};
    const auto predictor=riva_ib_host_stress_ratio_predictor(&p,&m,&initial,deps);
    auto boundary=initial;
    boundary.base.alpha0=riva_add(predictor,{.1,-.1,0,0,0,0});
    boundary.base.alpha=riva_sub(predictor,{.1,-.1,0,0,0,0});
    boundary.base.last_reversal_deviator={.01,-.01,0,0,0,0};
    const double threshold=riva_norm(boundary.base.last_reversal_deviator)/100.0;
    require(riva_ib_host_stress_ratio_reversal(&p,&m,&boundary,deps)==1,"literal boundary event");
    require(riva_ib_host_stress_ratio_reversal(&p,&m,&boundary,deps,threshold*.999)==1,"above excursion threshold");
    require(riva_ib_host_stress_ratio_reversal(&p,&m,&boundary,deps,threshold*1.001)==0,"below excursion threshold");
    require(riva_ib_host_stress_ratio_reversal(&p,&m,&boundary,riva_zero(),threshold)==0,"zero increment");
    require(riva_ib_host_stress_ratio_reversal(&p,&m,&boundary,{0,0,0,0,0,1e-14},threshold*.999)==0,"strain deadband");
    auto inactive=boundary; inactive.base.cyclic_phase_active=0;
    require(riva_ib_host_stress_ratio_reversal(&p,&m,&inactive,deps,threshold*.999)==0,"inactive guard event");

    for (int nsub: {1,4,16,40}) {
        for (double amplitude: {1e-10,.003}) {
            auto off=initial, explicitOff=initial, on=initial;
            double previous=0;
            for (int i=1;i<=256;++i) {
                const double gamma=amplitude*std::sin(2*std::acos(-1.0)*i/32);
                const tensor_t increment={0,0,0,0,0,.5*(gamma-previous)};
                require(riva_ib_update_material_reversal_ex(&p,&m,increment,nsub,&off,nullptr,nullptr,-1,3),"off integration");
                require(riva_ib_update_material_reversal_ex(&p,&m,increment,nsub,&explicitOff,nullptr,nullptr,-1,3,0),"zero integration");
                require(same(off,explicitOff),"explicit zero changes state");
                auto rejected=on, accepted=on;
                require(riva_ib_update_material_reversal_ex(&p,&m,{0,0,0,0,0,-.001},nsub,&rejected,nullptr,nullptr,-1,3,1e-4),"speculative trial");
                const int before=on.base.reversals;
                riva_update_info_t info={};
                require(riva_ib_update_material_reversal_ex(&p,&m,increment,nsub,&accepted,nullptr,&info,-1,3,1e-4),"guard integration");
                require(info.accepted_substeps==nsub && accepted.base.reversals-before>=0 && accepted.base.reversals-before<=1,"host event count");
                auto repeated=on;
                require(riva_ib_update_material_reversal_ex(&p,&m,increment,nsub,&repeated,nullptr,nullptr,-1,3,1e-4),"repeat");
                require(same(accepted,repeated),"trial contamination");
                on=accepted; previous=gamma;
            }
            if (amplitude<1e-8) require(on.base.reversals==0 && off.base.reversals>0,"tiny cycle suppression");
            else require(on.base.reversals>0,"physical cycles suppressed");
        }
    }
    for (int type: {1,2,3}) for (double guard: {-1.0,1e-4,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
        if (type==3 && guard==1e-4) continue;
        auto state=initial;
        require(!riva_ib_update_material_reversal_ex(&p,&m,deps,4,&state,nullptr,nullptr,-1,type,guard),"invalid guard accepted");
        require(same(state,initial),"invalid guard mutated state");
    }
    std::cout << "REVERSAL_GUARD_KERNEL_COMPLETE\n";
}
