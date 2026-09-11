// Read-only, frozen-state continuum backbone operator. No stress integration.
#ifndef RIVASandContinuumTangent_h
#define RIVASandContinuumTangent_h

#include <cmath>

namespace riva_continuum {

template<class Tensor>
inline void components(const Tensor &t, double v[6])
{
    v[0]=t.xx; v[1]=t.yy; v[2]=t.zz;
    v[3]=t.xy; v[4]=t.yz; v[5]=t.xz;
}

// OpenSees uses engineering shear strain and physical shear stress. Hence
// BOTH outer-product vectors use physical tensor components (no factor 2
// on shear entries). The returned operator is generally nonsymmetric.
template<class Tensor>
inline bool build(double shear, double bulk, const Tensor &normal,
                  const Tensor &flow, double alphaNormal, double dilatancy,
                  double hardeningTerm, double denominatorFloorRatio,
                  double result[6][6])
{
    double n[6], f[6]; components(normal,n); components(flow,f);
    double nf=0.0, nn=0.0, ff=0.0;
    for (int i=0;i<6;++i) {
        const double w=i<3?1.0:2.0;
        nf+=w*n[i]*f[i]; nn+=w*n[i]*n[i]; ff+=w*f[i]*f[i];
    }
    const double denominator=2.0*shear*nf+hardeningTerm-
        bulk*dilatancy*alphaNormal;
    if (!std::isfinite(shear) || !std::isfinite(bulk) || shear<=0 || bulk<=0 ||
        !std::isfinite(denominator) ||
        denominator<=denominatorFloorRatio*2.0*shear ||
        std::fabs(nn-1.0)>1.e-6 || std::fabs(ff-1.0)>1.e-6 || nf<=0.05)
        return false;
    double candidate[6][6];
    const double lame=bulk-2.0*shear/3.0;
    for (int i=0;i<6;++i) {
        const double a=2.0*shear*f[i]-(i<3?bulk*dilatancy:0.0);
        for (int j=0;j<6;++j) {
            const double b=2.0*shear*n[j]+(j<3?bulk*alphaNormal:0.0);
            const double elastic=(i<3 && j<3?lame:0.0)+
                (i==j?(i<3?2.0*shear:shear):0.0);
            candidate[i][j]=elastic-a*(b/denominator);
            if (!std::isfinite(candidate[i][j])) return false;
        }
    }
    for (int i=0;i<6;++i)
        for (int j=0;j<6;++j) result[i][j]=candidate[i][j];
    return true;
}

// Use committed/trial plastic activity, not a new reversal decision. Avoid
// differentiating steps containing a discrete event, clamp or projection.
template<class State>
inline bool smoothPlasticStep(const State &old, const State &trial)
{
    return trial.lambda_total>old.lambda_total &&
        trial.reversals==old.reversals &&
        trial.pressure_floor_hits==old.pressure_floor_hits &&
        trial.denominator_floor_hits==old.denominator_floor_hits &&
        trial.beta_fallbacks==old.beta_fallbacks;
}

} // namespace riva_continuum
#endif
