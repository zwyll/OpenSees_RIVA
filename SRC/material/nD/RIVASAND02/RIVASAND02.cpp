/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
** ****************************************************************** */

#include "RIVASAND02.h"
#include "../RIVASand/RIVASandContinuumTangent.h"

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <Parameter.h>
#include <classTags.h>
#include <elementAPI.h>

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>

using namespace riva_ib_native;

namespace {

const int RIVASerializedSize = 183;
const int RIVAAdapterRevision = 2; // tangent selector/activity in configuration bits

enum RIVAAdapterConfigurationFlag {
    RIVAGeostaticAdmissionFlag = 1,
    RIVAReversalLatchFlag = 2,
    RIVAFieldBiasMeanCorrectionFlag = 4,
    RIVANoBiasVolumeFlag = 8,
    RIVAContinuumTangentFlag = 16,
    RIVAPlasticLoadingFlag = 32
};

bool finiteVector(const Vector &value)
{
    for (int i = 0; i < value.Size(); ++i)
        if (!std::isfinite(value(i))) return false;
    return true;
}

} // namespace

void *
OPS_RIVASAND02Material(void)
{
    const int requiredValues = 12;
    if (OPS_GetNumRemainingInputArgs() < requiredValues + 1) {
        opserr << "Want: nDMaterial RIVASAND02 tag Dr G0 M kd h m zeta "
               << "eMax eMin Q R nG <-rho value> <-nSub value> "
               << "<-stressScale value> <-pMin value> "
               << "<-tangentPMin value> <-TanType 0|1> <-pResidual value> "
               << "<-geostaticAdmission> <-reversalLatch> "
               << "<-BiasVolume 0|1|2> <-stage 0|1|2> "
               << "<-initialStress sxx syy szz sxy syz sxz>" << endln;
        return 0;
    }

    int tag = 0;
    int count = 1;
    if (OPS_GetIntInput(&count, &tag) < 0) {
        opserr << "WARNING invalid RIVASAND02 tag" << endln;
        return 0;
    }

    double values[requiredValues];
    count = requiredValues;
    if (OPS_GetDoubleInput(&count, values) < 0) {
        opserr << "WARNING invalid RIVASAND02 material values for tag "
               << tag << endln;
        return 0;
    }

    double rho = 0.0;
    double stressScale = 1.0;
    double pMin = -1.0;
    double tangentPressureFloor = -1.0;
    double residualPressure = 0.0;
    bool geostaticAdmission = false;
    int tangentType = 0;
    bool tangentTypeSpecified = false;
    bool reversalLatch = false;
    int biasVolumeMode = -1; // Unspecified; use mode 0 after parsing.
    int fixedSubsteps = 1;
    int stage = 0;
    bool stageSpecified = false;
    bool initialStressSpecified = false;
    Vector initialStress(6);
    initialStress.Zero();

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *option = OPS_GetString();
        if (std::strcmp(option, "-rho") == 0) {
            count = 1;
            if (OPS_GetDoubleInput(&count, &rho) < 0) {
                opserr << "WARNING invalid -rho for RIVASAND02 tag "
                       << tag << endln;
                return 0;
            }
        } else if (std::strcmp(option, "-nSub") == 0 ||
                   std::strcmp(option, "-noSubsteps") == 0) {
            count = 1;
            if (OPS_GetIntInput(&count, &fixedSubsteps) < 0) {
                opserr << "WARNING invalid -nSub for RIVASAND02 tag "
                       << tag << endln;
                return 0;
            }
        } else if (std::strcmp(option, "-stressScale") == 0) {
            count = 1;
            if (OPS_GetDoubleInput(&count, &stressScale) < 0) {
                opserr << "WARNING invalid -stressScale for RIVASAND02 tag "
                       << tag << endln;
                return 0;
            }
        } else if (std::strcmp(option, "-pMin") == 0) {
            count = 1;
            if (OPS_GetDoubleInput(&count, &pMin) < 0 ||
                !std::isfinite(pMin) || !(pMin > 0.0)) {
                opserr << "WARNING invalid -pMin for RIVASAND02 tag "
                       << tag << "; value must be positive" << endln;
                return 0;
            }
        } else if (std::strcmp(option, "-tangentPMin") == 0) {
            count = 1;
            if (OPS_GetDoubleInput(&count, &tangentPressureFloor) < 0 ||
                !std::isfinite(tangentPressureFloor) ||
                !(tangentPressureFloor > 0.0)) {
                opserr << "WARNING invalid -tangentPMin for "
                       << "RIVASAND02 tag " << tag
                       << "; value must be positive" << endln;
                return 0;
            }
        } else if (std::strcmp(option, "-TanType") == 0) {
            int requested = -1;
            count = 1;
            if (OPS_GetIntInput(&count, &requested) < 0 ||
                (requested != 0 && requested != 1) ||
                (tangentTypeSpecified && requested != tangentType)) {
                opserr << "WARNING RIVASAND02 -TanType requires 0 or 1; "
                       << "conflicting duplicate values are not allowed" << endln;
                return 0;
            }
            tangentType = requested;
            tangentTypeSpecified = true;
        } else if (std::strcmp(option, "-pResidual") == 0) {
            count = 1;
            if (OPS_GetDoubleInput(&count, &residualPressure) < 0 ||
                !std::isfinite(residualPressure) || residualPressure < 0.0) {
                opserr << "WARNING invalid -pResidual for RIVASAND02 tag "
                       << tag << "; value must be nonnegative" << endln;
                return 0;
            }
        } else if (std::strcmp(option, "-geostaticAdmission") == 0) {
            geostaticAdmission = true;
        } else if (std::strcmp(option, "-reversalLatch") == 0) {
            reversalLatch = true;
        } else if (std::strcmp(option, "-BiasVolume") == 0 ||
                   std::strcmp(option, "-fieldBiasVolume") == 0 ||
                   std::strcmp(option, "-noBiasVolume") == 0) {
            int requestedMode = 0;
            if (std::strcmp(option, "-BiasVolume") == 0) {
                count = 1;
                if (OPS_GetIntInput(&count, &requestedMode) < 0 ||
                    requestedMode < 0 || requestedMode > 2) {
                    opserr << "WARNING invalid -BiasVolume for "
                           << "RIVASAND02 tag " << tag
                           << "; expected integer 0 (default), 1 (no bias volume), "
                           << "or 2 (field bias correction)" << endln;
                    return 0;
                }
            } else {
                requestedMode = std::strcmp(option, "-noBiasVolume") == 0 ? 1 : 2;
            }
            if (biasVolumeMode >= 0 && biasVolumeMode != requestedMode) {
                opserr << "WARNING conflicting bias-volume options for "
                       << "RIVASAND02 tag " << tag
                       << "; select one -BiasVolume mode (legacy aliases: "
                       << "-noBiasVolume=1, -fieldBiasVolume=2)" << endln;
                return 0;
            }
            biasVolumeMode = requestedMode;
        } else if (std::strcmp(option, "-stage") == 0) {
            count = 1;
            if (OPS_GetIntInput(&count, &stage) < 0) {
                opserr << "WARNING invalid -stage for RIVASAND02 tag "
                       << tag << endln;
                return 0;
            }
            stageSpecified = true;
        } else if (std::strcmp(option, "-initialStress") == 0) {
            double stress[6];
            count = 6;
            if (OPS_GetDoubleInput(&count, stress) < 0) {
                opserr << "WARNING invalid -initialStress for RIVASAND02 tag "
                       << tag << endln;
                return 0;
            }
            for (int i = 0; i < 6; ++i) initialStress(i) = stress[i];
            initialStressSpecified = true;
        } else {
            opserr << "WARNING unknown RIVASAND02 option '" << option
                   << "' for tag " << tag << endln;
            return 0;
        }
    }

    if (initialStressSpecified && !stageSpecified) stage = 1;
    if (biasVolumeMode < 0) biasVolumeMode = 0;
    if (stage != 0 && !initialStressSpecified) {
        opserr << "WARNING RIVASAND02 -stage 1 or 2 requires a compressive "
               << "-initialStress; otherwise create at stage 0, establish "
               << "geostatic stress, and use updateMaterialStage" << endln;
        return 0;
    }

    RIVASAND02 *material = new RIVASAND02(
        tag, values[0], values[1], values[2], values[3], values[4],
        values[5], values[6], values[7], values[8], values[9], values[10], values[11],
        rho, fixedSubsteps, stressScale, pMin, tangentPressureFloor,
        residualPressure, geostaticAdmission, stage, initialStress);
    if (material == 0 || !material->isValid()) {
        opserr << "WARNING invalid RIVASAND02 material with tag "
               << tag << endln;
        delete material;
        return 0;
    }
    material->setReversalLatch(reversalLatch);
    material->setFieldBiasMeanCorrection(biasVolumeMode == 2);
    material->setBiasReversibleVolumeEnabled(biasVolumeMode != 1);
    material->setTangentType(tangentType);
    return material;
}

RIVASAND02::RIVASAND02(
    int tag, double Dr, double G0, double M, double kd, double h, double m,
    double zeta, double eMax, double eMin, double Q, double R, double nG,
    double rho, int fixedSubsteps, double stressScale, double pMin,
    double tangentPressureFloor, double residualPressure,
    bool geostaticAdmission, int initialStage, const Vector &initialStress)
    : NDMaterial(tag, ND_TAG_RIVASAND02),
      mDr(Dr), mG0(G0), mRho(rho), mStressScale(stressScale),
      mTangentPressureFloor(0.0),
      mFixedSubsteps(fixedSubsteps), mStage(initialStage),
      mInitialStage(initialStage), mValid(true),
      mGeostaticAdmission(geostaticAdmission),
      mReversalLatch(false), mLatchValid(false), mLatchedReversal(0),
      mInitialStress(6), mCommittedStrain(6), mTrialStrain(6),
      mCommittedStress(6), mTrialStress(6), mTangent(6, 6),
      mInitialTangent(6, 6), mStateOutput(RIVA_IB_STATE_VALUE_COUNT),
      mScalarOutput(1)
{
    mInitialStress.Zero();
    if (initialStress.Size() == 6) mInitialStress = initialStress;
    setReferenceParameters();
    if (pMin > 0.0) mParameters.base.p_min = pMin;
    mParameters.base.p_residual = residualPressure;
    mTangentPressureFloor = tangentPressureFloor > 0.0 ?
        tangentPressureFloor : mParameters.base.p_ref/200.0;
    mTangentPressureFloor = riva_max(
        mTangentPressureFloor, mParameters.base.p_min);
    setMaterialParameters(G0, M, kd, h, m, zeta, eMax, eMin, Q, R, nG);

    if (!(mStressScale > 0.0) || !(mRho >= 0.0) ||
        mFixedSubsteps < 1 || (mStage < 0 || mStage > 2) ||
        !std::isfinite(mDr) || mDr < 0.0 || mDr > 1.0 ||
        !std::isfinite(mParameters.base.p_min) || !(mParameters.base.p_min > 0.0) ||
        !std::isfinite(mParameters.base.p_residual) ||
        mParameters.base.p_residual < 0.0 ||
        !std::isfinite(mTangentPressureFloor) ||
        !(mTangentPressureFloor > 0.0) ||
        !finiteVector(mInitialStress) ||
        !riva_material_parameters_valid(&mParameters.base, &mMaterial)) {
        mValid = false;
    }

    const double shear = mMaterial.E_ref/(2.0*(1.0+mParameters.base.nu));
    const double bulk = mMaterial.E_ref/(3.0*(1.0-2.0*mParameters.base.nu));
    buildTangent(bulk, shear, mInitialTangent);
    revertToStart();
}

RIVASAND02::RIVASAND02()
    : NDMaterial(0, ND_TAG_RIVASAND02),
      mDr(0.0), mG0(RIVA_REFERENCE_G0), mRho(0.0), mStressScale(1.0),
      mTangentPressureFloor(0.0), mFixedSubsteps(1),
      mStage(0), mInitialStage(0), mValid(false),
      mGeostaticAdmission(false),
      mReversalLatch(false), mLatchValid(false), mLatchedReversal(0),
      mInitialStress(6), mCommittedStrain(6), mTrialStrain(6),
      mCommittedStress(6), mTrialStress(6), mTangent(6, 6),
      mInitialTangent(6, 6), mStateOutput(RIVA_IB_STATE_VALUE_COUNT),
      mScalarOutput(1)
{
    mInitialStress.Zero();
    mCommittedStrain.Zero();
    mTrialStrain.Zero();
    mCommittedStress.Zero();
    mTrialStress.Zero();
    mCommittedState = riva_ib_state_t{};
    mTrialState = riva_ib_state_t{};
    setReferenceParameters();
    mTangentPressureFloor = riva_max(
        mParameters.base.p_ref/200.0, mParameters.base.p_min);
    mMaterial = riva_reference_material_parameters(&mParameters.base);
    const double shear = mMaterial.E_ref/(2.0*(1.0+mParameters.base.nu));
    const double bulk = mMaterial.E_ref/(3.0*(1.0-2.0*mParameters.base.nu));
    buildTangent(bulk, shear, mInitialTangent);
    mTangent = mInitialTangent;
}

RIVASAND02::~RIVASAND02()
{
}

void
RIVASAND02::setReferenceParameters(void)
{
    mParameters = riva_ib_reference_parameters(mStressScale);
}

void
RIVASAND02::setMaterialParameters(double G0, double M, double kd, double h,
    double m, double zeta, double eMax, double eMin, double Q, double R,
    double nG)
{
    mMaterial = riva_reference_material_parameters(&mParameters.base);
    mG0 = G0;
    riva_material_set_G0(&mParameters.base, &mMaterial, G0);
    mMaterial.M = M;
    mMaterial.kd = kd;
    mMaterial.h = h;
    mMaterial.m = m;
    mMaterial.zeta = zeta;
    mMaterial.e_max = eMax;
    mMaterial.e_min = eMin;
    mMaterial.Q = Q;
    mMaterial.R = R;
    mMaterial.n_G = nG;
    mParameters.reference_relative_density_value =
        (eMax-mParameters.base.state_shakedown_reference_void_ratio)/(eMax-eMin);
}

double
RIVASAND02::initialVoidRatio(void) const
{
    return riva_void_ratio_from_material_relative_density(&mMaterial, mDr);
}

int
RIVASAND02::activateFromCommittedStress(void)
{
    mTrialPlasticLoading = mCommittedPlasticLoading = false;
    tensor_t stress = stressToTensor(mCommittedStress);
    const double physicalPressure = riva_pressure(stress);
    const double conePressure = riva_cone_pressure(&mParameters.base, stress);
    if (!(conePressure > mParameters.base.p_min)) {
        if (!mGeostaticAdmission ||
            !std::isfinite(physicalPressure) ||
            !std::isfinite(conePressure)) {
            opserr << "RIVASAND02 tag " << this->getTag()
                   << " cannot enter a nonlinear stage: translated cone pressure "
                   << "p'+pResidual must exceed pMin"
                   << endln;
            return -1;
        }
        /* Keep the physical skeleton stress unchanged whenever pResidual
         * already places it inside the translated cone.  Only a state at or
         * beyond that translated apex needs a mean-stress projection. */
        const double admittedConePressure = std::nextafter(
            mParameters.base.p_min, std::numeric_limits<double>::infinity());
        const double admittedPhysicalPressure = riva_physical_pressure(
            &mParameters.base, admittedConePressure);
        stress = riva_sub(riva_dev(stress),
                          riva_iso(admittedPhysicalPressure));
    }
    riva_ib_state_t state = {};
    if (!riva_ib_initialize_material(&mParameters, &mMaterial, stress,
                                 initialVoidRatio(), &state)) return -1;
    if (mGeostaticAdmission) {
        int32_t admitted = 0;
        if (!riva_ib_admit_geostatic_state(
                &mParameters, &mMaterial, &state, &admitted)) return -1;
    }
    if (!mGeostaticAdmission || mStage == 2) {
        tensor_t reference = stress;
        reference.xy = reference.yz = reference.xz = 0.0;
        if (!riva_ib_begin_dynamic_phase(
                &mParameters, &mMaterial, &reference, &state)) return -1;
    }
    mCommittedState = state;
    mTrialState = state;
    tensorToStress(state.base.stress, mCommittedStress);
    mTrialStress = mCommittedStress;
    return 0;
}

int
RIVASAND02::beginDynamicFromCommittedState(void)
{
    if (!mCommittedState.base.initialized) return -1;
    riva_ib_state_t state = mCommittedState;
    tensor_t reference = state.base.stress;
    reference.xy = reference.yz = reference.xz = 0.0;
    if (!riva_ib_begin_dynamic_phase(
            &mParameters, &mMaterial, &reference, &state)) return -1;
    mTrialPlasticLoading = mCommittedPlasticLoading = false;
    mCommittedState = state;
    mTrialState = state;
    tensorToStress(state.base.stress, mCommittedStress);
    mTrialStress = mCommittedStress;
    return 0;
}

tensor_t
RIVASAND02::strainIncrementToTensor(const Vector &increment)
{
    tensor_t result = {increment(0), increment(1), increment(2),
                           0.5*increment(3), 0.5*increment(4),
                           0.5*increment(5)};
    return result;
}

tensor_t
RIVASAND02::stressToTensor(const Vector &stress)
{
    tensor_t result = {stress(0), stress(1), stress(2),
                           stress(3), stress(4), stress(5)};
    return result;
}

void
RIVASAND02::tensorToStress(tensor_t tensor, Vector &stress)
{
    stress(0) = tensor.xx;
    stress(1) = tensor.yy;
    stress(2) = tensor.zz;
    stress(3) = tensor.xy;
    stress(4) = tensor.yz;
    stress(5) = tensor.xz;
}

void
RIVASAND02::buildTangent(double bulk, double shear, Matrix &matrix) const
{
    matrix.Zero();
    const double lambda = bulk - 2.0*shear/3.0;
    for (int i = 0; i < 3; ++i) {
        matrix(i, i) = lambda + 2.0*shear;
        for (int j = 0; j < 3; ++j)
            if (i != j) matrix(i, j) = lambda;
    }
    matrix(3, 3) = shear;
    matrix(4, 4) = shear;
    matrix(5, 5) = shear;
}

bool
RIVASAND02::buildContinuumTangent(double result[6][6]) const
{
    const riva_ib_state_t &s=mTrialState;
    const riva_parameters_t &b=mParameters.base;
    const double pressure=riva_cone_pressure(&b,s.base.stress);
    if (!s.base.initialized || pressure<=mTangentPressureFloor ||
        pressure<=riva_cone_pressure_floor(&b) || s.base.beta<=1.e-6 ||
        (s.base.D_re>0.0 && s.base.eps_v_reversible<=0.0)) return false;
    double mb,md,xi;
    riva_surfaces(&b,&mMaterial,pressure,s.base.void_ratio,&mb,&md,&xi);
    if (riva_norm(s.base.alpha)>=sqrt(2.0/3.0)*mb*(1.0-1.e-10)) return false;
    double shear,bulk,hardening;
    double beta=s.base.beta, dilatancy=s.base.D, kinematic=0.0;
    tensor_t normal=s.base.n,flow=s.base.n;
    const double gate=riva_ib_mapping_gate(&mParameters,&s);
    if (gate>1.e-14 && s.base.cyclic_phase_active) {
        riva_moduli(&b,&mMaterial,pressure,&shear,&bulk);
        const double capacity=riva_ib_mapping_capacity(&mParameters,&s,mb);
        shear*=riva_ib_mapping_shear_ratio(&mParameters,s.mapping_backstress,capacity);
        const tensor_t center=riva_add(riva_add(riva_scale(s.base.alpha0,1.0-gate),
            riva_scale(s.mapping_anchor,gate)),riva_scale(s.mapping_backstress,gate));
        tensor_t ray; int failed=0;
        riva_ib_mapping_intersection(&mParameters,s.base.alpha,center,mb,
            s.base.last_host_deviatoric_strain_direction,&beta,&normal,&ray,&failed);
        if (failed) return false;
        flow=riva_ib_mapping_flow(&mParameters,normal,ray,s.mapping_directional_fabric,
            &s.base.static_bias_tensor);
        const tensor_t fabric=riva_add(s.base.fabric,riva_scale(
            s.mapping_directional_fabric,gate*mParameters.mapping_fabric_dilatancy_weight));
        double dir,dre;
        riva_ib_dilatancy(&b,&mMaterial,s.base.alpha,flow,beta,fabric,pressure,
            s.base.void_ratio,s.base.eps_v_irreversible,s.base.eps_v_reversible,&dir,&dre);
        dir*=riva_ib_irreversible_factor(&mParameters,&s)*
            riva_ib_directional_phase_scale(&mParameters,&s,s.mapping_backstress,capacity);
        dilatancy=dir+dre;
        tensor_t rate;
        riva_ib_backstress_update(&mParameters,s.mapping_backstress,flow,capacity,0.0,&rate);
        kinematic=gate*pressure*riva_ddot(normal,rate);
        hardening=pressure*mMaterial.h*pow(riva_max(pressure,b.p_min)/b.p_ref,-b.q_H)*
            pow(riva_max(beta,1.e-12),mMaterial.m);
    } else {
        riva_ib_moduli_for_state(&mParameters,&mMaterial,&s,pressure,&shear,&bulk);
        hardening=pressure*riva_ib_hardening_for_state(&mParameters,&mMaterial,&s,pressure)*
            pow(riva_max(beta,1.e-12),mMaterial.m);
    }
    return riva_continuum::build(shear,bulk,normal,flow,riva_ddot(s.base.alpha,normal),
        dilatancy,(2.0/3.0)*hardening+kinematic,b.denominator_floor_ratio,result);
}

void
RIVASAND02::updateTrialTangent(void)
{
    mTangentStatus = 0;
    if (mStage != 0 && mTrialState.base.initialized) {
        double shear = 0.0;
        double bulk = 0.0;
        // Regularize only the tangent advertised to OpenSees here.  The
        // explicit research stress update continues to use the constitutive
        // pressure floor, so
        // the default p_ref/200 tangent floor does not silently alter the
        // calibrated constitutive pressure path.
        const double pressure = riva_max(
            riva_cone_pressure(&mParameters.base, mTrialState.base.stress),
            mTangentPressureFloor);
        riva_ib_moduli_for_state(&mParameters, &mMaterial, &mTrialState,
                                 pressure, &shear, &bulk);
        // A rejected global Newton trial can drive auxiliary state variables
        // far outside the committed neighborhood even when the kernel later
        // reverts. Keep the fallback elastic operator finite and positive.
        // The opt-in continuum operator below is generally nonsymmetric;
        // unlike the elastic fallback it is not guaranteed positive definite.
        if (std::isfinite(shear) && std::isfinite(bulk) &&
            shear > 0.0 && bulk > 0.0)
            buildTangent(bulk, shear, mTangent);
        else
            mTangent = mInitialTangent;
    } else {
        mTangent = mInitialTangent;
    }
    if (mTangentType == 1 && mStage != 0) {
        mTangentStatus = 2; // elastic: inactive or discrete integration event
        if (mTrialPlasticLoading) {
            double continuum[6][6];
            mTangentStatus = 3; // elastic safeguard at a nonsmooth/invalid state
            if (buildContinuumTangent(continuum)) {
                for (int i=0;i<6;++i)
                    for (int j=0;j<6;++j) mTangent(i,j)=continuum[i][j];
                mTangentStatus = 1;
            }
        }
    }
}

int
RIVASAND02::setTangentType(int type)
{
    if (type != 0 && type != 1) return -1;
    if (type == mTangentType) return 0;
    mTangentType = type;
    if (type == 0) mTrialPlasticLoading = mCommittedPlasticLoading = false;
    updateTrialTangent();
    return 0;
}

int
RIVASAND02::setTrialStrain(const Vector &strain)
{
    if (!mValid || strain.Size() != 6 || !finiteVector(strain)) return -1;
    mTrialStrain = strain;
    Vector increment = mTrialStrain - mCommittedStrain;
    mTrialPlasticLoading = mCommittedPlasticLoading;

    if (mStage == 0) {
        mTrialPlasticLoading = false;
        mTangentStatus = 0;
        mTrialStress = mCommittedStress;
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j)
                mTrialStress(i) += mInitialTangent(i, j)*increment(j);
        mTrialState = mCommittedState;
        mTangent = mInitialTangent;
        return 0;
    }

    if (!mCommittedState.base.initialized && activateFromCommittedStress() != 0)
        return -1;
    mTrialState = mCommittedState;
    double incrementNorm2 = 0.0;
    for (int i = 0; i < 6; ++i)
        incrementNorm2 += increment(i)*increment(i);
    if (incrementNorm2 == 0.0) {
        // OpenSees may repeatedly query an unchanged trial strain while
        // assembling an equilibrium iteration.  Do not advance a
        // rate-independent material state for such host-only evaluations.
        mTrialStress = mCommittedStress;
        updateTrialTangent();
        return 0;
    }
    tensor_t stress = mTrialState.base.stress;
    riva_update_info_t information = {};
    const int reversalOverride =
        (mReversalLatch && mLatchValid) ? mLatchedReversal : -1;
    if (!riva_ib_update_material_ex(
            &mParameters, &mMaterial, strainIncrementToTensor(increment),
            mFixedSubsteps, &mTrialState, &stress, &information,
            reversalOverride)) {
        mTrialState = mCommittedState;
        mTrialStress = mCommittedStress;
        mTrialStrain = mCommittedStrain;
        if (mTangentType == 1) updateTrialTangent();
        return -1;
    }
    if (mReversalLatch && !mLatchValid && information.accepted_substeps > 0) {
        mLatchedReversal = information.reversal_registered ? 1 : 0;
        mLatchValid = true;
    }
    tensorToStress(stress, mTrialStress);
    if (!finiteVector(mTrialStress)) {
        mTrialState = mCommittedState;
        mTrialStress = mCommittedStress;
        mTrialStrain = mCommittedStrain;
        updateTrialTangent();
        return -1;
    }
    if (mTangentType == 1)
        mTrialPlasticLoading=riva_continuum::smoothPlasticStep(
            mCommittedState.base,mTrialState.base) &&
            mTrialState.mapping_stress_corrections==mCommittedState.mapping_stress_corrections &&
            mTrialState.mapping_monotone_caps==mCommittedState.mapping_monotone_caps;
    updateTrialTangent();
    return 0;
}

int
RIVASAND02::setTrialStrain(const Vector &strain, const Vector &rate)
{
    return setTrialStrain(strain);
}

int
RIVASAND02::setTrialStrainIncr(const Vector &increment)
{
    if (increment.Size() != 6) return -1;
    Vector target = mTrialStrain;
    target += increment;
    return setTrialStrain(target);
}

int
RIVASAND02::setTrialStrainIncr(const Vector &increment,
                                  const Vector &rate)
{
    return setTrialStrainIncr(increment);
}

const Vector &
RIVASAND02::getStress(void)
{
    return mTrialStress;
}

const Vector &
RIVASAND02::getStrain(void)
{
    return mTrialStrain;
}

const Matrix &
RIVASAND02::getTangent(void)
{
    return mTangent;
}

const Matrix &
RIVASAND02::getInitialTangent(void)
{
    return mInitialTangent;
}

double
RIVASAND02::getRho(void)
{
    return mRho;
}

int
RIVASAND02::commitState(void)
{
    mCommittedPlasticLoading = mTrialPlasticLoading;
    mCommittedStrain = mTrialStrain;
    mCommittedStress = mTrialStress;
    mCommittedState = mTrialState;
    mLatchValid = false;
    mLatchedReversal = 0;
    return 0;
}

int
RIVASAND02::revertToLastCommit(void)
{
    mTrialPlasticLoading = mCommittedPlasticLoading;
    mTrialStrain = mCommittedStrain;
    mTrialStress = mCommittedStress;
    mTrialState = mCommittedState;
    mLatchValid = false;
    mLatchedReversal = 0;
    updateTrialTangent();
    return 0;
}

int
RIVASAND02::revertToStart(void)
{
    mTrialPlasticLoading = mCommittedPlasticLoading = false;
    mTangentStatus = 0;
    mStage = mInitialStage;
    mCommittedStrain.Zero();
    mTrialStrain.Zero();
    mCommittedStress = mInitialStress;
    mTrialStress = mInitialStress;
    mCommittedState = riva_ib_state_t{};
    mTrialState = riva_ib_state_t{};
    mLatchValid = false;
    mLatchedReversal = 0;
    mTangent = mInitialTangent;
    if (mValid && mStage != 0 && activateFromCommittedStress() != 0)
        mValid = false;
    return mValid ? 0 : -1;
}

NDMaterial *
RIVASAND02::getCopy(void)
{
    return new RIVASAND02(*this);
}

NDMaterial *
RIVASAND02::getCopy(const char *code)
{
    if (std::strcmp(code, "ThreeDimensional") == 0 ||
        std::strcmp(code, "3D") == 0 ||
        std::strcmp(code, "RIVASAND02") == 0)
        return getCopy();
    return 0;
}

const char *
RIVASAND02::getType(void) const
{
    return "ThreeDimensional";
}

int
RIVASAND02::getOrder(void) const
{
    return 6;
}

bool
RIVASAND02::isValid(void) const
{
    return mValid;
}

int
RIVASAND02::sendSelf(int commitTag, Channel &theChannel)
{
    Vector data(RIVASerializedSize);
    data.Zero();
    data(0) = this->getTag();
    data(1) = mDr;
    data(2) = mStressScale;
    data(3) = mRho;
    data(4) = mFixedSubsteps;
    data(5) = mStage;
    data(6) = mInitialStage;
    for (int i = 0; i < 6; ++i) data(7+i) = mInitialStress(i);
    data(13) = mMaterial.M;
    data(14) = mMaterial.kd;
    data(15) = mMaterial.h;
    data(16) = mMaterial.m;
    data(17) = mMaterial.zeta;
    data(18) = mMaterial.e_max;
    data(19) = mMaterial.e_min;
    data(20) = mMaterial.Q;
    data(21) = mMaterial.R;
    data(22) = mMaterial.n_G;
    for (int i = 0; i < 6; ++i) data(23+i) = mCommittedStrain(i);
    for (int i = 0; i < 6; ++i) data(29+i) = mCommittedStress(i);
    double stateValues[RIVA_IB_STATE_VALUE_COUNT] = {};
    riva_ib_state_values(&mCommittedState, stateValues);
    for (int i = 0; i < RIVA_IB_STATE_VALUE_COUNT; ++i)
        data(35+i) = stateValues[i];
    data(174) = mCommittedState.base.initialized;
    data(175) = mParameters.base.p_min;
    data(176) = mTangentPressureFloor;
    data(177) = mParameters.base.p_residual;
    const int configurationFlags =
        (mGeostaticAdmission ? RIVAGeostaticAdmissionFlag : 0) |
        (mReversalLatch ? RIVAReversalLatchFlag : 0) |
        (mParameters.field_bias_mean_correction_enabled ?
            RIVAFieldBiasMeanCorrectionFlag : 0) |
        (!mParameters.base.bias_reversible_volume_enabled ?
            RIVANoBiasVolumeFlag : 0) |
        (mTangentType == 1 ? RIVAContinuumTangentFlag : 0) |
        (mCommittedPlasticLoading ? RIVAPlasticLoadingFlag : 0);
    data(178) = configurationFlags;
    data(179) = mCommittedState.base.geostatic_admitted;
    data(180) = RIVA_IB_KERNEL_REVISION;

    data(181) = mG0;
    data(182) = RIVAAdapterRevision;
    if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
        opserr << "RIVASAND02::sendSelf failed for tag "
               << this->getTag() << endln;
        return -1;
    }
    return 0;
}

int
RIVASAND02::restoreState(
    const double values[RIVA_IB_STATE_VALUE_COUNT], int initialized,
    int geostaticAdmitted, riva_ib_state_t &state)
{
    int i = 0;
    state = riva_ib_state_t{};
#define RIVA_RESTORE_TENSOR(t) do { \
    (t).xx=values[i++]; (t).yy=values[i++]; (t).zz=values[i++]; \
    (t).xy=values[i++]; (t).yz=values[i++]; (t).xz=values[i++]; \
} while (0)
    RIVA_RESTORE_TENSOR(state.base.stress);
    RIVA_RESTORE_TENSOR(state.base.alpha);
    RIVA_RESTORE_TENSOR(state.base.alpha0);
    RIVA_RESTORE_TENSOR(state.base.alpha01);
    RIVA_RESTORE_TENSOR(state.base.n);
    RIVA_RESTORE_TENSOR(state.base.fabric);
    state.base.D=values[i++]; state.base.beta=values[i++];
    state.base.lambda_total=values[i++];
    state.base.ep_eq_since_reversal=values[i++];
    state.base.void_ratio=values[i++];
    state.base.reversals=(int64_t)std::llround(values[i++]);
    state.base.pressure_floor_hits=(int64_t)std::llround(values[i++]);
    state.base.denominator_floor_hits=(int64_t)std::llround(values[i++]);
    state.base.beta_fallbacks=(int64_t)std::llround(values[i++]);
    state.base.eps_v_total=values[i++];
    state.base.eps_v_confining=values[i++];
    state.base.eps_v_irreversible=values[i++];
    state.base.eps_v_reversible=values[i++];
    state.base.pressure_anchor=values[i++]; state.base.D_ir=values[i++];
    state.base.D_re=values[i++];
    RIVA_RESTORE_TENSOR(state.base.last_reversal_deviator);
    state.base.cyclic_amplitude=values[i++];
    state.base.amplitude_factor=values[i++];
    state.base.amplitude_reversals=(int64_t)std::llround(values[i++]);
    state.base.initial_relative_state=values[i++];
    state.base.state_contraction_factor=values[i++];
    state.base.effective_knee_ratio=values[i++];
    RIVA_RESTORE_TENSOR(state.base.geostatic_deviator);
    RIVA_RESTORE_TENSOR(state.base.static_bias_tensor);
    RIVA_RESTORE_TENSOR(state.base.cyclic_direction);
    state.base.static_bias_index=values[i++];
    state.base.cyclic_phase_active=(int32_t)std::llround(values[i++]);
    state.base.bias_ratchet_strain=values[i++];
    state.base.physical_eps_v_total=values[i++];
    state.base.bias_reversible_volume=values[i++];
    RIVA_RESTORE_TENSOR(state.base.last_host_deviatoric_strain_direction);

    state.phase_irreversible_volume=values[i++];
    state.phase_reversible_volume=values[i++];
    state.phase_potential_anchor=values[i++];
    state.phase_accumulation_lambda_anchor=values[i++];
    state.phase_accumulation_hardening_state=values[i++];
    RIVA_RESTORE_TENSOR(state.unbiased_phase_direction);
    state.loose_shear_lambda_anchor=values[i++];
    state.loose_shear_hardening_state=values[i++];
    state.loose_shear_gate_value=values[i++];
    RIVA_RESTORE_TENSOR(state.mapping_anchor);
    RIVA_RESTORE_TENSOR(state.mapping_backstress);
    RIVA_RESTORE_TENSOR(state.mapping_directional_fabric);
    state.mapping_gate_value=values[i++];
    state.mapping_capacity=values[i++];
    state.mapping_kinematic_denominator=values[i++];
    state.mapping_shear_modulus_ratio=values[i++];
    state.mapping_phase_contraction_scale=values[i++];
    state.mapping_outer_residual=values[i++];
    state.mapping_stress_corrections=(int64_t)std::llround(values[i++]);
    state.mapping_corrector_passes=(int64_t)std::llround(values[i++]);
    state.mapping_monotone_caps=(int64_t)std::llround(values[i++]);
    state.initial_relative_density_value=values[i++];
    state.intermediate_low_gate_value=values[i++];
    state.intermediate_high_gate_base=values[i++];
    state.ep_half_last=values[i++];
    state.field_bias_mean_activity=values[i++];
#undef RIVA_RESTORE_TENSOR
    state.base.initialized = initialized;
    state.base.geostatic_admitted = geostaticAdmitted;
    return i == RIVA_IB_STATE_VALUE_COUNT ? 0 : -1;
}

int
RIVASAND02::recvSelf(int commitTag, Channel &theChannel,
                        FEM_ObjectBroker &theBroker)
{
    Vector data(RIVASerializedSize);
    if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
        opserr << "RIVASAND02::recvSelf failed" << endln;
        return -1;
    }
    if (data.Size() != RIVASerializedSize ||
        (data(182) != 1 && data(182) != RIVAAdapterRevision)) {
        opserr << "RIVASAND02::recvSelf incompatible adapter checkpoint" << endln;
        return -1;
    }
    if (!std::isfinite(data(178)) || data(178)!=std::floor(data(178)) ||
        data(178)<0 || data(178)>(data(182)==1?15:63)) return -1;
    if (data(180) != RIVA_IB_KERNEL_REVISION) {
        opserr << "RIVASAND02::recvSelf incompatible "
               << "kernel revision " << data(180) << endln;
        return -1;
    }
    this->setTag((int)data(0));
    mDr = data(1);
    mStressScale = data(2);
    mRho = data(3);
    mFixedSubsteps = (int)std::llround(data(4));
    mStage = (int)std::llround(data(5));
    mInitialStage = (int)std::llround(data(6));
    for (int i = 0; i < 6; ++i) mInitialStress(i) = data(7+i);
    setReferenceParameters();
    mParameters.base.p_min = data(175);
    mParameters.base.p_residual = data(177);
    const int configurationFlags = (int)std::llround(data(178));
    mTangentType=(configurationFlags & RIVAContinuumTangentFlag)?1:0;
    mTrialPlasticLoading=mCommittedPlasticLoading=
        (configurationFlags & RIVAPlasticLoadingFlag)!=0;
    if (mCommittedPlasticLoading && mTangentType==0) return -1;
    mGeostaticAdmission =
        (configurationFlags & RIVAGeostaticAdmissionFlag) != 0;
    mReversalLatch = (configurationFlags & RIVAReversalLatchFlag) != 0;
    mParameters.field_bias_mean_correction_enabled =
        (configurationFlags & RIVAFieldBiasMeanCorrectionFlag) != 0;
    mParameters.base.bias_reversible_volume_enabled =
        (configurationFlags & RIVANoBiasVolumeFlag) == 0;
    mLatchValid = false;
    mLatchedReversal = 0;
    mTangentPressureFloor = riva_max(data(176), mParameters.base.p_min);
    setMaterialParameters(data(181), data(13), data(14), data(15), data(16), data(17),
                          data(18), data(19), data(20), data(21), data(22));
    for (int i = 0; i < 6; ++i) mCommittedStrain(i) = data(23+i);
    for (int i = 0; i < 6; ++i) mCommittedStress(i) = data(29+i);
    double stateValues[RIVA_IB_STATE_VALUE_COUNT];
    for (int i = 0; i < RIVA_IB_STATE_VALUE_COUNT; ++i)
        stateValues[i] = data(35+i);
    if (restoreState(stateValues, (int)std::llround(data(174)),
                     (int)std::llround(data(179)),
                     mCommittedState) != 0) return -1;

    mValid = mStressScale > 0.0 && mRho >= 0.0 && mFixedSubsteps >= 1 &&
        (mStage >= 0 && mStage <= 2) && mDr >= 0.0 && mDr <= 1.0 &&
        std::isfinite(mParameters.base.p_min) && mParameters.base.p_min > 0.0 &&
        std::isfinite(mParameters.base.p_residual) &&
        mParameters.base.p_residual >= 0.0 &&
        std::isfinite(mTangentPressureFloor) &&
        mTangentPressureFloor > 0.0 &&
        finiteVector(mInitialStress) &&
        riva_material_parameters_valid(&mParameters.base, &mMaterial) &&
        (mStage == 0 || mCommittedState.base.initialized);
    const double shear = mMaterial.E_ref/
        (2.0*(1.0+mParameters.base.nu));
    const double bulk = mMaterial.E_ref/
        (3.0*(1.0-2.0*mParameters.base.nu));
    buildTangent(bulk, shear, mInitialTangent);
    mTrialStrain = mCommittedStrain;
    mTrialStress = mCommittedStress;
    mTrialState = mCommittedState;
    updateTrialTangent();
    return mValid ? 0 : -1;
}

const Vector &
RIVASAND02::getStateVector(void)
{
    double values[RIVA_IB_STATE_VALUE_COUNT] = {};
    riva_ib_state_values(&mTrialState, values);
    for (int i = 0; i < RIVA_IB_STATE_VALUE_COUNT; ++i)
        mStateOutput(i) = values[i];
    return mStateOutput;
}

const Vector &
RIVASAND02::getScalarResponse(int responseID)
{
    if (responseID == 18) { mScalarOutput(0)=mTangentType; return mScalarOutput; }
    if (responseID == 19) { mScalarOutput(0)=mTangentStatus; return mScalarOutput; }
    if (responseID == 4) {
        mScalarOutput(0) = mTrialState.base.initialized ?
            mTrialState.base.void_ratio : initialVoidRatio();
    } else if (responseID == 5) {
        const double anchor = mTrialState.base.pressure_anchor-
            mParameters.base.p_residual;
        mScalarOutput(0) = mTrialState.base.initialized && anchor > 0.0 ?
            1.0-riva_pressure(mTrialState.base.stress)/anchor : 0.0;
    } else if (responseID == 6) {
        mScalarOutput(0) = (double)mTrialState.base.reversals;
    } else if (responseID == 7) {
        mScalarOutput(0) = mTrialState.base.initialized ?
            riva_compatibility_residual(&mTrialState.base) : 0.0;
    } else if (responseID == 8) {
        mScalarOutput(0) = mParameters.base.p_min;
    } else if (responseID == 9) {
        mScalarOutput(0) = mTangentPressureFloor;
    } else if (responseID == 10) {
        mScalarOutput(0) = (double)mStage;
    } else if (responseID == 11) {
        mScalarOutput(0) = mReversalLatch ? 1.0 : 0.0;
    } else if (responseID == 12) {
        mScalarOutput(0) =
            mParameters.field_bias_mean_correction_enabled ? 1.0 : 0.0;
    } else if (responseID == 13) {
        mScalarOutput(0) =
            mParameters.base.bias_reversible_volume_enabled ? 0.0 : 1.0;
    } else if (responseID == 14) {
        mScalarOutput(0) = getBiasVolumeMode();
    } else if (responseID == 15) {
        mScalarOutput(0) = mG0;
    } else if (responseID == 16) {
        mScalarOutput(0) = mMaterial.E_ref/(2.0*(1.0+mParameters.base.nu));
    } else if (responseID == 17) {
        mScalarOutput(0) = mParameters.reference_relative_density_value;
    }
    return mScalarOutput;
}

Response *
RIVASAND02::setResponse(const char **argv, int argc, OPS_Stream &output)
{
    if (argc < 1) return 0;
    if (std::strcmp(argv[0], "TanType") == 0)
        return new MaterialResponse(this, 18, getScalarResponse(18));
    if (std::strcmp(argv[0], "tangentStatus") == 0)
        return new MaterialResponse(this, 19, getScalarResponse(19));
    if (std::strcmp(argv[0], "tangent") == 0)
        return new MaterialResponse(this, 20, getTangent());
    if (std::strcmp(argv[0], "G0") == 0)
        return new MaterialResponse(this, 15, getScalarResponse(15));
    if (std::strcmp(argv[0], "Gref") == 0)
        return new MaterialResponse(this, 16, getScalarResponse(16));
    if (std::strcmp(argv[0], "referenceRelativeDensity") == 0)
        return new MaterialResponse(this, 17, getScalarResponse(17));
    if (std::strcmp(argv[0], "stress") == 0 ||
        std::strcmp(argv[0], "stresses") == 0)
        return new MaterialResponse(this, 1, getStress());
    if (std::strcmp(argv[0], "strain") == 0 ||
        std::strcmp(argv[0], "strains") == 0)
        return new MaterialResponse(this, 2, getStrain());
    if (std::strcmp(argv[0], "state") == 0)
        return new MaterialResponse(this, 3, getStateVector());
    if (std::strcmp(argv[0], "voidRatio") == 0)
        return new MaterialResponse(this, 4, getScalarResponse(4));
    if (std::strcmp(argv[0], "effectivePressureRatio") == 0)
        return new MaterialResponse(this, 5, getScalarResponse(5));
    if (std::strcmp(argv[0], "reversals") == 0)
        return new MaterialResponse(this, 6, getScalarResponse(6));
    if (std::strcmp(argv[0], "compatibilityResidual") == 0)
        return new MaterialResponse(this, 7, getScalarResponse(7));
    if (std::strcmp(argv[0], "pressureFloor") == 0)
        return new MaterialResponse(this, 8, getScalarResponse(8));
    if (std::strcmp(argv[0], "tangentPressureFloor") == 0)
        return new MaterialResponse(this, 9, getScalarResponse(9));
    if (std::strcmp(argv[0], "stage") == 0)
        return new MaterialResponse(this, 10, getScalarResponse(10));
    if (std::strcmp(argv[0], "reversalLatch") == 0)
        return new MaterialResponse(this, 11, getScalarResponse(11));
    if (std::strcmp(argv[0], "fieldBiasVolume") == 0)
        return new MaterialResponse(this, 12, getScalarResponse(12));
    if (std::strcmp(argv[0], "noBiasVolume") == 0)
        return new MaterialResponse(this, 13, getScalarResponse(13));
    if (std::strcmp(argv[0], "BiasVolume") == 0 ||
        std::strcmp(argv[0], "biasVolume") == 0)
        return new MaterialResponse(this, 14, getScalarResponse(14));
    return NDMaterial::setResponse(argv, argc, output);
}

int
RIVASAND02::getResponse(int responseID, Information &materialInfo)
{
    if (responseID == 18 || responseID == 19)
        return materialInfo.setVector(getScalarResponse(responseID));
    if (responseID == 20) return materialInfo.setMatrix(getTangent());
    if (responseID == 1) return materialInfo.setVector(getStress());
    if (responseID == 2) return materialInfo.setVector(getStrain());
    if (responseID == 3) return materialInfo.setVector(getStateVector());
    if (responseID >= 4 && responseID <= 17)
        return materialInfo.setVector(getScalarResponse(responseID));
    return NDMaterial::getResponse(responseID, materialInfo);
}

int
RIVASAND02::setParameter(const char **argv, int argc,
                            Parameter &parameter)
{
    if (argc < 2 || std::atoi(argv[1]) != this->getTag()) return -1;
    if (std::strcmp(argv[0], "updateMaterialStage") == 0 ||
        std::strcmp(argv[0], "materialState") == 0)
        return parameter.addObject(StageParameter, this);
    if (std::strcmp(argv[0], "noSubsteps") == 0 ||
        std::strcmp(argv[0], "nSub") == 0)
        return parameter.addObject(SubstepParameter, this);
    return -1;
}

int
RIVASAND02::updateParameter(int responseID, Information &information)
{
    if (responseID == StageParameter) {
        const int requested = information.theInt;
        if (requested < 0 || requested > 2) return -1;
        if (requested == mStage) return 0;
        if (requested == 1) {
            if (mStage != 0) return -1;
            if (activateFromCommittedStress() != 0) return -1;
            mStage = 1;
        } else if (requested == 2) {
            if (mStage == 0) {
                mStage = 2;
                if (activateFromCommittedStress() != 0) {
                    mStage = 0;
                    return -1;
                }
            } else if (mStage == 1) {
                if (beginDynamicFromCommittedState() != 0) return -1;
                mStage = 2;
            } else return -1;
        } else {
            mStage = 0;
            mCommittedState = riva_ib_state_t{};
            mTrialState = riva_ib_state_t{};
        }
        updateTrialTangent();
        return 0;
    }
    if (responseID == SubstepParameter) {
        const int requested = information.theInt;
        if (requested < 1) return -1;
        mFixedSubsteps = requested;
        return 0;
    }
    return -1;
}

void
RIVASAND02::Print(OPS_Stream &output, int flag)
{
    output << "RIVASAND02, tag: " << this->getTag() << endln;
    output << "  Dr=" << mDr << " G0=" << mG0 << " M=" << mMaterial.M
           << " kd=" << mMaterial.kd << " h=" << mMaterial.h
           << " m=" << mMaterial.m << " zeta=" << mMaterial.zeta
           << " eMax=" << mMaterial.e_max << " eMin=" << mMaterial.e_min
           << " Q=" << mMaterial.Q << " R=" << mMaterial.R
           << " nG=" << mMaterial.n_G << endln;
    output << "  stage=" << mStage << " nSub=" << mFixedSubsteps
           << " stressScale=" << mStressScale
           << " pMin=" << mParameters.base.p_min
           << " pResidual=" << mParameters.base.p_residual
           << " geostaticAdmission="
           << (mGeostaticAdmission ? 1 : 0)
           << " reversalLatch=" << (mReversalLatch ? 1 : 0)
           << " BiasVolume=" << getBiasVolumeMode()
           << " fieldBiasVolume="
           << (mParameters.field_bias_mean_correction_enabled ? 1 : 0)
           << " noBiasVolume="
           << (mParameters.base.bias_reversible_volume_enabled ? 0 : 1)
           << " tangentPMin=" << mTangentPressureFloor
           << " TanType=" << mTangentType
           << " parameterSHA=" << RIVA_IB_PARAMETER_SHA256 << endln;
}
