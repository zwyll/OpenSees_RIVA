/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
** ****************************************************************** */

#include "RIVASandIntermediateBiasResearch.h"

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

const int RIVASerializedSize = 180;

bool finiteVector(const Vector &value)
{
    for (int i = 0; i < value.Size(); ++i)
        if (!std::isfinite(value(i))) return false;
    return true;
}

} // namespace

void *
OPS_RIVASandIntermediateBiasResearchMaterial(void)
{
    const int requiredValues = 11;
    if (OPS_GetNumRemainingInputArgs() < requiredValues + 1) {
        opserr << "Want: nDMaterial RIVASandIntermediateBiasResearch tag Dr M kd h m zeta "
               << "eMax eMin Q R nG <-rho value> <-nSub value> "
               << "<-stressScale value> <-pMin value> "
               << "<-tangentPMin value> <-pResidual value> "
               << "<-geostaticAdmission> <-stage 0|1|2> "
               << "<-initialStress sxx syy szz sxy syz sxz>" << endln;
        return 0;
    }

    int tag = 0;
    int count = 1;
    if (OPS_GetIntInput(&count, &tag) < 0) {
        opserr << "WARNING invalid RIVASandIntermediateBiasResearch tag" << endln;
        return 0;
    }

    double values[requiredValues];
    count = requiredValues;
    if (OPS_GetDoubleInput(&count, values) < 0) {
        opserr << "WARNING invalid RIVASandIntermediateBiasResearch material values for tag "
               << tag << endln;
        return 0;
    }

    double rho = 0.0;
    double stressScale = 1.0;
    double pMin = -1.0;
    double tangentPressureFloor = -1.0;
    double residualPressure = 0.0;
    bool geostaticAdmission = false;
    bool reversalLatch = false;
    bool noBiasVolume = false;
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
                opserr << "WARNING invalid -rho for RIVASandIntermediateBiasResearch tag "
                       << tag << endln;
                return 0;
            }
        } else if (std::strcmp(option, "-nSub") == 0 ||
                   std::strcmp(option, "-noSubsteps") == 0) {
            count = 1;
            if (OPS_GetIntInput(&count, &fixedSubsteps) < 0) {
                opserr << "WARNING invalid -nSub for RIVASandIntermediateBiasResearch tag "
                       << tag << endln;
                return 0;
            }
        } else if (std::strcmp(option, "-stressScale") == 0) {
            count = 1;
            if (OPS_GetDoubleInput(&count, &stressScale) < 0) {
                opserr << "WARNING invalid -stressScale for RIVASandIntermediateBiasResearch tag "
                       << tag << endln;
                return 0;
            }
        } else if (std::strcmp(option, "-pMin") == 0) {
            count = 1;
            if (OPS_GetDoubleInput(&count, &pMin) < 0 ||
                !std::isfinite(pMin) || !(pMin > 0.0)) {
                opserr << "WARNING invalid -pMin for RIVASandIntermediateBiasResearch tag "
                       << tag << "; value must be positive" << endln;
                return 0;
            }
        } else if (std::strcmp(option, "-tangentPMin") == 0) {
            count = 1;
            if (OPS_GetDoubleInput(&count, &tangentPressureFloor) < 0 ||
                !std::isfinite(tangentPressureFloor) ||
                !(tangentPressureFloor > 0.0)) {
                opserr << "WARNING invalid -tangentPMin for "
                       << "RIVASandIntermediateBiasResearch tag " << tag
                       << "; value must be positive" << endln;
                return 0;
            }
        } else if (std::strcmp(option, "-pResidual") == 0) {
            count = 1;
            if (OPS_GetDoubleInput(&count, &residualPressure) < 0 ||
                !std::isfinite(residualPressure) || residualPressure < 0.0) {
                opserr << "WARNING invalid -pResidual for RIVASandIntermediateBiasResearch tag "
                       << tag << "; value must be nonnegative" << endln;
                return 0;
            }
        } else if (std::strcmp(option, "-geostaticAdmission") == 0) {
            geostaticAdmission = true;
        } else if (std::strcmp(option, "-reversalLatch") == 0) {
            reversalLatch = true;
        } else if (std::strcmp(option, "-noBiasVolume") == 0) {
            noBiasVolume = true;
        } else if (std::strcmp(option, "-stage") == 0) {
            count = 1;
            if (OPS_GetIntInput(&count, &stage) < 0) {
                opserr << "WARNING invalid -stage for RIVASandIntermediateBiasResearch tag "
                       << tag << endln;
                return 0;
            }
            stageSpecified = true;
        } else if (std::strcmp(option, "-initialStress") == 0) {
            double stress[6];
            count = 6;
            if (OPS_GetDoubleInput(&count, stress) < 0) {
                opserr << "WARNING invalid -initialStress for RIVASandIntermediateBiasResearch tag "
                       << tag << endln;
                return 0;
            }
            for (int i = 0; i < 6; ++i) initialStress(i) = stress[i];
            initialStressSpecified = true;
        } else {
            opserr << "WARNING unknown RIVASandIntermediateBiasResearch option '" << option
                   << "' for tag " << tag << endln;
            return 0;
        }
    }

    if (initialStressSpecified && !stageSpecified) stage = 1;
    if (stage != 0 && !initialStressSpecified) {
        opserr << "WARNING RIVASandIntermediateBiasResearch -stage 1 or 2 requires a compressive "
               << "-initialStress; otherwise create at stage 0, establish "
               << "geostatic stress, and use updateMaterialStage" << endln;
        return 0;
    }

    RIVASandIntermediateBiasResearch *material = new RIVASandIntermediateBiasResearch(
        tag, values[0], values[1], values[2], values[3], values[4],
        values[5], values[6], values[7], values[8], values[9], values[10],
        rho, fixedSubsteps, stressScale, pMin, tangentPressureFloor,
        residualPressure, geostaticAdmission, stage, initialStress);
    if (material == 0 || !material->isValid()) {
        opserr << "WARNING invalid RIVASandIntermediateBiasResearch material with tag "
               << tag << endln;
        delete material;
        return 0;
    }
    material->setReversalLatch(reversalLatch);
    if (noBiasVolume) material->disableBiasReversibleVolume();
    return material;
}

RIVASandIntermediateBiasResearch::RIVASandIntermediateBiasResearch(
    int tag, double Dr, double M, double kd, double h, double m,
    double zeta, double eMax, double eMin, double Q, double R, double nG,
    double rho, int fixedSubsteps, double stressScale, double pMin,
    double tangentPressureFloor, double residualPressure,
    bool geostaticAdmission, int initialStage, const Vector &initialStress)
    : NDMaterial(tag, ND_TAG_RIVASandIntermediateBiasResearch),
      mDr(Dr), mRho(rho), mStressScale(stressScale),
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
    setMaterialParameters(M, kd, h, m, zeta, eMax, eMin, Q, R, nG);

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

    const double shear = mParameters.base.E_ref/(2.0*(1.0+mParameters.base.nu));
    const double bulk = mParameters.base.E_ref/(3.0*(1.0-2.0*mParameters.base.nu));
    buildTangent(bulk, shear, mInitialTangent);
    revertToStart();
}

RIVASandIntermediateBiasResearch::RIVASandIntermediateBiasResearch()
    : NDMaterial(0, ND_TAG_RIVASandIntermediateBiasResearch),
      mDr(0.0), mRho(0.0), mStressScale(1.0),
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
    mLatchValid = false;
    setReferenceParameters();
    mTangentPressureFloor = riva_max(
        mParameters.base.p_ref/200.0, mParameters.base.p_min);
    mMaterial = riva_reference_material_parameters(&mParameters.base);
    const double shear = mParameters.base.E_ref/(2.0*(1.0+mParameters.base.nu));
    const double bulk = mParameters.base.E_ref/(3.0*(1.0-2.0*mParameters.base.nu));
    buildTangent(bulk, shear, mInitialTangent);
    mTangent = mInitialTangent;
}

RIVASandIntermediateBiasResearch::~RIVASandIntermediateBiasResearch()
{
}

void
RIVASandIntermediateBiasResearch::setReferenceParameters(void)
{
    mParameters = riva_ib_reference_parameters(mStressScale);
}

void
RIVASandIntermediateBiasResearch::setMaterialParameters(double M, double kd, double h,
    double m, double zeta, double eMax, double eMin, double Q, double R,
    double nG)
{
    mMaterial = riva_reference_material_parameters(&mParameters.base);
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
}

double
RIVASandIntermediateBiasResearch::initialVoidRatio(void) const
{
    return riva_void_ratio_from_material_relative_density(&mMaterial, mDr);
}

int
RIVASandIntermediateBiasResearch::activateFromCommittedStress(void)
{
    tensor_t stress = stressToTensor(mCommittedStress);
    const double physicalPressure = riva_pressure(stress);
    const double conePressure = riva_cone_pressure(&mParameters.base, stress);
    if (!(conePressure > mParameters.base.p_min)) {
        if (!mGeostaticAdmission ||
            !std::isfinite(physicalPressure) ||
            !std::isfinite(conePressure)) {
            opserr << "RIVASandIntermediateBiasResearch tag " << this->getTag()
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
RIVASandIntermediateBiasResearch::beginDynamicFromCommittedState(void)
{
    if (!mCommittedState.base.initialized) return -1;
    riva_ib_state_t state = mCommittedState;
    tensor_t reference = state.base.stress;
    reference.xy = reference.yz = reference.xz = 0.0;
    if (!riva_ib_begin_dynamic_phase(
            &mParameters, &mMaterial, &reference, &state)) return -1;
    mCommittedState = state;
    mTrialState = state;
    tensorToStress(state.base.stress, mCommittedStress);
    mTrialStress = mCommittedStress;
    return 0;
}

tensor_t
RIVASandIntermediateBiasResearch::strainIncrementToTensor(const Vector &increment)
{
    tensor_t result = {increment(0), increment(1), increment(2),
                           0.5*increment(3), 0.5*increment(4),
                           0.5*increment(5)};
    return result;
}

tensor_t
RIVASandIntermediateBiasResearch::stressToTensor(const Vector &stress)
{
    tensor_t result = {stress(0), stress(1), stress(2),
                           stress(3), stress(4), stress(5)};
    return result;
}

void
RIVASandIntermediateBiasResearch::tensorToStress(tensor_t tensor, Vector &stress)
{
    stress(0) = tensor.xx;
    stress(1) = tensor.yy;
    stress(2) = tensor.zz;
    stress(3) = tensor.xy;
    stress(4) = tensor.yz;
    stress(5) = tensor.xz;
}

void
RIVASandIntermediateBiasResearch::buildTangent(double bulk, double shear, Matrix &matrix) const
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

void
RIVASandIntermediateBiasResearch::updateTrialTangent(void)
{
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
        // reverts. Never expose a non-finite/indefinite material tangent to
        // the assembled u-p system; use the finite reference tangent for that
        // rejected trial and let the global algorithm reduce the step.
        if (std::isfinite(shear) && std::isfinite(bulk) &&
            shear > 0.0 && bulk > 0.0)
            buildTangent(bulk, shear, mTangent);
        else
            mTangent = mInitialTangent;
    } else {
        mTangent = mInitialTangent;
    }
}

int
RIVASandIntermediateBiasResearch::setTrialStrain(const Vector &strain)
{
    if (!mValid || strain.Size() != 6 || !finiteVector(strain)) return -1;
    mTrialStrain = strain;
    Vector increment = mTrialStrain - mCommittedStrain;

    if (mStage == 0) {
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
    updateTrialTangent();
    return 0;
}

int
RIVASandIntermediateBiasResearch::setTrialStrain(const Vector &strain, const Vector &rate)
{
    return setTrialStrain(strain);
}

int
RIVASandIntermediateBiasResearch::setTrialStrainIncr(const Vector &increment)
{
    if (increment.Size() != 6) return -1;
    Vector target = mTrialStrain;
    target += increment;
    return setTrialStrain(target);
}

int
RIVASandIntermediateBiasResearch::setTrialStrainIncr(const Vector &increment,
                                  const Vector &rate)
{
    return setTrialStrainIncr(increment);
}

const Vector &
RIVASandIntermediateBiasResearch::getStress(void)
{
    return mTrialStress;
}

const Vector &
RIVASandIntermediateBiasResearch::getStrain(void)
{
    return mTrialStrain;
}

const Matrix &
RIVASandIntermediateBiasResearch::getTangent(void)
{
    return mTangent;
}

const Matrix &
RIVASandIntermediateBiasResearch::getInitialTangent(void)
{
    return mInitialTangent;
}

double
RIVASandIntermediateBiasResearch::getRho(void)
{
    return mRho;
}

int
RIVASandIntermediateBiasResearch::commitState(void)
{
    mCommittedStrain = mTrialStrain;
    mCommittedStress = mTrialStress;
    mCommittedState = mTrialState;
    mLatchValid = false;
    return 0;
}

int
RIVASandIntermediateBiasResearch::revertToLastCommit(void)
{
    mTrialStrain = mCommittedStrain;
    mTrialStress = mCommittedStress;
    mTrialState = mCommittedState;
    mLatchValid = false;
    updateTrialTangent();
    return 0;
}

int
RIVASandIntermediateBiasResearch::revertToStart(void)
{
    mStage = mInitialStage;
    mCommittedStrain.Zero();
    mTrialStrain.Zero();
    mCommittedStress = mInitialStress;
    mTrialStress = mInitialStress;
    mCommittedState = riva_ib_state_t{};
    mTrialState = riva_ib_state_t{};
    mLatchValid = false;
    mTangent = mInitialTangent;
    if (mValid && mStage != 0 && activateFromCommittedStress() != 0)
        mValid = false;
    return mValid ? 0 : -1;
}

NDMaterial *
RIVASandIntermediateBiasResearch::getCopy(void)
{
    return new RIVASandIntermediateBiasResearch(*this);
}

NDMaterial *
RIVASandIntermediateBiasResearch::getCopy(const char *code)
{
    if (std::strcmp(code, "ThreeDimensional") == 0 ||
        std::strcmp(code, "3D") == 0 ||
        std::strcmp(code, "RIVASandIntermediateBiasResearch") == 0)
        return getCopy();
    return 0;
}

const char *
RIVASandIntermediateBiasResearch::getType(void) const
{
    return "ThreeDimensional";
}

int
RIVASandIntermediateBiasResearch::getOrder(void) const
{
    return 6;
}

bool
RIVASandIntermediateBiasResearch::isValid(void) const
{
    return mValid;
}

int
RIVASandIntermediateBiasResearch::sendSelf(int commitTag, Channel &theChannel)
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
    data(173) = mCommittedState.base.initialized;
    data(174) = mParameters.base.p_min;
    data(175) = mTangentPressureFloor;
    data(176) = mParameters.base.p_residual;
    data(177) = mGeostaticAdmission ? 1.0 : 0.0;
    data(178) = mCommittedState.base.geostatic_admitted;
    data(179) = RIVA_IB_KERNEL_REVISION;

    if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
        opserr << "RIVASandIntermediateBiasResearch::sendSelf failed for tag "
               << this->getTag() << endln;
        return -1;
    }
    return 0;
}

int
RIVASandIntermediateBiasResearch::restoreState(
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
#undef RIVA_RESTORE_TENSOR
    state.base.initialized = initialized;
    state.base.geostatic_admitted = geostaticAdmitted;
    return i == RIVA_IB_STATE_VALUE_COUNT ? 0 : -1;
}

int
RIVASandIntermediateBiasResearch::recvSelf(int commitTag, Channel &theChannel,
                        FEM_ObjectBroker &theBroker)
{
    Vector data(RIVASerializedSize);
    if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
        opserr << "RIVASandIntermediateBiasResearch::recvSelf failed" << endln;
        return -1;
    }
    if ((uint32_t)std::llround(data(179)) != RIVA_IB_KERNEL_REVISION) {
        opserr << "RIVASandIntermediateBiasResearch::recvSelf incompatible "
               << "kernel revision " << data(179) << endln;
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
    mParameters.base.p_min = data(174);
    mParameters.base.p_residual = data(176);
    mGeostaticAdmission = std::llround(data(177)) != 0;
    mTangentPressureFloor = riva_max(data(175), mParameters.base.p_min);
    setMaterialParameters(data(13), data(14), data(15), data(16), data(17),
                          data(18), data(19), data(20), data(21), data(22));
    for (int i = 0; i < 6; ++i) mCommittedStrain(i) = data(23+i);
    for (int i = 0; i < 6; ++i) mCommittedStress(i) = data(29+i);
    double stateValues[RIVA_IB_STATE_VALUE_COUNT];
    for (int i = 0; i < RIVA_IB_STATE_VALUE_COUNT; ++i)
        stateValues[i] = data(35+i);
    if (restoreState(stateValues, (int)std::llround(data(173)),
                     (int)std::llround(data(178)),
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
    const double shear = mParameters.base.E_ref/
        (2.0*(1.0+mParameters.base.nu));
    const double bulk = mParameters.base.E_ref/
        (3.0*(1.0-2.0*mParameters.base.nu));
    buildTangent(bulk, shear, mInitialTangent);
    mTrialStrain = mCommittedStrain;
    mTrialStress = mCommittedStress;
    mTrialState = mCommittedState;
    updateTrialTangent();
    return mValid ? 0 : -1;
}

const Vector &
RIVASandIntermediateBiasResearch::getStateVector(void)
{
    double values[RIVA_IB_STATE_VALUE_COUNT] = {};
    riva_ib_state_values(&mTrialState, values);
    for (int i = 0; i < RIVA_IB_STATE_VALUE_COUNT; ++i)
        mStateOutput(i) = values[i];
    return mStateOutput;
}

const Vector &
RIVASandIntermediateBiasResearch::getScalarResponse(int responseID)
{
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
    }
    return mScalarOutput;
}

Response *
RIVASandIntermediateBiasResearch::setResponse(const char **argv, int argc, OPS_Stream &output)
{
    if (argc < 1) return 0;
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
    return NDMaterial::setResponse(argv, argc, output);
}

int
RIVASandIntermediateBiasResearch::getResponse(int responseID, Information &materialInfo)
{
    if (responseID == 1) return materialInfo.setVector(getStress());
    if (responseID == 2) return materialInfo.setVector(getStrain());
    if (responseID == 3) return materialInfo.setVector(getStateVector());
    if (responseID >= 4 && responseID <= 10)
        return materialInfo.setVector(getScalarResponse(responseID));
    return NDMaterial::getResponse(responseID, materialInfo);
}

int
RIVASandIntermediateBiasResearch::setParameter(const char **argv, int argc,
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
RIVASandIntermediateBiasResearch::updateParameter(int responseID, Information &information)
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
RIVASandIntermediateBiasResearch::Print(OPS_Stream &output, int flag)
{
    output << "RIVASandIntermediateBiasResearch, tag: " << this->getTag() << endln;
    output << "  Dr=" << mDr << " M=" << mMaterial.M
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
           << " tangentPMin=" << mTangentPressureFloor
           << " parameterSHA=" << RIVA_IB_PARAMETER_SHA256 << endln;
}
