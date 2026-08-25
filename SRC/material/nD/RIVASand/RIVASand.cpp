/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
** ****************************************************************** */

#include "RIVASand.h"

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

namespace {

const int RIVASerializedSize = 154;

bool finiteVector(const Vector &value)
{
    for (int i = 0; i < value.Size(); ++i)
        if (!std::isfinite(value(i))) return false;
    return true;
}

} // namespace

void *
OPS_RIVASandMaterial(void)
{
    const int requiredValues = 11;
    if (OPS_GetNumRemainingInputArgs() < requiredValues + 1) {
        opserr << "Want: nDMaterial RIVASand tag Dr M kd h m zeta "
               << "eMax eMin Q R nG <-rho value> <-nSub value> "
               << "<-stressScale value> <-pMin value> "
               << "<-tangentPMin value> <-pResidual value> "
               << "<-geostaticAdmission> <-stage 0|1> "
               << "<-initialStress sxx syy szz sxy syz sxz>" << endln;
        return 0;
    }

    int tag = 0;
    int count = 1;
    if (OPS_GetIntInput(&count, &tag) < 0) {
        opserr << "WARNING invalid RIVASand tag" << endln;
        return 0;
    }

    double values[requiredValues];
    count = requiredValues;
    if (OPS_GetDoubleInput(&count, values) < 0) {
        opserr << "WARNING invalid RIVASand material values for tag "
               << tag << endln;
        return 0;
    }

    double rho = 0.0;
    double stressScale = 1.0;
    double pMin = -1.0;
    double tangentPressureFloor = -1.0;
    double residualPressure = 0.0;
    bool geostaticAdmission = false;
    int fixedSubsteps = 1;
    int stage = 0;
    bool stageSpecified = false;
    bool initialStressSpecified = false;
    bool reversalLatch = false;
    bool admitOverbound = false;
    bool numericalTangent = false;
    bool noBiasVolume = false;
    bool recenterActivation = false;
    double recenterThreshold = 0.9;
    double betaFloor = 0.0;
    double betaCap = 0.0;
    bool noBiasHardening = false;
    bool noCyclicFlow = false;
    double l3[3] = {0.0, 0.0, 0.0};
    bool v11Enabled = false;
    double v11[3] = {2.0, 0.35, 1.0};
    double betaReserve = 0.0;
    double pMax = 0.0;
    double projectActivation = 0.0;
    Vector initialStress(6);
    initialStress.Zero();

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *option = OPS_GetString();
        if (std::strcmp(option, "-rho") == 0) {
            count = 1;
            if (OPS_GetDoubleInput(&count, &rho) < 0) {
                opserr << "WARNING invalid -rho for RIVASand tag "
                       << tag << endln;
                return 0;
            }
        } else if (std::strcmp(option, "-nSub") == 0 ||
                   std::strcmp(option, "-noSubsteps") == 0) {
            count = 1;
            if (OPS_GetIntInput(&count, &fixedSubsteps) < 0) {
                opserr << "WARNING invalid -nSub for RIVASand tag "
                       << tag << endln;
                return 0;
            }
        } else if (std::strcmp(option, "-stressScale") == 0) {
            count = 1;
            if (OPS_GetDoubleInput(&count, &stressScale) < 0) {
                opserr << "WARNING invalid -stressScale for RIVASand tag "
                       << tag << endln;
                return 0;
            }
        } else if (std::strcmp(option, "-pMin") == 0) {
            count = 1;
            if (OPS_GetDoubleInput(&count, &pMin) < 0 ||
                !std::isfinite(pMin) || !(pMin > 0.0)) {
                opserr << "WARNING invalid -pMin for RIVASand tag "
                       << tag << "; value must be positive" << endln;
                return 0;
            }
        } else if (std::strcmp(option, "-tangentPMin") == 0) {
            count = 1;
            if (OPS_GetDoubleInput(&count, &tangentPressureFloor) < 0 ||
                !std::isfinite(tangentPressureFloor) ||
                !(tangentPressureFloor > 0.0)) {
                opserr << "WARNING invalid -tangentPMin for "
                       << "RIVASand tag " << tag
                       << "; value must be positive" << endln;
                return 0;
            }
        } else if (std::strcmp(option, "-reversalLatch") == 0) {
            reversalLatch = true;
        } else if (std::strcmp(option, "-admitOverbound") == 0) {
            admitOverbound = true;
        } else if (std::strcmp(option, "-numericalTangent") == 0) {
            numericalTangent = true;
        } else if (std::strcmp(option, "-noBiasVolume") == 0) {
            // diagnostic switch: disable the bias reversible-volume block
            noBiasVolume = true;
        } else if (std::strcmp(option, "-recenterActivation") == 0) {
            recenterActivation = true;
        } else if (std::strcmp(option, "-recenterThreshold") == 0) {
            count = 1;
            if (OPS_GetDoubleInput(&count, &recenterThreshold) < 0 ||
                !std::isfinite(recenterThreshold) ||
                recenterThreshold < 0.0) {
                opserr << "WARNING invalid -recenterThreshold for RIVASand "
                       << "tag " << tag << "; value must be >= 0" << endln;
                return 0;
            }
        } else if (std::strcmp(option, "-betaFloor") == 0) {
            count = 1;
            if (OPS_GetDoubleInput(&count, &betaFloor) < 0 ||
                !std::isfinite(betaFloor) || betaFloor < 0.0) {
                opserr << "WARNING invalid -betaFloor for RIVASand tag "
                       << tag << "; value must be >= 0" << endln;
                return 0;
            }
        } else if (std::strcmp(option, "-betaCap") == 0) {
            count = 1;
            if (OPS_GetDoubleInput(&count, &betaCap) < 0 ||
                !std::isfinite(betaCap) || betaCap < 0.0) {
                opserr << "WARNING invalid -betaCap for RIVASand tag "
                       << tag << "; value must be >= 0" << endln;
                return 0;
            }
        } else if (std::strcmp(option, "-projectActivation") == 0) {
            count = 1;
            if (OPS_GetDoubleInput(&count, &projectActivation) < 0 ||
                !std::isfinite(projectActivation) ||
                projectActivation < 0.0 || projectActivation >= 1.0) {
                opserr << "WARNING invalid -projectActivation for RIVASand "
                       << "tag " << tag << "; margin must be in [0,1)"
                       << endln;
                return 0;
            }
        } else if (std::strcmp(option, "-pMax") == 0) {
            count = 1;
            if (OPS_GetDoubleInput(&count, &pMax) < 0 ||
                !std::isfinite(pMax) || pMax < 0.0) {
                opserr << "WARNING invalid -pMax for RIVASand tag "
                       << tag << "; value must be >= 0" << endln;
                return 0;
            }
        } else if (std::strcmp(option, "-betaReserve") == 0) {
            count = 1;
            if (OPS_GetDoubleInput(&count, &betaReserve) < 0 ||
                !std::isfinite(betaReserve) || betaReserve < 0.0) {
                opserr << "WARNING invalid -betaReserve for RIVASand tag "
                       << tag << "; value must be >= 0" << endln;
                return 0;
            }
        } else if (std::strcmp(option, "-pResidual") == 0) {
            count = 1;
            if (OPS_GetDoubleInput(&count, &residualPressure) < 0 ||
                !std::isfinite(residualPressure) || residualPressure < 0.0) {
                opserr << "WARNING invalid -pResidual for RIVASand tag "
                       << tag << "; value must be nonnegative" << endln;
                return 0;
            }
        } else if (std::strcmp(option, "-geostaticAdmission") == 0) {
            geostaticAdmission = true;
        } else if (std::strcmp(option, "-noBiasHardening") == 0) {
            noBiasHardening = true;   // research diagnostic
        } else if (std::strcmp(option, "-noCyclicFlow") == 0) {
            noCyclicFlow = true;      // research diagnostic
        } else if (std::strcmp(option, "-v11") == 0) {
            count = 3;
            if (OPS_GetDoubleInput(&count, v11) < 0 ||
                !std::isfinite(v11[0]) || v11[0] < 0.0 ||
                !std::isfinite(v11[1]) || v11[1] <= 0.0 || v11[1] > 1.0 ||
                !std::isfinite(v11[2]) || v11[2] <= 0.0 || v11[2] > 1.0) {
                opserr << "WARNING invalid -v11 (k etaFloor dirFloor) for "
                       << "RIVASand tag " << tag << endln;
                return 0;
            }
            v11Enabled = true;
        } else if (std::strcmp(option, "-l3") == 0) {
            count = 3;
            if (OPS_GetDoubleInput(&count, l3) < 0 ||
                !std::isfinite(l3[0]) || !std::isfinite(l3[1]) ||
                !std::isfinite(l3[2]) || l3[1] < 0.0 || l3[1] > 1.0 ||
                l3[2] < 0.0 || l3[2] > 1.0) {
                opserr << "WARNING invalid -l3 (epRef boostFloor cutFloor) "
                       << "for RIVASand tag " << tag << endln;
                return 0;
            }
        } else if (std::strcmp(option, "-stage") == 0) {
            count = 1;
            if (OPS_GetIntInput(&count, &stage) < 0) {
                opserr << "WARNING invalid -stage for RIVASand tag "
                       << tag << endln;
                return 0;
            }
            stageSpecified = true;
        } else if (std::strcmp(option, "-initialStress") == 0) {
            double stress[6];
            count = 6;
            if (OPS_GetDoubleInput(&count, stress) < 0) {
                opserr << "WARNING invalid -initialStress for RIVASand tag "
                       << tag << endln;
                return 0;
            }
            for (int i = 0; i < 6; ++i) initialStress(i) = stress[i];
            initialStressSpecified = true;
        } else {
            opserr << "WARNING unknown RIVASand option '" << option
                   << "' for tag " << tag << endln;
            return 0;
        }
    }

    if (initialStressSpecified && !stageSpecified) stage = 1;
    if (stage == 1 && !initialStressSpecified) {
        opserr << "WARNING RIVASand -stage 1 requires a compressive "
               << "-initialStress; otherwise create at stage 0, establish "
               << "geostatic stress, and use updateMaterialStage" << endln;
        return 0;
    }

    RIVASand *material = new RIVASand(
        tag, values[0], values[1], values[2], values[3], values[4],
        values[5], values[6], values[7], values[8], values[9], values[10],
        rho, fixedSubsteps, stressScale, pMin, tangentPressureFloor,
        residualPressure, geostaticAdmission, stage, initialStress);
    if (material != 0) material->setReversalLatch(reversalLatch);
    if (material != 0) material->setAdmitOverbound(admitOverbound);
    if (material != 0 && noBiasVolume) material->setBiasVolumeEnabled(false);
    if (material != 0) material->setNumericalTangent(numericalTangent);
    if (material != 0) material->setBetaFloor(betaFloor);
    if (material != 0) material->setBetaCap(betaCap);
    if (material != 0) material->setDiagNoBiasHardening(noBiasHardening);
    if (material != 0) material->setDiagNoCyclicFlow(noCyclicFlow);
    if (material != 0) material->setL3(l3[0], l3[1], l3[2]);
    if (material != 0) material->setV11(v11Enabled, v11[0], v11[1], v11[2]);
    if (material != 0) material->setBetaReserve(betaReserve);
    if (material != 0) material->setPMax(pMax);
    if (material != 0) material->setProjectActivation(projectActivation);
    if (material != 0) material->setRecenterActivation(recenterActivation);
    if (material != 0) material->setRecenterThreshold(recenterThreshold);
    if (material == 0 || !material->isValid()) {
        opserr << "WARNING invalid RIVASand material with tag "
               << tag << endln;
        delete material;
        return 0;
    }
    return material;
}

RIVASand::RIVASand(
    int tag, double Dr, double M, double kd, double h, double m,
    double zeta, double eMax, double eMin, double Q, double R, double nG,
    double rho, int fixedSubsteps, double stressScale, double pMin,
    double tangentPressureFloor, double residualPressure,
    bool geostaticAdmission, int initialStage, const Vector &initialStress)
    : NDMaterial(tag, ND_TAG_RIVASand),
      mDr(Dr), mRho(rho), mStressScale(stressScale),
      mTangentPressureFloor(0.0),
      mFixedSubsteps(fixedSubsteps), mStage(initialStage),
      mInitialStage(initialStage), mValid(true),
      mReversalLatch(false), mLatchValid(false), mLatchedReversal(0),
      mNumericalTangent(false), mFDPending(false), mFDReversal(0),
      mRecenterActivation(false), mRecenterThreshold(0.9), mBetaReserve(0.0), mProjectActivation(0.0),
      mInitialStress(6), mCommittedStrain(6), mTrialStrain(6),
      mCommittedStress(6), mTrialStress(6), mTangent(6, 6),
      mInitialTangent(6, 6), mStateOutput(RIVA_STATE_VALUE_COUNT),
      mScalarOutput(1)
{
    mInitialStress.Zero();
    if (initialStress.Size() == 6) mInitialStress = initialStress;
    setReferenceParameters();
    if (pMin > 0.0) mParameters.p_min = pMin;
    mParameters.p_residual = residualPressure;
    mParameters.geostatic_admission_enabled = geostaticAdmission ? 1 : 0;
    mTangentPressureFloor = tangentPressureFloor > 0.0 ?
        tangentPressureFloor : mParameters.p_ref/200.0;
    mTangentPressureFloor = riva_max(
        mTangentPressureFloor, mParameters.p_min);
    setMaterialParameters(M, kd, h, m, zeta, eMax, eMin, Q, R, nG);

    if (!(mStressScale > 0.0) || !(mRho >= 0.0) ||
        mFixedSubsteps < 1 || (mStage != 0 && mStage != 1) ||
        !std::isfinite(mDr) || mDr < 0.0 || mDr > 1.0 ||
        !std::isfinite(mParameters.p_min) || !(mParameters.p_min > 0.0) ||
        !std::isfinite(mParameters.p_residual) ||
        mParameters.p_residual < 0.0 ||
        !std::isfinite(mTangentPressureFloor) ||
        !(mTangentPressureFloor > 0.0) ||
        !finiteVector(mInitialStress) ||
        !riva_material_parameters_valid(&mParameters, &mMaterial)) {
        mValid = false;
    }

    const double shear = mParameters.E_ref/(2.0*(1.0+mParameters.nu));
    const double bulk = mParameters.E_ref/(3.0*(1.0-2.0*mParameters.nu));
    buildTangent(bulk, shear, mInitialTangent);
    revertToStart();
}

RIVASand::RIVASand()
    : NDMaterial(0, ND_TAG_RIVASand),
      mDr(0.0), mRho(0.0), mStressScale(1.0),
      mTangentPressureFloor(0.0), mFixedSubsteps(1),
      mStage(0), mInitialStage(0), mValid(false),
      mReversalLatch(false), mLatchValid(false), mLatchedReversal(0),
      mNumericalTangent(false), mFDPending(false), mFDReversal(0),
      mRecenterActivation(false), mRecenterThreshold(0.9), mBetaReserve(0.0), mProjectActivation(0.0),
      mInitialStress(6), mCommittedStrain(6), mTrialStrain(6),
      mCommittedStress(6), mTrialStress(6), mTangent(6, 6),
      mInitialTangent(6, 6), mStateOutput(RIVA_STATE_VALUE_COUNT),
      mScalarOutput(1)
{
    mInitialStress.Zero();
    mCommittedStrain.Zero();
    mTrialStrain.Zero();
    mCommittedStress.Zero();
    mTrialStress.Zero();
    mCommittedState = riva_state_t{};
    mTrialState = riva_state_t{};
    setReferenceParameters();
    mTangentPressureFloor = riva_max(
        mParameters.p_ref/200.0, mParameters.p_min);
    mMaterial = riva_reference_material_parameters(&mParameters);
    const double shear = mParameters.E_ref/(2.0*(1.0+mParameters.nu));
    const double bulk = mParameters.E_ref/(3.0*(1.0-2.0*mParameters.nu));
    buildTangent(bulk, shear, mInitialTangent);
    mTangent = mInitialTangent;
}

RIVASand::~RIVASand()
{
}

void
RIVASand::setReferenceParameters(void)
{
    mParameters = riva_reference_parameters(mStressScale);
}

void
RIVASand::setMaterialParameters(double M, double kd, double h,
    double m, double zeta, double eMax, double eMin, double Q, double R,
    double nG)
{
    mMaterial = riva_reference_material_parameters(&mParameters);
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
RIVASand::initialVoidRatio(void) const
{
    return riva_void_ratio_from_material_relative_density(&mMaterial, mDr);
}

int
RIVASand::activateFromCommittedStress(void)
{
    riva_tensor_t stress = stressToTensor(mCommittedStress);
    const double physicalPressure = riva_pressure(stress);
    if (!(physicalPressure > mParameters.p_min)) {
        if (!mParameters.geostatic_admission_enabled ||
            !std::isfinite(physicalPressure)) {
            opserr << "RIVASand tag " << this->getTag()
                   << " cannot enter stage 1: physical p' must exceed pMin"
                   << endln;
            return -1;
        }
        /* A no-tension cone cannot admit a slightly tensile stress left by
         * the host gravity solve.  Project only its mean component to pMin;
         * the following plastic-equilibration stage restores global balance. */
        stress = riva_sub(riva_dev(stress), riva_iso(mParameters.p_min));
    }
    riva_tensor_t stressUse = stress;
    if (mProjectActivation > 0.0) {
        // PM4Sand-style one-time initial stress reset: elastic gravity can
        // leave stress ratios beyond anything the plastic model can hold
        // statically (real soil redistributes DURING gravity ramp-up; the
        // elastic+switch shortcut cannot). Scale the deviatoric part once,
        // at activation, so every point starts inside its bounding surface
        // with margin; the equilibrium deficit is finite and the transient
        // redistributes it through a well-posed system. Response stays
        // continuous from the corrected state ??unlike per-call projection.
        const double p0 = riva_max(riva_pressure(stressUse),
                                   mParameters.p_min);
        double mb = 0.0, md = 0.0, xi = 0.0;
        riva_surfaces(&mParameters, &mMaterial, p0, initialVoidRatio(),
                      &mb, &md, &xi);
        const double bound = sqrt(2.0/3.0)*mb;
        const double target = (1.0 - mProjectActivation)*bound;
        riva_tensor_t dev = riva_dev(stressUse);
        const double alphaNorm = riva_norm(dev)/p0;
        if (bound > 0.0 && alphaNorm > target) {
            dev = riva_scale(dev, target/alphaNorm);
            stressUse = riva_sub(dev, riva_iso(p0));
        }
    }
    riva_state_t state = {};
    if (!riva_initialize_material(&mParameters, &mMaterial, stressUse,
                                 initialVoidRatio(), &state)) return -1;
    if (mParameters.geostatic_admission_enabled) {
        double mb = 0.0, md = 0.0, xi = 0.0;
        riva_surfaces(&mParameters, &mMaterial, state.pressure_anchor,
                      state.void_ratio, &mb, &md, &xi);
        const double bound = sqrt(2.0/3.0)*mb;
        if (bound > 0.0 && riva_norm(state.alpha) > bound) {
            /* Preserve the equilibrated stress exactly.  Recenter only the
             * bounding-surface mapping origin so the inherited over-bound
             * state does not create a finite correction at zero strain. */
            state.alpha0 = state.alpha;
            state.alpha01 = state.alpha;
        }
    }
    riva_tensor_t reference = stressUse;
    reference.xy = reference.yz = reference.xz = 0.0;
    if (!riva_begin_dynamic_phase(&mParameters, reference, &state)) return -1;
    // Activation severity: how close this point's stage-in stress ratio sits
    // to its bounding surface. Level-ground K0 states (the calibration
    // condition) sit near 0.7; SSI foundation-edge and low-p' surface states
    // sit at or beyond 1.0. Both stage-in aids below apply ONLY to severe
    // points, so every calibrated loading path is untouched bit-for-bit.
    double severity = 0.0;
    {
        double mb = 0.0, md = 0.0, xi = 0.0;
        riva_surfaces(&mParameters, &mMaterial,
                      riva_max(riva_pressure(state.stress), mParameters.p_min),
                      state.void_ratio, &mb, &md, &xi);
        const double bound = sqrt(2.0/3.0)*mb;
        if (bound > 0.0) severity = riva_norm(state.alpha)/bound;
    }
    const bool severe = severity > mRecenterThreshold;
    if (mRecenterActivation && severe) {
        // PM4Sand FirstCall semantics: re-center the bounding-surface mapping
        // reference at the activation stress ratio, as if a reversal had just
        // occurred there. The kernel's default (alpha0 = 0) leaves such
        // states with little or no bounding distance ??over-bound foundation
        // states get beta <= 0 permanently until a real reversal.
        state.alpha0 = state.alpha;
        state.alpha01 = state.alpha;
    }
    // The hardening reserve is spatially UNGATED (the redistribution wave
    // reaches activation-benign points too) but stage-scoped in time: the
    // host zeroes it via the betaReserve parameter before the shaking stage,
    // so the frozen cyclic response applies where it was calibrated.
    mParameters.beta_reserve = mBetaReserve;
    mCommittedState = state;
    mTrialState = state;
    tensorToStress(state.stress, mCommittedStress);
    mTrialStress = mCommittedStress;
    return 0;
}

riva_tensor_t
RIVASand::strainIncrementToTensor(const Vector &increment)
{
    riva_tensor_t result = {increment(0), increment(1), increment(2),
                           0.5*increment(3), 0.5*increment(4),
                           0.5*increment(5)};
    return result;
}

riva_tensor_t
RIVASand::stressToTensor(const Vector &stress)
{
    riva_tensor_t result = {stress(0), stress(1), stress(2),
                           stress(3), stress(4), stress(5)};
    return result;
}

void
RIVASand::tensorToStress(riva_tensor_t tensor, Vector &stress)
{
    stress(0) = tensor.xx;
    stress(1) = tensor.yy;
    stress(2) = tensor.zz;
    stress(3) = tensor.xy;
    stress(4) = tensor.yz;
    stress(5) = tensor.xz;
}

void
RIVASand::buildTangent(double bulk, double shear, Matrix &matrix) const
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
RIVASand::updateTrialTangent(void)
{
    if (mStage == 1 && mTrialState.initialized) {
        double shear = 0.0;
        double bulk = 0.0;
        // Regularize only the tangent advertised to OpenSees here.  The
        // explicit V8 stress update continues to use mParameters.p_min, so
        // the default p_ref/200 tangent floor does not silently alter the
        // calibrated constitutive pressure path.
        const double pressure = riva_max(
            riva_cone_pressure(&mParameters, mTrialState.stress),
            mTangentPressureFloor);
        riva_moduli_for_state(&mParameters, &mMaterial, pressure,
                             &mTrialState, &shear, &bulk);
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
RIVASand::setTrialStrain(const Vector &strain)
{
    if (!mValid || strain.Size() != 6) return -1;
    if (!finiteVector(strain)) return -1;
    mTrialStrain = strain;
    mFDPending = false;
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

    if (!mCommittedState.initialized && activateFromCommittedStress() != 0)
        return -1;
    mTrialState = mCommittedState;
    riva_tensor_t stress = mTrialState.stress;
    riva_update_info_t information = {};
    const int reversalOverride =
        (mReversalLatch && mLatchValid) ? mLatchedReversal : -1;
    if (!riva_update_material_ex(
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
    if (mNumericalTangent) {
        mFDPending = true;
        mFDReversal = information.reversal_registered ? 1 : 0;
    }
    return 0;
}

int
RIVASand::setTrialStrain(const Vector &strain, const Vector &rate)
{
    return setTrialStrain(strain);
}

int
RIVASand::setTrialStrainIncr(const Vector &increment)
{
    if (increment.Size() != 6) return -1;
    Vector target = mTrialStrain;
    target += increment;
    return setTrialStrain(target);
}

int
RIVASand::setTrialStrainIncr(const Vector &increment,
                                  const Vector &rate)
{
    return setTrialStrainIncr(increment);
}

const Vector &
RIVASand::getStress(void)
{
    return mTrialStress;
}

const Vector &
RIVASand::getStrain(void)
{
    return mTrialStrain;
}

const Matrix &
RIVASand::getTangent(void)
{
    if (mFDPending) {
        mFDPending = false;
        computeNumericalTangent();
    }
    return mTangent;
}

void
RIVASand::computeNumericalTangent(void)
{
    // SANISAND-style consistent tangent, obtained by forward differences:
    // six extra kernel evaluations from the SAME committed state with the
    // SAME reversal decision as the accepted trial. Keeps the elastic
    // tangent (already in mTangent) on any FD failure, on implausibly large
    // trial increments (doomed Newton probes), or on badly scaled columns
    // that would poison the factorization.
    const double delta = 1.0e-8;
    const double maxIncrement = 1.0e-2;
    Vector increment = mTrialStrain - mCommittedStrain;
    for (int i = 0; i < 6; ++i)
        if (std::fabs(increment(i)) > maxIncrement) return;
    double elasticScale = 0.0;
    for (int i = 0; i < 6; ++i)
        elasticScale = riva_max(elasticScale, std::fabs(mTangent(i, i)));
    // Legitimate elastoplastic columns never exceed the elastic scale; a
    // column of 1e7+ means the +delta perturbation crossed a kernel branch
    // (projection cap, p_min floor, dilatancy switch) and the divided
    // difference measured the kink, not the stiffness. Reject those.
    const double entryLimit = 10.0*elasticScale;
    Matrix fd(6, 6);
    for (int j = 0; j < 6; ++j) {
        Vector strainP = mTrialStrain;
        strainP(j) += delta;
        Vector incP = strainP - mCommittedStrain;
        riva_state_t stP = mCommittedState;
        riva_tensor_t sP = stP.stress;
        riva_update_info_t infP = {};
        if (!riva_update_material_ex(
                &mParameters, &mMaterial, strainIncrementToTensor(incP),
                mFixedSubsteps, &stP, &sP, &infP, mFDReversal)) return;
        Vector sigP(6);
        tensorToStress(sP, sigP);
        for (int i = 0; i < 6; ++i) {
            const double value = (sigP(i) - mTrialStress(i))/delta;
            if (!std::isfinite(value) || std::fabs(value) > entryLimit)
                return;
            fd(i, j) = value;
        }
        // A near-zero (or negative) diagonal means this point is in a
        // mechanism direction; assembling it makes the global matrix
        // near-singular and the linear solve returns garbage corrections
        // (observed as sudden 1e8-1e68 residual explosions mid-descent).
        // Keep the elastic tangent for such points instead.
        if (fd(j, j) < 1.0e-4*mTangent(j, j)) return;
    }
    mTangent = fd;
}

const Matrix &
RIVASand::getInitialTangent(void)
{
    return mInitialTangent;
}

double
RIVASand::getRho(void)
{
    return mRho;
}

void
RIVASand::setReversalLatch(bool enabled)
{
    mReversalLatch = enabled;
    mLatchValid = false;
    mLatchedReversal = 0;
}

void
RIVASand::setAdmitOverbound(bool enabled)
{
    mParameters.admit_inherited_overbound = enabled ? 1 : 0;
}

void
RIVASand::setBiasVolumeEnabled(bool enabled)
{
    mParameters.bias_reversible_volume_enabled = enabled ? 1 : 0;
}

void
RIVASand::setNumericalTangent(bool enabled)
{
    mNumericalTangent = enabled;
}

void
RIVASand::setBetaFloor(double value)
{
    mParameters.beta_floor =
        (std::isfinite(value) && value > 0.0) ? value : 0.0;
}

void
RIVASand::setBetaCap(double value)
{
    mParameters.beta_cap =
        (std::isfinite(value) && value > 0.0) ? value : 0.0;
}

void
RIVASand::setDiagNoBiasHardening(bool enabled)
{
    mParameters.diag_no_bias_hardening = enabled ? 1 : 0;
}

void
RIVASand::setDiagNoCyclicFlow(bool enabled)
{
    mParameters.diag_no_cyclic_flow = enabled ? 1 : 0;
}

void
RIVASand::setL3(double epRef, double boostFloor, double cutFloor)
{
    mParameters.l3_ep_ref =
        (std::isfinite(epRef) && epRef > 0.0) ? epRef : 0.0;
    mParameters.l3_boost_floor = riva_clip(boostFloor, 0.0, 1.0);
    mParameters.l3_cut_floor = riva_clip(cutFloor, 0.0, 1.0);
}

void
RIVASand::setV11(bool enabled, double k, double etaFloor, double dirFloor)
{
    mParameters.v11_enabled = enabled ? 1 : 0;
    mParameters.v11_k =
        (std::isfinite(k) && k >= 0.0) ? k : 2.0;
    mParameters.v11_eta_floor = riva_clip(etaFloor, 1.0e-3, 1.0);
    mParameters.v11_dir_floor = riva_clip(dirFloor, 1.0e-3, 1.0);
}

void
RIVASand::setRecenterActivation(bool enabled)
{
    mRecenterActivation = enabled;
}

void
RIVASand::setProjectActivation(double value)
{
    mProjectActivation =
        (std::isfinite(value) && value > 0.0 && value < 1.0) ? value : 0.0;
}

void
RIVASand::setPMax(double value)
{
    mParameters.p_max =
        (std::isfinite(value) && value > 0.0) ? value : 0.0;
}

void
RIVASand::setRecenterThreshold(double value)
{
    mRecenterThreshold =
        (std::isfinite(value) && value >= 0.0) ? value : 0.9;
}

void
RIVASand::setBetaReserve(double value)
{
    mBetaReserve = (std::isfinite(value) && value > 0.0) ? value : 0.0;
    mParameters.beta_reserve = mBetaReserve;
}

int
RIVASand::commitState(void)
{
    mCommittedStrain = mTrialStrain;
    mCommittedStress = mTrialStress;
    mCommittedState = mTrialState;
    mLatchValid = false;
    mFDPending = false;
    return 0;
}

int
RIVASand::revertToLastCommit(void)
{
    mTrialStrain = mCommittedStrain;
    mTrialStress = mCommittedStress;
    mTrialState = mCommittedState;
    mLatchValid = false;
    mFDPending = false;
    updateTrialTangent();
    return 0;
}

int
RIVASand::revertToStart(void)
{
    mLatchValid = false;
    mFDPending = false;
    mStage = mInitialStage;
    mCommittedStrain.Zero();
    mTrialStrain.Zero();
    mCommittedStress = mInitialStress;
    mTrialStress = mInitialStress;
    mCommittedState = riva_state_t{};
    mTrialState = riva_state_t{};
    mTangent = mInitialTangent;
    if (mValid && mStage == 1 && activateFromCommittedStress() != 0)
        mValid = false;
    return mValid ? 0 : -1;
}

NDMaterial *
RIVASand::getCopy(void)
{
    return new RIVASand(*this);
}

NDMaterial *
RIVASand::getCopy(const char *code)
{
    if (std::strcmp(code, "ThreeDimensional") == 0 ||
        std::strcmp(code, "3D") == 0 ||
        std::strcmp(code, "RIVASand") == 0)
        return getCopy();
    return 0;
}

const char *
RIVASand::getType(void) const
{
    return "ThreeDimensional";
}

int
RIVASand::getOrder(void) const
{
    return 6;
}

bool
RIVASand::isValid(void) const
{
    return mValid;
}

int
RIVASand::sendSelf(int commitTag, Channel &theChannel)
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
    double stateValues[RIVA_STATE_VALUE_COUNT] = {};
    riva_state_values(&mCommittedState, stateValues);
    for (int i = 0; i < RIVA_STATE_VALUE_COUNT; ++i)
        data(35+i) = stateValues[i];
    // state block: data(35 .. 35+RIVA_STATE_VALUE_COUNT-1) = 35..129
    data(130) = mCommittedState.initialized;
    data(131) = mParameters.p_min;
    data(132) = mTangentPressureFloor;
    data(133) = mReversalLatch ? 1.0 : 0.0;
    data(134) = mParameters.admit_inherited_overbound;
    data(135) = mNumericalTangent ? 1.0 : 0.0;
    data(136) = mParameters.beta_floor;
    data(137) = mRecenterActivation ? 1.0 : 0.0;
    data(138) = mBetaReserve;
    data(139) = mRecenterThreshold;
    data(140) = mParameters.p_max;
    data(141) = mProjectActivation;
    data(142) = mParameters.p_residual;
    data(143) = mParameters.geostatic_admission_enabled;
    data(144) = mParameters.beta_cap;
    data(145) = mParameters.diag_no_bias_hardening;
    data(146) = mParameters.diag_no_cyclic_flow;
    data(147) = mParameters.l3_ep_ref;
    data(148) = mParameters.l3_boost_floor;
    data(149) = mParameters.l3_cut_floor;
    data(150) = mParameters.v11_enabled;
    data(151) = mParameters.v11_k;
    data(152) = mParameters.v11_eta_floor;
    data(153) = mParameters.v11_dir_floor;

    if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
        opserr << "RIVASand::sendSelf failed for tag "
               << this->getTag() << endln;
        return -1;
    }
    return 0;
}

int
RIVASand::restoreState(
    const double values[RIVA_STATE_VALUE_COUNT], int initialized,
    riva_state_t &state)
{
    int i = 0;
    state = riva_state_t{};
#define RIVA_RESTORE_TENSOR(t) do { \
    (t).xx=values[i++]; (t).yy=values[i++]; (t).zz=values[i++]; \
    (t).xy=values[i++]; (t).yz=values[i++]; (t).xz=values[i++]; \
} while (0)
    RIVA_RESTORE_TENSOR(state.stress);
    RIVA_RESTORE_TENSOR(state.alpha);
    RIVA_RESTORE_TENSOR(state.alpha0);
    RIVA_RESTORE_TENSOR(state.alpha01);
    RIVA_RESTORE_TENSOR(state.n);
    RIVA_RESTORE_TENSOR(state.fabric);
    state.D=values[i++]; state.beta=values[i++];
    state.lambda_total=values[i++]; state.ep_eq_since_reversal=values[i++];
    state.void_ratio=values[i++];
    state.reversals=(int64_t)std::llround(values[i++]);
    state.pressure_floor_hits=(int64_t)std::llround(values[i++]);
    state.denominator_floor_hits=(int64_t)std::llround(values[i++]);
    state.beta_fallbacks=(int64_t)std::llround(values[i++]);
    state.eps_v_total=values[i++]; state.eps_v_confining=values[i++];
    state.eps_v_irreversible=values[i++]; state.eps_v_reversible=values[i++];
    state.pressure_anchor=values[i++]; state.D_ir=values[i++];
    state.D_re=values[i++];
    RIVA_RESTORE_TENSOR(state.last_reversal_deviator);
    state.cyclic_amplitude=values[i++]; state.amplitude_factor=values[i++];
    state.amplitude_reversals=(int64_t)std::llround(values[i++]);
    state.initial_relative_state=values[i++];
    state.state_contraction_factor=values[i++];
    state.effective_knee_ratio=values[i++];
    RIVA_RESTORE_TENSOR(state.geostatic_deviator);
    RIVA_RESTORE_TENSOR(state.static_bias_tensor);
    RIVA_RESTORE_TENSOR(state.cyclic_direction);
    state.static_bias_index=values[i++];
    state.cyclic_phase_active=(int32_t)std::llround(values[i++]);
    state.bias_ratchet_strain=values[i++];
    state.physical_eps_v_total=values[i++];
    state.bias_reversible_volume=values[i++];
    RIVA_RESTORE_TENSOR(state.last_host_deviatoric_strain_direction);
    state.ep_half_last=values[i++];
    state.eta_m=values[i++];
#undef RIVA_RESTORE_TENSOR
    state.initialized = initialized;
    return i == RIVA_STATE_VALUE_COUNT ? 0 : -1;
}

int
RIVASand::recvSelf(int commitTag, Channel &theChannel,
                        FEM_ObjectBroker &theBroker)
{
    Vector data(RIVASerializedSize);
    if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
        opserr << "RIVASand::recvSelf failed" << endln;
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
    mParameters.p_min = data(131);
    mParameters.p_residual = data(142);
    mParameters.geostatic_admission_enabled =
        (int32_t)std::llround(data(143));
    mTangentPressureFloor = riva_max(data(132), mParameters.p_min);
    setMaterialParameters(data(13), data(14), data(15), data(16), data(17),
                          data(18), data(19), data(20), data(21), data(22));
    for (int i = 0; i < 6; ++i) mCommittedStrain(i) = data(23+i);
    for (int i = 0; i < 6; ++i) mCommittedStress(i) = data(29+i);
    double stateValues[RIVA_STATE_VALUE_COUNT];
    for (int i = 0; i < RIVA_STATE_VALUE_COUNT; ++i)
        stateValues[i] = data(35+i);
    if (restoreState(stateValues, (int)std::llround(data(130)),
                     mCommittedState) != 0) return -1;
    setReversalLatch(std::llround(data(133)) != 0);
    setAdmitOverbound(std::llround(data(134)) != 0);
    setNumericalTangent(std::llround(data(135)) != 0);
    setBetaFloor(data(136));
    setRecenterActivation(std::llround(data(137)) != 0);
    setBetaReserve(data(138));
    setRecenterThreshold(data(139));
    setPMax(data(140));
    setProjectActivation(data(141));
    setBetaCap(data(144));
    setDiagNoBiasHardening(std::llround(data(145)) != 0);
    setDiagNoCyclicFlow(std::llround(data(146)) != 0);
    setL3(data(147), data(148), data(149));
    setV11(std::llround(data(150)) != 0, data(151), data(152), data(153));

    mValid = mStressScale > 0.0 && mRho >= 0.0 && mFixedSubsteps >= 1 &&
        (mStage == 0 || mStage == 1) && mDr >= 0.0 && mDr <= 1.0 &&
        std::isfinite(mParameters.p_min) && mParameters.p_min > 0.0 &&
        std::isfinite(mParameters.p_residual) &&
        mParameters.p_residual >= 0.0 &&
        (mParameters.geostatic_admission_enabled == 0 ||
         mParameters.geostatic_admission_enabled == 1) &&
        std::isfinite(mTangentPressureFloor) &&
        mTangentPressureFloor > 0.0 &&
        finiteVector(mInitialStress) &&
        riva_material_parameters_valid(&mParameters, &mMaterial) &&
        (mStage == 0 || mCommittedState.initialized);
    const double shear = mParameters.E_ref/(2.0*(1.0+mParameters.nu));
    const double bulk = mParameters.E_ref/(3.0*(1.0-2.0*mParameters.nu));
    buildTangent(bulk, shear, mInitialTangent);
    mTrialStrain = mCommittedStrain;
    mTrialStress = mCommittedStress;
    mTrialState = mCommittedState;
    updateTrialTangent();
    return mValid ? 0 : -1;
}

const Vector &
RIVASand::getStateVector(void)
{
    double values[RIVA_STATE_VALUE_COUNT] = {};
    riva_state_values(&mTrialState, values);
    for (int i = 0; i < RIVA_STATE_VALUE_COUNT; ++i)
        mStateOutput(i) = values[i];
    return mStateOutput;
}

const Vector &
RIVASand::getScalarResponse(int responseID)
{
    if (responseID == 4) {
        mScalarOutput(0) = mTrialState.initialized ?
            mTrialState.void_ratio : initialVoidRatio();
    } else if (responseID == 5) {
        const double anchor = mTrialState.pressure_anchor-
            mParameters.p_residual;
        mScalarOutput(0) = mTrialState.initialized && anchor > 0.0 ?
            1.0-riva_pressure(mTrialState.stress)/anchor : 0.0;
    } else if (responseID == 6) {
        mScalarOutput(0) = (double)mTrialState.reversals;
    } else if (responseID == 7) {
        mScalarOutput(0) = mTrialState.initialized ?
            riva_compatibility_residual(&mTrialState) : 0.0;
    } else if (responseID == 8) {
        mScalarOutput(0) = mParameters.p_min;
    } else if (responseID == 9) {
        mScalarOutput(0) = mTangentPressureFloor;
    } else if (responseID == 10) {
        mScalarOutput(0) = (double)mStage;
    }
    return mScalarOutput;
}

Response *
RIVASand::setResponse(const char **argv, int argc, OPS_Stream &output)
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
RIVASand::getResponse(int responseID, Information &materialInfo)
{
    if (responseID == 1) return materialInfo.setVector(getStress());
    if (responseID == 2) return materialInfo.setVector(getStrain());
    if (responseID == 3) return materialInfo.setVector(getStateVector());
    if (responseID >= 4 && responseID <= 10)
        return materialInfo.setVector(getScalarResponse(responseID));
    return NDMaterial::getResponse(responseID, materialInfo);
}

int
RIVASand::setParameter(const char **argv, int argc,
                            Parameter &parameter)
{
    if (argc >= 1 && std::strcmp(argv[0], "numericalTangent") == 0) {
        // Stage-scoped like betaReserve: FD consistent tangent for the
        // implicit consolidation solve, elastic tangent for shaking.
        if (argc >= 2 && std::atoi(argv[1]) != this->getTag()) return -1;
        return parameter.addObject(NumericalTangentParameter, this);
    }
    if (argc >= 1 && std::strcmp(argv[0], "betaReserve") == 0) {
        // Stage-scoped hardening reserve: enabled for consolidation, then
        // zeroed via `setParameter -val 0 ... betaReserve` before shaking so
        // the frozen cyclic response applies. Accepts the tag-less
        // element-forwarded form and the tagged form alike.
        if (argc >= 2 && std::atoi(argv[1]) != this->getTag()) return -1;
        return parameter.addObject(BetaReserveParameter, this);
    }
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
RIVASand::updateParameter(int responseID, Information &information)
{
    if (responseID == StageParameter) {
        const int requested = information.theInt;
        if (requested != 0 && requested != 1) return -1;
        if (requested == mStage) return 0;
        mLatchValid = false;
        mFDPending = false;
        if (requested == 1) {
            if (activateFromCommittedStress() != 0) return -1;
            mStage = 1;
        } else {
            mStage = 0;
            mCommittedState = riva_state_t{};
            mTrialState = riva_state_t{};
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
    if (responseID == BetaReserveParameter) {
        const double requested = information.theDouble;
        if (!std::isfinite(requested) || requested < 0.0) return -1;
        setBetaReserve(requested);
        return 0;
    }
    if (responseID == NumericalTangentParameter) {
        const double requested = information.theDouble;
        if (!std::isfinite(requested)) return -1;
        setNumericalTangent(requested != 0.0);
        if (!mNumericalTangent) mFDPending = false;
        return 0;
    }
    return -1;
}

void
RIVASand::Print(OPS_Stream &output, int flag)
{
    output << "RIVASand, tag: " << this->getTag() << endln;
    output << "  Dr=" << mDr << " M=" << mMaterial.M
           << " kd=" << mMaterial.kd << " h=" << mMaterial.h
           << " m=" << mMaterial.m << " zeta=" << mMaterial.zeta
           << " eMax=" << mMaterial.e_max << " eMin=" << mMaterial.e_min
           << " Q=" << mMaterial.Q << " R=" << mMaterial.R
           << " nG=" << mMaterial.n_G << endln;
    output << "  stage=" << mStage << " nSub=" << mFixedSubsteps
           << " stressScale=" << mStressScale
           << " pMin=" << mParameters.p_min
           << " pResidual=" << mParameters.p_residual
           << " geostaticAdmission="
           << mParameters.geostatic_admission_enabled
           << " tangentPMin=" << mTangentPressureFloor
           << " parameterSHA=" << RIVA_PARAMETER_SHA256 << endln;
}
