/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
** ****************************************************************** */

#ifndef RIVASandIntermediateBiasResearch_h
#define RIVASandIntermediateBiasResearch_h

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>

#include "RIVASandIntermediateBiasResearchKernel.h"

class Channel;
class FEM_ObjectBroker;
class Information;
class Parameter;
class Response;

class RIVASandIntermediateBiasResearch : public NDMaterial
{
public:
    RIVASandIntermediateBiasResearch(int tag, double Dr, double M, double kd, double h,
                  double m, double zeta, double eMax, double eMin,
                  double Q, double R, double nG, double rho,
                  int fixedSubsteps, double stressScale, double pMin,
                  double tangentPressureFloor, double residualPressure,
                  bool geostaticAdmission, int initialStage,
                  const Vector &initialStress);
    RIVASandIntermediateBiasResearch();
    virtual ~RIVASandIntermediateBiasResearch();

    int setTrialStrain(const Vector &strain);
    int setTrialStrain(const Vector &strain, const Vector &rate);
    int setTrialStrainIncr(const Vector &strain);
    int setTrialStrainIncr(const Vector &strain, const Vector &rate);

    const Vector &getStress(void);
    const Vector &getStrain(void);
    const Matrix &getTangent(void);
    const Matrix &getInitialTangent(void);
    double getRho(void);

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    NDMaterial *getCopy(void);
    NDMaterial *getCopy(const char *code);
    const char *getType(void) const;
    int getOrder(void) const;

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker);

    Response *setResponse(const char **argv, int argc, OPS_Stream &output);
    int getResponse(int responseID, Information &materialInfo);
    int setParameter(const char **argv, int argc, Parameter &parameter);
    int updateParameter(int responseID, Information &information);
    void Print(OPS_Stream &output, int flag = 0);

    bool isValid(void) const;
    void setReversalLatch(bool on) { mReversalLatch = on; }
    void setFieldBiasMeanCorrection(bool on) {
        mParameters.field_bias_mean_correction_enabled = on ? 1 : 0;
    }

private:
    enum {
        StageParameter = 1,
        SubstepParameter = 2
    };

    void setReferenceParameters(void);
    void setMaterialParameters(double M, double kd, double h, double m,
                               double zeta, double eMax, double eMin,
                               double Q, double R, double nG);
    int activateFromCommittedStress(void);
    int beginDynamicFromCommittedState(void);
    void buildTangent(double bulk, double shear, Matrix &matrix) const;
    void updateTrialTangent(void);
    double initialVoidRatio(void) const;
    const Vector &getStateVector(void);
    const Vector &getScalarResponse(int responseID);

    static riva_ib_native::tensor_t strainIncrementToTensor(const Vector &increment);
    static riva_ib_native::tensor_t stressToTensor(const Vector &stress);
    static void tensorToStress(riva_ib_native::tensor_t tensor, Vector &stress);
    static int restoreState(const double values[RIVA_IB_STATE_VALUE_COUNT],
                            int initialized, int geostaticAdmitted,
                            riva_ib_native::riva_ib_state_t &state);

    double mDr;
    double mRho;
    double mStressScale;
    double mTangentPressureFloor;
    int mFixedSubsteps;
    int mStage;
    int mInitialStage;
    bool mValid;
    bool mGeostaticAdmission;
    /* Freeze one host-level reversal decision across the trial evaluations
     * belonging to a single OpenSees load step. The enabled setting is
     * serialized; the transient decision is cleared at every committed or
     * restored state. */
    bool mReversalLatch;
    bool mLatchValid;
    int mLatchedReversal;

    riva_ib_native::riva_ib_parameters_t mParameters;
    riva_ib_native::riva_material_parameters_t mMaterial;
    riva_ib_native::riva_ib_state_t mCommittedState;
    riva_ib_native::riva_ib_state_t mTrialState;

    Vector mInitialStress;
    Vector mCommittedStrain;
    Vector mTrialStrain;
    Vector mCommittedStress;
    Vector mTrialStress;
    Matrix mTangent;
    Matrix mInitialTangent;
    Vector mStateOutput;
    Vector mScalarOutput;
};

void *OPS_RIVASandIntermediateBiasResearchMaterial(void);

#endif
