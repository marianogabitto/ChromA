#include <iostream>
#include "Eigen/Dense"

using namespace Eigen;
using namespace std;

extern "C" {
  void FwdBwdAlg(
    double * initPiIN,
    double * transPiIN,
    double * SoftEvIN,
    double * respOUT,
    double * TransStateCountOUT,
    double * MargSeqProbOUT,
    int K,
    int T);
}

// ======================================================== Custom Type Defs
// ========================================================
typedef Array<double, Dynamic, Dynamic, RowMajor> Arr2D;
typedef Array<double, 1, Dynamic, RowMajor> Arr1D;
typedef Map<Arr2D> ExtArr2D;
typedef Map<Arr1D> ExtArr1D;


void FwdBwdAlg(
    double * initPiIN,
    double * transPiIN,
    double * SoftEvIN,
    double * respOUT,
    double * TransStateCountOUT,
    double * MargSeqProbOUT,
    int K,
    int T)
{
    // Prep input
    ExtArr1D initPi (initPiIN, K);
    ExtArr2D transPi (transPiIN, K, K);
    ExtArr2D SoftEv (SoftEvIN, T, K);

    // Prep output
    ExtArr2D resp (respOUT, T, K);
    ExtArr2D TransStateCount (TransStateCountOUT, K, K);
    ExtArr1D lsum (MargSeqProbOUT, 1);

    // Temporary Arrays
    Arr2D fwdMsg = ArrayXXd::Zero(T, K);
    Arr2D bwdMsg = ArrayXXd::Zero(T, K);
    Arr1D margPrObs = ArrayXd::Zero(T);
    Arr2D respPair_t = ArrayXXd::Zero(K, K);

    // Base case update for first time-step
    fwdMsg.row(0) = initPi * SoftEv.row(0);
    margPrObs(0) = fwdMsg.row(0).sum();
    fwdMsg.row(0) /= margPrObs(0);

    // Recursive update of timesteps 1, 2, ... T-1
    // Note: fwdMsg.row(t) is a *row vector*
    //       so needs to be left-multiplied to square matrix transPi
    for (int t = 1; t < T; t++) {
        fwdMsg.row(t) = fwdMsg.row(t-1).matrix() * transPi.matrix();
        // fwdMsg.row(t) = transPi.matrix() * fwdMsg.row(t-1);
        fwdMsg.row(t) *= SoftEv.row(t);
        margPrObs(t) = fwdMsg.row(t).sum();
        fwdMsg.row(t) /= margPrObs(t);
    }

    // Base case update for last time-step
    bwdMsg.row(T-1).fill(1.0);

    // Recursive update of timesteps T-2, T-3, ... 3, 2, 1, 0
    // Note: bMsg.row(t) is a *row vector*
    //       so needs to be left-multiplied to square matrix transPi.T
    for (int t = T-2; t >= 0; t--) {
        bwdMsg.row(t) = (bwdMsg.row(t+1) * SoftEv.row(t+1)).matrix() \
                      * transPi.transpose().matrix();
        bwdMsg.row(t) /= margPrObs(t+1);
    }

    for (int t = 0; t < T; t++) {
        resp.row(t) = fwdMsg.row(t) * bwdMsg.row(t);
    }

    for (int t = 1; t < T; t++) {
        // In Python, we want:
        // >>> respPair[t] = np.outer(fmsg[t-1], bmsg[t] * SoftEv[t])
        // >>> respPair[t] *= PiMat / margPrObs[t]
        respPair_t = fwdMsg.row(t-1).transpose().matrix() \
                      * (bwdMsg.row(t) * SoftEv.row(t)).matrix();
        respPair_t *= transPi;
        respPair_t /= margPrObs(t);

        // Aggregate pairwise transition counts
        TransStateCount += respPair_t;
    }

    lsum = margPrObs.log().sum();
    // std::cout << "Sum n:" << lsum << std::endl;

}
