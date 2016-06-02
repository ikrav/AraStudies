#ifndef PRECISERADIALRAYTRACER_H
#define PRECISERADIALRAYTRACER_H

#include "iostream"
#include "assert.h"
#define _USE_MATH_DEFINES
#include "math.h"
#include "float.h"
#include <string>
#include <vector>
#include "gmp.h"
#include "assert.h"

#include "TMath.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH3.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TArc.h"
#include "TStopwatch.h"

#include "ASlibEdit.h"

namespace PRECISERADIALRAYTRACER_NAMESPACE {
  const int    MaxSoln = 3; 
  const double C = .299792458;
  const double AngStep = 1e-60;
  
  const double DistStep = 1.0;
  const double DS = 1.0;
  
  const int    maxSurfaceFindIter = 50;
  const double earthRadius = 6371000.0;
  const double r2d = (180.0/TMath::Pi());

  //const double fDefaultTolerance = 0.00001;
  //const double fDefaultMaxAngChange = 0.0001;

  const double fDefaultTolerance = 0.1; // How close to target final ray must be
  const double fDefaultMaxAngChange = 0.0001; // Maximum possible angular change (excluding reflections) for a single step

};

class PreciseRadialRayTracer {

 private:

  enum verbosity {NONE, LOW, MEDIUM, HIGH};

  //debugging

  double fEntryOrExitGlobalAngle;
  double fEntryOrExitXYDist;
  double fEntryOrExitTime;

  //

  double fOrigin[3];
  double fDestination[3];
  int    fMaxSolution;
  int    fNCuts;
  int    iCase;
  int    fVerbose;
  int    fNSolution;
  int    fSuccess;

  double fMaxAngChange;
  double fTolerance;


  //  PRT Tracking

  double fGlobalAng;
  double fLocalAng;
  double fTempRecieveAngle;

  double fGlobalZ;
  double fGlobalR;

  double fXYDist;
  double fDistance;
  double fTime;

  double fXYLimit;
  double fDistLimit;
  bool fReflection;
  bool fTopReflection;

  double fMin;
  double fAbsMin;

  bool fDraw;
  bool fDeBug;
  int  fIVolley;

  //index of ref
  
  double fIoRA;
  double fIoRB;
  double fIoRC;
  
  std::vector<double> fFoundAngleVector;
  std::vector<double> fTrajectoryDistance; 
  std::vector<double> fTravelTime;
  std::vector<double> fRecieveAngle;
  std::vector<bool>   fReflectionVector;

  std::vector<double> fMinAngleVector;
  std::vector<double> fMaxAngleVector;
  std::vector<int> fMode;

  std::vector<double> fMinAry; 
  std::vector<double> fAngAry; 
  std::vector<bool> fTopRefAry;
  std::vector<bool> fReflectionAry;

 public:

  PreciseRadialRayTracer();
  
  inline void    ChangeIceModel(double A, double B, double C){fIoRA=A;fIoRB=B;fIoRC=C;};
  inline int     GetSuccess() {return fSuccess;};
  inline int     GetNSolutions() {return fNSolution;};
  inline void    SetVerbosity(int verb) {fVerbose=verb;};
  inline double  GetRecieveAngle(int iSol) {return fRecieveAngle[iSol];};
  inline double  GetRecieveTime(int iSol) {return fTravelTime[iSol];};
  inline bool    GetReflectionOccur(int iSol) {return fReflectionVector[iSol];};
  inline double  GetEntryOrExitXYDistance(){ return fEntryOrExitXYDist;}
  inline double  GetEntryOrExitTime(){ return fEntryOrExitTime;}

  void    SetDestination(double x,double y,double z);
  void    SetOrigin( double x, double y, double z);

  int            AdvanceRay(int stopCond, double stoptime=0, bool debug=false);
  double         BinarySearch(double minTheta, double maxTheta, int nSolution, int trialNum);
  int            FindPointsOfInterest();
  double         GetLowestTime();
  double         GetLowestDistance();
  void           GetCase(double &minTheta, double maxTheta);
  double         GetDrop(double minTheta);
  double         NofZ2(double z);
  bool           TraceRay(bool draw=false);
  int            Volley(double TrajectoryAngle, bool debug, bool draw, int color);
  void           TimeVolley(double TrajectoryAngle, double stoptime, int color);
  void           FindDownCoordsFromUpRay(double &xystep, double &zstep, double &diststep, bool debug);
  void           FindDownCoords(double &xystep, double &zstep, double &diststep, bool debug);
  void           FindUpCoords(double &xystep, double &zstep, double &diststep, bool debug);
  double         GetAngularChange(double xcompearth, double ycompearth, double xcompray, double ycompray, double IncidentIndex, double RefractionIndex, bool &reflection);
  void IncreaseAnglePrecision();
  void Default();

};

#endif
