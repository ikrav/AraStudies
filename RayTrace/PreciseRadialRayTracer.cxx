#include "PreciseRadialRayTracer.h"

using namespace PRECISERADIALRAYTRACER_NAMESPACE;

PreciseRadialRayTracer::PreciseRadialRayTracer() {

  fOrigin[0]= 999.;
  fOrigin[1]= 999.;
  fOrigin[2]= 999.;

  fDestination[0] = 999.;
  fDestination[1] = 999.;
  fDestination[2] = 999.;

  fNCuts = 10;
  fMaxSolution = 8;

  fFoundAngleVector.resize(fMaxSolution);
  fMinAngleVector.resize(fMaxSolution);
  fMode.resize(fMaxSolution);
  fMaxAngleVector.resize(fMaxSolution);
  fTrajectoryDistance.resize(fMaxSolution);
  fTravelTime.resize(fMaxSolution);
  fRecieveAngle.resize(fMaxSolution);
  fReflectionVector.resize(fMaxSolution);

  fMinAry.resize(fNCuts+3);
  fAngAry.resize(fNCuts+3);
  fTopRefAry.resize(fNCuts+3);
  fReflectionAry.resize(fNCuts+3);

  fVerbose = NONE;
  
  fIoRA = 1.78;
  fIoRB = 0.4272;
  fIoRC = 0.016;

  fMaxAngChange = fDefaultMaxAngChange;
  fTolerance = fDefaultTolerance;
  
}

void PreciseRadialRayTracer::Default() {
  fMaxAngChange = fDefaultMaxAngChange;
  fTolerance = fDefaultTolerance;
}
void PreciseRadialRayTracer::IncreaseAnglePrecision() {
  fMaxAngChange = fMaxAngChange*0.5;
}
//
//
//
//
//
//
void    PreciseRadialRayTracer::SetDestination(double x,double y,double z){
  fDestination[0]=x;
  fDestination[1]=y;
  fDestination[2]=z;
  //printf("Setting Destination %.3f %.3f %.3f -> %.3f %.3f %.3f \n",x,y,z,fDestination[0],fDestination[1],fDestination[2]);
}
//
//
//
//
//
//
void    PreciseRadialRayTracer::SetOrigin( double x, double y, double z) { 
  fOrigin[0]=x;
  fOrigin[1]=y;
  fOrigin[2]=z;
  //printf("Setting Origin %.3f %.3f %.3f -> %.3f %.3f %.3f \n",x,y,z,fOrigin[0],fOrigin[1],fOrigin[2]);
}
//
//
//
//
//
//
double PreciseRadialRayTracer::BinarySearch(double minTheta, double maxTheta, int nSolution, int trialNum) {

  double tempMinTheta = minTheta;
  double tempMaxTheta = maxTheta;

  double halfAngle = (tempMinTheta+tempMaxTheta)/2;
  
  if(fVerbose>=MEDIUM)printf("Mode = %i \n",fMode[trialNum]);

  for(int i=0; i<50; i++) {
    
    halfAngle = (tempMinTheta+tempMaxTheta)/2;
    /*if(i==0) Volley(halfAngle,_xyDist,RS,true);
      else Volley(halfAngle,false,fDraw,kBlue);*/
    
    //if(i==29) Volley(halfAngle,true,fDraw,kBlue);
    //else Volley(halfAngle,false,fDraw,kBlue);

    /*if(i==29) Volley(halfAngle,true,false,kBlue);
      else*/ 
    int failCheck = Volley(halfAngle,false,false,kBlue);

    
    if(fVerbose>=LOW) printf("Trial: Angles: %.3f - %.3f - %.3f  Min: %.9f Time: %.3f XY: %.4f Distance: %.8f TopRef: %d Ref: %d \n",tempMinTheta*(180/TMath::Pi()),halfAngle*(180/TMath::Pi()),tempMaxTheta*(180/TMath::Pi()),fMin,fTime,fXYDist,fDistance,fTopReflection,fReflection);
    
    if( fabs(tempMaxTheta-tempMinTheta)<AngStep || failCheck==99 ) { 
      if(fVerbose>=LOW){printf("Unable to resolve angle \n");} 
      return 999;
    }
    
    if(fAbsMin<fTolerance) {
      /*if(iCase==1 || iCase==2) {
	fTime += fMin*(NofZ2(fDestination[2])/C);
	fDistance += fMin;
      }
      if(iCase==3 || iCase==4) {
	fTime -= fMin*(NofZ2(fDestination[2])/C);
	fDistance -= fMin;
	}*/

      fTrajectoryDistance[nSolution] = fDistance;
      fTravelTime[nSolution] = fTime;
      fRecieveAngle[nSolution] = fTempRecieveAngle;
      if(fTopReflection || fReflection) fReflectionVector[nSolution] = true;
      else fReflectionVector[nSolution] = false;

      // DEBUG
      //Volley(halfAngle,true,false,kBlue);
      // *****
      if(fDraw) Volley(halfAngle,false,fDraw,kRed);
      return halfAngle;
    }
    
    if(fMode[trialNum]==1) {
      if(fMin<0) tempMaxTheta = halfAngle;
      else if(fMin>0) tempMinTheta = halfAngle;
      else { printf("Unresolved situation,  Min = %.9f AbsMin = %.9f \n",fMin,fAbsMin); return 999;; }
    }
    
    else if(fMode[trialNum]==2) {
      if(fMin>0) tempMaxTheta = halfAngle;
      else if(fMin<0) tempMinTheta = halfAngle;
      else { printf("Unresolved situation 2,  \n"); return 999;; }
    }

    else if(fMode[trialNum]==3) {
      if(fMin<0 && fTopReflection) tempMinTheta = halfAngle;
      else if(fMin<0) tempMaxTheta = halfAngle;
      else if(fMin>0) tempMinTheta = halfAngle;
      else { printf("Unresolved situation 4,  \n"); return 999;; }
    }

    else if(fMode[trialNum]==4) {
      if(fMin>0 && !fTopReflection) tempMinTheta = halfAngle;
      else if(fMin>0) tempMaxTheta = halfAngle;
      else if(fMin<0) tempMinTheta = halfAngle;
      else { printf("Unresolved situation 5,  \n"); return 999;; }
    }

    else if(fMode[trialNum]==5) {
      if(fMin>0 && fReflection) tempMinTheta = halfAngle;
      else if(fMin<0 && fReflection) tempMinTheta = halfAngle;
      else if(fMin>0 && !fReflection) tempMinTheta = halfAngle;
      else if(fMin<0 && !fReflection) tempMaxTheta = halfAngle;
      else { printf("Unresolved situation 5,  \n"); return 999; }
    }

    else { printf("Unresolved situation 3,  \n"); return 999;; }

    //sleep(1);

  }

  return 999;

}
//
//
//
//
//
//
int PreciseRadialRayTracer::FindPointsOfInterest() {

  int nPoints = 0;
  
  for(int i=0; i<fNCuts; i++) {

    double first = fMinAry[i];
    double second = fMinAry[i+1];
    
    bool mode1check = ((first>0) && (second<0));
    bool mode2check = ((first<0) && (second>0));
    bool mode3check = ( ((first<0) && (second<0)) && (fTopRefAry[i] && !fTopRefAry[i+1]) );
    bool mode4check = ( ((first>0) && (second>0)) && (!fTopRefAry[i] && fTopRefAry[i+1]) );
    bool mode5check = ( (!fTopRefAry[i] && !fTopRefAry[i+1]) && (fReflectionAry[i] && !fReflectionAry[i+1]) );

    if( mode1check || mode2check || mode3check || mode4check || mode5check) {
  
      if(mode1check) fMode[nPoints] = 1;
      else if(mode2check) fMode[nPoints] = 2;
      else if(mode3check) fMode[nPoints] = 3;
      else if(mode4check) fMode[nPoints] = 4;
      else if(mode5check) fMode[nPoints] = 5;
      else fMode[nPoints] = 0;

      fMinAngleVector[nPoints] = fAngAry[i];
      fMaxAngleVector[nPoints] = fAngAry[i+1];

      //printf("First: %.3f  Second: %.3f  Theta: %.3f-%.3f",first,second,fAngAry[i],fAngAry[i+1]);cout << "FirRef: " << fTopRefAry[i] << " SecRef: " << fTopRefAry[i+1] << endl;


      nPoints++;  
    }
    
  }
  
  return nPoints;

}
//
//
//
//
//
//
double PreciseRadialRayTracer::GetLowestTime() {

  if(fNSolution==0) return 0;

  double min = fTravelTime[0];

  for(int i=0; i<fNSolution; i++) if(fTravelTime[i]<min) min = fTravelTime[i];
  
  return min;

}
//
//
//
//
//
//
double PreciseRadialRayTracer::GetLowestDistance() {

  if(fNSolution==0) return 0;

  double min = fTrajectoryDistance[0];

  for(int i=0; i<fNSolution; i++) if(fTrajectoryDistance[i]<min) min = fTrajectoryDistance[i];
  
  return min;

}
//
//
//
//
//
//
double PreciseRadialRayTracer::NofZ2(double z) {

  if(z>=0) return 1.0;
  
  double n = fIoRA - fIoRB * exp(fIoRC*z);
  
  return n;
}
//
//
//
//
//
//  
double PreciseRadialRayTracer::GetDrop(double minTheta) {
  double xyshift = fabs(fOrigin[2]*tan(minTheta));
  double theta, phi, r;
  ConvertCartesian2Spherical(0.0,0.0,0.0,0.0,xyshift,earthRadius,theta,phi,r);
  double xshift, yshift, zshift;
  ConvertSpherical2Cartesian(0.0,0.0,0.0,theta,0.0,earthRadius,xshift,yshift,zshift);
  return (earthRadius-zshift);
}
//
//
//
//
//
//  
void PreciseRadialRayTracer::GetCase(double &minTheta, double maxTheta) {

  double theta,phi,r;

  //Case 1 (Source above ice, destination below ice)

  if( (fOrigin[2]>fDestination[2]) && ( (fOrigin[2]>0) && (fDestination[2]<0) ) ) {
    ConvertCartesian2Spherical(fOrigin[0],fOrigin[1],fOrigin[2],fDestination[0],fDestination[1],0.0,theta,phi,r);
    minTheta = theta;
    if(fVerbose>=MEDIUM) printf("Case 1 Used. Min/Max Theta = %.3f/%.3f \n",minTheta*(180/(TMath::Pi())),maxTheta*(180/(TMath::Pi())));
    iCase = 1;
  }

  //Case 2 (Both below ice, Source above destination)

  else if( ((fOrigin[2]<0) && (fDestination[2]<0)) && (fOrigin[2]>fDestination[2])) {
    minTheta = asin(1/NofZ2(fOrigin[2]));
    //printf("Index = %.9f \n",NofZ2(fOrigin[2]));
    //double zshift = GetDrop(minTheta);
    //minTheta = asin(1/NofZ2(fOrigin[2]+zshift));
    if(fVerbose>=MEDIUM) printf("Case 2 Used. Min/Max Theta = %.3f/%.3f \n",minTheta*(180/(TMath::Pi())),maxTheta*(180/(TMath::Pi())));
    iCase = 2;
  }

  //Case 3 (Source below ice, destination above ice)

  else if( (fOrigin[2]<fDestination[2]) && ( (fOrigin[2]<0) && (fDestination[2]>0) ) ) {
    //minTheta = (sin(maxTheta))/NofZ2(fOrigin[2]);
    minTheta = asin((sin(maxTheta))/NofZ2(fOrigin[2]));
    if(fVerbose>=MEDIUM) printf("Case 3 Used. Min/Max Theta = %.3f/%.3f \n",minTheta*(180/(TMath::Pi())),maxTheta*(180/(TMath::Pi())));
    iCase = 3;
  }

  //Case 4 (Both below ice, Source below destination)

  else if( ((fOrigin[2]<0) && (fDestination[2]<0)) && (fOrigin[2]<fDestination[2])) {
    minTheta = 0.0;/*asin(1/NofZ2(fOrigin[2]));
    //double zshift = GetDrop(minTheta);
    //minTheta = asin(1/NofZ2(fOrigin[2]+zshift));
    //minTheta = asin((sin(maxTheta))/NofZ2(fOrigin[2]));
    if(minTheta>maxTheta) minTheta = 0.0;*/
    if(fVerbose>=MEDIUM) printf("Case 4 Used. Min/Max Theta = %.3f/%.3f \n",minTheta*(180/(TMath::Pi())),maxTheta*(180/(TMath::Pi())));
    iCase = 4;
  }

  //Case 5 (Both above ice)

  else {
    minTheta = maxTheta;
    if(fVerbose>=MEDIUM) printf("Case 5 Used. Min/Max Theta = %.3f/%.3f \n",minTheta*(180/(TMath::Pi())),maxTheta*(180/(TMath::Pi())));iCase=5;
  }

}
//
//
//
//
//
//
bool PreciseRadialRayTracer::TraceRay(bool draw, &n) {

  fIVolley = 0;
  fDraw = draw;
  //if(fDraw) fVerbose = HIGH;
  fSuccess = 0;

  fXYLimit = Dist2d(fOrigin[0],fOrigin[1],fDestination[0],fDestination[1]);
  fDistLimit = Dist3d(fOrigin[0],fOrigin[1],fOrigin[2],fDestination[0],fDestination[1],fDestination[2]);

  double theta,phi,r;

  ConvertCartesian2Spherical(fOrigin[0],fOrigin[1],fOrigin[2],fDestination[0],fDestination[1],fDestination[2],theta,phi,r);
  //printf("Origin: %.3f %.3f %.3f Destination: %.3f %.3f %.3f \n",fOrigin[0],fOrigin[1],fOrigin[2],fDestination[0],fDestination[1],fDestination[2]);
  double maxTheta = theta+(1e-4);
  double minTheta = 0;

  GetCase(minTheta,maxTheta);
  double cutStep = (maxTheta-minTheta)/fNCuts;

  if(draw) printf("Case: %i  MinTheta = %.4f  MaxTheta = %.4f  Step = %.4f \n",iCase,minTheta,maxTheta,cutStep);

  for(int iCut=0; iCut<=fNCuts && iCase!=5; iCut++) {

    double testTheta = minTheta + (iCut*cutStep);

    //if(iCut==10) { fDeBug = true; Volley(testTheta,true,fDraw,kBlack); fDeBug = false; }
    //else Volley(testTheta,false,fDraw,kBlack);

    Volley(testTheta,false,false,kBlack);
    
    fMinAry[iCut] = fMin;
    fAngAry[iCut] = testTheta;
    fTopRefAry[iCut] = fTopReflection;
    fReflectionAry[iCut] = fReflection;

    if(iCut==fNCuts && fMinAry[iCut]==0.0) fMinAry[iCut] = -(1e-30);

    if(fVerbose>=MEDIUM) {std::cout << "Top Reflection = " << fTopReflection << " Reflection = " << fReflection; printf(" |   Angle = %.3f  Min = %.3f \n",testTheta*(180/(TMath::Pi())),fMin); }

    }

  //fVerbose = NONE;
  //fDraw = false;

  int nPointsFound = 0;
  if(iCase!=5) nPointsFound = FindPointsOfInterest();
  else if(iCase==5) nPointsFound = 1;
  int nSol = 0;
  
  if(fVerbose>=MEDIUM) printf("Points Found: %i \n",nPointsFound);
  
  for(int i=0; i<nPointsFound && iCase!=5; i++) {
  
    //if(nPointsFound>1 && i==0 && iCase==2) i++;
  
    fFoundAngleVector[nSol] = BinarySearch(fMinAngleVector[i],fMaxAngleVector[i],nSol,i);

    if(fFoundAngleVector[nSol]!=999) { 
      fTravelTime[nSol] += fMin*(NofZ2(fDestination[2])/C);
      if(fVerbose>=LOW) printf("Search Angle %i: %.5f - %.5f Found = %.5f Travel Time = %.5f Dist = %.5f RecieveAng: %.4f min = %.8f \n",i, fMinAngleVector[i]*(180/(TMath::Pi())), fMaxAngleVector[i]*(180/(TMath::Pi())), fFoundAngleVector[nSol]*(180/(TMath::Pi())),fTravelTime[nSol],fTrajectoryDistance[nSol],fRecieveAngle[nSol]*(180/(TMath::Pi())),fMin);
      fSuccess = 1;
      nSol++;
    }
  }
  
  fNSolution = nSol;

  if(iCase==5) {
    fNSolution = 1;
    Volley(maxTheta,false,false,kBlack);
  }

	std::cout << fNSolution << std::endl;
	n = fNSolution;

  if(fSuccess==1) return true;
  else return false;
  
}
//
//
//
//
//
//
int PreciseRadialRayTracer::Volley(double TrajectoryAngle, bool debug, bool draw, int color)  {
  
  //if(fVerbose==LOW){ fIVolley++; printf("Volley(%i)\n",fIVolley);}
  //if(debug) fVerbose = HIGH;
  //  good = (fXYDist<fXYLimit) && (fDistance<(2*fXYLimit));

  vector<double> xy;
  vector<double> z;

  fGlobalAng = 0.0;
  fLocalAng = TrajectoryAngle;

  fGlobalZ = fOrigin[2] + earthRadius;
  fGlobalR = fGlobalZ;

  fXYDist = 0.0;
  fDistance = 0.0;
  fTime = 0.0;

  fReflection = false;
  fTopReflection = false;

  int good = 0;
  
  for(int i=0; good==0; i++) { 
    if(draw) {
      xy.push_back(fXYDist);
      z.push_back(fGlobalZ-earthRadius);
    }
    good = AdvanceRay(1,0.0,debug);
  }
  
  if(draw) {
    TGraph *raygra = new TGraph(xy.size());
    raygra->SetLineColor(color);
    raygra->SetLineWidth(3);
    //raygra->SetTitle("Effect of Depth Dependant Refractive Index");
    for(int i=0; i<xy.size(); i++) raygra->SetPoint(i,xy[i],z[i]);
    if(gPad) {
      if(gPad->GetListOfPrimitives()->GetEntries()==0) { raygra->DrawClone("l"); printf("Drawing first trace at TrajAng = %.4f \n",TrajectoryAngle); }
      else { raygra->DrawClone("same,l"); printf("Drawing another trace at TrajAng = %.4f \n",TrajectoryAngle); }
    }
    delete raygra;
  }

  double finalZ = fGlobalZ-earthRadius;
  
  fMin = (finalZ-fDestination[2]);
  fAbsMin = fabs(finalZ-fDestination[2]);
  //RS->RecieveAngle = L_O_C_Angle(0.0,0.0,is,js,il,jl);

  if(good==99) return 99;

  //if(debug) fVerbose = NONE;
  
  if(fReflection && finalZ<fDestination[2]) return 1;

  else if(fReflection && finalZ>fDestination[2]) return 2;

  else if(finalZ>fDestination[2]) return 1;
  
  else if(finalZ<fDestination[2]) return 2;
  
  else return 0;
  
}
//
//
//
//
//
//
void PreciseRadialRayTracer::TimeVolley(double TrajectoryAngle, double stoptime, int color)  {

  if(gPad) {  
    
    vector<double> xy;
    vector<double> z;
    
    fGlobalAng = 0.0;
    fLocalAng = TrajectoryAngle;
        
    fGlobalZ = fOrigin[2] + earthRadius;
    fGlobalR = fGlobalZ;
    
    fXYDist = 0.0;
    fDistance = 0.0;
    fTime = 0.0;
    
    fReflection = false;
    fTopReflection = false;
    
    int good = 0;
    
    for(int i=0; good==0; i++) { 
      xy.push_back(fXYDist);
      z.push_back(fGlobalZ-earthRadius);
      good = AdvanceRay(2,stoptime);
    }
    
    TGraph *raygra = new TGraph(xy.size());
    raygra->SetLineColor(color);
    raygra->SetLineWidth(3);
    //raygra->SetTitle("Effect of Depth Dependant Refractive Index");
    for(int i=0; i<xy.size(); i++) raygra->SetPoint(i,xy[i],z[i]);
    if(gPad->GetListOfPrimitives()->GetEntries()==0) { raygra->DrawClone("l"); printf("Drawing first trace at TrajAng = %.4f \n",TrajectoryAngle); }
    else { raygra->DrawClone("same,l"); printf("Drawing another trace at TrajAng = %.4f \n",TrajectoryAngle);}
    delete raygra;
  
  }

}
//
//
//
//
//
//
int PreciseRadialRayTracer::AdvanceRay(int stopCond, double stoptime, bool debug) {

  double inBetweenAng  = 0;
  double angchange     = 0.0;
  bool   refref        = false;
  double nextGlobalAng = 0;
  double nextGlobalR   = 0;
  int    retInt       = false;
  //bool   debug         = false;

  double diststep = DistStep;
  double xyStep = diststep*sin(fLocalAng);
  double zStep  = diststep*cos(fLocalAng);

  if(isnan(zStep)) zStep = 0.0;
  
  double fake1, fake2;

  ConvertCartesian2Spherical(0.0,0.0,0.0,0.0,fXYDist+xyStep,fGlobalZ+zStep,nextGlobalAng,fake1,nextGlobalR);

  double il = xyStep/diststep;
  double jl = zStep/diststep;
  double ig = sin(nextGlobalAng);
  double jg = cos(nextGlobalAng);
      
  bool upAndOut;
  bool downAndIn = (fGlobalR>earthRadius) && (nextGlobalR<earthRadius);
    
  if(downAndIn) {
    FindDownCoords(xyStep,zStep,diststep,debug);
    ConvertCartesian2Spherical(0.0,0.0,0.0,0.0,fXYDist+xyStep,fGlobalZ+zStep,nextGlobalAng,fake1,nextGlobalR);
    fEntryOrExitGlobalAngle = nextGlobalAng;
    fEntryOrExitXYDist = fXYDist+xyStep;
    fEntryOrExitTime = fTime + diststep/(C/NofZ2(fGlobalR-earthRadius));
    //if(debug) printf("Down in Current Depth = %.9f \n",nextGlobalR-earthRadius);
    il = xyStep/diststep;
    jl = zStep/diststep;
    ig = sin(nextGlobalAng);
    jg = cos(nextGlobalAng);
    angchange = GetAngularChange(ig,jg,il,jl,NofZ2(fGlobalR-earthRadius),NofZ2(nextGlobalR-earthRadius),refref);
  }
  else { 
    angchange = GetAngularChange(ig,jg,il,jl,NofZ2(fGlobalR-earthRadius),NofZ2(nextGlobalR-earthRadius),refref);
    bool anggood = angchange>fMaxAngChange;
    while(anggood) {
      xyStep /= 2;
      zStep /= 2;
      diststep /= 2;
      if(diststep<1e-15) return 99;
      ConvertCartesian2Spherical(0.0,0.0,0.0,0.0,fXYDist+xyStep,fGlobalZ+zStep,nextGlobalAng,fake1,nextGlobalR);
      upAndOut = (fGlobalR<earthRadius) && (nextGlobalR>earthRadius);
      angchange = GetAngularChange(ig,jg,il,jl,NofZ2(fGlobalR-earthRadius),NofZ2(nextGlobalR-earthRadius),refref);
      //if(fDeBug) printf("AngChange = %.15f  Ref: %d  Anggood: %d \n",angchange,refref,angchange>fMaxAngChange);
      //fReflection = refref;
      //if(diststep==0) {printf("DistStep equals zero!!! \n"); assert(0);}
      if(refref && diststep<(fTolerance) ) {fReflection = true; /*FindDownCoordsFromUpRay(xyStep,zStep,diststep,debug);*/ anggood = false;}
      else if(upAndOut && diststep<(fTolerance)) {
	FindUpCoords(xyStep,zStep,diststep,debug);
	ConvertCartesian2Spherical(0.0,0.0,0.0,0.0,fXYDist+xyStep,fGlobalZ+zStep,nextGlobalAng,fake1,nextGlobalR);
	fEntryOrExitGlobalAngle = nextGlobalAng;
	fEntryOrExitXYDist = fXYDist+xyStep;
	fEntryOrExitTime = fTime + diststep/(C/NofZ2(fGlobalR-earthRadius));
	il = xyStep/diststep;
	jl = zStep/diststep;
	ig = sin(nextGlobalAng);
	jg = cos(nextGlobalAng);
	angchange = GetAngularChange(ig,jg,il,jl,NofZ2(fGlobalR-earthRadius),NofZ2(nextGlobalR-earthRadius),refref);
	anggood = false;
      }
      else anggood = angchange>fMaxAngChange;
    }
  }

  if(refref) fReflection = true;
  
  if(refref && (fGlobalR<earthRadius) && (nextGlobalR>earthRadius)) { 
    FindDownCoordsFromUpRay(xyStep,zStep,diststep,debug);
    ConvertCartesian2Spherical(0.0,0.0,0.0,0.0,fXYDist+xyStep,fGlobalZ+zStep,nextGlobalAng,fake1,nextGlobalR);
    fTopReflection = true;
  }

  if((fXYDist+xyStep)>fXYLimit && stopCond==1) {
    double mult = (fXYLimit-fXYDist)/xyStep;
    xyStep *= mult;
    zStep *= mult;
    diststep *= mult;
    ConvertCartesian2Spherical(0.0,0.0,0.0,0.0,fXYDist+xyStep,fGlobalZ+zStep,nextGlobalAng,fake1,nextGlobalR);
    angchange = GetAngularChange(ig,jg,il,jl,NofZ2(fGlobalR-earthRadius),NofZ2(nextGlobalR-earthRadius),refref);
    retInt = 1;
    //printf("Hit xy dist \n");
  }

  fGlobalAng = nextGlobalAng;
  fLocalAng += angchange;

  fGlobalZ  += zStep;    
  fGlobalR = nextGlobalR;

  fXYDist  += xyStep;
  fDistance += diststep; 
  if(downAndIn) fTime += diststep/C;
  else fTime += diststep/(C/NofZ2(fGlobalR-earthRadius));

  /*if(fVerbose>=HIGH)*/ if(debug) { printf("xyDist: %.3f zDist: %.3f totDist: %.3f fLocalAng: %.6f Depth: %.3f -> %.3f steps: (%.4f,%.8f,%.8f) Time: %.6f n: %.3f \n",fXYDist,fOrigin[2]-(fGlobalZ-earthRadius),fDistance,fLocalAng*180./TMath::Pi(),fGlobalR-earthRadius,nextGlobalR-earthRadius,diststep,xyStep,zStep,fTime,NofZ2(fGlobalR-earthRadius)); /*usleep(10000);*/ }

  if(stopCond==2) { if(fTime>stoptime) return true; }
  if(fDistance>2*fDistLimit && stopCond==1) return true;

  //if(fVerbose==LOW) cout << "*";
  
  //if(fLocalAng==TMath::Pi()/2.0) printf("90 degrees hit \n");
    
  return retInt;
        
}
//
//
//
//
//
//
void PreciseRadialRayTracer::FindDownCoords(double &xystep, double &zstep, double &diststep, bool debug) {

  double tempXY, tempZ;
  double tempang, tempRup, tempRmid, tempRdown;
  double fake1;

  double xystepup,xystepmid,xystepdown;
  double zstepup,zstepmid,zstepdown;

  double upfactor = 1.0;
  double downfactor = 0.0;
  double midfactor = (upfactor+downfactor)/2.0;

  int iter = 0;

  if(fVerbose>=HIGH) printf("Searching for entry... \n");
  
  bool good1 = (tempRmid-earthRadius)>0.0;
  bool good2 = fabs(tempRmid-earthRadius)>(1e-9);
  
  while( (good1 || good2) && iter<maxSurfaceFindIter) {

    xystepup = (xystep*upfactor);
    xystepmid = (xystep*midfactor);
    xystepdown = (xystep*downfactor);
    zstepup = (zstep*upfactor);
    zstepmid = (zstep*midfactor);
    zstepdown = (zstep*downfactor);

    ConvertCartesian2Spherical(0.0,0.0,0.0,0.0,fXYDist+xystepdown,fGlobalZ+zstepdown,tempang,fake1,tempRdown);
    ConvertCartesian2Spherical(0.0,0.0,0.0,0.0,fXYDist+xystepmid,fGlobalZ+zstepmid,tempang,fake1,tempRmid);
    ConvertCartesian2Spherical(0.0,0.0,0.0,0.0,fXYDist+xystepup,fGlobalZ+zstepup,tempang,fake1,tempRup);
    
    if((tempRmid-earthRadius)==0.0) {
      //printf("Equals Zero !!!! \n");
      downfactor = midfactor;
      midfactor = (upfactor+downfactor)/2.0;
    }
    else if( (tempRdown-earthRadius)>0.0 && (tempRmid-earthRadius)>=0.0) {
      downfactor = midfactor;
      midfactor = (upfactor+downfactor)/2.0;
    }
    else if( (tempRup-earthRadius)<0.0 && (tempRmid-earthRadius)<=0.0) {
      upfactor = midfactor;
      midfactor = (upfactor+downfactor)/2.0;
    }
    
    if(fVerbose>=HIGH)/* if(debug)*/ printf("Iteration: %i up/mid/down: (%.12f,%.12f,%.12f) Guess resultant depth = %.12f-%.12f-%.12f \n",iter,upfactor,midfactor,downfactor,(tempRdown-earthRadius),(tempRmid-earthRadius),(tempRup-earthRadius));
    good1 = (tempRmid-earthRadius)>=0.0;
    good2 = fabs(tempRmid-earthRadius)>(1e-9);
    iter++;    
  }
  
  xystep = xystepmid;
  zstep = zstepmid;
  diststep = fabs(zstep/cos(fLocalAng));

}
//
//
//
//
//
//
void PreciseRadialRayTracer::FindUpCoords(double &xystep, double &zstep, double &diststep, bool debug) {

  double tempXY, tempZ;
  double tempang, tempRup, tempRmid, tempRdown;
  double fake1;

  double xystepup,xystepmid,xystepdown;
  double zstepup,zstepmid,zstepdown;

  double upfactor = 1.0;
  double downfactor = 0.0;
  double midfactor = (upfactor+downfactor)/2.0;

  int iter = 0;

  if(fVerbose>=HIGH) printf("Searching for exit... \n");
  
  bool good1 = (tempRmid-earthRadius)<0.0;
  bool good2 = fabs(tempRmid-earthRadius)>(1e-9);
  
  while( (good1 || good2) && iter<maxSurfaceFindIter) {

    xystepup = (xystep*upfactor);
    xystepmid = (xystep*midfactor);
    xystepdown = (xystep*downfactor);
    zstepup = (zstep*upfactor);
    zstepmid = (zstep*midfactor);
    zstepdown = (zstep*downfactor);

    ConvertCartesian2Spherical(0.0,0.0,0.0,0.0,fXYDist+xystepdown,fGlobalZ+zstepdown,tempang,fake1,tempRdown);
    ConvertCartesian2Spherical(0.0,0.0,0.0,0.0,fXYDist+xystepmid,fGlobalZ+zstepmid,tempang,fake1,tempRmid);
    ConvertCartesian2Spherical(0.0,0.0,0.0,0.0,fXYDist+xystepup,fGlobalZ+zstepup,tempang,fake1,tempRup);
    
    if( (tempRdown-earthRadius)<0.0 && (tempRmid-earthRadius)<=0.0) {
      downfactor = midfactor;
      midfactor = (upfactor+downfactor)/2.0;
    }
    if( (tempRup-earthRadius)>0.0 && (tempRmid-earthRadius)>=0.0) {
      upfactor = midfactor;
      midfactor = (upfactor+downfactor)/2.0;
    }
    
    if(fVerbose>=HIGH) /*if(debug)*/ printf("Iteration: %i up/mid/down: (%.6f,%.6f,%.6f) Guess resultant depth = %.12f \n",iter,upfactor,midfactor,downfactor,(tempRmid-earthRadius));
    good1 = (tempRmid-earthRadius)<=0.0;
    good2 = fabs(tempRmid-earthRadius)>(1e-9);
    iter++;    
  }
  
  xystep = xystepmid;
  zstep = zstepmid;
  diststep = fabs(zstep/cos(fLocalAng));

}
//
//
//
//
//
//
void PreciseRadialRayTracer::FindDownCoordsFromUpRay(double &xystep, double &zstep, double &diststep, bool debug) {

  double tempXY, tempZ;
  double tempang, tempRup, tempRmid, tempRdown;
  double fake1;

  double xystepup,xystepmid,xystepdown;
  double zstepup,zstepmid,zstepdown;

  double upfactor = 1.0;
  double downfactor = 0.0;
  double midfactor = (upfactor+downfactor)/2.0;

  int iter = 0;

  if(fVerbose>=HIGH) printf("Searching for entry... \n");
  
  bool good1 = (tempRmid-earthRadius)>0.0;
  bool good2 = fabs(tempRmid-earthRadius)>(1e-9);
  
  while( (good1 || good2) && iter<maxSurfaceFindIter) {

    xystepup = (xystep*upfactor);
    xystepmid = (xystep*midfactor);
    xystepdown = (xystep*downfactor);
    zstepup = (zstep*upfactor);
    zstepmid = (zstep*midfactor);
    zstepdown = (zstep*downfactor);

    ConvertCartesian2Spherical(0.0,0.0,0.0,0.0,fXYDist+xystepdown,fGlobalZ+zstepdown,tempang,fake1,tempRdown);
    ConvertCartesian2Spherical(0.0,0.0,0.0,0.0,fXYDist+xystepmid,fGlobalZ+zstepmid,tempang,fake1,tempRmid);
    ConvertCartesian2Spherical(0.0,0.0,0.0,0.0,fXYDist+xystepup,fGlobalZ+zstepup,tempang,fake1,tempRup);
    
    if( (tempRdown-earthRadius)<0.0 && (tempRmid-earthRadius)>=0.0) {
      upfactor = midfactor;
      midfactor = (upfactor+downfactor)/2.0;
    }
    else if( (tempRup-earthRadius)>0.0 && (tempRmid-earthRadius)<=0.0) {
      downfactor = midfactor;
      midfactor = (upfactor+downfactor)/2.0;
    }
    
    if(fVerbose>=HIGH) /*if(debug)*/ printf("Iteration: %i up/mid/down: (%.6f,%.6f,%.6f) Guess resultant depth = %.12f \n",iter,upfactor,midfactor,downfactor,(tempRmid-earthRadius));
    good1 = (tempRmid-earthRadius)>=0.0;
    good2 = fabs(tempRmid-earthRadius)>(1e-9);
    iter++;    
  }
  
  xystep = xystepmid;
  zstep = zstepmid;
  diststep = fabs(zstep/cos(fLocalAng));

}
//
//
//
//
//
//
double PreciseRadialRayTracer::GetAngularChange(double xcompearth, double ycompearth, double xcompray, double ycompray, double IncidentIndex, double RefractionIndex, bool &reflection) {

  reflection = false;

  double CriticalAngle = asin(RefractionIndex/IncidentIndex);
  double inBetweenAng = L_O_C_Angle(0.0,0.0,xcompearth,ycompearth,xcompray,ycompray);
  fTempRecieveAngle = inBetweenAng;
  //if(inBetweenAng==TMath::Pi()/2.0) printf("90 degrees hit \n");

  double theta = inBetweenAng;
  if(inBetweenAng>(M_PI/2.0)) theta = M_PI-inBetweenAng;
  double opptheta = (M_PI/2.0)-theta;
  
  if(theta>CriticalAngle && !isnan(CriticalAngle)) { reflection = true; /*printf("In angle = %.9f  Out angle = %.9f \n",inBetweenAng*r2d,2*opptheta*r2d) ;*/return (2.0*opptheta);}
  double nextAngle = asin( (IncidentIndex/RefractionIndex)*sin(theta));
  if(theta<M_PI/2.0 && (theta+fabs(nextAngle-theta))>M_PI/2.0) { reflection = true; }
  return fabs(nextAngle-theta);

}
