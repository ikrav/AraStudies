#include "ASlibEdit.h"

//using namespace ASLIBEDIT_NSPACE;

void Orderer(double *values, int *order, int n) {

  double valord[n];
  double high=0;

  for (int i=0; i<n; i++) {

    valord[i] = values[i];
    high += fabs(values[i]);

    }

  for(int i = 0;i<n;i++) {
    double checker = high;

    for (int j = 0;j<n;j++) {

      if (valord[j]<=checker) {

	checker = values[j];
	order[i] = j; 

      }
    }

    valord[order[i]] = high+1;

  }

}

void Orderer(vector<double> &values, vector<int> &order, int n) {
  double valord[n];
  double high=0;
  order.resize(n);
  for (int i=0; i<n; i++) {
    valord[i] = values[i];
    high += fabs(values[i]);
  }
  for(int i = 0;i<n;i++) {
    double checker = high;
    for (int j = 0;j<n;j++) {
      if (valord[j]<=checker) {
	checker = values[j];
	order[i] = j; 
      }
    }
    valord[order[i]] = high+1;
  }
  //for(int i=0; i<order.size(); i++) printf("Order[%i] = %i \n",i,order[i]);
}

void Orderer(vector<int> &values, vector<int> &order) {
  int n = values.size();
  int valord[n];
  int high=0;
  vector<int> order2;
  order.resize(n);
  order2.resize(n);
  for (int i=0; i<n; i++) {
    valord[i] = values[i];
    high += fabs(values[i]);
  }
  for(int i = 0;i<n;i++) {
    int checker = high;
    for (int j = 0;j<n;j++) {
      if (valord[j]<=checker) {
	checker = values[j];
	order[i] = checker; 
	order2[i] = j;
      }
    }
    valord[order2[i]] = high+1;
  }
  //for(int i=0; i<order.size(); i++) printf("Order[%i] = %i \n",i,order[i]);
}

void Orderer(vector<double> &values, vector<int> &order) {
  int n = values.size();
  double valord[n];
  double high=0;
  vector<double> order2;
  order.resize(n);
  order2.resize(n);
  for (int i=0; i<n; i++) {
    valord[i] = values[i];
    high += fabs(values[i]);
  }
  for(int i = 0;i<n;i++) {
    double checker = high;
    for (int j = 0;j<n;j++) {
      if (valord[j]<=checker) {
	checker = values[j];
	order2[i] = checker; 
	order[i] = j;
      }
    }
    valord[order[i]] = high+1;
  }
  //for(int i=0; i<order.size(); i++) printf("Order[%i] = %i \n",i,order[i]);
}

void StatFind(double *values, Int_t numval, double &mean, double &stdev) {

  mean = 0;
  stdev = 0;

  for (int i=0; i<numval; i++) mean += values[i];

  mean /= numval;

  for (int i=0; i<numval; i++) stdev += SQ(values[i] - mean);

  stdev = sqrt(stdev/numval);

}

float ntimexp(float z1, float z2, float dtot){

  const float km2m=1e3;
  const float cVac=2.998e-1;

  z1 /= 1000;
  z2 /= 1000;
  dtot /= 1000;

  float A=1.78;
  float B=A*-0.24;
  const float C=16.; //was 16. but must be shifted since units are km
  //float B = 0.498;
  //const float C = 0.0225;
  float zmax=z2;
  float zmin=z1;
  if(zmin>zmax){
    zmin=z2;
    zmax=z1;
    }
  float ctheta;
  if(dtot>0)
    ctheta=fabs(z2-z1)/dtot;
  else
    return 0;
  float ftimexp=dtot/cVac; //was =0
  if(ctheta!=0){
    if(zmax<0. && zmin<0.)
      ftimexp=(A*(zmax-zmin)+(B/C)*(exp(C*zmax)-exp(C*zmin)))/(cVac * ctheta);
    if(zmax>=0. && zmin<0.)
      ftimexp=(A*(0.-zmin)+(B/C)*(exp(C*0.)-exp(C*zmin))+zmax)/(cVac * ctheta);
    if(zmax>=0. && zmin>=0.)
      ftimexp=(zmax-zmin)/(cVac * ctheta);
    if(zmax<0. && zmin>0.)
      cout<<" CONDITION NOT COVERED IN NTIMEZ!!!"<<endl;
  }
  return ftimexp*km2m; //km->meters
}

double InvTimeExp(double time, double z1, double z2) {

  const float cVac=2.998e-1;

  double zmax = z1;
  double zmin = z2;

  if (zmin>zmax) {

    zmin = z1;
    zmax = z2;

  }
  
  double A = 1.78;
  double B = A * (-.24);
  double C = .016;

  double dist = (time * cVac * fabs(zmax-zmin) ) /
    ( A * (zmax-zmin) + (B/C)* ( exp(C*zmax) - exp(C*zmin)) );

  return dist;

}

double fact(double n) {

  double val = 1;
  for (int i=1; i<=n; i++)val*=i;
  return val;
}
  
double DP(double *a, double *b, Int_t NDim){

  double sp = 0;

  for(int i=0; i<NDim; i++) {

    sp += a[i]*b[i];
  }

  return sp;
}

void QUAD(double A, double B, double C, double &Pos, double &Neg) {
  
  double determinant = B*B - 4*A*C;
  Pos = (-B + sqrt(determinant))/(2*A);
  Neg = (-B - sqrt(determinant))/(2*A);
}

void QUAD(long double A, long double B, long double C, long double &Pos, long double &Neg) {
  
  long double determinant = B*B - 4*A*C;
  Pos = (-B + sqrt(determinant))/(2*A);
  Neg = (-B - sqrt(determinant))/(2*A);
}


double NofZ(double z) {

  if(z>=0) return 1.0;

  z /= 1000;
  
  const double A = 1.78;
  const double B = -.24 * A;
  const double C = 16;
  
  double n = A + B * exp(C*z);
  
  return n;
}

bool TimeCheck(TStopwatch *SW, double UpdateTimeInterval, double TVal, double BVal) {

  double time = SW->RealTime();

  if(time<UpdateTimeInterval) {

    SW->Start(false);
    return false;
  }

  else {

    SW->Start(true);
    printf("Percent Complete = %f \n",(TVal/BVal)*100.0);
    return true;
  }

}

void ConvertCartesian2Spherical(double origx, double origy, double origz,
				double vecx, double vecy, double vecz,
				double &Theta, double &Phi, double &R) {

  R = Dist3d(origx,origy,origz,vecx,vecy,vecz);

  if(R==0) {
    Phi   = nan("");
    Theta = nan("");
    return;
      }
  
  else {
    Theta = acos((vecz-origz)/R);

    double top = (vecy - origy);
    double bottom = (vecx - origx);
    
    Phi = atan2(top,bottom);
    
    if(Phi<0)
      Phi += 2*TMath::Pi();
    
  }

}

void ConvertSpherical2Cartesian(double origx, double origy, double origz,
				  double Theta, double Phi, double R,
				double &vecx, double &vecy, double &vecz) {
  
  if(R==0) {
    vecx = origx;
    vecy = origy;
    vecz = origz;  
    return;
  }
  
  else {
    vecx = R*sin(Theta)*cos(Phi)+origx;
    vecy = R*sin(Theta)*sin(Phi)+origy;
    vecz = R*cos(Theta)+origz;
  }
  
}

void ConvertCartesian2Spherical(long double origx, long double origy, long double origz,
				long double vecx, long double vecy, long double vecz,
				long double &Theta, long double &Phi, long double &R) {

  R = Dist3d(origx,origy,origz,vecx,vecy,vecz);

  if(R==0) {
    Phi   = nan("");
    Theta = nan("");
    return;
      }
  
  else {
    Theta = acos((vecz-origz)/R);

    long double top = (vecy - origy);
    long double bottom = (vecx - origx);
    
    Phi = atan2(top,bottom);
    
    if(Phi<0)
      Phi += 2*TMath::Pi();
    
  }

}

void ConvertSpherical2Cartesian(long double origx, long double origy, long double origz,
				  long double Theta, long double Phi, long double R,
				long double &vecx, long double &vecy, long double &vecz) {
  
  if(R==0) {
    vecx = origx;
    vecy = origy;
    vecz = origz;  
    return;
  }
  
  else {
    vecx = R*sin(Theta)*cos(Phi)+origx;
    vecy = R*sin(Theta)*sin(Phi)+origy;
    vecz = R*cos(Theta)+origz;
  }
  
}


double Rounder(double value, double increment) {

  double a = value/increment;
  double b = (int)a;
  int c = (int)a;

  b+=.5;

  if(b<=a) { c++; return (c*increment);}

  else return (c*increment);

}

double EffectiveIoR(double z1, double z2) {

  double A = 1.78;
  double B = 0.4272;
  double C = .016;

  double NFrac = ( (A * fabs(z2-z1)) - (B/C)*fabs( exp(C*z2) - exp(C*z1) ) ) / fabs(z2-z1);

  return NFrac;

}

int digitCheck(int number) {
  int iter = 1;
  int check = 0;
  while(iter<100) {
    for(int i=0; i<iter; i++) check += 9*pow(10,i);
    if( (number-check)<0 ) return iter;
    else {check = 0;iter++;}
  }
  return 0;
}

bool nextFile(TString dirname, TString findfile, TString &ext) {

  TSystemDirectory dir(dirname, dirname);
  TList *files = dir.GetListOfFiles();

  bool good = false;

  if (files && !good) {

    TSystemFile *file;
    TString fname;

    TIter next(files);

    while ((file=(TSystemFile*)next()) && !good) {
      fname = file->GetName();
      if (file->IsDirectory() && !(fname=="." || fname==".." || fname.BeginsWith("."))) {
        TString dirstr = "";
        dirstr += dirname + "/" + fname;

        good = nextFile(dirstr,findfile,ext);
        if(good) return true;
      }

      else if(fname==findfile) {ext = dirname+"/"+fname ;/* cout << dirname << "/" << fname << endl ;*/ good = true; return true;}

    }
  }
  delete files;

  return false;

}

double L_O_C_Angle(double centerX, double centerY, double x1, double y1, double x2, double y2) {

  double distA = Dist2d(centerX,centerY,x1,y1);
  double distB = Dist2d(centerX,centerY,x2,y2);
  double distC = Dist2d(x1,y1,x2,y2);

  double returnVal = acos(( SQ(distA)+SQ(distB)-SQ(distC) ) / ( 2*distA*distB));

  return returnVal;

}

bool checkFileExistance(const std::string& name) {
  ifstream f(name.c_str());
  if (f.good()) {
    f.close();
    return true;
  } else {
    f.close();
    return false;
  }   
}

double getIoRFromParameters(double A, double B, double C, double depth) {
  if(depth>=0) return 1.0;
  double n = A - B * exp(C*depth);
  return n;
}

void SaveToPdf(vector<TCanvas*> can, TString pdfname) {
  TString beg = pdfname;
  beg += "[";
  TString end = pdfname;
  end += "]";
  can[0]->Print(beg);
  for(int i=0; i<can.size(); i++) can[i]->Print(pdfname);
  can[can.size()-1]->Print(end);
}

void SaveToPdf(vector < vector<TCanvas*> > can, TString pdfname) {
  TString beg = pdfname;
  beg += "[";
  TString end = pdfname;
  end += "]";
  can[0][0]->Print(beg);
  for(int i=0; i<can.size(); i++) {
    for(int j=0; j<can[i].size(); j++) {
      can[i][j]->Print(pdfname);
    }
  }
  can[can.size()-1][can[can.size()-1].size()-1]->Print(end);
}

void ConvertUTCunixtime2DateTimeUTC(double unixtime, int &sec, int &min, int &hour, int &day, int &month, int &year) {
  time_t  time = unixtime;
  struct tm *timeStruct;
  timeStruct = gmtime(&time);
  sec = timeStruct->tm_sec;
  min = timeStruct->tm_min;
  hour = timeStruct->tm_hour;
  day = timeStruct->tm_mday;
  month = timeStruct->tm_mon;
  year  = timeStruct->tm_year;
}
