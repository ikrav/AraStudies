#ifndef ASlibEdit_h
#define ASlibEdit_h

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include "time.h"

#include <TVector3.h>
#include <TBenchmark.h>
#include <TH3.h>
#include <TCanvas.h>
#include <TString.h>
#include <TRandom3.h>
#include <TH1.h>
#include <TTree.h>
#include <TFile.h>
#include <TMinuit.h>
#include <TGraph.h>
#include "TMath.h"
#include <TSystemDirectory.h>
#include <TList.h>
#include <TSystem.h>
using namespace std;

// Version 1.6
// Date: 05/28/15
/*
    *******************************************************************
    *******************************************************************
    ************************** CONTENTS *******************************
    *******************************************************************
    *******************************************************************
	     

1) inline double SQ(double a);
Info: Square equation x^2.

2) void Orderer(double *values, int *order, int n) 
Info: Outputs the "order" of an array from least to greatest

3) void StatFind( double *values, Int_t numval, double &mean, double &stdev);
Info: input array of *values and number (numval) of values in the array;
output the mean and standard deviation (stdev) of the *values.

4) float ntimexp(float z1, float z2, float dtot)
Info: Calculates time of wave propagation through SP ice

5) double InvTimeExp(double time, double z1, double z2)
Info: Inverse of ntimexp, has different ice models within

6) inline double Dist3d(double x1, double y1, double z1,double x2, double y2, double z2)
Info: Distance between two 3d points

7) inline double Dist2d(double x1, double y1,double x2, double y2)
Info: Distance between two 2d points

8) double fact(Int_t n)
Info: Factorial Function

9) inline double nCr(Int_t n, Int_t r)
Info: n choose r, combination finder

10) double DP(double *a, double *b, Int_t NDim)
Info: Vector Scalar Product

11) void QUAD(double A, double B, double C, double &Pos, double &Neg)
Info: Lightweight quadratic equation solver

12) double NofZ(double z)
Info: Index of refraction at south pole as a function of depth (must be negative)

13) bool TimeCheck(TStopwatch *SW, double UpdateTimeInterval)
Info: Useful for completion percent printing 

14) void ConvertCartesian2Spherical(double origx, double origy, double origz,
double vecx, double vecy, double vecz,
double &Theta, double &Phi, double &R)
Info: Convert vector orig(x,y,z) to vec(x,y,z) into Theta, Phi, and R

15) void ConvertSpherical2Cartesian(double origx, double origy, double origz,
double Theta, double Phi, double R,
double &vecx, double &vecy, double &vecz);
Info: Convert spherical theta, phi, and r into cartesian vec(x,y,z) from center orig(x,y,z)

16) double Rounder(double value, double increment)
Info: Prototype function, do not use, further rewriting and testing needed

17) double EffectiveIoR(double z1, double z2)
Info: Calculates effective index of refraction between to points

18) int digitCheck(int number)
Info: Returns number of digits in a int type number

19) bool nextFile(TString dirname, TString findfile, TString &ext)
Info: Search the directories in a directory for a file, can only handle
one directory above the starting directory for now.

20) double L_O_C_Angle(double x1, double y1, double x1, double y1)
Info: Law Of Cosines, returns angle between point 1 and point 2

21) bool checkFileExistance(const std::string& name)

22) double getIoRFromParameters(double A, double B, double C, double depth);

    *******************************************************************
    *******************************************************************
    ************************ END CONTENTS *****************************
    *******************************************************************
    ******************************************************************/


  
inline double SQ(double a){return(a*a);}

void Orderer(double *values, int *order, int n);
void Orderer(vector<double> &values, vector<int> &order, int n);
void Orderer(vector<int> &values, vector<int> &order);
void Orderer(vector<double> &values, vector<int> &order);

void StatFind(double *values, Int_t numval, double &mean, double &stdev);

float ntimexp(float z1, float z2, float dtot);

double InvTimeExp(double time, double z1, double z2);

inline double Dist3d(double x1, double y1, double z1,double x2, double y2, double z2)
{return sqrt(SQ(x2-x1) + SQ(y2-y1) + SQ(z2-z1));}

inline double Dist2d(double x1, double y1,double x2, double y2)
{return sqrt(SQ(x2-x1) + SQ(y2-y1));}

double fact(double n);

inline double nCr(double n, double r) {return ( fact(n)/( fact(r)*fact(n-r)));}

double DP(double *a, double *b, Int_t NDim);

void QUAD(double A, double B, double C, double &Pos, double &Neg);
void QUAD(long double A, long double B, long double C, long double &Pos, long double &Neg);

double NofZ(double z);

bool TimeCheck(TStopwatch *SW, double UpdateTimeInterval, double TVal, double BVal);

void ConvertCartesian2Spherical(double origx, double origy, double origz,
				  double vecx, double vecy, double vecz,
				  double &Theta, double &Phi, double &R);

void ConvertSpherical2Cartesian(double origx, double origy, double origz,
				  double Theta, double Phi, double R,
				  double &vecx, double &vecy, double &vecz);

void ConvertCartesian2Spherical(long double origx, long double origy, long double origz,
				  long double vecx, long double vecy, long double vecz,
				  long double &Theta, long double &Phi, long double &R);

void ConvertSpherical2Cartesian(long double origx, long double origy, long double origz,
				  long double Theta, long double Phi, long double R,
				  long double &vecx, long double &vecy, long double &vecz);

double Rounder(double value, double increment);

double EffectiveIoR(double z1, double z2);
  
int digitCheck(int number);

bool nextFile(TString dirname, TString findfile, TString &ext);

double L_O_C_Angle(double centerX, double centerY, double x1, double y1, double x2, double y2);

bool checkFileExistance(const std::string& name);

double getIoRFromParameters(double A, double B, double C, double depth);

void SaveToPdf(vector<TCanvas*> can, TString pdfname);
void SaveToPdf(vector < vector<TCanvas*> > can, TString pdfname);

void ConvertUTCunixtime2DateTimeUTC(double unixtime, int &sec, int &min, int &hour, int &day, int &month, int &year);

#endif

