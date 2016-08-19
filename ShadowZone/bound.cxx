#include "../RayTrace/PreciseRadialRayTracer.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TMarker.h"
#include "TMath.h"
#include <iostream>
#include <vector>

TGraph * bound(double xGlobal, double yGlobal);
double thetaNext(double z);

// iA, iB, and iC are index of refraction constants 
const double iA = 1.77, iB = 0.4272, iC = 0.016;
const double FLOOR = -2400;
const double ANTZ = -170;
const int INC = 10;

TGraph * bound(double xGlobal, double yGlobal) {
	
	PreciseRadialRayTracer *tracer = new PreciseRadialRayTracer();
	tracer->SetDestination(0, 0, -170);
	
	TCanvas *c = new TCanvas("zones", "# of Solutions", 0, 0, 5000, 1750);
	c->Draw();
	
	double z = ANTZ - 10;
	double r = 0;
	double theta = thetaNext(z);
	
	// find number of indices needed to hold results for the graph
	int n = ((FLOOR - ANTZ) / INC) * -1;
	double zPoints[n + 1];
	double rPoints[n + 1];
	
	zPoints[0] = xGlobal;
	rPoints[0] = yGlobal;
	
	for(int i = 1; i < n; i++) {
		
		zPoints[i] = z;
		rPoints[i] = r;
	
		// update var values
		z -= INC;
		r += INC * tan(theta);
		theta = thetaNext(z);
	}
	
	TGraph *g = new TGraph(n, rPoints, zPoints);
	//g->SetMarkerStyle(3);
	g->Draw();
	
	c->Update();
	c->Show();
	
	return g;
}

// calculate theta using n() = A - B * e^(Cz)
double thetaNext(double z) {
	double n1 = iA - iB * exp(iC*ANTZ);
	double n2 = iA - iB * exp(iC*z);
	return asin(n1 / n2);
}
