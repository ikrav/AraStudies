#include "../RayTrace/PreciseRadialRayTracer.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TMarker.h"
#include "TMath.h"
#include <iostream>

int volume(TGraph *g);

int volume(TGraph *g) {
	
	int n = g->GetN();
	double r, z;
	
	if(n <= 1) {
		return 1;
	}
	
	g->GetPoint(0, r, z);
	double z1 = z;
	g->GetPoint(1, r, z);
	double z2 = z;
	
	double inc = z1 - z2;
		
	if(inc <= 0) {
		return 1;
	}
	
	double volume = 0;
	for(int i = 0; i < n; i++) {
	
		g->GetPoint(i, r, z);
		volume += (TMath::Pi() * (r * r)) * inc;
	}
	
	std::cout << volume << std::endl;
	return 0;
}
