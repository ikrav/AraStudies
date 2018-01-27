#include "../RayTrace/PreciseRadialRayTracer.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TMarker.h"
#include <iostream>

void solutionPlot();

const int VMIN = -3000;
const int VMAX = -49;
const int VINC = 50;
const int HINC = 100;

const int RED = 2;
const int BLUE = 4;

void solutionPlot() {
	
	PreciseRadialRayTracer *tracer = new PreciseRadialRayTracer();
	tracer->SetDestination(0, 0, -170);
			
	ofstream file;
	file.open("solutionPlot_out");
	
	TCanvas *c = new TCanvas("zones", "# of Solutions", 0, 0, 5000, 1750);
	TH2F *graph = new TH2F("points", "# of Solutions", 100, 0, 20000, 100, -3000, 0);
	c->Draw();
	graph->Draw();
	
	// iterate over horizontal slices
	for(int z = VMIN; z <= VMAX; z += VINC) {
		
		bool success = true;
		int r = 1, numSolutions;
		
		// keep plotting points until reaching edge of SZ
		while(success) {
		
			tracer->SetOrigin(r, 0, z);
			success = tracer->TraceRay(numSolutions, false);
			// uncomment below line for file output
			//file << numSolutions << " " << r << " " << z << std::endl;
			
			TMarker *m = new TMarker(r, z, 8); // 8 is a scalable circle
			m->SetMarkerSize(1);
			
			if(numSolutions == 1) {
				m->SetMarkerColor(BLUE); 
			} else if(numSolutions == 2) {
				m->SetMarkerColor(RED);
			}
			
			m->Draw("Same");
			r += HINC;
		}
	}
	
	c->Update();
	c->Show();
}
