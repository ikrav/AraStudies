#include "../RayTrace/PreciseRadialRayTracer.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TMarker.h"
#include "TMath.h"
#include "TGraph2D.h"
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int volume(TGraph *graphs[], int numGraphs);
int FindInc(TGraph *graph, double &inc);
int Configure(TGraph *graphs[], int numGraphs, long &xBound, long &yBound);
bool isHit(TGraph *graphs[], int numGraphs, int layer, int xMark, int yMark);
int FindBounds(TGraph *graph, double &floor, double &ceil);
double distance(double x1, double y1, double x2, double y2);

const int SAMPLE_SIZE = 100000;

int volume(TGraph *graphs[], int numGraphs) {
	
	if(numGraphs <= 0) {
		return 1;
	}
	
	// Monte Carlo requires that all graphs have the same cardinality
	// Use 0th index as benchmark
	cardinality = graphs[0]->GetN();
	for(int i = 1; i < numGraphs; i++) {
		if(graphs[i]->GetN() != cardinality) {
			return 1;
		}
	}

	//TODO: check layer alignment in the z-axis
		
	// These values are vulnerable to user error
	// until the above validation is implemented
	double inc, z, r;
	long xBound, yBound;
	FindInc(graphs[0], inc);
	Configure(graphs, numGraphs, xBound, yBound);
	
	//cout << xBound << "×" << yBound << endl;
	
	double volume = 0;
	int xMark, yMark;
	int hits = 0;
	int misses = 0;
	
	for(int layer = 1; layer < graphs[0]->GetN(); layer++) {
	
		//cout << layer << endl;
		
		srand(time(NULL));
		
		while(hits + misses <= SAMPLE_SIZE) {
			
			xMark = std::rand() % xBound;
			yMark = std::rand() % yBound;
			
			if(isHit(graphs, numGraphs, layer, xMark, yMark)) {
				hits++;
			} else {
				misses++;
			}
		}
		
		//cout << hits << "/" << SAMPLE_SIZE << endl;
		volume += (((double)hits / SAMPLE_SIZE) * (xBound * yBound)) * inc;
		//cout << (double)hits / SAMPLE_SIZE << endl;
		//cout << xBound * yBound << endl;
		//cout << endl;
		hits = 0;
		misses = 0;
	}
	
	std::cout << volume << std::endl;
	return 0;
}

bool isHit(TGraph *graphs[], int numGraphs, int layer, int xMark, int yMark) {
	
	double x, y, r, z;
	for(int i = 0; i < numGraphs; i++) {
		
		graphs[i]->GetPoint(0, x, y);
		graphs[i]->GetPoint(layer, r, z);
		//TODO:RM: cout << x << ":" << y << ":" << r << endl;

		//cout << distance(xMark, yMark, x, y) << " — " << r << endl;
		if(distance(xMark, yMark, x, y) <= r) {			
			return true;
		}
	}
	
	return false;
}

int FindBounds(TGraph *graph, double &floor, double &ceil) {
	
	double r, z;
	
	graph->GetPoint(1, r, z);
	floor = z;
	ceil = z;
	
	for(int i = 2; i < graph->GetN(); i++) {
	
		graph->GetPoint(i, r, z);
		
		if(z < floor) {
			floor = z;
		} else if(z > ceil) {
			ceil = z;
		}
	}

	return 0;
}

int FindInc(TGraph *graph, double &inc) {

	double r, z1, z2;
	graph->GetPoint(1, r, z1);
	graph->GetPoint(2, r, z2);	
	
	inc = z1 - z2;
	return 0;	
}

//7.86266e+11

int Configure(TGraph *graphs[], int numGraphs, long &xBound, long &yBound) {

	double x, y, r, z;
	
	// xy-coords of every graph should be translated
	// to keep everything in quadrant I
	graphs[0]->GetPoint(0, x, y);
	double xOffset = 0 - x;
	double yOffset = 0 - y;
	
	for(int i = 1; i < numGraphs; i++) {
		graphs[i]->GetPoint(0, x, y);
		xOffset = x > xOffset ? x : xOffset;
		yOffset = y > yOffset ? y : yOffset;
	}
	
	xOffset = xOffset > 0 ? xOffset : 0;
	yOffset = yOffset > 0 ? yOffset : 0;
	
	// find the max radius to add to the offsets
	graphs[0]->GetPoint(graphs[0]->GetN() - 1, r, z);
	double rMax;
	
	// assume that the last point in the graph has the greatest radius
	for(int i = 1; i < numGraphs; i++) {	
		graphs[i]->GetPoint(graphs[i]->GetN() - 1, r, z);
		rMax = r > rMax ? r : rMax;
	}
	
	// original offset + max radius = offset needed to place entire area in Q1
	xOffset += rMax;
	yOffset += rMax;
	
	for(int i = 0; i < numGraphs; i++) {
		graphs[i]->GetPoint(0, x, y);
		graphs[i]->SetPoint(0, x + xOffset, y + yOffset);
	}
	
	// Determine xy-boundary with offset applied
	// Used to bound rand() output
	xBound = 0;
	yBound = 0;
	
	for(int i = 0; i < numGraphs; i++) {
		graphs[i]->GetPoint(0, x, y);
		xBound = (x + rMax) > xBound ? ceil(x + rMax) : xBound;
		yBound = (y + rMax) > yBound ? ceil(y + rMax) : yBound;
	}
	
	return 0;
}

double distance(double x1, double y1, double x2, double y2) {
	return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}















