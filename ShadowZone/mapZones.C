#include "../RayTrace/PreciseRadialRayTracer.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TF1.h"

bool findRayTrace(PreciseRadialRayTracer *tracer, 
		  double R1, double Z1, double R2, double Z2);
void fillMapOfZones(TCanvas *canvas, TH2D *hist, int nPrimitives0);
TGraphErrors *makeBoundaryGraph(TH2D *hist);

void mapZones(){

  //gSystem->Load("PreciseRadialRayTracer_cxx.so");

  PreciseRadialRayTracer *tracer = new PreciseRadialRayTracer();

  TCanvas *canRays = new TCanvas("canRays","",0,0,600,600);
  canRays->Draw();

  const double Zmin = -1500; 
  const double Zmax = -10;
  const double dZ = 100;
  const double Rmin = 10;
  const double Rmax = 5000;
  const double dR = 500;
  
  TH2D *hist = new TH2D("hist","",100,0,Rmax,100, Zmin, 20);
  hist->SetTitle("Ray Trace; XY Distance [m]; Depth [m]");
  hist->SetStats(0);
  hist->Draw();
//   hist->DrawClone();
//   delete hist;
  canRays->Update();
  // Save the number of primitives on the canvas before any
  // rays are drawn
  int nPrimitives0 = canRays->GetListOfPrimitives()->GetEntries();

  TH2D *hZones = new TH2D("hZones","",1000, 0, Rmax, 1000, Zmin, 20);
  hZones->GetXaxis()->SetTitle("XY distance [m]");
  hZones->GetYaxis()->SetTitle("Z [m]");

  double originZ = -170.0;
  int nZ = (Zmax-Zmin)/dZ + 1;
  int nR = (Rmax-Rmin)/dR;
  for(int iZ = 0; iZ < nZ; iZ++){
    double destinationZ = Zmin + dZ * iZ;

    // Binary search for R follows    
    double Rlow = 0;
    double Rhigh = Rmax;
    const double Rtolerance = 1; //1 meter
    const int nIterationsMax = 100;
    int nIterations = 0;
    bool result = false;
    double R;
    while( (Rhigh - Rlow) > Rtolerance && nIterations < nIterationsMax){

      nIterations++;
      R = (Rlow + Rhigh)/2.0;
      result = findRayTrace(tracer, R, originZ, 0.0, destinationZ);
      canRays->Update();
      if( result ){
	Rlow = R;
      } else {
	Rhigh = R;
      }
      
    }// end while over R
    printf("For Z= %f   radius is %.0f   iterations %d\n", destinationZ, R, nIterations);
  }
  
  // Fill the map of zones and draw it
  fillMapOfZones(canRays, hZones, nPrimitives0);
  
  TCanvas *canMap = new TCanvas("canMap","canMap",100,10,600,600);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  hZones->GetZaxis()->SetRangeUser(0,100);
  hZones->Draw("col");

  // Prepare a graph of the boundary, fit it, and draw it
  TGraphErrors *grBoundary = makeBoundaryGraph(hZones);
  TCanvas *canBoundary = new TCanvas("canBoundary","canBoundary",200,10,600,600);

  TH2D *dummy = (TH2D*)hZones->Clone("dummy");
  dummy->Reset();
  dummy->Draw();
  grBoundary->Draw("PE,same");

  // In fitting, avoid area close to surface
  float lowX = dummy->GetXaxis()->GetBinCenter(1);
  float highX = dummy->GetXaxis()->GetBinCenter(dummy->GetNbinsX());
  float lowFitX = lowX + (highX-lowX)/3.0;
  float highFitX = highX;
  TF1 *func = new TF1("func","pol1",lowFitX, highFitX);
  func->SetLineColor(kRed);

  grBoundary->Fit("func","R");
  float thetaRad = atan( fabs( func->GetParameter(1)));
  float thetaDeg = thetaRad * 180 / TMath::Pi();
  printf("The angle of the shadow zone boundary is %4.1f degrees\n", thetaDeg);

  TString str = TString::Format("#theta = %4.1f degrees", thetaDeg);
  TLatex *lat = new TLatex(0.5, 0.7, str.Data());
  lat->SetNDC(kTRUE);
  lat->Draw("same");

}


bool findRayTrace(PreciseRadialRayTracer *tracer, 
		  double R1, double Z1, double R2, double Z2){

  tracer->SetOrigin     (0.0, R1, Z1); //set location of event
  tracer->SetDestination(0.0, R2, Z2); // set location of signal observation
  tracer->SetVerbosity(0); // medium=2
  
  // true: draw rays onto TCanvas, false: do not, if no argument supplied defaults to false
  int nSolutions = 0;
  bool success = tracer->TraceRay(nSolutions, true); 
  
//   if(success) {
//     printf("For Z= %f R= %f solution ok\n", Z2, R1);
//     //can->SaveAs("PreciseRadialRayTracer.pdf");
//   }else printf("--For Z= %f R= %f solution not possible ... \n", Z2, R1);
  
  return success;

}

void fillMapOfZones(TCanvas *canvas, TH2D *hist, int nPrimitives0){
  
  int nPrimitives = canvas->GetListOfPrimitives()->GetEntries();
  printf("Total %d rays, transfer them to a map\n", nPrimitives-nPrimitives0);
  for( int i = nPrimitives0; i < nPrimitives; i++){

    // Get the next ray graph
//     printf("  process ray %d\n", i-nPrimitives0);
    TObject *obj = (TObject*)canvas->GetListOfPrimitives()->At(i);
    if( obj->GetTitle() != TString("Graph") ) continue;
    TGraph *gr = (TGraph*)obj;

    // Loop over points on the graph and transfer them to points on the map
    for( int ipoint = 0; ipoint < gr->GetN(); ipoint++){
      double x,y;
      gr->GetPoint(ipoint, x,y);
      int binx = hist->GetXaxis()->FindBin(x);
      int biny = hist->GetYaxis()->FindBin(y);
      hist->SetBinContent(binx, biny, 1);
    }

  } 

  // Apply a uniform fill to the map
  for(int ix = 1; ix <= hist->GetNbinsX(); ix++){
    bool blank = true;
    for(int iy = hist->GetNbinsY(); iy>0; iy--){

      // Go from the surface down, and once the first
      // ray is encountered, fill everything below it
      // in the vertical column
      if( blank == true && hist->GetBinContent(ix,iy) > 0 )
	blank = false;
      
      if( !blank )
	hist->SetBinContent(ix,iy,1);
    } // end loop over y
  } // end loop over x

}

TGraphErrors *makeBoundaryGraph(TH2D *hist){

  TGraphErrors *gr = new TGraphErrors();
  for( int ix = 1; ix <= hist->GetNbinsX(); ix++){

    float x = hist->GetXaxis()->GetBinCenter(ix);
    float dx = hist->GetXaxis()->GetBinWidth(ix)/2.0;
    float y=0, dy=0;
    // Loop over the vertical dimension until a non-zero entry is found
    for(int iy = hist->GetNbinsY(); iy>0; iy--){
      if( hist->GetBinContent(ix,iy) > 0 ){
	y  = hist->GetYaxis()->GetBinCenter(iy);
	dy = hist->GetYaxis()->GetBinWidth(iy)/2.0;
	break;
      }
    } // end loop over Y
    int nextPoint = gr->GetN();
    gr->SetPoint( nextPoint, x, y);
    gr->SetPointError(nextPoint, dx, dy);
  } // end loop over X

  return gr;
}



