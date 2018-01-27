{

  gSystem->Load("PreciseRadialRayTracer_cxx.so");

  PreciseRadialRayTracer *rayTracer = new PreciseRadialRayTracer();

  TCanvas *cv = new TCanvas("cv","",0,0,600,600);

  double originZ = -500.0;
  double destinationXYDistance = 500.0;
  double destinationZ = -180.0;

  // Find limits for the figure of the rays
  float margin = 20.0; //meters
  // The top is just above the surface if source and destination are in-ice. Otherwise,
  // the top of the figure is just above the max(source, destination)
  float zMax = std::max( (float)0.0, std::max((float)originZ, (float)destinationZ) ) + margin;
  // The bottom of the plot is similar
  float zMin = std::min( (float)0.0, std::min((float)originZ, (float)destinationZ) ) - margin;
  // The left and right are a bit wider than the XY distance
  float xyMin = -1*margin;
  float xyMax = destinationXYDistance + margin;
  TH2D *hist = new TH2D("hist", "", 100, xyMin, xyMax, 100, zMin, zMax);   
  hist->SetTitle("Ray Trace; XY Distance [m]; Depth [m]");
  hist->SetStats(0);
  hist->DrawClone();
  delete hist;
  // Draw a line at the surface of the ice
  TLine *line = new TLine(xyMin, 0.0, xyMax, 0.0);
  line->Draw("same");


  rayTracer->SetOrigin(0.0,0.0,originZ); //set location of event
  rayTracer->SetDestination(0.0,destinationXYDistance,destinationZ); // set location of signal observation
  //rayTracer->SetVerbosity(2);

  int numSolutions = 0;
  bool success = rayTracer->TraceRay(numSolutions, true); // true: draw rays onto TCanvas, false: do not, if no argument supplied defaults to false
  printf("Number of solutions: %d\n", numSolutions);
  
  if(success) cv->SaveAs("PreciseRadialRayTracer.pdf");
  else printf("Solution not possible ... \n");
}
