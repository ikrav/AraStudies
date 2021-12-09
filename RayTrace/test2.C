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
  rayTracer->SetVerbosity(2);

  // 
  // Run the ray tracer and find the solutions
  // 

  int numSolutions = 0;
  bool success = rayTracer->TraceRay(numSolutions, true); // true: draw rays onto TCanvas, false: do not, if no argument supplied defaults to false
  printf("Number of solutions: %d\n", numSolutions);


  // 
  // Now that we have run the ray tracer, retrieve the results and work with them
  // 

  if( success ){

    // This code below retrieves the propagation time and the ray trajectory length
    // as computed by the ray tracer code
    std::vector<double> distanceVector = rayTracer->GetTrajectoryDistanceVector();
    std::vector<double> timeVector     = rayTracer->GetTravelTimeVector();
    for( int i=0; i < numSolutions; i++){
      printf("Solution %2i   distance= %f    time= %f\n", i, distanceVector.at(i), timeVector.at(i));
    }
    
    // Retrieve trajectories 
    std::vector<TGraph*> rayTrajectories = rayTracer->GetRayTrajectoryVector();
    printf("Number of retrieved trajectories %ld\n", rayTrajectories.size());

    // Prepare variables for the rays
    std::vector<int> rayNumberOfPointsVector;
    std::vector<double*> rayXYValuesVector;
    std::vector<double*> rayZValuesVector;

    // Draw them on another canvas
    TCanvas *cv2 = new TCanvas("cv2","",50,50,600,600);
    TH2D *hist2 = new TH2D("hist22", "", 100, xyMin, xyMax, 100, zMin, zMax);   
    hist2->SetTitle("Ray Trace; XY Distance [m]; Depth [m]");
    hist2->SetStats(0);
    hist2->DrawClone();
    line->Draw();
    for(uint iGraph=0; iGraph < rayTrajectories.size(); iGraph++){
      TGraph *thisGraph = rayTrajectories.at(iGraph);
      thisGraph->Draw("L,same");
      // Save the array of ray points
      int numberOfPoints =  thisGraph->GetN();
      double *rayXYValues = thisGraph->GetX(); // horizontal coordinate of ray points
      double *rayZValues  = thisGraph->GetY(); // vertical coordinate of ray points
      rayNumberOfPointsVector.push_back( numberOfPoints );
      rayXYValuesVector.push_back( rayXYValues );
      rayZValuesVector .push_back( rayZValues  );
    }

    // As an example, compute the curve length of each ray from the stored values
    for(uint iRay = 0; iRay < rayTrajectories.size(); iRay++){
      // For this ray, compute the curve length
      double length = 0;
      int nPoints = rayNumberOfPointsVector.at(iRay);
      double *xyValues = rayXYValuesVector.at(iRay);
      double *zValues = rayZValuesVector.at(iRay);
      // Loop over points of this ray
      for(uint iPoint = 1; iPoint < nPoints; iPoint++){
	double deltaXY = xyValues[iPoint] - xyValues[iPoint-1];
	double deltaZ  = zValues [iPoint] - zValues [iPoint-1];
	double deltaLength = sqrt(deltaXY*deltaXY + deltaZ*deltaZ);
	length += deltaLength;
      }
      printf("Manual calculation of the ray length yielded L = %f\n", length);
    }

  }
  
  if(success) cv->SaveAs("PreciseRadialRayTracer.pdf");
  else printf("Solution not possible ... \n");
}
