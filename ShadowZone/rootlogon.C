{

  printf("Executing rootlogon.C\n");

  printf("  load ASlibEdit.cxx\n");
  gROOT->ProcessLine(".L ../RayTrace/ASlibEdit.cxx+");

  printf("  load PreciseRadialRayTracer.cxx\n");
  gROOT->ProcessLine(".L ../RayTrace/PreciseRadialRayTracer.cxx+");

}
