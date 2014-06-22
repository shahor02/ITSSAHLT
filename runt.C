{
  //
  gROOT->ProcessLine(".L AliXXXAux.cxx+");
  gROOT->ProcessLine(".L AliXXXLayer.cxx+");
  gROOT->ProcessLine(".L AliXXXITSTracker.cxx+");
  gROOT->ProcessLine(".L Process.C+");
  Process("~/ppbench/AliESDs.root"); 
  //
}
