void runt(){
  //
  gROOT->ProcessLine(".L AliXXXAux.cxx+");
  gROOT->ProcessLine(".L AliXXXLayer.cxx+");
  gROOT->ProcessLine(".L AliXXXITSTracker.cxx+");
  gROOT->ProcessLine(".L Process.C+");
  //  Process("~/ppbench/AliESDs.root"); 
  Process("/data1/LHC14c2_p2/195483/001/AliESDs.root");
  //
}
