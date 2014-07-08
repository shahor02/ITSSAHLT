void runt(){
  //
  TString dtPath = "/data1/LHC14c2_p2/195483/001";
  //  TString dtPath = "~/ppbench";
  //
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local:///home/shahoian/ALICE/Aliroot/OCDB");
  //  man->SetSpecificStorage("ITS/Align/Data","alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
  //  man->SetSpecificStorage("ITS/Align/Data","alien://folder=/alice/simulation/2008/v4-15-Release/Ideal");
  man->SetSpecificStorage("ITS/Align/Data","local:///alice/simulation/2008/v4-15-Release/Ideal");

  gROOT->ProcessLine(".L AliXXXAux.cxx+");
  gROOT->ProcessLine(".L AliXXXLayer.cxx+");
  gROOT->ProcessLine(".L AliXXXITSTracker.cxx+");
  gROOT->ProcessLine(".L Process.C+");
  gSystem->Exec(Form("ln -s -f %s/geometry.root ./",dtPath.Data()));
  //
  TString inpData = Form("%s/AliESDs.root",dtPath.Data());
  printf("InputData : %s\n",inpData.Data());
  Process(inpData.Data());
  //
}
