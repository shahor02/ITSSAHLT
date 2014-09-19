#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH1F.h>
#include <TH1I.h>
#include <TH2D.h>
#include <TFile.h>
#include <Riostream.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <TArrayF.h>
#include <TStyle.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TParticlePDG.h>
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliMultiplicity.h"
#include "AliESDtrackCuts.h"
#include "AliStack.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliITSRecPointContainer.h"
#include "AliITSgeomTGeo.h"
#include "AliITSRecPoint.h"
#include "AliCDBManager.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
//
#include "HistoManager.h"

#include "AliXXXAux.h"
#include "AliXXXLayer.h"
#include "AliXXXITSTracker.h"
#include "TProfile.h"
#include "TStopwatch.h"
#endif


AliXXXITSTracker* tracker=0;

TTree* treeInp = 0;
AliRunLoader* runLoader = 0;
AliESDEvent* esd=0;
typedef struct {
  Bool_t  validSPD;
  Bool_t  validTrc;
  Bool_t  isPrim;
  Char_t  ntlC;
  Char_t  ntlF;
  Char_t  mtlC;
  Char_t  mtlF;
  Char_t  ntrC;
  Char_t  ntrF;
  Int_t   pdg;
  Int_t   evID;
  Int_t   mult;
  Float_t pt;
  Float_t eta;
  Float_t phi;
  Float_t zv;
  //
  Float_t spdDPhi;
  Float_t spdDTht;
  Float_t spdChi2;
} tres_t;

tres_t tres;
TFile* flOut = 0;
TTree* trOut = 0;


void ProcessEvent(int iev);
void Process(const char* path);
void ProcChunk(const char* path);
void TestTracker(TTree* tRP, const AliESDVertex* vtx);
void LoadClusters(TTree* tRP);
Double_t* DefLogAx(double xMn,double xMx, int nbin);
void CheckRecStatus();
void InitOutTree(const char* fname="rcInfo.root");
void SaveOutTree();



#ifdef _TIMING_
TProfile* hTiming[AliXXXITSTracker::kNSW];
#endif

void Process(const char* inpData)
{
  //
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetRun(0);
  tracker = new AliXXXITSTracker();
  tracker->Init();
  tracker->SetMaxMissedLayers(1);
  //
  InitOutTree();
  //
#ifdef _TIMING_
  int nbMlt = 100;
  double *mltAx = DefLogAx(1,6000,nbMlt);
  for (int i=0;i<AliXXXITSTracker::kNSW;i++) {
    hTiming[i] = new TProfile(tracker->GetStopwatchName(i),tracker->GetStopwatchName(i),nbMlt,mltAx);
  }
#endif
  //
  TString inpDtStr = inpData;
  if (inpDtStr.EndsWith(".root") || inpDtStr.EndsWith(".zip#")) {
    ProcChunk(inpDtStr.Data());
  }
  else {
    ifstream inpf(inpData);
    if (!inpf.good()) {
      printf("Failed on input filename %s\n",inpData);
      return;
    }
    //
    inpDtStr.ReadLine(inpf);
    while ( !inpDtStr.IsNull() ) {
      inpDtStr = inpDtStr.Strip(TString::kBoth,' ');
      if (inpDtStr.BeginsWith("//") || inpDtStr.BeginsWith("#")) {inpDtStr.ReadLine(inpf); continue;}
      inpDtStr = inpDtStr.Strip(TString::kBoth,',');
      inpDtStr = inpDtStr.Strip(TString::kBoth,'"');
      ProcChunk(inpDtStr.Data());
      inpDtStr.ReadLine(inpf);
    }
  }
#ifdef _CONTROLH_
  tracker->SaveHistos();
#endif
  SaveOutTree();
}

void ProcChunk(const char* path)
{
  //
  TString dir=path;
  //
  printf("Processing %s\n",dir.Data());
  //
  if (dir.EndsWith("AliESDs.root")) {
    dir.ReplaceAll("AliESDs.root","");
  }
  //
  esd = new AliESDEvent();
  //
  runLoader = AliRunLoader::Open(Form("%sgalice.root",dir.Data()));
  if (!runLoader) {
    printf("galice not found\n");
    return;
  }
  runLoader->LoadgAlice();
  gAlice = runLoader->GetAliRun();
  runLoader->LoadHeader();
  runLoader->LoadKinematics();
  runLoader->LoadRecPoints("ITS");
  // ESD
  TFile* esdFile = TFile::Open(Form("%sAliESDs.root",dir.Data()));
  if (!esdFile || !esdFile->IsOpen()) {
    printf("Error in opening ESD file\n");
    runLoader->UnloadKinematics();
    runLoader->UnloadHeader();
    runLoader->UnloadgAlice();
    delete runLoader;
    return;
  }
  treeInp = (TTree*) esdFile->Get("esdTree");
  if (!treeInp) {
    printf("Error: no ESD tree found\n");
    runLoader->UnloadKinematics();
    runLoader->UnloadHeader();
    runLoader->UnloadgAlice();
    delete runLoader;
    return;
  }
  esd->ReadFromTree(treeInp);
  //
  //  for(Int_t iEv=0; iEv<100;iEv++) {
  for (int iEv=0;iEv<runLoader->GetNumberOfEvents(); iEv++) {
    //    for(Int_t iEv=93; iEv<=95; iEv++){
    printf("ev %d\n",iEv);
    ProcessEvent(iEv);
  }
  runLoader->UnloadRecPoints("all");
  runLoader->UnloadKinematics();
  runLoader->UnloadHeader();
  runLoader->UnloadgAlice();
  delete runLoader; 
  runLoader = 0;
  esdFile->Close();
  delete esdFile;
}

//_______________________________________________
void ProcessEvent(int iEv)
{
  runLoader->GetEvent(iEv);
  treeInp->GetEvent(iEv);
  AliStack *stack = runLoader->Stack(); 
  Int_t nParticles = stack->GetNtrack();
  Int_t nTracks= esd->GetNumberOfTracks();
  printf("Ntr: %d NPart: %d\n",nTracks,nParticles);
  const AliESDVertex* vtx = esd->GetPrimaryVertexSPD();
  if (!vtx || !vtx->GetStatus()) return;
  //
  //  if (esd->GetMultiplicity()->GetNumberOfTracklets()<1000) return;
  //if (esd->GetMultiplicity()->GetNumberOfTracklets()>300||esd->GetMultiplicity()->GetNumberOfTracklets()<30) return;
  TTree* treeITSRP = runLoader->GetTreeR("ITS",kFALSE);
  if (!treeITSRP) {
    printf("Failed to fetch ITS recpoints\n");
    exit(1);
  }
  //
#ifdef _TIMING_
  static Double_t cpuTime[AliXXXITSTracker::kNSW] = {0};
#endif
  printf("\n\n\nEvent: %d %s/%s\n",iEv,treeITSRP->GetDirectory()->GetFile()->GetName(),treeITSRP->GetDirectory()->GetName());
  //  treeITSRP->Print();
  vtx->Print();
  //
  // MC Vertex
  AliGenEventHeader *hdmc = runLoader->GetHeader()->GenEventHeader();
  TArrayF vtxF(3);
  hdmc->PrimaryVertex(vtxF);
  printf("\nMCVertex : %+.4f %+.4f %+.4f -> %+e %+e %+e\n",vtxF[0],vtxF[1],vtxF[2], vtxF[0]-vtx->GetX(),vtxF[1]-vtx->GetY(),vtxF[2]-vtx->GetZ());
  Double_t xyz[]={vtxF[0],vtxF[1],vtxF[2]};
  Double_t ers[]={1E-4,1E-4,1E-4};
  AliESDVertex vertexMC(xyz,ers);
  //TestTracker(treeITSRP, &vertexMC);
  TestTracker(treeITSRP, vtx);
  //
#ifdef _TIMING_
  if (vtx->GetStatus()==1) {
    double mlt = esd->GetMultiplicity()->GetNumberOfTracklets();
    for (int i=0;i<AliXXXITSTracker::kNSW;i++) {
      TStopwatch& sw = (TStopwatch&)tracker->GetStopwatch(i);
      double tm = sw.CpuTime();
      hTiming[i]->Fill(mlt, tm - cpuTime[i]);
      cpuTime[i] = tm;
    }
  }
#endif
  //
  CheckRecStatus();
  //esd->Reset();
  //
}

//_________________________________________________
void TestTracker(TTree* tRP, const AliESDVertex* vtx)
{
  tracker->Clear(); 
  LoadClusters(tRP); 
  tracker->SetSPDVertex(vtx);
  tracker->SetBz(esd->GetMagneticField());
  tracker->ProcessEvent();
  //  tracker->PrintTracklets();
  //  tracker->PrintTracks();  
  //
  //  vtx->Print();
  //  esd->GetMultiplicity()->Print("t");
  //
  /*
  AliHeader* hd = runLoader->GetHeader();
  AliGenEventHeader* hdmc;
  if (tracker->GetTrackVertex().GetStatus()==1 && hd && (hdmc=hd->GenEventHeader()) ) {
    TArrayF vtxMC;
    hdmc->PrimaryVertex(vtxMC);
    printf("MCvtx: %f %f %f\n",vtxMC[0],vtxMC[1],vtxMC[2]);
  }
  */
}

//_________________________________________________
void LoadClusters(TTree* tRP)
{
  //  AliITSRecPointContainer* rpcont = AliITSRecPointContainer::Instance();
  //  TClonesArray* itsClusters = rpcont->FetchClusters(0,tRP);
  //  if(!rpcont->IsSPDActive()){
  //    printf("No SPD rec points found, multiplicity not calculated\n");
  //    tRP->Print();
  //    return;
  //  }
  //  else printf("NSP0: %d\n",itsClusters->GetEntriesFast());
  static TClonesArray* itsModAdd[2198] = {0};
  static Bool_t first = kTRUE;
  if (first) {
    first = 0;
    for (int i=0;i<2198;i++) itsModAdd[i] = new TClonesArray("AliITSRecPoint");
  }
  int nMod = AliITSgeomTGeo::GetNModules();

  for (int imd=0;imd<nMod;imd++) {
    TClonesArray* itsClusters = itsModAdd[imd];
    tRP->SetBranchAddress("ITSRecPoints",&itsClusters);    
    tRP->GetEntry(imd);
    //    itsClusters = rpcont->UncheckedGetClusters(imd);
    
    int nClusters = itsClusters->GetEntriesFast();
    if (!nClusters) continue;
    while(nClusters--) {
      AliITSRecPoint* cluster = (AliITSRecPoint*)itsClusters->UncheckedAt(nClusters);
      if (!cluster) continue;
      tracker->AddCluster(cluster);
    }
  }
  //
}

//______________________________________________
Double_t* DefLogAx(double xMn,double xMx, int nbin)
{
  // get array for log axis
  if (xMn<=0 || xMx<=xMn || nbin<2) {
    printf("Wrong axis request: %f %f %d\n",xMn,xMx,nbin);
    return 0;
  }
  double dx = log(xMx/xMn)/nbin;
  double *xax = new Double_t[nbin+1];
  for (int i=0;i<=nbin;i++) xax[i]= xMn*exp(dx*i);
  return xax;
}

//______________________________________________
void CheckRecStatus()
{
  //
  static int nev = -1;
  static TBits patternMC;
  static vector<char> ntrCorr;
  static vector<char> ntrFake;
  static vector<char> mtlCorr;
  static vector<char> mtlFake;
  static vector<char> ntlCorr;
  static vector<char> ntlFake;
  static vector<float> spdDPhi;
  static vector<float> spdDTht;
  static vector<float> spdChi2;
  //
  AliStack* stack = 0;
  AliRunLoader* rl = AliRunLoader::Instance();
  if (!rl || !(stack=rl->Stack())) return;
  nev++;
  //
  enum {kIsPrim=AliXXXITSTracker::kNLrActive,kValidTracklet,kValidTrack,kRecDone,kBitPerTrack};
  int nTrkMC = stack->GetNtrack();
  patternMC.SetBitNumber(nTrkMC*kBitPerTrack,0);
  patternMC.ResetAllBits();
  //
  ntrCorr.clear();
  ntrFake.clear();
  ntlCorr.clear();
  ntlFake.clear();
  mtlCorr.clear();
  mtlFake.clear();
  //
  spdDPhi.clear();
  spdDTht.clear();
  spdChi2.clear();
  //  
  ntrCorr.resize(nTrkMC,0);
  ntrFake.resize(nTrkMC,0);
  ntlCorr.resize(nTrkMC,0);
  ntlFake.resize(nTrkMC,0);
  mtlCorr.resize(nTrkMC,0);
  mtlFake.resize(nTrkMC,0);
  spdDPhi.resize(nTrkMC,0);
  spdDTht.resize(nTrkMC,0);
  spdChi2.resize(nTrkMC,0);
  
  //
  // fill MC track patterns
  for (int ilr=6;ilr--;) {
    AliXXXLayer *lr = tracker->GetLayer(ilr);
    int ncl = lr->GetNClusters();
    for (int icl=ncl;icl--;) {
      AliITSRecPoint* cl = lr->GetClusterUnSorted(icl);
      for (int j=0;j<3;j++) {
	int lb = cl->GetLabel(j);
	if (lb<0 || lb>=nTrkMC) break;
	patternMC.SetBitNumber(lb*kBitPerTrack+ilr,kTRUE);
      }
    }
  }
  //
  int nTrk = tracker->GetNTracklets();
  if (!nTrk) return;
  for (int itr=0;itr<nTrk;itr++) {
    const AliXXXITSTracker::SPDtracklet_t& trlet = tracker->GetTracklet(itr);
    //    PrintTracklet(itr);
    //
    int lbl = trlet.label;
    if (lbl==-3141593) continue;
    int lblA = TMath::Abs(lbl);
    if (lblA==lbl) ntlCorr[lblA]++;
    else           ntlFake[lblA]++;
    //
    if (spdChi2[lblA]==0 || spdChi2[lblA]<trlet.chi2) {
      spdChi2[lblA] = trlet.chi2;
      spdDPhi[lblA] = trlet.dphi;
      spdDTht[lblA] = trlet.dtht;
    }
  }
  //
  AliMultiplicity* mltESD = esd->GetMultiplicity();
  nTrk = mltESD->GetNumberOfTracklets();
  for (int itr=0;itr<nTrk;itr++) {
    int lb0 = mltESD->GetLabel(itr,0);
    int lb1 = mltESD->GetLabel(itr,1);    
    if (lb0==lb1) mtlCorr[lb1]++;
    else          mtlFake[lb1]++;
  }
  //
  nTrk = tracker->GetNTracks();
  for (int itr=0;itr<nTrk;itr++) {
    const AliXXXITSTracker::ITStrack_t &track = tracker->GetTrack(itr);
    //
    int lbl = track.label;
    if (lbl==-3141593) continue;
    int lblA = TMath::Abs(lbl);
    if (lblA==lbl) ntrCorr[lblA]++;
    else           ntrFake[lblA]++;
  }
  //
  // set reconstructability
  for (int itr=nTrkMC;itr--;) {
    int bitoffs = itr*kBitPerTrack;
    //
    tres.validSPD = patternMC.TestBitNumber(bitoffs+AliXXXITSTracker::kALrSPD1) && 
      patternMC.TestBitNumber(bitoffs+AliXXXITSTracker::kALrSPD2);
    int nmiss = 0;
    for (int il=AliXXXITSTracker::kALrSDD1;il<=AliXXXITSTracker::kALrSSD2;il++) 
      if (tracker->IsObligatoryLayer(il) && !patternMC.TestBitNumber(bitoffs+il)) nmiss++;
    tres.validTrc = tres.validSPD && (nmiss<=tracker->GetMaxMissedLayers());
    //
    if ( !tres.validSPD && !tres.validTrc && 
	 (ntlCorr[itr]+ntlFake[itr])==0 && 
	 (ntrCorr[itr]+ntrFake[itr])==0 &&
	 (mtlCorr[itr]+mtlFake[itr])==0 ) continue;
    //
    TParticle* mctr = stack->Particle(itr);
    tres.isPrim = stack->IsPhysicalPrimary(itr);
    tres.pdg = mctr->GetPdgCode();
    tres.pt =  mctr->Pt();
    tres.eta = mctr->Eta();
    tres.phi = mctr->Phi();
    //
    tres.ntlC = ntlCorr[itr];
    tres.ntlF = ntlFake[itr];
    tres.ntrC = ntrCorr[itr];
    tres.ntrF = ntrFake[itr];
    //
    tres.mtlC = mtlCorr[itr];
    tres.mtlF = mtlFake[itr];
    //
    tres.evID = nev;
    tres.mult = tracker->GetNTracklets();
    tres.zv = tracker->GetSPDVertex()->GetZ();
    //
    tres.spdDPhi = spdDPhi[itr];
    tres.spdDTht = spdDTht[itr];
    tres.spdChi2 = spdChi2[itr];
    //
    trOut->Fill();
  }
  //
}

//___________________________________________________
void InitOutTree(const char* fname)
{
  // output tree
  flOut = TFile::Open(fname,"recreate");
  trOut = new TTree("rcinfo","rcinfo");
  trOut->Branch("res",&tres);
}

//___________________________________________________
void SaveOutTree()
{
  if (!trOut) return;
  flOut->cd();
  trOut->Write();
  delete trOut;
  flOut->Close();
  delete flOut;
  trOut = 0;
  flOut = 0;
}
