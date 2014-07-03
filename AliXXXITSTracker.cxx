#include "AliXXXITSTracker.h"
#include "AliXXXLayer.h"
#include "AliITSRecPoint.h"
#include "AliESDVertex.h"
#include "AliGeomManager.h"
#include "AliVParticle.h"



ClassImp(AliXXXITSTracker)

const Float_t AliXXXITSTracker::fgkZSpanITS[AliXXXITSTracker::kMaxLrITS] = 
{ 36. ,14.1,14.1,  38., 22.2,29.7, 51.   ,43.1,48.9};

const Float_t AliXXXITSTracker::fgkRLayITS[AliXXXITSTracker::kMaxLrITS] = 
  { 2.94, 3.9,7.6, 11.04, 15.0,23.9, 29.44 ,38.0,43.0};

const Float_t AliXXXITSTracker::fgkRSpanITS[AliXXXITSTracker::kMaxLrITS] = // half span in R
  { 0.04, 0.5,0.5, 0.5, 0.8, 0.8, 0.5 ,0.6,0.6};

const Int_t    AliXXXITSTracker::fgkPassivLrITS[AliXXXITSTracker::kNLrPassive] = 
  {AliXXXITSTracker::kLrBeamPime,AliXXXITSTracker::kLrShield1,AliXXXITSTracker::kLrShield2};

const Int_t    AliXXXITSTracker::fgkActiveLrITS[AliXXXITSTracker::kNLrActive] = 
  {AliXXXITSTracker::kLrSPD1,AliXXXITSTracker::kLrSPD2,
   AliXXXITSTracker::kLrSDD1,AliXXXITSTracker::kLrSDD2,
   AliXXXITSTracker::kLrSSD1,AliXXXITSTracker::kLrSSD2};

const Int_t    AliXXXITSTracker::fgkLr2Active[AliXXXITSTracker::kMaxLrITS] = // conversion to active lr.
  {-1, 0, 1, -1, 2, 3, -1, 4, 5};

const Float_t AliXXXITSTracker::fgkRhoLITS[AliXXXITSTracker::kMaxLrITS] = {
  0.162802, 0.321960,0.354588, 0.274995, 0.193789,0.198168, 0.435372, 0.195828,0.226940};

const Float_t AliXXXITSTracker::fgkX2X0ITS[AliXXXITSTracker::kMaxLrITS] = {
  0.002757, 0.011660,0.012614, 0.006488, 0.007714,0.007916, 0.012689, 0.007849,0.009128};


const Double_t AliXXXITSTracker::fgkClSystYErr2[AliXXXITSTracker::kNLrActive] = 
  {0.0010*0.0010, 0.0030*0.0030, 0.0500*0.0500, 0.0500*0.0500, 0.0020*0.0020, 0.0020*0.0020};

const Double_t AliXXXITSTracker::fgkClSystZErr2[AliXXXITSTracker::kNLrActive] = 
  {0.0050*0.0050, 0.0050*0.0050, 0.0050*0.0050, 0.0050*0.0050, 0.1000*0.1000, 0.1000*0.1000};


const Int_t    AliXXXITSTracker::fgkLrDefBins[AliXXXITSTracker::kNLrActive][2] = // n bins in z, phi
  { {20,20}, {20,20}, {20,20}, {20,20}, {20,20}, {20,20} };

const Float_t AliXXXITSTracker::fgkDefMass = 0.14;
const Int_t   AliXXXITSTracker::fgkDummyLabel = -3141593;

#ifdef _TIMING_
const char* AliXXXITSTracker::fgkSWNames[AliXXXITSTracker::kNSW] = {
  "Total:     "
  ,"Tracklets: "
  ,"Tracks:    "
  ,"Vertex:    "
};
#endif


//______________________________________________
AliXXXITSTracker::AliXXXITSTracker() :
  fBlacklist(0)
  ,fPhiShift(0.0045)
  ,fSigThetaTracklet(0.025)
  ,fSigPhiTracklet(0.08)
  ,fChi2CutTracklet(1.5)
  ,fPhiShiftSc(0.)
  ,fDThetaTrackletSc(0)
  ,fDPhiTrackletSc(0)
  ,fBz(5.0)
  ,fDPhiTol(0.)
  ,fDThSig2Inv(0.)
  ,fDPhSig2Inv(0.)
  //
  ,fMinPt(0.3)
  ,fCurvMax(0)
  ,fZSPD2CutMin(1e9)
  ,fZSPD2CutMax(-1e9)
  ,fMaxChi2Tr2Cl(35)
  ,fAddErr2YspdVtx(0.02*0.02)
  ,fAddErr2ZspdVtx(0.04*0.04)
  //
  ,fMissChi2Penalty(3)
  ,fMaxMissedLayers(1)
  ,fNTracks(0)
  //
  ,fSPDVertex(0)
{
  // def. c-tor
  for (int i=kNLrActive;i--;) fLayers[i] = 0;
}

//______________________________________________
AliXXXITSTracker::~AliXXXITSTracker()
{
  // d-tor
  for (int i=0;i<kNLrActive;i++) delete fLayers[i];
}

//______________________________________________
void AliXXXITSTracker::Init()
{
  // init tracker
  //
  if (!AliGeomManager::GetGeometry()) {
    AliGeomManager::LoadGeometry("geometry.root");
    AliGeomManager::ApplyAlignObjsFromCDB("ITS");
  }
  //
  for (int i=0;i<kNLrActive;i++) {
    int iAct = fgkActiveLrITS[i];
    fLayers[i] = new AliXXXLayer(i,fgkZSpanITS[iAct]+1,fgkLrDefBins[i][0],fgkLrDefBins[i][1]);
    fSkipLayer[i] = kFALSE;
    fNSigma2[i] = 5*5;
    fYToler2[i] = 10;//0.2*0.2;
    fZToler2[i] = 10;//0.2*0.2;
    fChi2TotCut[i] = 0;
  }  
  fChi2TotCut[1] = 10; // 2 cl+vtx -> NDF=1
  fChi2TotCut[2] = 15; 
  fChi2TotCut[3] = 20; 
  fChi2TotCut[4] = 24; 
  fChi2TotCut[5] = 26; 
  //
  fMissChi2Penalty = 3;
  fMaxMissedLayers = 1;
  //
  // auxialary precalculated variables
  if (fChi2CutTracklet<0.1) fChi2CutTracklet = 0.1;
  double scl = TMath::Sqrt(fChi2CutTracklet);
  fDThetaTrackletSc = fSigThetaTracklet*scl;
  fDPhiTrackletSc   = fSigPhiTracklet*scl;
  //
  fDThSig2Inv = 1./(fSigThetaTracklet*fSigThetaTracklet);
  fDPhSig2Inv = 1./(fSigPhiTracklet*fSigPhiTracklet);
  //
  fBlacklist = new TBits(100*100);
  //
#ifdef _TIMING_
  for (int i=kNSW;i--;) {
    fSW[i].Stop();
    fSW[i].Reset();
  }
#endif
}

//______________________________________________
void AliXXXITSTracker::ProcessEvent()
{
  // do full reconstruction
#ifdef _TIMING_
  fSW[kSWTotal].Start(0);
  fSW[kSWTracklets].Start(0);
#endif
  //
  FindTracklets();
  //
#ifdef _TIMING_
  fSW[kSWTracklets].Stop();
  fSW[kSWTracks].Start(0);
#endif
  //
  Tracklets2Tracks();
  RefitInward();
  //
#ifdef _TIMING_
  fSW[kSWTracks].Stop();
#endif
  //
#ifdef _TIMING_
  fSW[kSWTotal].Stop();
  PrintTiming();
#endif
  //
}


//______________________________________________
void AliXXXITSTracker::Clear(Option_t*)
{
  // reset event info
  ClearTracklets();
  ClearTracks();
  for (int i=kNLrActive;i--;) {
    fNClusters[i] = 0;
    if (fLayers[i]) fLayers[i]->Clear();
  }
}

//______________________________________________
void AliXXXITSTracker::ClearTracklets()
{
  // reset tracklets info
  fSPD2Discard.clear();
  fTracklets.clear();
  fSPD1Tracklet.clear();
  if (fBlacklist) fBlacklist->ResetAllBits();
}


//______________________________________________
void AliXXXITSTracker::AddCluster(AliITSRecPoint* cl)
{
  // add cluster to corresponding layer
  if (!cl->Misalign()) AliWarning("Can't misalign this cluster !"); 
  fLayers[cl->GetLayer()]->AddCluster(cl); 
}

//______________________________________________
Bool_t AliXXXITSTracker::FindTracklets()
{
  // find SPD tracklets
  //
  if (!fSPDVertex) {
    AliInfo("No SPD vertex set");
    return kFALSE;
  }
  float rv2 = fSPDVertex->GetX()*fSPDVertex->GetX()+fSPDVertex->GetY()*fSPDVertex->GetY();
  if (rv2>0.25*fgkRLayITS[kLrBeamPime]*fgkRLayITS[kLrBeamPime]) {
    AliInfo("SPD vertex is too far from beam line");
    fSPDVertex->Print();
    return kFALSE;    
  } 
  fPhiShiftSc = fPhiShift*TMath::Abs(fBz/5.0);
  fDPhiTol = fDPhiTrackletSc + fPhiShiftSc;
  //
  AliXXXLayer &spdL1 = *fLayers[kALrSPD1];
  AliXXXLayer &spdL2 = *fLayers[kALrSPD2];
  spdL1.SortClusters(fSPDVertex);
  spdL2.SortClusters(fSPDVertex);
  fNClusters[0] = spdL1.GetNClusters();
  fNClusters[1] = spdL2.GetNClusters();
  //
  fSPD2Discard.resize(fNClusters[1]);
  fSPD1Tracklet.resize(fNClusters[0]);
  //
  fBlacklist->SetBitNumber(TMath::Max(fNClusters[0]*fNClusters[1],10000),kFALSE); // to reserve the space
  //
  int nfound;
  do {
    nfound = 0;
    for (int icl2=fNClusters[1];icl2--;) if (!fSPD2Discard[icl2]) nfound += AssociateClusterOfL2(icl2);
  } while(nfound);
  //
  return kTRUE;
}

//______________________________________________
Int_t AliXXXITSTracker::AssociateClusterOfL2(int icl2)
{
  // find SPD1 cluster matching to SPD2 cluster icl2
  AliXXXLayer &spdL1 = *fLayers[kALrSPD1];
  AliXXXLayer &spdL2 = *fLayers[kALrSPD2];
  AliXXXLayer::ClsInfo* cli2 = spdL2.GetClusterInfo(icl2);
  // expected z at SPD1
  float zV = fSPDVertex->GetZ();
  float z2 = cli2->z - zV;
  float tg2Inv = z2/cli2->r;
  float dzt = (1.+tg2Inv*tg2Inv)*fDThetaTrackletSc;
  float dz = dzt*fgkRLayITS[kLrSPD1] + TMath::Abs(tg2Inv)*fgkRSpanITS[kLrSPD1]; // uncertainty from dTheta and from Layer1 R spread
  float zL1 = zV + tg2Inv*fgkRLayITS[kLrSPD1]; // center of expected Z1
  int nsel1 = spdL1.SelectClusters(zL1-dz,zL1+dz, cli2->phi-fDPhiTol,cli2->phi+fDPhiTol);
  if (!nsel1) {
    fSPD2Discard[icl2] = true;
    return 0; // no candidates
  }
  float chiBest = 9999;
  SPDtracklet_t trk;
  trk.id1 = -1;
  int icl1,nCand=0;
  while ( (icl1=spdL1.GetNextClusterInfoID())!=-1) {  // loop over matching clusters of lr1
    if (IsBlacklisted(icl1,icl2)) continue;
    AliXXXLayer::ClsInfo* cli1 = spdL1.GetClusterInfo(icl1);
    float z1 = cli1->z - zV;
    float tg1Inv = z1/cli1->r;
    //
    float dTheta = (tg2Inv-tg1Inv)/(1.+tg1Inv*tg1Inv);        // fast check on theta
    if (TMath::Abs(dTheta)>fDThetaTrackletSc) continue;
    //
    float dPhi = cli1->phi - cli2->phi;                       // fast check on phi
    if (dPhi>TMath::Pi()) dPhi = TMath::TwoPi()-dPhi;
    double dPhiS = TMath::Abs(dPhi)-fPhiShiftSc;
    if (TMath::Abs(dPhiS)>fDPhiTrackletSc) continue;
    //
    float chi2 = dTheta*dTheta*fDThSig2Inv + dPhiS*dPhiS*fDPhSig2Inv; // check final chi2
    if (chi2>1.) continue;
    nCand++;
    if (chi2>chiBest) continue;
    // check if cl1 is already associated with better 
    trk.id1 = icl1;
    trk.id2 = icl2;
    trk.dtht = dTheta;
    trk.dphi = dPhi;
    trk.chi2 = chiBest = chi2;
  }
  //
  if (trk.id1!=-1) { // check if there is no better icl1 candidate for icl2
    int oldId = fSPD1Tracklet[trk.id1];
    if (!oldId) { // store new tracklet
      fTracklets.push_back(trk);
      fSPD1Tracklet[trk.id1] = fTracklets.size(); // refer from clusters to tracklet (id+1)
      Blacklist(trk.id1,trk.id2);
      return 1;
    }
    SPDtracklet_t& oldTrk = (SPDtracklet_t&)fSPD1Tracklet[--oldId];
    if (oldTrk.chi2 < trk.chi2) { // previous is better 
      Blacklist(trk.id1,trk.id2);  // shall we blacklist new combination?
      if (nCand==1)  fSPD2Discard[icl2] = true; // there was just 1 candidate and it is discarded
      return 0;
    }
    oldTrk = trk;         // new combination is better, overwrite it with new one
    Blacklist(trk.id1,trk.id2);
    return 1;
  }
  //
  fSPD2Discard[icl2] = true; // no chance to find partner for this cluster
  return 0;
  //
}


//______________________________________________
void AliXXXITSTracker::Tracklets2Tracks()
{
  // try to extend tracklets to outer layers
  int nTrk = GetNTracklets();
  if (!nTrk) return;
  //
  CalcAuxTracking(); // RS??? do we need to repeat this?
  //
  for (int ila=kALrSDD1;ila<kNLrActive;ila++) {
    if (fSkipLayer[ila]) continue;
    fLayers[ila]->SortClusters(0);
    fNClusters[ila] = fLayers[ila]->GetNClusters();
  }
  //
  fTracks.resize(nTrk);
  fNTracks = 0;
  //
  for (int itr=0;itr<nTrk;itr++) {
    SPDtracklet_t& trlet = fTracklets[itr];
    //
    printf("TestTracklet %d\t|",itr);
    int stat = GetTrackletMCTruth(trlet);
    for (int i=0;i<kNLrActive;i++) printf("%c", (stat&(0x1<<i)) ? '*':'-'); printf("|\n");
    //
    PrintTracklet(itr);
    //
    float zspd2 = fLayers[kALrSPD2]->GetClusterInfo(trlet.id2)->z;
    if (zspd2<fZSPD2CutMin || zspd2>fZSPD2CutMax) continue;
    ITStrack_t &track = fTracks[fNTracks];
    if (!CreateTrack(track, trlet)) continue;
    track.trackletID = itr;
    Bool_t res;
    printf("process track pt:%f\n",track.paramOut.Pt());
    for (int lrID=kLrShield1;lrID<kMaxLrITS;lrID++) {
      res = FollowToLayer(track,lrID) && IsAcceptableTrack(track);
      if (!res) break;
    }
    printf("%s\n",res ? "OK" : "Fail");
    if (!res) continue;
    track.paramOut.ResetBit(kInvalidBit); // flag that outward fit succeeded
    CookLabel(track);
    fNTracks++;
    //
  }  
}

//______________________________________________
Bool_t AliXXXITSTracker::IsAcceptableTrack(const AliXXXITSTracker::ITStrack_t& track) const
{
  // check if the track is acceptable
  return kTRUE;
}

//______________________________________________
void AliXXXITSTracker::PrintTrack(const AliXXXITSTracker::ITStrack_t& track) const
{
  // print track info
  printf("Chi2 = %f for %d clusters. Tracklet %d\n",track.chi2,track.ncl,track.trackletID);
  for (int ilr=0;ilr<kNLrActive;ilr++) {
    if (track.clID[ilr]<0) continue;
    AliITSRecPoint* cl = fLayers[ilr]->GetClusterSorted(track.clID[ilr]);
    printf("L%d #%4d ",ilr,track.clID[ilr]);
    for (int i=0;i<3;i++) printf("%d ",cl->GetLabel(i)); printf("\n");
  }
  track.paramOut.Print();
  track.paramInw.Print();  
}

//______________________________________________
void AliXXXITSTracker::PrintTracklets() const
{
  // print traklets info
  int ntr = fTracklets.size();
  printf("NTracklets: %d\n",ntr);
  printf("Nspd1: %4d Nspd2: %4d, Ntracklets: %d\n",fNClusters[0],fNClusters[1],ntr);
  for (int itr=0;itr<ntr;itr++) PrintTracklet(itr);
  //
}

//______________________________________________
void AliXXXITSTracker::PrintTracklet(Int_t itr) const
{
  // print single tracklet
  const SPDtracklet_t* trk = &fTracklets[itr];
  AliITSRecPoint* cl1 = fLayers[kALrSPD1]->GetClusterSorted(trk->id1);
  AliITSRecPoint* cl2 = fLayers[kALrSPD2]->GetClusterSorted(trk->id2);
  AliXXXLayer::ClsInfo_t* cli0 = fLayers[kALrSPD1]->GetClusterInfo(trk->id1);
  printf("#%3d Phi:%+.3f Eta:%+.3f Dphi:%+.3f Dtht:%+.3f Chi2:%.3f | Lbl:",
	 itr,cli0->phi,
	 -TMath::Log(TMath::Tan(TMath::ATan2(cli0->r,cli0->z-fSPDVertex->GetZ())/2.)),
	 trk->dphi,trk->dtht,trk->chi2);
  int lb = 0;
  for (int i=0;i<3;i++) if ( (lb=cl1->GetLabel(i))>=0 ) printf(" %5d",lb); printf("|");
  for (int i=0;i<3;i++) if ( (lb=cl2->GetLabel(i))>=0 ) printf(" %5d",lb); 
  printf("\n");
}


//______________________________________________
void AliXXXITSTracker::PrintTracks() const
{
  // print tracks info
  printf("NTracks: %d\n",fNTracks);
  for (int itr=0;itr<fNTracks;itr++) PrintTrack(fTracks[itr]);
  //
}


//______________________________________________
void AliXXXITSTracker::CalcAuxTracking()
{
  // precalculate auxilarry variables for tracking
  //
  // largest track curvature to search
  const double ztolerEdge = 1.0;
  fCurvMax = TMath::Abs(fBz*kB2C/fMinPt);
  double thMin =-1e9;
  double thMax = 1e9;
  for (int ilA=kNLrActive-1;ilA>kALrSPD2;ilA--) {
    if (!IsObligatoryLayer(ilA)) continue;
    int ilr=fgkActiveLrITS[ilA];
    double r   = fgkRLayITS[ilr] - fgkRSpanITS[ilr];
    double dz = fgkZSpanITS[ilr]+ztolerEdge+fDThetaTrackletSc*r;
    double ri  = 1./r;
    double tmin= (-dz-fSPDVertex->GetZ())*ri;
    double tmax= ( dz-fSPDVertex->GetZ())*ri;
    if (tmin>thMin) thMin = tmin;
    if (tmax<thMax) thMax = tmax;
  }
  double r = fgkRLayITS[kLrSPD2] + fgkRSpanITS[kLrSPD2];
  fZSPD2CutMin = fSPDVertex->GetZ()+thMin*r; // min Z of SPD2 in tracklet to consider tracking
  fZSPD2CutMax = fSPDVertex->GetZ()+thMax*r; // max Z of SPD2 in tracklet to consider tracking
  //
}

//______________________________________________
Bool_t AliXXXITSTracker::CreateTrack(AliXXXITSTracker::ITStrack_t& track, 
				     AliXXXITSTracker::SPDtracklet_t& trlet)
{
  // create track seed from tracklet
  // init track
  AliXXXLayer::ClsInfo_t *cli1=fLayers[kALrSPD1]->GetClusterInfo(trlet.id1);
  AliXXXLayer::ClsInfo_t *cli2=fLayers[kALrSPD2]->GetClusterInfo(trlet.id2);
  AliITSRecPoint *cl1=fLayers[kALrSPD1]->GetClusterUnSorted(cli1->index);
  AliITSRecPoint *cl2=fLayers[kALrSPD2]->GetClusterUnSorted(cli2->index);
  int det1 = cl1->GetVolumeId()-fLayers[kALrSPD1]->GetVIDOffset();
  int det2 = cl2->GetVolumeId()-fLayers[kALrSPD2]->GetVIDOffset();
  AliXXXLayer::ITSDetInfo_t& detInfo1 = fLayers[kALrSPD1]->GetDetInfo(det1);
  AliXXXLayer::ITSDetInfo_t& detInfo2 = fLayers[kALrSPD2]->GetDetInfo(det2);
  //
  // crude momentun estimate
  float dx=cli1->x-cli2->x,dy=cli1->y-cli2->y,d=TMath::Sqrt(dx*dx+dy*dy);
  float qptInv = fBz ? 2*TMath::Sin(cli2->phi-cli1->phi)/d/fBz/kB2C : 0; // positive particle goes anticlockwise in B+
  //
  // we initialize the seed in the tracking frame of 1st detector
  float xv= fSPDVertex->GetX()*detInfo1.cosTF + fSPDVertex->GetY()*detInfo1.sinTF;
  float yv=-fSPDVertex->GetX()*detInfo1.sinTF + fSPDVertex->GetY()*detInfo1.cosTF;
  float zv= fSPDVertex->GetZ();
  float par[5] = {yv, zv, (float)TMath::Sin(cli1->phi-detInfo1.phiTF), (cli1->z-zv)/cli1->r, qptInv};
  double covVtx[6]; 
  fSPDVertex->GetCovarianceMatrix(covVtx);
  float cov[15] = {float(covVtx[0]+covVtx[2] + fAddErr2YspdVtx),
		   0, float(covVtx[5] + fAddErr2ZspdVtx),
		   0,0,1,
		   0,0,0,1,
		   0,0,0,0,100*100};
  AliExternalTrackParam& param = track.paramOut;
  param.Set(xv, detInfo1.phiTF, par, cov);
  track.chi2 = 0;   // chi2 at 1st two point is 0
  // go to 1st layer, ignoring the MS (errors are anyway not defined)
  if (!param.PropagateTo(detInfo1.xTF+cl1->GetX(), fBz)) return kFALSE;
  Double_t cpar0[2]={ cl1->GetY(), cl1->GetZ()};
  Double_t ccov0[3]={ cl1->GetSigmaY2() + GetClSystYErr2(kALrSPD1), 0., cl1->GetSigmaZ2() + GetClSystZErr2(kALrSPD1)};
  if (!param.Update(cpar0,ccov0)) return kFALSE;
  if (!param.CorrectForMeanMaterial(fgkX2X0ITS[kLrSPD1],-fgkRhoLITS[kLrSPD1],fgkDefMass)) return kFALSE;
  // go to 2nd layer
  if (!param.Rotate(detInfo2.phiTF)) return kFALSE;
  if (!param.PropagateTo(detInfo2.xTF+cl2->GetX(), fBz)) return kFALSE;
  Double_t cpar1[2]={ cl2->GetY(), cl2->GetZ()};
  Double_t ccov1[3]={ cl2->GetSigmaY2() + GetClSystYErr2(kALrSPD2), 0., cl2->GetSigmaZ2() + GetClSystZErr2(kALrSPD2)};
  track.chi2 += param.GetPredictedChi2(cpar1,ccov1);
  if (!param.Update(cpar1,ccov1)) return kFALSE;
  //
  track.clID[0] = trlet.id1;
  track.clID[1] = trlet.id2;
  track.clID[2] = track.clID[3] = track.clID[4] = track.clID[5] = -1;
  track.ncl = 2;
  track.nmiss=0;
  track.label = fgkDummyLabel;
  //
  param.SetBit(kInvalidBit); // flag that track is not yer refitted outward 
  track.paramOut.SetBit(kInvalidBit); // flag that track was not refitter inward
  return kTRUE;
}

//______________________________________________
Bool_t AliXXXITSTracker::CrossPassiveLayer(AliExternalTrackParam& param, Int_t lrID)
{
  // cross the layer, applying mat. corrections
  double xStart=param.GetX();
  double xToGo = GetXatLabRLin(param,fgkRLayITS[lrID]);
  if (xToGo<0 || !param.PropagateTo(xToGo,fBz)) return kFALSE;
  double x2x0=fgkX2X0ITS[lrID],xrho=fgkRhoLITS[lrID];
  if (xStart<xToGo) xrho = -xrho; // inward propagation
  return param.CorrectForMeanMaterial(x2x0,xrho,fgkDefMass,kFALSE);
//
}

//______________________________________________
Bool_t AliXXXITSTracker::FollowToLayer(AliXXXITSTracker::ITStrack_t& track, Int_t lrID)
{
  // take track to given layer, searching hits if needed and applying mat. corrections
  int lrIDA = fgkLr2Active[lrID]; // active layer ID
  if (lrIDA<0 || fSkipLayer[lrIDA]) return CrossPassiveLayer(track.paramOut,lrID);
  //
  AliExternalTrackParam trCopy(track.paramOut);
  double xToGo = GetXatLabRLin(trCopy,fgkRLayITS[lrID]); // aproximate X at lrID
  if (!trCopy.PropagateTo(xToGo,fBz)) return kFALSE;
  double xyz[3];
  trCopy.GetXYZ(xyz);
  double phi=TMath::ATan2(xyz[1],xyz[0]),z=trCopy.GetZ();
  // we need track errors in the plane nearly tangential to crossing point
  if (!trCopy.Rotate(phi)) return kFALSE;
  double dphi = TMath::Sqrt(trCopy.GetSigmaY2()*fNSigma2[lrIDA]+fYToler2[lrIDA])/fgkRLayITS[lrID];
  double dz   = TMath::Sqrt(trCopy.GetSigmaZ2()*fNSigma2[lrIDA]+fZToler2[lrIDA]);
  AliXXXLayer* lrA = fLayers[lrIDA];
  int nCl = lrA->SelectClusters(z-dz,z+dz,phi-dphi,phi+dphi);
  Bool_t updDone = kFALSE;
  //
  if (nCl) {
    int icl,iclBest=-1;
    double chi2Best = fMaxChi2Tr2Cl;
    AliITSRecPoint* bestCl = 0;
    AliExternalTrackParam bestTr;
    //
    while ( (icl=lrA->GetNextClusterInfoID())!=-1) {
      AliXXXLayer::ClsInfo_t *cli = lrA->GetClusterInfo(icl);
      AliITSRecPoint *cl=lrA->GetClusterUnSorted(cli->index);
      int detId = cl->GetVolumeId()-lrA->GetVIDOffset();
      AliXXXLayer::ITSDetInfo_t& detInfo = lrA->GetDetInfo(detId);
      trCopy = track.paramOut;
      if (!trCopy.Propagate(detInfo.phiTF, detInfo.xTF+cl->GetX(), fBz)) continue;
      double cpar[2]={ cl->GetY(), cl->GetZ()};
      double ccov[3]={ cl->GetSigmaY2() + GetClSystYErr2(lrIDA) , 0., cl->GetSigmaZ2() + GetClSystZErr2(lrIDA)};
      double chi2cl = trCopy.GetPredictedChi2(cpar,ccov);
      if (chi2cl>fMaxChi2Tr2Cl) continue;
      //    SaveCandidate(lrIDA,trCopy,chi2cl,icl);  // RS: do we need this?
      if (chi2cl>chi2Best) continue;
      chi2Best = chi2cl;
      iclBest = icl;
      bestCl = cl;
      bestTr = trCopy;
      if (nCl==1) { // in absence of competitors, do the fit on spot
	if (!bestTr.Update(cpar,ccov)) return kFALSE;
	updDone = kTRUE;
      }
    }
    printf("Lr%d -> %f\n",lrIDA,chi2Best);
    //
    if (bestCl) {
      if (!updDone) {
	double cpar[2]={ bestCl->GetY(), bestCl->GetZ()};
	double ccov[3]={ bestCl->GetSigmaY2(), 0., bestCl->GetSigmaZ2()}; // RS: add syst errors    
	if (!bestTr.Update(cpar,ccov)) return kFALSE;
	updDone = kTRUE;
      }
      track.paramOut = bestTr;
      track.clID[lrIDA] = iclBest;      
      track.ncl++;
      track.chi2 += chi2Best;      
    }
  }
  //
  if (!updDone) {
    if (++track.nmiss > fMaxMissedLayers)  return kFALSE;
    track.paramOut = trCopy;
    track.chi2 += fMissChi2Penalty;
  }
  //
  if (track.chi2 > GetChi2TotCut(track.ncl+1)) return kFALSE;
  //
  return track.paramOut.CorrectForMeanMaterial(fgkX2X0ITS[lrID],-fgkRhoLITS[lrID],fgkDefMass,kFALSE);
  //
}

//______________________________________________
void AliXXXITSTracker::CookLabel(AliXXXITSTracker::ITStrack_t track)
{
  // cook mc label for the track
  track.label = fgkDummyLabel;
  if (!track.ncl) return;
  const int kMaxLbPerCl = 3;
  int lbID[kNLrActive*6],lbStat[kNLrActive*6];
  Int_t nLab=0;
  for (int i=kNLrActive;i--;) {
    int clid = track.clID[i];
    if (clid<0) continue;
    AliITSRecPoint* cl = fLayers[i]->GetClusterSorted(clid);
    for (int imc=0;imc<kMaxLbPerCl;imc++) { // labels within single cluster
      int trLb = cl->GetLabel(imc);
      if (trLb<0) break;
      // search this mc track in already accounted ones
      int iLab;
      for (iLab=0;iLab<nLab;iLab++) if (lbID[iLab]==trLb) break;
      if (iLab<nLab) lbStat[iLab]++;
      else {
	lbID[nLab] = trLb;
	lbStat[nLab++] = 1;
      }
    } // loop over given cluster's labels
  } // loop over all clusters
  //
  if (nLab) {
    int maxLab=0;
    for (int ilb=nLab;ilb--;) if (lbStat[maxLab]<lbStat[ilb]) maxLab=ilb;
    track.label = maxLab==track.ncl ? lbID[maxLab] : -lbID[maxLab];
  }
  //
}

//______________________________________________
Double_t AliXXXITSTracker::GetXatLabRLin(AliExternalTrackParam& track, double r)
{
  // X of track circle intersection in current tracking frame, neglecting the curvature
  // Solution of equation (x+d)^2+(y+b*d)^2 - r^2, where x,y are current coordinates of 
  // track and d=X-x0. b = tg(phi)
  //double sn=tr.GetSnp();
  double sn=track.GetSnp();
  if (TMath::Abs(sn)>kAlmost1) return -999;
  double x=track.GetX(), y=track.GetY();
  double cs2=(1.-sn)*(1.+sn), tg=sn/TMath::Sqrt(cs2);
  double t0=x+tg*y, t1=x*x+y*y-r*r, det=t0*t0-t1/cs2;
  if (det<0) return -999; // does not touch circle
  det = TMath::Sqrt(det);
  return x+(det-t0)*cs2;
  //
}

//______________________________________________
Int_t AliXXXITSTracker::GetTrackletMCTruth(AliXXXITSTracker::SPDtracklet_t& trlet) const
{
  int status = 0;
  AliXXXLayer::ClsInfo_t *cli1=fLayers[kALrSPD1]->GetClusterInfo(trlet.id1);
  AliXXXLayer::ClsInfo_t *cli2=fLayers[kALrSPD2]->GetClusterInfo(trlet.id2);
  AliITSRecPoint *cl1=fLayers[kALrSPD1]->GetClusterUnSorted(cli1->index);
  AliITSRecPoint *cl2=fLayers[kALrSPD2]->GetClusterUnSorted(cli2->index);
  //
  int lab = -1;
  //
  for (int i=0;i<3;i++) {
    int lb1 = cl1->GetLabel(i); 
    if (lb1<0) continue;
    for (int j=0;j<3;j++) {
      int lb2 = cl2->GetLabel(i); 
      if (lb2<0) break;
      if (lb1==lb2) {lab = lb1; break;}
    }
    if (lab>=0) break;
  }
  if (lab<0) return 0;
  //
  for (int ila=kALrSDD1;ila<kNLrActive;ila++) {
    for (int icl=fNClusters[ila];icl--;) {
      AliITSRecPoint *cl=fLayers[ila]->GetClusterUnSorted(icl);
      for (int i=0;i<3;i++) {
	if (cl->GetLabel(i)<0) break;
	if (cl->GetLabel(i)==lab) {status |= 0x1<<ila; break;}
      }
      if (status & (0x1<<ila)) break;
    }
  }
  return status;
}

//______________________________________________
Bool_t AliXXXITSTracker::RefitInward(int itr)
{
  // refit track inward with material correction
  ITStrack_t &track = fTracks[itr];
  AliExternalTrackParam &trout = track.paramOut;
  if (trout.TestBit(kInvalidBit)) return kFALSE;
  AliExternalTrackParam &trin  = track.paramInw;
  trin = trout;
  int ilA = kNLrActive;
  for (;ilA--;) {                    // find outermost layer with cluster
    if (track.clID[ilA]<0) continue;
    break;
  }
  int ilStart = fgkActiveLrITS[ilA]; // corresponding total lr id 
  AliXXXLayer* lrA = fLayers[ilA];
  AliITSRecPoint *cl=lrA->GetClusterSorted(track.clID[ilA]);
  AliXXXLayer::ITSDetInfo_t& detInfo = lrA->GetDetInfo(cl->GetVolumeId()-lrA->GetVIDOffset());
  if (!trin.RotateParamOnly(detInfo.phiTF)) return kFALSE;
  if (!trin.PropagateParamOnlyTo(detInfo.xTF+cl->GetX(), fBz)) return kFALSE;
 // init with outer cluster y,z and slopes, q/pt of outward track
  double par[5] = {cl->GetY(), cl->GetZ(), trin.GetSnp(), trin.GetTgl(), trin.GetSigned1Pt()}; 
  double cov[15] = {cl->GetSigmaY2() + GetClSystYErr2(kALrSPD1), 
		   0., cl->GetSigmaZ2() + GetClSystZErr2(kALrSPD1),
		   0,0,1,
		   0,0,0,1,
		   0,0,0,0,100*100};
  trin.Set(double(detInfo.xTF+cl->GetX()),double(detInfo.phiTF), par, cov);
  // !!! no material correction is needed: errors are not defined yer
  //
  for (int ilr=ilStart;ilr--;) {
    //
    if ( (ilA=fgkLr2Active[ilr])<0 || track.clID[ilA]<0) { // either passive layer or no cluster
      if (CrossPassiveLayer(trin,ilr)) continue;
      else return kFALSE;
    }
    // there is a cluster, need to update
    lrA = fLayers[ilA];
    cl = lrA->GetClusterSorted(track.clID[ilA]);
    detInfo = lrA->GetDetInfo(cl->GetVolumeId()-lrA->GetVIDOffset());
    if (!trin.Propagate(detInfo.phiTF, detInfo.xTF+cl->GetX(), fBz)) return kFALSE;
    double cpar[2]={ cl->GetY(), cl->GetZ()};
    double ccov[3]={ cl->GetSigmaY2() + GetClSystYErr2(ilA) , 0., cl->GetSigmaZ2() + GetClSystZErr2(ilA)};
    if (!trin.Update(cpar,ccov)) return kFALSE;
    //
    // correct for layer materials
    if (!trin.CorrectForMeanMaterial(fgkX2X0ITS[ilr],fgkRhoLITS[ilr],fgkDefMass,kFALSE)) return kFALSE;
    //
  }
  //
  // now go to PCA to vertex
  //double dca[2],dcaCov[3];
  if (!trin.PropagateToDCA(fSPDVertex,fBz,fgkRLayITS[kLrBeamPime])) return kFALSE; //,dca,dcaCov);
  //
  trin.ResetBit(kInvalidBit); // flag that inward fit succeeded
  return kTRUE;
  //
}

//______________________________________________
void AliXXXITSTracker::RefitInward()
{
  // refit tracks inward with material correction
  for (int itr=fNTracks;itr--;) {
    if (!RefitInward(itr)) {
      printf("RefitInward failed for track %d\n",itr);
      PrintTrack(fTracks[itr]);
    }
  }
  //
}

/*
//______________________________________________
Bool_t AliXXXITSTracker::FitTrackVertex()
{
  // refit tracks inward with material correction
  for (int itr=fNTracks;itr--;) {
    //
    AliExternalTrackParam& trc = fTracks[itr].paramInw;
    if (trc.TestBit(kInvalidBit)) continue; // the track is invalidated, skip
    //
    if (!RefitInward(itr)) {
      printf("RefitInward failed for track %d\n",itr);
      PrintTrack(fTracks[itr]);
    }
  }
  //
}
*/
 
#ifdef _TIMING_
//______________________________________________
void AliXXXITSTracker::PrintTiming()
{
  // print timing info
  for (int i=0;i<kNSW;i++) {printf("%s\t",fgkSWNames[i]); fSW[i].Print();}
}
#endif
