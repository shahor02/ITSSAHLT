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
  1.48e-01, 2.48e-01,2.57e-01, 1.34e-01, 3.34e-01,3.50e-01, 2.22e-01, 2.38e-01,2.25e-01};

const Float_t AliXXXITSTracker::fgkX2X0ITS[AliXXXITSTracker::kMaxLrITS] = {
  0.5e-2, 1.e-2,1.e-2,  1.e-2,  1.e-2,1.e-2,  1.e-2,  1.e-2,1e-2};


const Int_t    AliXXXITSTracker::fgkLrDefBins[AliXXXITSTracker::kNLrActive][2] = // n bins in z, phi
  { {20,20}, {20,20}, {20,20}, {20,20}, {20,20}, {20,20} };

const Float_t AliXXXITSTracker::fgkDefMass = 0.14;
const Int_t   AliXXXITSTracker::fgkDummyLabel = -3141593;

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
  ,fMaxChi2Tr2Cl(15)
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
    fNSigma2[i] = 4*4;
    fYToler2[i] = 0.1*0.1;
    fZToler2[i] = 0.1*0.1;
    fChi2TotCut[i] = 0;
  }  
  fChi2TotCut[1] = 10; // 2 cl+vtx -> NDF=1
  fChi2TotCut[2] = 15; 
  fChi2TotCut[3] = 20; 
  fChi2TotCut[4] = 24; 
  fChi2TotCut[5] = 26; 

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
  AliXXXLayer &spdL1 = *fLayers[0];
  AliXXXLayer &spdL2 = *fLayers[1];
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
  AliXXXLayer &spdL1 = *fLayers[0];
  AliXXXLayer &spdL2 = *fLayers[1];
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
  for (int ila=2;ila<kNLrActive;ila++) {
    if (fSkipLayer[ila]) continue;
    fLayers[ila]->SortClusters(fSPDVertex);
    fNClusters[ila] = fLayers[ila]->GetNClusters();
  }
  //
  fTracks.resize(nTrk);
  fNTracks = 0;
  //
  for (int itr=0;itr<nTrk;itr++) {
    ITStrack_t &track = fTracks[fNTracks];
    if (!CreateTrack(track, fTracklets[itr])) continue;
    for (int lrID=kLrShield1;lrID<kMaxLrITS;lrID++) if (!FollowToLayer(track,lrID)) break;
    CookLabel(track);
    fNTracks++;
    //
  }  
}

//______________________________________________
void AliXXXITSTracker::PrintTracklets() const
{
  // print traklets info
  int ntr = fTracklets.size();
  printf("Nspd1: %4d Nspd2: %4d, Ntracklets: %d\n",fNClusters[0],fNClusters[1],ntr);
  for (int itr=0;itr<ntr;itr++) {
    const SPDtracklet_t* trk = &fTracklets[itr];
    AliITSRecPoint* cl1 = fLayers[0]->GetClusterSorted(trk->id1);
    AliITSRecPoint* cl2 = fLayers[1]->GetClusterSorted(trk->id2);
    AliXXXLayer::ClsInfo_t* cli0 = fLayers[0]->GetClusterInfo(trk->id1);
    printf("#%3d Phi:%+.3f Tht:%+.3f Dphi:%+.3f Dtht:%+.3f Chi2:%.3f | Lbl:",
	   itr,cli0->phi,TMath::ATan2(cli0->r,cli0->z-fSPDVertex->GetZ()), trk->dphi,trk->dtht,trk->chi2);
    int lb = 0;
    for (int i=0;i<3;i++) if ( (lb=cl1->GetLabel(i))>=0 ) printf(" %5d",lb); printf("|");
    for (int i=0;i<3;i++) if ( (lb=cl2->GetLabel(i))>=0 ) printf(" %5d",lb); 
    printf("\n");
  }
}

//______________________________________________
void AliXXXITSTracker::CalcAuxTracking()
{
  // precalculate auxilarry variables for tracking
  //
  // largest track curvature to search
  fCurvMax = TMath::Abs(fBz*kB2C/fMinPt);
  //
  for (int ilr=0;ilr<kMaxLrITS;ilr++) {
    double sgms2 = 0.014*0.014*fgkX2X0ITS[ilr];
    for (int jlr=ilr+1;jlr<kMaxLrITS;jlr++) {
      int jlA = GetActiveLayerID(jlr);
      if (jlA<0) continue;
      double dr = fgkRLayITS[jlr] - fgkRLayITS[ilr];
      fMSDist[jlA] += sgms2*dr*dr;  // expected position uncertainty squared for 1 GeV particle
    }
  }
  //
  for (int ila=0;ila<kNLrActive;ila++) {
    fMSDist[ila] = fMSDist[ila]>0 ? TMath::Sqrt(fMSDist[ila]) : 0; // tolerance in coordinate
    fMSPhi[ila]  = fMSDist[ila]/fgkRLayITS[ila]; // tolerance in az. angle
    //
    fTolPhiCrude[ila] = TMath::ASin(0.5*fgkRLayITS[ila]*fCurvMax)  // max deflection in phi
      + 0.5*fCurvMax*fgkRSpanITS[ila]            // account for the spread in cluster R
      + fMSPhi[ila]/fMinPt;                      // MS contribution
    //
    // contribution to deviation from straight line
    fTolZCrude[ila] = TMath::Abs(2*TMath::ASin(0.5*fgkRLayITS[ila]*fCurvMax)-fgkRLayITS[ila])  // arc vs segment
      + fMSDist[ila]/fMinPt;                      // MS contribution
    
  }
  //
}

//______________________________________________
Bool_t AliXXXITSTracker::CreateTrack(AliXXXITSTracker::ITStrack_t& track, 
				     AliXXXITSTracker::SPDtracklet_t& trlet)
{
  // create track seed from tracklet
  double zv = fSPDVertex->GetZ();
  // init track
  AliXXXLayer::ClsInfo_t *cli1=fLayers[0]->GetClusterInfo(trlet.id1);
  AliXXXLayer::ClsInfo_t *cli2=fLayers[1]->GetClusterInfo(trlet.id2);
  AliITSRecPoint *cl1=fLayers[0]->GetClusterUnSorted(cli1->index);
  AliITSRecPoint *cl2=fLayers[1]->GetClusterUnSorted(cli2->index);
  int det1 = cl1->GetVolumeId()-fLayers[0]->GetVIDOffset();
  int det2 = cl2->GetVolumeId()-fLayers[1]->GetVIDOffset();
  AliXXXLayer::ITSDetInfo_t& detInfo1 = fLayers[0]->GetDetInfo(det1);
  AliXXXLayer::ITSDetInfo_t& detInfo2 = fLayers[1]->GetDetInfo(det2);
  //
  // crude momentun estimate
  double dx=cli1->x-cli2->x,dy=cli1->y-cli2->y,d=TMath::Sqrt(dx*dx+dy*dy);
  double qptInv = fBz ? 2*TMath::Sin(cli2->phi-cli1->phi)/d/fBz/kB2C : 0; // positive particle goes anticlockwise in B+
  //
  // we initialize the seed in the tracking frame of 1st detector
  float xv= fSPDVertex->GetX()*detInfo1.cosTF + fSPDVertex->GetY()*detInfo1.sinTF;
  float yv=-fSPDVertex->GetX()*detInfo1.sinTF + fSPDVertex->GetY()*detInfo1.cosTF;
  float par[5] = {yv, zv, TMath::Sin(cli1->phi-detInfo1.phiTF), (cli1->z-zv)/cli1->r, qptInv};
  Double_t covVtx[6]; 
  fSPDVertex->GetCovarianceMatrix(covVtx);
  float cov[15] = {covVtx[0]+covVtx[2],
		   0, covVtx[5],
		   0,0,1,
		   0,0,0,1,
		   0,0,0,0,100*100};
  AliExternalTrackParam& param = track.param;
  param.Set(xv, detInfo1.phiTF, par, cov);
  track.chi2 = 0;   // chi2 at 1st two point is 0
  // go to 1st layer, ignoring the MS (errors are anyway not defined)
  if (!param.PropagateTo(detInfo1.xTF+cl1->GetX(), fBz)) return kFALSE;
  Double_t cpar0[2]={ cl1->GetY(), cl1->GetZ()};
  Double_t ccov0[3]={ cl1->GetSigmaY2(), 0., cl1->GetSigmaZ2()};
  if (!param.Update(cpar0,ccov0)) return kFALSE;
  if (!param.CorrectForMeanMaterial( fgkX2X0ITS[kLrSPD1], fgkRhoLITS[kLrSPD1], 
					   fgkDefMass)) return kFALSE;
  // go to 2nd layer
  if (!param.Rotate(detInfo2.phiTF)) return kFALSE;
  if (!param.PropagateTo(detInfo2.xTF+cl2->GetX(), fBz)) return kFALSE;
  Double_t cpar1[2]={ cl2->GetY(), cl2->GetZ()};
  Double_t ccov1[3]={ cl2->GetSigmaY2(), 0., cl2->GetSigmaZ2()};
  track.chi2 += param.GetPredictedChi2(cpar1,ccov1);
  if (!param.Update(cpar1,ccov1)) return kFALSE;
  //
  track.clID[0] = trlet.id1;
  track.clID[1] = trlet.id2;
  track.clID[2] = track.clID[3] = track.clID[4] = track.clID[5] = -1;
  track.ncl = 2;
  track.label = fgkDummyLabel;
  return kTRUE;
}

//______________________________________________
Bool_t AliXXXITSTracker::CrossPassiveLayer(AliXXXITSTracker::ITStrack_t& track, Int_t lrID)
{
  // cross the layer, applying mat. corrections
  AliExternalTrackParam& param = track.param;
  double xStart=param.GetX();
  double xToGo = GetXatLabRLin(param,fgkRLayITS[lrID]);
  if (xToGo<0 || !param.PropagateTo(xToGo,fBz)) return kFALSE;
  double x2x0=fgkX2X0ITS[lrID],xrho=fgkRhoLITS[lrID];
  if (xStart>xToGo) { // inward propagation
    x2x0 = -x2x0;
    xrho = -xrho;
  }
  return param.CorrectForMeanMaterial(x2x0,xrho,fgkDefMass,kFALSE);
//
}

//______________________________________________
Bool_t AliXXXITSTracker::FollowToLayer(AliXXXITSTracker::ITStrack_t& track, Int_t lrID)
{
  // take track to given layer, searching hits if needed and applying mat. corrections
  int lrIDA = fgkLr2Active[lrID]; // active layer ID
  if (lrIDA<0 || fSkipLayer[lrIDA]) return CrossPassiveLayer(track,lrID);
  //
  AliExternalTrackParam trCopy(track.param);
  double xToGo = GetXatLabRLin(trCopy,fgkRLayITS[lrID]); // aproximate X at lrID
  if (!trCopy.PropagateTo(xToGo,fBz)) return kFALSE;
  double xyz[3];
  trCopy.GetXYZ(xyz);
  double phi=TMath::ATan2(xyz[1],xyz[2]),z=trCopy.GetZ();
  // we need track errors in the plane nearly tangential to crossing point
  if (!trCopy.Rotate(phi)) return kFALSE;
  double dphi = TMath::Sqrt(trCopy.GetSigmaY2()*fNSigma2[lrIDA]+fYToler2[lrIDA])/fgkRLayITS[lrID];
  double dz   = TMath::Sqrt(trCopy.GetSigmaZ2()*fNSigma2[lrIDA]+fZToler2[lrIDA]);
  AliXXXLayer* lrA = fLayers[lrIDA];
  int nCl = lrA->SelectClusters(z-dz,z+dz,phi-dphi,phi+dphi);
  if (!nCl) return kFALSE;
  //
  int icl,iclBest=-1;
  double chi2Best = fMaxChi2Tr2Cl;
  AliITSRecPoint* bestCl = 0;
  AliExternalTrackParam bestTr;
  Bool_t updDone = kFALSE;
  while ( (icl=lrA->GetNextClusterInfoID())!=-1) {
    AliXXXLayer::ClsInfo_t *cli = lrA->GetClusterInfo(icl);
    AliITSRecPoint *cl=lrA->GetClusterUnSorted(cli->index);
    int detId = cl->GetVolumeId()-lrA->GetVIDOffset();
    AliXXXLayer::ITSDetInfo_t& detInfo = lrA->GetDetInfo(detId);
    trCopy = track.param;
    if (!trCopy.Propagate(detInfo.phiTF, detInfo.xTF+cl->GetX(), fBz)) continue;
    double cpar[2]={ cl->GetY(), cl->GetZ()};
    double ccov[3]={ cl->GetSigmaY2(), 0., cl->GetSigmaZ2()}; // RS: add syst errors    
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
  if (!bestCl) return kFALSE;
  if (!updDone) {
    double cpar[2]={ bestCl->GetY(), bestCl->GetZ()};
    double ccov[3]={ bestCl->GetSigmaY2(), 0., bestCl->GetSigmaZ2()}; // RS: add syst errors    
    if (!bestTr.Update(cpar,ccov)) return kFALSE;
  }
  track.ncl++;
  track.chi2 += chi2Best;
  if (track.chi2 > GetChi2TotCut(track.ncl+1)) return kFALSE;
  track.param = bestTr;
  track.clID[lrIDA] = iclBest;
  if (bestTr.CorrectForMeanMaterial(fgkX2X0ITS[lrID],fgkRhoLITS[lrID],fgkDefMass,kFALSE)) return kFALSE;
  return kTRUE;
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
  return x+(t0+det)*cs2;
  //
}
