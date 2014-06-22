/*

Assume helix passing through the vertex:
x = r*[cos(phi0) - cos(phi0+t)]
y = r*[sin(phi0) - sin(phi0+t)]
z = zv - r*t*tgl

*/

using namespace TMath;


namespace AliAlgAux {
  const double kAlmostZeroD = 1e-15;
  const double kAlmostZeroF = 1e-11;
  const double kAlmostOneD = 1.-kAlmostZeroD;
  const double kAlmostOneF = 1.-kAlmostZeroF;
  //
  void   PrintBits(ULong64_t patt, Int_t maxBits);
  Bool_t SmallerAbs(double d, double tolD)    {return Abs(d)<tolD;}
  Bool_t SmallerAbs(float  f, double tolF)    {return Abs(f)<tolF;}
  Bool_t Smaller(double d, double tolD)       {return d<tolD;}
  Bool_t Smaller(float  f, double tolF)       {return f<tolF;}
  Bool_t IsZeroAbs(double d)                  {return SmallerAbs(d,kAlmostZeroD);}
  Bool_t IsZeroAbs(float  f)                  {return SmallerAbs(f,kAlmostZeroF);}
  Bool_t IsZeroPos(double d)                  {return Smaller(d,kAlmostZeroD);}
  Bool_t IsZeroPos(float  f)                  {return Smaller(f,kAlmostZeroF);}
}


Bool_t calcToler(double minPt, double bz, double rDest)
{
  const kb2c = -0.299792458e-3;
  double cMax = TMath::Abs(bz*kb2c/minPt);
  //
  double cRMax = cMax*rDest; // ratio destination R to track radius
  if (cRmax>2*kAlmost1) return kFALSE;
  double cRMax2 = cRMax*cRMax;
  double tMaxCos = 1.-0.5*cRMax2; // cos of max allowed helix param at rDest
  double tMax = TMath::ACos(tMaxCos); 
  //  double tMaxSin = cRMax*TMath::Sqrt((1-0.5*cRMax)*(1+0.5*cRMax)); // sin...  
  //
  double zTol2Tgl = rDest - tMax/cMax; // difference between straight line and helix segment
}


void ProcessTracklet()
{
  
}

//______________________________________________________________________________
Bool_t InitTrackParams(trackC &track)
{
  // Set the initial guess on track kinematics for propagation using 3 points
  // Assume at least 3 points available
  int lrOcc[AliITSUAux::kMaxLayers], nCl=0;
  //
  // we will need endpoints and middle layer
  for (int i=fITS->GetNLayersActive();i--;) if (track.fPoints[i<<0x1]>-1) lrOcc[nCl++] = i;
  if (nCl<3) {
    AliError(Form("Cannot estimate momentum of tracks with %d clusters",nCl));
    cout << track << endl;
    return kFALSE;
  }
  //
  int lr0   = lrOcc[0];
  int lr1   = lrOcc[nCl/2];
  int lr2   = lrOcc[nCl-1];
  //
  const itsCluster& cl0 = fClusters[lr0][ track.fPoints[lr0<<0x1] ];
  const itsCluster& cl1 = fClusters[lr1][ track.fPoints[lr1<<0x1] ];
  const itsCluster& cl2 = fClusters[lr2][ track.fPoints[lr2<<0x1] ];
  double cv = Curvature(cl0.x,cl0.y, cl1.x,cl1.y, cl2.x,cl2.y);
  double tgl = (cl2.z-cl0.z)/TMath::Sqrt((cl2.x-cl0.x)*(cl2.x-cl0.x)+(cl2.y-cl0.y)*(cl2.y-cl0.y));
  //  double phi = TMath::ATan2((cl2.y-cl1.y),(cl2.x-cl1.x));
  //
  AliITSUClusterPix* clus = (AliITSUClusterPix*)fClustersTC[ lr0 ]->At( track.fPoints[lr0<<0x1] );
  AliITSURecoLayer* lr = fITS->GetLayerActive(lr0);
  AliITSURecoSens* sens = lr->GetSensorFromID(clus->GetVolumeId());
  double x = sens->GetXTF() + clus->GetX();
  double alp = sens->GetPhiTF();
  //  printf("Alp: %f phi: %f\n",alp,phi);
  double par[5] = {clus->GetY(),clus->GetZ(),0,tgl,cv};
  double cov[15] = {
    5*5,
    0, 5*5,
    0, 0, 0.7*0.7,
    0,0,0,0.7*0.7,
    0,0,0,0,10
  };
  track.Set(x,alp,par,cov);
  return kTRUE;
}

//____________________________________________________________________
Double_t Curvature(Double_t x1,Double_t y1,Double_t x2,Double_t y2,Double_t x3,Double_t y3)
{
  //calculates the curvature of track using 3 points 
  double dy31=y3-y1,dy21=y2-y1;
  double dx31=x3-x1,dx21=x2-x1;
  double den = dx31*dy21-dx21*dy31;
  //
  if (IsZeroAbs(den)) return 0;
  //
  Double_t a = (dy31*( (x2+x1)*dx21 + (y2+y1)*dy21 ) -
		dy21*( (x3+x1)*dx31 + (y3+y1)*dy31 ))/den;
  Double_t b = -( (x2+x1)*dx21+(y2+y1)*dy21+a*dx21 )/dy21;
  Double_t c = -(x1*x1+y1*y1+a*x1+b*y1);
  Double_t xc=-a/2.;
  //
  if((a*a+b*b-4*c)<0) return 0;
  Double_t rad = TMath::Sqrt(a*a+b*b-4*c)/2.;
  if(IsZeroPos(red)) return 0;
  
  if((x1>0 && y1>0 && x1<xc)) rad*=-1;
  if((x1<0 && y1>0 && x1<xc)) rad*=-1;
  //  if((x1<0 && y1<0 && x1<xc)) rad*=-1;
  // if((x1>0 && y1<0 && x1<xc)) rad*=-1;
  
  return 1/rad;
 
}






