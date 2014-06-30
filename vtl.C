{
  double sn=0.84, cs=TMath::Sqrt((1-sn)*(1+sn)), vX=0.14, vY=-0.13,vZ=3.14, x0i=10.23, y0i=11.2, 
    z0i=14.3, tgp=1.5, tgl=-0.55, syy=1.2, szz=0.35, syz=-0.2;

  double tmpSP = sn*tgp;
  double tmpCP = cs*tgp;
  double tmpSC = sn+tmpCP;
  double tmpCS =-cs+tmpSP;
  double tmpCL = cs*tgl;
  double tmpSL = sn*tgl;
  double tmpYXP = y0i-tgp*x0i;
  double tmpZXL = z0i-tgl*x0i;
  //
  double tmpCLzz = tmpCL*szz;
  double tmpSLzz = tmpSL*szz;
  double tmpSCyz = tmpSC*syz;
  double tmpCSyz = tmpCS*syz;
  double tmpCSyy = tmpCS*syy;
  double tmpSCyy = tmpSC*syy;
  double tmpSLyz = tmpSL*syz;
  double tmpCLyz = tmpCL*syz;
  //
  double cxx = tmpCL*(tmpCLzz+tmpSCyz+tmpSCyz)+tmpSC*tmpSCyy;
  double cxy = tmpCL*(tmpSLzz+tmpCSyz)+tmpSL*tmpSCyz+tmpSC*tmpCSyy;
  double cxz = -sn*syz-tmpCLzz-tmpCP*syz;
  double cx0 = -(tmpCLyz+tmpSCyy)*tmpYXP-(tmpCLzz+tmpSCyz)*tmpZXL;
  //
  //double cyx
  double cyy = tmpSL*(tmpSLzz+tmpCSyz+tmpCSyz)+tmpCS*tmpCSyy;
  double cyz = -(tmpCSyz+tmpSLzz);
  double cy0 = -tmpYXP*(tmpCSyy+tmpSLyz)-tmpZXL*(tmpCSyz+tmpSLzz);
  //
  //double czx
  //double czy
  double czz = szz;
  double cz0 = tmpZXL*szz+tmpYXP*syz;
  printf("%.5f %.5f %.5f  %.5f %.5f %.5f %.5f  %.5f %.5f\n",cxx,cxy,cxz,cx0,cyy,cyz,cy0,czz,cz0);
 
}


