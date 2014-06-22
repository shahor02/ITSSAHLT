//______________________________________________
Double_t GetXatLabRLin(double sn, double x, double y, double r)
{
  // X of track circle intersection in current tracking frame, neglecting the curvature
  // Solution of equation (x+d)^2+(y+b*d)^2 - r^2, where x,y are current coordinates of 
  // track and d=X-x0. b = tg(phi)
  //double sn=tr.GetSnp();
  if (TMath::Abs(sn)>1.-1e-5) return -999;
  //double x=tr.GetX(), y=tr.GetY();
  double cs2=(1.-sn)*(1.+sn), tg=sn/TMath::Sqrt(cs2);
  double t0=x+tg*y, t1=x*x+y*y-r*r, det=t0*t0-t1/cs2;
  if (det<0) return -999; // does not touch circle
  det = TMath::Sqrt(det);
  return x+(t0+det)*cs2;
  

}
