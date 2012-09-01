#ifndef commonTools_C
#define commonTools_C

// Adds overflow of a TH1F in the final bin of the histogram
TH1F * overflow(TH1F *h)
{
  int  N= h->GetNbinsX();
  h->SetBinContent(N,h->GetBinContent(N)+h->GetBinContent(N+1));
  return h;
}

// Creates a TLorentzVector from pt, eta, phi and mass
TLorentzVector convert(const double & pt, const double & eta, const double & phi, const double & mass)
{
  double px = pt*cos(phi);
  double py = pt*sin(phi);
  double tmp = 2*atan(exp(-eta));
  double pz = pt*cos(tmp)/sin(tmp);
  double E  = sqrt(px*px+py*py+pz*pz+mass*mass);

  return TLorentzVector(px,py,pz,E);
}

#endif
