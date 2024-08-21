#include <iostream>
#include <math.h>
#include "eNCrossSection.hh"

eNCrossSection::eNCrossSection()
{
  // Set defaults
  myModel=kelly;
  myMethod=cc1;
}

eNCrossSection::eNCrossSection(csMethod thisMeth, ffModel thisMod)
{
  std::cerr << "eNCrossSection: you have selected configuration: " << thisMeth << " " << thisMod <<"\n";
  myModel=thisMod;
  myMethod=thisMeth;
}

eNCrossSection::~eNCrossSection(){}

double eNCrossSection::sigma_eN(double Ebeam,TVector3 k, TVector3 p, bool isProton)
{
  switch (myMethod)
    {
    case onshell:
      return sigma_onShell_by_Etheta(Ebeam,k,isProton);
    case cc1:
      return sigmacc1(Ebeam,k,p,isProton);
    case cc2:
      return sigmacc2(Ebeam,k,p,isProton);
    default:
      {
	std::cerr << "Invalid cross section method! Double check and fix!\n";
        exit(-1);
      }
    }
  return 0;
}

double eNCrossSection::sigma_CC(double Ebeam,TVector3 k, TVector3 p, bool isProton)
{
  TVector3 q = TVector3(0.,0.,Ebeam) - k;
  double omega = Ebeam - k.Mag();
  double QSq = q.Mag2() - sq(omega);

  return sq(QSq/sq(mW))*sq(Vud)*sigma_eN(Ebeam, k, p, isProton);
  
}

// Because DeForest uses the opposite Lorentz convention
double dot4(double x0, TVector3 x, double y0, TVector3 y)
{
  return ((x0*y0)-(x*y));
}

double eNCrossSection::sigmaccn(double Ebeam, TVector3 k, TVector3 p, bool isProton, int n)

{
  TVector3 q = TVector3(0.,0.,Ebeam) - k;
  TVector3 pM = p-q;
  
  double omega = Ebeam - k.Mag();
  double QSq = q.Mag2() - sq(omega);
  double E = sqrt(p.Mag2() + sq(mN));
  double Ebar = sqrt(pM.Mag2() + sq(mN));
  double omegabar = E-Ebar;
  double QSqbar = q.Mag2() - sq(omegabar);

  // Calculate form factors
  double GE = (isProton)? GEp(QSq) : GEn(QSq);
  double GM = (isProton)? GMp(QSq) : GMn(QSq);

  double F1 = (GE + GM * QSq/(4.*sq(mN)))/(1. + QSq/(4.*sq(mN)));
  double kF2 = (GM - GE)/(1. + QSq/(4.*sq(mN)));

  double wC;
  double wT;
  double wS;
  double wI;
  
  if (n==1)
    {
      wC = (sq(E+Ebar)*(sq(F1) + QSqbar/(4.*mN*mN) * sq(kF2)) - q.Mag2()*sq(F1 + kF2))/(4.*E*Ebar);
      wT = QSqbar*sq(F1 + kF2)/(2.*Ebar*E);
      wS = p.Mag2() * sq(sin(p.Angle(q))) * (sq(F1) + QSqbar/(4.*mN*mN) * sq(kF2))/(E*Ebar);
      wI = -p.Mag()*sin(p.Angle(q))*(Ebar + E)*(sq(F1) + QSqbar/(4.*mN*mN) * sq(kF2))/(E*Ebar);
    }
  else if (n==2)
    {  
      double pbarp = dot4(Ebar,pM,E,p);
      double pbarq = dot4(Ebar,pM,omega,q);
      double pq = dot4(E,p,omega,q);
      double qbarq = dot4(omegabar,q,omega,q);
      double sumq = dot4((Ebar+E),(pM+p),omega,q);

      wC = (E*Ebar
		   + 0.5 * (pbarp + sq(mN)) * sq(F1)
		   - 0.5 * q.Mag2() * F1 * kF2
		   - ((pbarq*E + pq*Ebar)*omega
		      - Ebar * E * QSq
		      + pbarq * pq
		      - 0.5 * (pbarp - sq(mN)) * q.Mag2())
		   * sq(kF2)/(4*sq(mN))
		   )/(E*Ebar);
      wT = (-(pbarp + sq(mN)) * sq(F1)
		   + qbarq * F1 * kF2
		   + (2*pbarq*pq
		      - (pbarp - sq(mN))*QSq)
		   * sq(kF2)/(4*sq(mN))
		   )/(Ebar*E);
      wS = p.Mag2() * sq(sin(p.Angle(q))) * (sq(F1)
						    + QSq/(4.*mN*mN) * sq(kF2))/(E*Ebar);
      wI = p.Mag()*sin(p.Angle(q))*(-(Ebar + E) * sq(F1)
					   + (sumq * omega
					      - (Ebar + E) * QSq)
					   * sq(kF2)/(4*sq(mN))
					   )/(E*Ebar);
    }
  else
    {
      std::cerr << "Invalid cross section designation. Check and fix. Exiting\n\n\n";
    }
      
  double sigmaMott = nbGeVSq * 4. * sq(alpha) * k.Mag2() * sq(cos(k.Theta()/2.)) / sq(QSq);

  double phi = q.Cross(k).Angle( q.Cross(p) );
  return sigmaMott * ( sq(QSq)/q.Mag2() * wC +
                       (QSq/(2.*q.Mag2()) + sq(tan(k.Theta()/2.))) * wT +
                       QSq/q.Mag2() * sqrt(QSq/q.Mag2() + sq(tan(k.Theta()/2.))) * wI * cos(phi) +
                       (QSq/q.Mag2() * sq(cos(phi)) + sq(tan(k.Theta()/2.))) * wS
                       );
}

double eNCrossSection::sigmacc1(double Ebeam, TVector3 k, TVector3 p, bool isProton)
{
  return sigmaccn(Ebeam, k, p, isProton, 1);
}

double eNCrossSection::sigmacc2(double Ebeam, TVector3 k, TVector3 p, bool isProton)
{
  return sigmaccn(Ebeam, k, p, isProton, 2);
}

double eNCrossSection::sigma_onShell_by_Etheta(double Ebeam, TVector3 k, bool isProton)
{
  double theta=k.Theta();
  double E3 = Ebeam * mN/ (mN + Ebeam*(1.-k.CosTheta()));
  double QSq = 2. * Ebeam * E3 * (1.-k.CosTheta());
  double tau = QSq/(4.*mN*mN);
  double GE = isProton ? GEp(QSq) : GEn(QSq);
  double GM = isProton ? GMp(QSq) : GMn(QSq);
  double epsilon = 1./(1.+2.*(1.+tau)*sq(tan(theta/2.)));

  double sigmaMott = nbGeVSq * sq(2.*alpha*E3 * cos(theta/2.)/QSq) * (E3/Ebeam);

  double sigmaRosenbluth = sigmaMott * (sq(GE) + tau/epsilon * sq(GM))/(1. + tau);
  return sigmaRosenbluth * Ebeam / (E3 * (2.*tau + 1.));
}

double eNCrossSection::GEp(double QSq)
{
  switch (myModel)
    {
    case dipole:
      return Gdipole(QSq);
    case kelly:
      return Gkelly(QSq,-0.24,10.98,12.82,21.97);
    default:
      std::cerr << "Error in GEp: invalid form factor model!\n";
      exit(-1);
    }
  return 0.;
}

double eNCrossSection::GEn(double QSq) // This will use the Galster parameterization
{
  double tau = QSq/(4.*mN*mN);
  return 1.70 * tau / (1. + 3.3 * tau) * Gdipole(QSq); // params from Kelly paper
  //return mu_n * tau / (1. + 5.6 * tau) * Gdipole(QSq); // the original Galster numbers
}

double eNCrossSection::GMp(double QSq)
{
  switch (myModel)
    {
    case dipole:
      return mu_p * Gdipole(QSq);
    case kelly:
      return mu_p * Gkelly(QSq,0.12,10.97,18.86,6.55);
    default:
      std::cerr << "Error in GMp: invalid form factor model!\n";
      exit(-1);
    }
  return 0.;
}

double eNCrossSection::GMn(double QSq)
{
  switch (myModel)
    {
    case dipole:
      return mu_n * Gdipole(QSq);
    case kelly:
      return mu_n * Gkelly(QSq,2.33,14.72,24.20,84.1);
    default:
      std::cerr << "Error in GMn: invalid form factor model!\n";
      exit(-1);
    }
  return 0.;
}

double eNCrossSection::Gdipole(double QSq){ return 1. / sq(1 + QSq/0.71); }

double eNCrossSection::Gkelly(double QSq,double a1, double b1, double b2, double b3)
{
  double tau = QSq/(4.*mN*mN);
  double denom = 1. + b1*tau + b2*tau*tau + b3*tau*tau*tau;
  double numer = 1. + a1*tau;
  return numer/denom;
}

