#ifndef __EN_CROSS_SECTIONS_H__
#define __EN_CROSS_SECTIONS_H__

#include "TVector3.h"
#include "gcfSRC.hh"
#include "functions.h"

enum ffModel {dipole, kelly};
enum csMethod {onshell, cc1, cc2};
inline double sq(double x){
  return x*x;
}


class eNCrossSection
{
 public:
  eNCrossSection();
  eNCrossSection(csMethod thisMeth, ffModel thisMod);
  ~eNCrossSection();
  double sigma_CC(double Ebeam, TVector3 k, TVector3 p, bool isProton);
  double sigma_eN(double Ebeam, TVector3 k, TVector3 p, bool isProton);
  double sigmaccn(double Ebeam, TVector3 k, TVector3 p, bool isProton, int n);
  double sigmacc1(double Ebeam, TVector3 k, TVector3 p, bool isProton);
  double sigmacc2(double Ebeam, TVector3 k, TVector3 p, bool isProton);
  double sigma_onShell_by_Etheta(double Ebeam, TVector3 k, bool isProton);
  double GEp(double QSq);
  double GEn(double QSq);
  double GMp(double QSq);
  double GMn(double QSq);

 private:
  ffModel myModel;
  csMethod myMethod;
  static double Gdipole(double QSq);
  static double Gkelly(double QSq,double a1, double b1, double b2, double b3);

};

#endif
