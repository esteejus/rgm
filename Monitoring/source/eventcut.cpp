#include "eventcut.h"

eventcut::eventcut(double E)
{
  Ebeam = E;
  vbeam.SetXYZ(0,0,Ebeam);
  cutInfo structcut;

  structcut.docut = true;
  structcut.min = 1;
  structcut.max = 100;
  cutmap[e_nphe]=structcut;
 
  structcut.docut = true;
  structcut.min = 14;
  structcut.max = 3000;
  cutmap[e_calv]=structcut;

  structcut.docut = true;
  structcut.min = 14;
  structcut.max = 3000;
  cutmap[e_calw]=structcut;

  structcut.docut = true;
  structcut.min = 0.18;
  structcut.max = 0.28;
  cutmap[e_SF]=structcut;

  structcut.docut = true;
  structcut.min = 1;
  structcut.max = 6.6;
  cutmap[e_mom]=structcut;


  structcut.docut = true;
  structcut.count = 2212;  
  cutmap[l_pid]=structcut;

  structcut.docut = true;
  structcut.count = 12;
  cutmap[l_scint]=structcut;

  structcut.docut = true;
  structcut.min = 0;
  structcut.max = 50;
  cutmap[l_theta]=structcut;

  structcut.docut = true;
  structcut.min = 0;
  structcut.max = 25;
  cutmap[l_thetalq]=structcut;

  structcut.docut = true;
  structcut.min = -3;
  structcut.max = 3;
  cutmap[l_chipid]=structcut;

  structcut.docut = true;
  structcut.min = 1.5;
  structcut.max = 100;
  cutmap[lsrc_Q2]=structcut;

  structcut.docut = true;
  structcut.min = 1.2;
  structcut.max = 2;
  cutmap[lsrc_xB]=structcut;

  structcut.docut = true;
  structcut.min = 0.3;
  structcut.max = 100;
  cutmap[lsrc_pmiss]=structcut;

  structcut.docut = true;
  structcut.min = 0.84;
  structcut.max = 1.04;
  cutmap[lsrc_mmiss]=structcut;

  structcut.docut = true;
  structcut.min = 0.62;
  structcut.max = 0.96;
  cutmap[lsrc_loq]=structcut;


  std::cout<<"\n\n\n\n\n\n\n\n\n Event selection class created \n\n\n\n\n\n\n";

}
eventcut::~eventcut()
{
}

bool eventcut::electroncut(const std::unique_ptr<clas12::clas12reader>& c12)
{
  auto electrons=c12->getByID(11);
  if(electrons.size()!=1){ return false;}
  if(!e_nphecut(c12)){ return false; }
  if(!e_calvcut(c12)){ return false; }
  if(!e_calwcut(c12)){ return false; }
  if(!e_SFcut(c12)){ return false; }
  if(!e_momcut(c12)){ return false; }
  return true;
}

int eventcut::leadnucleoncut(const std::unique_ptr<clas12::clas12reader>& c12)
{
  int pid = cutmap[l_pid].count;
  auto nucleons=c12->getByID(pid);
  int num_L = 0;
  int index_L = -1;
  for(int i = 0; i < nucleons.size(); i++){
    if(!l_scintcut(c12,i,pid)){ continue; }
    //if(!l_thetacut(c12,i,pid)){ continue; }
    if(!l_thetalqcut(c12,i,pid)){ continue; }
    if(!l_chipidcut(c12,i,pid)){ continue; }
    num_L++;
    index_L = i;
  }
  if(num_L != 1){return -1;}
  return index_L;
}

bool eventcut::leadSRCnucleoncut(const std::unique_ptr<clas12::clas12reader>& c12, int i)
{
  int pid = cutmap[l_pid].count;
  if(!lsrc_Q2cut(c12)){ return false; }
  if(!lsrc_xBcut(c12)){ return false; }
  if(!lsrc_pmisscut(c12,i,pid)){ return false; }
  if(!lsrc_mmisscut(c12,i,pid)){ return false; }
  if(!lsrc_loqcut(c12,i,pid)){ return false; }
  return true;
}



bool eventcut::e_nphecut(const std::unique_ptr<clas12::clas12reader>& c12)
{
  auto electrons=c12->getByID(11);
  if(electrons[0]->che(clas12::HTCC)->getNphe() <= cutmap[e_nphe].min){ return false; }  
  if(electrons[0]->che(clas12::HTCC)->getNphe() >= cutmap[e_nphe].max){ return false; }
  return true;
}
bool eventcut::e_calvcut(const std::unique_ptr<clas12::clas12reader>& c12)
{
  auto electrons=c12->getByID(11);
  
  if(electrons[0]->cal(clas12::PCAL)->getLv() <= cutmap[e_calv].min){ return false; }  
  if(electrons[0]->cal(clas12::PCAL)->getLv() >= cutmap[e_calv].max){ return false; }
  return true;
}
bool eventcut::e_calwcut(const std::unique_ptr<clas12::clas12reader>& c12)
{
  auto electrons=c12->getByID(11);
  if(electrons[0]->cal(clas12::PCAL)->getLw() <= cutmap[e_calw].min){ return false; }  
  if(electrons[0]->cal(clas12::PCAL)->getLw() >= cutmap[e_calw].max){ return false; }
  return true;
}
bool eventcut::e_SFcut(const std::unique_ptr<clas12::clas12reader>& c12)
{
  auto electrons=c12->getByID(11);

  TVector3 ve;
  ve.SetMagThetaPhi(electrons[0]->getP(),electrons[0]->getTheta(),electrons[0]->getPhi());
  double EoP_e =  (electrons[0]->cal(clas12::PCAL)->getEnergy() +  electrons[0]->cal(clas12::ECIN)->getEnergy() +  electrons[0]->cal(clas12::ECOUT)->getEnergy()) / ve.Mag();
  
  if(EoP_e <= cutmap[e_SF].min){ return false; }  
  if(EoP_e >= cutmap[e_SF].max){ return false; }
  return true;
}
bool eventcut::e_momcut(const std::unique_ptr<clas12::clas12reader>& c12)
{
  auto electrons=c12->getByID(11);

  if(electrons[0]->getP() <= cutmap[e_mom].min){ return false; }  
  if(electrons[0]->getP() >= cutmap[e_mom].max){ return false; }
  return true;
}




bool eventcut::l_scintcut(const std::unique_ptr<clas12::clas12reader>& c12, int i, int pid)
{
  auto leadnucleons=c12->getByID(pid);
  int FTOF1A = leadnucleons[i]->sci(clas12::FTOF1A)->getDetector();
  int FTOF1B = leadnucleons[i]->sci(clas12::FTOF1B)->getDetector();
  int FTOF2 = leadnucleons[i]->sci(clas12::FTOF2)->getDetector();
  int CTOF = leadnucleons[i]->sci(clas12::CTOF)->getDetector();
  if(FTOF1A == cutmap[l_scint].count){return true;}
  if(FTOF1B==cutmap[l_scint].count){return true;}
  if(FTOF2==cutmap[l_scint].count){return true;}
  if(CTOF==cutmap[l_scint].count){return true;}
  
  return false;
}

bool eventcut::l_thetacut(const std::unique_ptr<clas12::clas12reader>& c12, int i, int pid)
{
  auto leadnucleons=c12->getByID(pid);
  double theta = leadnucleons[i]->getTheta() * 180 / M_PI;
  if(theta < cutmap[l_theta].min){ return false; }
  if(theta > cutmap[l_theta].max){ return false; }
  return true;
}
bool eventcut::l_thetalqcut(const std::unique_ptr<clas12::clas12reader>& c12, int i, int pid)
{
  auto electrons=c12->getByID(11);
  TVector3 ve;
  ve.SetMagThetaPhi(electrons[0]->getP(),electrons[0]->getTheta(),electrons[0]->getPhi());

  auto leadnucleons=c12->getByID(pid);
  TVector3 vL;
  vL.SetMagThetaPhi(leadnucleons[i]->getP(),leadnucleons[i]->getTheta(),leadnucleons[i]->getPhi());
  
  TVector3 vq = vbeam - ve;
  
  double thetalq = vq.Angle(vL) * 180 / M_PI;
  if(thetalq < cutmap[l_thetalq].min){ return false; }
  if(thetalq > cutmap[l_thetalq].max){ return false; }
  return true;
}

bool eventcut::l_chipidcut(const std::unique_ptr<clas12::clas12reader>& c12, int i, int pid)
{
  auto leadnucleons=c12->getByID(pid);
  double Chi2Pid_L = leadnucleons[i]->par()->getChi2Pid();
  if(Chi2Pid_L < cutmap[l_chipid].min){ return false; }
  if(Chi2Pid_L > cutmap[l_chipid].max){ return false; }
  return true;
}




bool eventcut::lsrc_Q2cut(const std::unique_ptr<clas12::clas12reader>& c12)
{
  auto electrons=c12->getByID(11);
  TVector3 ve;
  ve.SetMagThetaPhi(electrons[0]->getP(),electrons[0]->getTheta(),electrons[0]->getPhi());
  TVector3 vq = vbeam - ve;
  double nu = Ebeam - ve.Mag();
  double Q2 = vq.Mag2() - (nu*nu);

  if(Q2 < cutmap[lsrc_Q2].min){ return false; }
  if(Q2 > cutmap[lsrc_Q2].max){ return false; }
  return true;
}

bool eventcut::lsrc_xBcut(const std::unique_ptr<clas12::clas12reader>& c12)
{
  auto electrons=c12->getByID(11);
  TVector3 ve;
  ve.SetMagThetaPhi(electrons[0]->getP(),electrons[0]->getTheta(),electrons[0]->getPhi());
  TVector3 vq = vbeam - ve;
  double nu = Ebeam - ve.Mag();
  double Q2 = vq.Mag2() - (nu*nu);
  double xB = Q2 / (2 * mN * nu);
  
  if(xB < cutmap[lsrc_xB].min){ return false; }
  if(xB > cutmap[lsrc_xB].max){ return false; }
  return true;
}

bool eventcut::lsrc_pmisscut(const std::unique_ptr<clas12::clas12reader>& c12, int i, int pid)
{
  auto electrons=c12->getByID(11);
  TVector3 ve;
  ve.SetMagThetaPhi(electrons[0]->getP(),electrons[0]->getTheta(),electrons[0]->getPhi());

  auto leadnucleons=c12->getByID(pid);
  TVector3 vL;
  vL.SetMagThetaPhi(leadnucleons[i]->getP(),leadnucleons[i]->getTheta(),leadnucleons[i]->getPhi());
  
  TVector3 vmiss = vbeam - ve - vL;

  if(vmiss.Mag() < cutmap[lsrc_pmiss].min){ return false; }
  if(vmiss.Mag() > cutmap[lsrc_pmiss].max){ return false; }
  return true;
}

bool eventcut::lsrc_mmisscut(const std::unique_ptr<clas12::clas12reader>& c12, int i, int pid)
{

  auto electrons=c12->getByID(11);
  TVector3 ve;
  ve.SetMagThetaPhi(electrons[0]->getP(),electrons[0]->getTheta(),electrons[0]->getPhi());

  auto leadnucleons=c12->getByID(pid);
  TVector3 vL;
  vL.SetMagThetaPhi(leadnucleons[i]->getP(),leadnucleons[i]->getTheta(),leadnucleons[i]->getPhi());
  
  TVector3 vmiss = vbeam - ve - vL;

  double Ee = ve.Mag();
  double Ep = sqrt((mN * mN) + vL.Mag2());
  double emiss = Ebeam + mD - Ee - Ep;
  double mmiss = sqrt((emiss * emiss) - vmiss.Mag2());

  if(mmiss < cutmap[lsrc_mmiss].min){ return false; }
  if(mmiss > cutmap[lsrc_mmiss].max){ return false; }
  return true;
  
}

bool eventcut::lsrc_loqcut(const std::unique_ptr<clas12::clas12reader>& c12, int i, int pid)
{

  auto electrons=c12->getByID(11);
  TVector3 ve;
  ve.SetMagThetaPhi(electrons[0]->getP(),electrons[0]->getTheta(),electrons[0]->getPhi());

  auto leadnucleons=c12->getByID(pid);
  TVector3 vL;
  vL.SetMagThetaPhi(leadnucleons[i]->getP(),leadnucleons[i]->getTheta(),leadnucleons[i]->getPhi());
  
  TVector3 vq = vbeam - ve;
  double Loq = vL.Mag()/vq.Mag();
  
  if(Loq < cutmap[lsrc_loq].min){ return false; }
  if(Loq > cutmap[lsrc_loq].max){ return false; }
  return true;

}

void eventcut::setl_scint(int i)
{
  cutInfo structcut;
  structcut.docut = true;
  structcut.count = i;
  cutmap[l_scint]=structcut;
}
