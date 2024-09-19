#include "reweighter.h"

reweighter::reweighter(double E, int Z, int N)
{
  Ebeam = E;
  
  Z_nuc = Z;
  N_nuc = N;

  uType_init="AV18";
  gcf_config_init = new gcfSRC(2,2,uType_init);
  sigma_cm_init = 0.2;
  CS_config_init = new eNCrossSection(cc1,kelly);

  uType_fin="AV18";
  gcf_config_fin = new gcfSRC(Z_nuc,N_nuc,uType_fin);
  sigma_cm_fin = 0.15;
  CS_config_fin = new eNCrossSection(cc1,kelly);
  /*
  double P_new[4][2] = {{10,10},
			{10,10},
			{10,10},
			{10,10}};
  */ 
 
  double P_new[4][2] = {{4.1,3.5},
			{4.8,4.1},
			{4.8,4.1},
			{4.1,3.5}};
  memcpy(P,P_new,sizeof(P));
  TN = 0.53;
  TNN = 0.44;  
}

reweighter::~reweighter()
{
}

double reweighter::get_weight_noT(clas12::mcparticle* mcInfo)
{
  double den = 1;
  double num = 1;

  //Grabe the momentum values
  TVector3 vbeam(0,0,Ebeam);
  TVector3 ve(mcInfo->getPx(0),mcInfo->getPy(0),mcInfo->getPz(0));
  TVector3 vlead(mcInfo->getPx(1),mcInfo->getPy(1),mcInfo->getPz(1));
  TVector3 vrec(mcInfo->getPx(2),mcInfo->getPy(2),mcInfo->getPz(2));

  TVector3 vq = vbeam - ve;
  TVector3 vmiss = vlead - vq;
  TVector3 vcm = vmiss + vrec;
  TVector3 vrel = 0.5 * (vmiss - vrec);

  //Get the PIDs and PIDs under Single Charge Exchange
  int leadCode = mcInfo->getPid(1);
  int leadCodeX = (leadCode==pCode)?nCode:pCode;
  int recCode = mcInfo->getPid(2);
  int recCodeX = (recCode==pCode)?nCode:pCode;


  //Grab the correct index of the 2d array
  //pp=0
  //pn=1
  //np=2
  //nn=3
  int indexP = 0;
  indexP += (leadCode==nCode)?2:0;
  indexP += (recCode==nCode)?1:0;
  double P_L_RX = P[indexP][0]/100;
  double P_LX_R = P[indexP][1]/100;
  double P_L_R = 1;
  if((leadCode==pCode)&&(recCode==pCode)){
    P_L_R -= (P[1][0]+P[2][1])/100;
  }
  else if((leadCode==pCode)&&(recCode==nCode)){
    P_L_R -= (P[0][0]+P[3][1])/100;
  }
  else if((leadCode==nCode)&&(recCode==pCode)){
    P_L_R -= (P[0][1]+P[3][0])/100;
  }
  else if((leadCode==nCode)&&(recCode==nCode)){
    P_L_R -= (P[1][1]+P[2][0])/100;
  }


  //Reweight for center of mass momentum
  den *= Gauss(vcm.X(),0,sigma_cm_init)*Gauss(vcm.Y(),0,sigma_cm_init)*Gauss(vcm.Z(),0,sigma_cm_init); 
  num *= Gauss(vcm.X(),0,sigma_cm_fin)*Gauss(vcm.Y(),0,sigma_cm_fin)*Gauss(vcm.Z(),0,sigma_cm_fin);
  
  //Reweight for Potential and Single Charge Exchange
  den *= CS_config_init->sigma_eN(Ebeam,ve,vlead,(leadCode==pCode))*gcf_config_init->get_S(vrel.Mag(),leadCode,recCode);

  //get sigma*S
  double sig_S = 0;
  sig_S += CS_config_fin->sigma_eN(Ebeam,ve,vlead,(leadCode==pCode))*gcf_config_fin->get_S(vrel.Mag(),leadCode,recCode)*P_L_R;
  sig_S += CS_config_fin->sigma_eN(Ebeam,ve,vlead,(leadCode==pCode))*gcf_config_fin->get_S(vrel.Mag(),leadCode,recCodeX)*P_L_RX;
  sig_S += CS_config_fin->sigma_eN(Ebeam,ve,vlead,(leadCodeX==pCode))*gcf_config_fin->get_S(vrel.Mag(),leadCodeX,recCode)*P_LX_R;
  num *= sig_S;

  return num/den;
}

double reweighter::get_weight_ep(clas12::mcparticle* mcInfo)
{
  return get_weight_noT(mcInfo)*TN;
}

double reweighter::get_weight_epp(clas12::mcparticle* mcInfo)
{
  return get_weight_noT(mcInfo)*TNN;
}

double reweighter::Gauss(double x, double mu, double sigma){
  return (1/(sigma*sqrt(2*M_PI))) * exp(-pow((x-mu)/sigma,2)/2);
}
