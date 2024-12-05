#ifndef CORRECTIONS_HH
#define CORRECTIONS_HH
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

#include <cstdlib>
#include <iostream>
#include <chrono>
#include <vector>
#include <typeinfo>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TDatabasePDG.h>
#include "HipoChain.h"
#include "clas12ana.h"
#include "TRandom3.h"

TRandom3 * thisRand = new TRandom3(0);;

double params_Theta_FD[6][4][4]={{{-0.191151,-0.669629,10,1.00023},
			    {0.113741,-0.0599251,19.0221,31.6623},
			    {350.715,460.172,17.1609,1},
			    {3.32022,2.7311,29.6054,25.1635}},
			   {{-0.0891744,-0.14685,7.50659,25.2638},
			    {0.0450352,-0.109521,6.0011,19.0888},
			    {350.357,140.077,6,15.8616},
			    {-3.25284,-1.24789,10.6854,24.956}},
			   {{-0.0849098,-0.0651043,11.7406,19.6928},
			    {0.0209555,-0.0273879,6.00308,20.1268},
			    {98.1241,20.8179,6,29.528},
			    {-5.06788,-4.9607,20.1145,1.00211}},
			   {{-0.0739777,0.0963433,6.00048,17.6548},
			    {0.014477,-0.156077,9.29576,12.1747},
			    {341.693,154.122,4.66104,22.2281},
			    {-6.38361,-1.95148,6,29.7998}},
			   {{-0.0528101,-0.100172,6.00982,12.4738},
			    {0.010489,-0.0339853,6.09838,16.0642},
			    {348.057,161.215,6,21.4269},
			    {-4.07125,-10,6,2.57348}},
			   {{-0.0623622,-0.179411,7.77701,9.66338},
			    {0.040617,-0.185013,6.00574,9.96181},
			    {199.345,74.073,6,28.6864},
			    {-3.83103,-1.41067,6,11.1274}}};
  
double params_Theta_CD[3][2]={{-0.19764,8.14947},
			      {-0.0207092,4.52151},
			      {-1.45427,0}};

double params_EnergyLoss_FD[3][3]={{-0.000695124,-0.000355869,0},
				   {0.00182181,7.77933e-05,0},
				   {0.000266541,0.424055,49.068}};

double params_EnergyLoss_CD[3][2]={{-0.00555986,-6.06201e-05},
				   {0.00695634,8.24535e-05},
				   {0.00155103,1.74283e-05}};

double params_Momentum_negparts_FD[6][4][4]={{{-0.00214337,-0.0202413,19.9966,19.9994},
					      {0.0673539,0.0273376,5.8074,24.9818},
					      {145.567,50.4663,4,16.4935},
					      {0.470742,-1.15944,4.00003,19.082}},
					     {{0.0376835,0.0166851,8.23076,24.9989},
					      {0.0235114,0.00170224,12.7113,15.3318},
					      {140.222,31.2228,1.5,13.9711},
					      {-0.0372922,-1.9027,1.5,20.4744}},
					     {{0.0205703,0.0146918,5.7492,19.9938},
					      {0.0711701,0.0328288,4.29665,24.9937},
					      {115.006,22.6111,1.5,23.0711},
					      {2.67023,1.55987,4,15.5237}},
					     {{-0.00204655,-0.0207699,19.9975,19.9986},
					      {0.0588647,0.0229838,5.62843,24.9975},
					      {142.111,66.05,8.64234,17.5045},
					      {2.25598,1.26862,5.72383,13.3042}},
					     {{0.0424364,0.0494246,19.2225,23.2837},
					      {0.0489891,0.0145193,8.31129,24.9937},
					      {148.98,65.4393,4,10},
					      {1.01252,-0.181366,4.00001,10.0004}},
					     {{0.0159458,0.0106683,1.51661,18.2011},
					      {0.0591658,0.0269474,4.17777,24.9879},
					      {125.981,50.2509,4,15.6732},
					      {0.97521,0.672086,4.00002,24.9999}}};

double params_Momentum_posparts_FD[6][4][4]={{{0.0994698,0.06203,4.00001,37.9995},
					      {0.083833,0.000338038,20.8983,30.0176},
					      {102.196,-14.058,4,27.4109},
					      {-2.70134,-0.335195,4.00079,37.9999}},
					     {{0.100609,0.0714968,4.01151,37.9997},
					      {0.0674467,0.0195951,4.53667,26.0002},
					      {101.988,-14.4089,4,27.1633},
					      {-2.66578,-0.241631,4,38}},
					     {{0.131204,0.0895802,3.53093,37.9672},
					      {0.0189443,0.0420609,19.8174,29.9951},
					      {107.198,26.2374,2,38},
					      {-2.92875,-0.0467252,2.00003,32.0315}},
					     {{0.0902545,0.0582306,2.02946,37.0886},
					      {0.0368395,0.00706995,13.4986,26.001},
					      {98.4986,-12.3151,2.97928,36.2819},
					      {-1.8963,-0.652783,5.99785,34.34}},
					     {{0.139861,0.104554,5.27712,37.924},
					      {-0.0329756,-0.0325642,4.00869,26.0012},
					      {106.329,16.3735,4,33.941},
					      {0.241212,-1.18799,23.8338,37.9821}},
					     {{0.146429,0.0995005,4.14491,36.9621},
					      {-0.0281189,-0.0121676,4.36001,30.4651},
					      {107.814,28.5356,4,33.207},
					      {-0.324309,-2.82037,31.3704,32.1691}}};

double params_Momentum_CD[2][4]={{0.0585217,0.00989627,5.00872,42.3376},
				 {-1.11617,0.0736438,5.0042,48.9754}};


double params_Smear_FD[6][6]={{1.27164,-0.103322,0.00374144,0.22998,0.00155623,0.000888604},
		      {0.716006,-0.0386559,0.00203629,0.557439,-0.0323154,0.00163883},
		      {1.10877,-0.0853149,0.00335627,-0.465281,0.0735502,-0.000816378},
		      {0.916793,-0.0645254,0.0026928,0.267525,-0.00692347,0.00114717},
		      {1.18955,-0.0993153,0.00359098,0.262575,-0.00198241,0.000950544},
		      {0.486033,-0.0162104,0.0015784,-0.143655,0.0362722,0.000199701}};

double params_Smear_CD[8]={10.049,3.07365,0.461349,1.47116,-0.832429,1.45479,1.67907,3.90999};

/*
//Vertex Smearing

Electrons=
0.309387&1.85845
0.27904&2.80253
0.336598&1.9909
0.310611&2.25277
0.245308&2.59497
0.232689&2.9493

Protons=
0.763939&16.1967
0.629776&18.9645
0.588898&19.4205
0.500838&21.3636
0.729782&15.4736
0.643492&16.808

Protons CD 
0.82702&-0.016491&9.69232e-05
*/
double Function_Erf(double x, double A, double B, double C, double D){
  return A - B*(1+erf(((-x+D)/C))); 
}

double Function_Trig(double x, double A, double B, double C, double D){
  return A + B*sin((x*2*M_PI/C)+D); 
}

double Function_TrigFixedPeriodFixedOffset(double x, double A, double B){
  return A*sin((x*2*M_PI/180)+B); 
}

double Function_TrigFixedPeriod(double x, double A, double B, double C){
  return A + B*sin((x*2*M_PI/180)+C); 
}

double Function_TrigTripleFixedPeriod(double x, double A, double B, double C, double D, double E, double F, double G){
  return A + B*sin((x*1*M_PI/180)+C) + D*sin((x*2*M_PI/180)+E) + F*sin((x*3*M_PI/180)+G); 
}

double Function_Algebraic1(double x, double A, double B){
  return A + B*x; 
}

double Function_Algebraic2(double x, double A, double B){
  return A + B/(x); 
}

double Function_Algebraic3(double x, double A, double B, double C){
  return A + B*x + (C/x); 
}

double Function_Algebraic4(double x, double A, double B, double C){
  return A - B/(x-C); 
}

double Function_Algebraic5(double x, double A, double B, double C){
  return A + B*x + C*x*x; 
}

double SubtractQuadriture_FD(double x, double A, double B, double C, double D, double E, double F){
  double X = Function_Algebraic5(x,A,B,C)*Function_Algebraic5(x,A,B,C) - Function_Algebraic5(x,D,E,F)*Function_Algebraic5(x,D,E,F);
  if(X>0.01){
    return sqrt(X); 
  }
  return 0.1;
}

double SubtractQuadriture_CD(double x, double A, double B, double C, double D, double E, double F, double G, double H){
  double X = Function_TrigTripleFixedPeriod(x,A,B,C,D,E,F,G)*Function_TrigTripleFixedPeriod(x,A,B,C,D,E,F,G) - H*H;
  if(X>0.01){
    return sqrt(X); 
  }
  return 0.1;
}

void GetLorentzVector_ReconVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),rp->par()->getPz(),p4.M());
}

void SetLorentzVector_ThetaCorrection(TLorentzVector &p4, clas12::region_particle* p){
  TVector3 v3 = p4.Vect();
  double mom = v3.Mag();
  double theta = v3.Theta()*180/M_PI;
  double phi = v3.Phi()*180/M_PI;
  double M = p4.M();        
  if(p->getRegion()==FD){
    int sector = p->getSector();
    double shift = 0;
    shift += (sector==1)?0:(sector==2)?60:(sector==3)?120:(sector==4 && phi>0)?180:(sector==4 && phi<0)?-180:(sector==5)?-120:(sector==6)?-60:0;
    phi-=shift;

    double params[4];
    for(int i = 0; i < 4; i++){
      params[i] = Function_Erf(theta,params_Theta_FD[sector-1][i][0],params_Theta_FD[sector-1][i][1],params_Theta_FD[sector-1][i][2],params_Theta_FD[sector-1][i][3]);
    }
    theta+=Function_Trig(phi,params[0],params[1],params[2],params[3]);
  }
  else if(p->getRegion()==CD){
    double params[3];
    for(int i = 0; i < 3; i++){
      params[i] = Function_Algebraic2(theta,params_Theta_CD[i][0],params_Theta_CD[i][1]);
    }
    theta+=Function_TrigFixedPeriod(phi,params[0],params[1],params[2]);        
  }
  else{
    cout<<"Problem\n\n\n\n\n\n\n\n\n";
    exit(-2);
  }
  v3.SetTheta(theta*M_PI/180);
  p4.SetXYZM(v3.X(),v3.Y(),v3.Z(),M);  
}

void SetLorentzVector_EnergyLossCorrection(TLorentzVector &p4, clas12::region_particle* p){
  TVector3 v3 = p4.Vect();
  double mom = v3.Mag();
  double theta = v3.Theta()*180/M_PI;
  double phi = v3.Phi()*180/M_PI;
  double M = p4.M();        
  if(p->par()->getCharge()<=0){return;}
  if(p->getRegion()==FD){
    double params0 = Function_Algebraic1(theta,params_EnergyLoss_FD[0][0],params_EnergyLoss_FD[0][1]);
    double params1 = Function_Algebraic1(theta,params_EnergyLoss_FD[1][0],params_EnergyLoss_FD[1][1]);
    double params2 = Function_Algebraic4(theta,params_EnergyLoss_FD[2][0],params_EnergyLoss_FD[2][1],params_EnergyLoss_FD[2][2]);
    mom+=Function_Algebraic3(mom,params0,params1,params2);
  }
  else if(p->getRegion()==CD){
    double params[3];
    for(int i = 0; i < 3; i++){
      params[i] = Function_Algebraic1(theta,params_EnergyLoss_CD[i][0],params_EnergyLoss_CD[i][1]);
    }
    mom+=Function_Algebraic3(mom,params[0],params[1],params[2]);    
      }
  else{
    cout<<"Problem\n\n\n\n\n\n\n\n\n";
    exit(-2);
  }
  v3.SetMag(mom);
  p4.SetXYZM(v3.X(),v3.Y(),v3.Z(),M);  
} 

void SetLorentzVector_MomentumCorrection(TLorentzVector &p4, clas12::region_particle* p){
  TVector3 v3 = p4.Vect();
  double mom = v3.Mag();
  double theta = v3.Theta()*180/M_PI;
  double phi = v3.Phi()*180/M_PI;
  double M = p4.M();        
  if(p->getRegion()==FD){
    int sector = p->getSector();
    double shift = 0;
    shift += (sector==1)?0:(sector==2)?60:(sector==3)?120:(sector==4 && phi>0)?180:(sector==4 && phi<0)?-180:(sector==5)?-120:(sector==6)?-60:0;
    phi-=shift;

    double params[4];
    for(int i = 0; i < 4; i++){
      if(p->par()->getCharge()<0){
	params[i] = Function_Erf(theta,params_Momentum_negparts_FD[sector-1][i][0],params_Momentum_negparts_FD[sector-1][i][1],params_Momentum_negparts_FD[sector-1][i][2],params_Momentum_negparts_FD[sector-1][i][3]);
      }
      else{
	params[i] = Function_Erf(theta,params_Momentum_posparts_FD[sector-1][i][0],params_Momentum_posparts_FD[sector-1][i][1],params_Momentum_posparts_FD[sector-1][i][2],params_Momentum_posparts_FD[sector-1][i][3]);
      }
    }    
    mom+=Function_Trig(phi,params[0],params[1],params[2],params[3]);
  }
  else if(p->getRegion()==CD){
    double params[2];
    for(int i = 0; i < 2; i++){
      params[i] = Function_Erf(theta,params_Momentum_CD[i][0],params_Momentum_CD[i][1],params_Momentum_CD[i][2],params_Momentum_CD[i][3]);
    }
    mom+=Function_TrigFixedPeriodFixedOffset(phi,params[0],params[1]);

    
  }
  else{
    cout<<"Problem\n\n\n\n\n\n\n\n\n";
    exit(-2);
  }
  v3.SetMag(mom);
  p4.SetXYZM(v3.X(),v3.Y(),v3.Z(),M);  
}

void SetLorentzVector_MomentumSimulationSmear(TLorentzVector &p4, clas12::region_particle* p){
  TVector3 v3 = p4.Vect();
  double mom = v3.Mag();
  double theta = v3.Theta()*180/M_PI;
  double phi = v3.Phi()*180/M_PI;
  double M = p4.M();        
  if(p->getRegion()==FD){
    int sector = p->getSector();
    double smear = 0.01*SubtractQuadriture_FD(theta,params_Smear_FD[sector-1][0],params_Smear_FD[sector-1][1],params_Smear_FD[sector-1][2],params_Smear_FD[sector-1][3],params_Smear_FD[sector-1][4],params_Smear_FD[sector-1][5]);
    mom*=thisRand->Gaus(1.0,smear);
  }
  else if(p->getRegion()==CD){
    double smear = 0.01*SubtractQuadriture_CD(phi,params_Smear_CD[0],params_Smear_CD[1],params_Smear_CD[2],params_Smear_CD[3],params_Smear_CD[4],params_Smear_CD[5],params_Smear_CD[6],params_Smear_CD[7]);
    mom*=thisRand->Gaus(1.0,smear);
  }
  else{
    cout<<"Problem\n\n\n\n\n\n\n\n\n";
    exit(-2);
  }
  v3.SetMag(mom);
  p4.SetXYZM(v3.X(),v3.Y(),v3.Z(),M);  
}


void GetLorentzVector_Corrected(TLorentzVector &p4, clas12::region_particle* p, bool isMC){
  GetLorentzVector_ReconVector(p4,p);
  if(!isMC){SetLorentzVector_ThetaCorrection(p4,p);}
  SetLorentzVector_EnergyLossCorrection(p4,p);  
  if(!isMC){SetLorentzVector_MomentumCorrection(p4,p);}
  if(isMC){SetLorentzVector_MomentumSimulationSmear(p4,p);}
}


#endif

