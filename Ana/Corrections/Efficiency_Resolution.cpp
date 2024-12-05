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
#include "reweighter.h"
#include "Corrections.h"

using namespace std;
using namespace clas12;

const double c = 29.9792458;
////////////////////////
///////////////////////
double getPhiDiff(TVector3 A, TVector3 B){
  TVector3 Aperp(A.X(),A.Y(),0);
  TVector3 Bperp(B.X(),B.Y(),0);
  double phidiff = (A.Angle(B))*180/M_PI;
  if(A.Cross(B).Z()<0){phidiff=phidiff*-1;}
  return phidiff;
}

//vector<double> bE_ThetaCD = {35,40,45,50,55,60,70};
vector<double> bE_MomCD = {0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,3.0};//{0.5,1.0,1.3,1.6,2.0,2.5,3.0};
vector<double> bE_Theta = {8,10,12,14,16,18,20,23,26,29,32,35,38,41,45};
vector<double> bE_ThetaCD = {38,41,44,47,50,53,56,59,62,65,68,71,74,77,80,83,86,89,92};//{40,45,50,55,60,70,90};
vector<double> bE_Phi = {-35,-15,-5,0,5,10,15,25,35};
vector<double> bE_ThetaE = {10,13,16,19,22,25,28,31,34,37};
vector<double> bE_ThetapFD = {19,22,25,28,31,34,37,40,43,46};

int binX(vector<double> XS, double X){
  for(int i = 0; i <= XS.size(); i++){if(X<XS[i]){return i-1;}}
  return -1;
}

auto db=TDatabasePDG::Instance();
double mass_p = db->GetParticle(2212)->Mass();
double mD = 1.8756;
double beam_E = 5.984792;
double beam_E_sigma = 0.00299;
const double me = 0.000511;
const double mU = 0.9314941024;
const double m_4He = 4.00260325415 * mU - 2*me;


double G(double x, double N, double mu, double sigma){
  return (N/(sigma*sqrt(2*M_PI))) * exp(-0.5 * sq((x-mu)/sigma)) ; 
}

double SQ(double x){ return x*x;}

void getGraph(TH2D * h_myhist, TGraphErrors * g_mygraph){
  int ctr = 0;
  char temp[100];
  //Now project the histogram    
  for(int j = 0; j < h_myhist->GetXaxis()->GetNbins(); j++){
    //Define x and y(1D histogram)
    double x = h_myhist->GetXaxis()->GetBinCenter(j+1);
    ctr++;
    sprintf(temp,"Proj_num%d",ctr);
    TH1D * proj = h_myhist->ProjectionY(temp,j+1,j+1);
    //proj->Rebin(2);
    //Now preform a guassian fit
    if(proj->GetEntries()<50){continue;}

    TF1 * gFit = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-1,1,3);
    gFit->SetParameter(0,proj->GetMaximum()/G(0,1,0,0.1));
    gFit->SetParameter(1,0.0);
    gFit->SetParLimits(1,-0.8,0.8);
    gFit->SetParameter(2,0.2);
    gFit->SetParLimits(2,0.001,1.0);
    
    TFitResultPtr gPoint = proj->Fit(gFit,"SrBeqn","",-1,1);
    if(gPoint == 0){
      g_mygraph->SetPoint(g_mygraph->GetN(),x,gPoint->Parameter(1));
      g_mygraph->SetPointError(g_mygraph->GetN()-1,0,gPoint->Parameter(2));
    }
    proj->Write();
  }
}

void Usage()
{
  std::cerr << "Usage: ./code isMC outputfile.root outputfile.pdf inputfiles.hipo \n\n\n";

}

int main(int argc, char ** argv)
{

  if(argc < 4)
    {
      Usage();
      return -1;
    }



  int isMC = atoi(argv[1]);
  TString outFile = argv[2];
  char * pdfFile = argv[3];

  cout<<"Ouput file "<< outFile <<endl;
  cout<<"Ouput PDF file "<< pdfFile <<endl;


  clas12ana clasAna;
  clasAna.printParams();
    
  clas12root::HipoChain chain;
  for(int k = 4; k < argc; k++){
    cout<<"Input file "<<argv[k]<<endl;
    chain.Add(argv[k]);
  }
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();
  auto config_c12=chain.GetC12Reader();
  auto &c12=chain.C12ref();
  
  double mN = db->GetParticle(2212)->Mass();
  //some particles
  TLorentzVector beam(0,0,beam_E,beam_E);
  TVector3 vbeam(0,0,beam_E);
  TLorentzVector nucleus_ptr(0,0,0,m_4He);
  TLorentzVector deut_ptr(0,0,0,mD);
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector el_corrected(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector proton_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector proton_ptr_corrected(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector lead_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector lead_pion_ptr(0,0,0,db->GetParticle(211)->Mass());
  TLorentzVector stationary_proton_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector ntr(0,0,0,db->GetParticle(2112)->Mass());
  reweighter newWeight(beam_E,2,2,kelly,"AV18");

  //////////////////////////////////
  char temp_name[100];
  char temp_title[100];

  vector<TH1*> hist_list;

  TH1D * h_Q2FD_e = new TH1D("Q2FD_e","Q^{2} (e,e');Q^{2};Counts",100,0,5);
  hist_list.push_back(h_Q2FD_e);
  TH1D * h_Q2FD_ep = new TH1D("Q2FD_ep","Q^{2} (e,e'p);Q^{2};Counts",100,0,5);
  hist_list.push_back(h_Q2FD_ep);

  TH1D * h_thetaFD_e = new TH1D("thetaFD_e","#theta_{q} (e,e');#theta_{q};Counts",100,5,55);
  hist_list.push_back(h_thetaFD_e);
  TH1D * h_thetaFD_ep = new TH1D("thetaFD_ep","#theta_{q} (e,e'p);#theta_{q};Counts",100,5,55);
  hist_list.push_back(h_thetaFD_ep);

  TH1D * h_DEpFD_ep = new TH1D("DEpFD_ep","#Delta E' (e,e'p);#Delta #E';Counts",100,-0.5,0.5);
  hist_list.push_back(h_DEpFD_ep);
  TH1D * h_phidiffFD_ep = new TH1D("phidiffFD_ep","#Delta #phi (e,e'p);#Delta #phi;Counts",100,-15,15);
  hist_list.push_back(h_phidiffFD_ep);
  TH1D * h_angleFD_ep = new TH1D("angleFD_ep","#theta (e,e'p);#theta;Counts",100,0,50);
  hist_list.push_back(h_angleFD_ep);

  
  
  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Sumw2();
    hist_list[i]->GetXaxis()->CenterTitle();
    hist_list[i]->GetYaxis()->CenterTitle();
  }

  int counter = 0;
  //while(chain.Next())
  while(chain.Next() && (counter<100000000))
    {
      double wep = 1;
      double wepp = 1;
      if(isMC==1){
	double original_weight = c12->mcevent()->getWeight(); //used if MC events have a weight
	wep = original_weight;
      }
      //Display completed  
      counter++;
      if((counter%1000000) == 0){
	cerr << "\n" <<counter/1000000 <<" million completed";
      }    
      if((counter%100000) == 0){
	cerr << ".";
      }    

      clasAna.Run(c12);
      auto electrons = clasAna.getByPid(11);
      auto protons = clasAna.getByPid(2212);
      auto particles = c12->getDetParticles(); //particles is now 
      if(electrons.size() == 1 && protons.size() >= 0)
	{

	  GetLorentzVector_ReconVector(el,electrons[0]);
	  if(isMC==0){
	    SetLorentzVector_ThetaCorrection(el,electrons[0]);
	    SetLorentzVector_MomentumCorrection(el,electrons[0]);
	  }

	  double Eprime = ( mN*beam_E )/( beam_E*(1-cos(el.Theta())) + mN );
	  TVector3 vel_fromAngle;
	  vel_fromAngle.SetMagThetaPhi(Eprime,el.Theta(),el.Phi());
	  TVector3 vp_fromAngle = vbeam-vel_fromAngle;
	  double Delta_Eprime = ( mN*beam_E )/( beam_E*(1-cos(el.Theta())) + mN ) - el.P();
	  ///
	  
	  int sector_e = electrons[0]->getSector();
	  TLorentzVector q = beam - el;
          double Q2        = -q.M2();
	  double omega = q.E();
          double xB       = Q2/(2 * mass_p * (beam.E() - el.E()) );

	  double vtz_e = electrons[0]->par()->getVz();
	  double WSq = (mN*mN) - Q2 + (2*omega*mN);
	  double W = sqrt(WSq);
	  double phi_e = el.Phi() * 180/M_PI;
	  double theta_e = el.Theta()*180/M_PI;
	  double theta_q = q.Theta()*180/M_PI;
	  double shift_e = 0;
	  shift_e += (sector_e==1)?0:(sector_e==2)?60:(sector_e==3)?120:(sector_e==4 && phi_e>0)?180:(sector_e==4 && phi_e<0)?-180:(sector_e==5)?-120:(sector_e==6)?-60:0;

	  if(fabs(Delta_Eprime)>0.1){continue;}
	  h_thetaFD_e->Fill(theta_q,wep);
	  h_Q2FD_e->Fill(Q2,wep);
	  
	  if(protons.size() <= 0){continue;}
	  if(protons[0]->getRegion()!=FD){continue;}
	  GetLorentzVector_ReconVector(proton_ptr,protons[0]);
	  if(isMC==0){
	    SetLorentzVector_ThetaCorrection(proton_ptr,protons[0]);
	  }
	  SetLorentzVector_EnergyLossCorrection(proton_ptr,protons[0]);
	  if(isMC==0){
	    SetLorentzVector_MomentumCorrection(proton_ptr,protons[0]);
	  }
	  
	  double mom_q = q.P();
	  double phi_q = q.Phi()*180/M_PI;
	  double mom_p = proton_ptr.P();
	  double Delta_Mom_FromAngle = vp_fromAngle.Mag()-mom_p;
	  double theta_p = proton_ptr.Theta()*180/M_PI;
	  double phi_p = proton_ptr.Phi()*180/M_PI;
	  double Delta_mom = (mom_p-mom_q)/mom_q;
	  double Delta_theta = theta_p-theta_q;
	  //double Delta_phi = getPhiDiff(q.Vect(),proton_ptr.Vect());;
	  double Delta_phi = (phi_p-phi_q);
	  if(Delta_phi<-180){Delta_phi+=360;}
	  else if(Delta_phi>180){Delta_phi-=360;}

	  h_DEpFD_ep->Fill(Delta_Eprime,wep);
	  h_phidiffFD_ep->Fill(Delta_phi,wep);
	  h_angleFD_ep->Fill(q.Vect().Angle(proton_ptr.Vect())*180/M_PI,wep);
	  if(fabs(Delta_phi)>3){continue;}
	  if((q.Vect().Angle(proton_ptr.Vect())*180/M_PI)>5){continue;}
	  h_Q2FD_ep->Fill(Q2,wep);	  
	  h_thetaFD_ep->Fill(theta_q,wep);

	}
    }
  
  /////////////////////////////////////////////////////
  //Now create the output PDFs
  /////////////////////////////////////////////////////
  TFile *f = new TFile(outFile,"RECREATE");
  f->cd();
  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Write();
  }
  //g_omega_diff_mu->Write();
  for(int j=1; j<=6; j++){
    //g_omega_diff_mu_sectors[j-1]->Write();
  }
  int pixelx = 1980;
  int pixely = 1530;
  TCanvas * myCanvas = new TCanvas("myPage","myPage",pixelx,pixely);
  TCanvas * myText = new TCanvas("myText","myText",pixelx,pixely);
  TLatex text;
  text.SetTextSize(0.05);
  
  char fileName[100];
  sprintf(fileName,"%s[",pdfFile);
  myText->SaveAs(fileName);
  sprintf(fileName,"%s",pdfFile);

  ///////////////////////////////////
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_DEpFD_ep->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_phidiffFD_ep->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_angleFD_ep->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_thetaFD_e->Draw();
  h_thetaFD_ep->SetLineColor(2);
  h_thetaFD_ep->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_thetaFD_ep->Divide(h_thetaFD_e);
  h_thetaFD_ep->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  f->Close();


  return 0;
}

  


/*
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

double params_Momentum_FD[6][4][4]={{{0.00999685,-0.0274693,97.8319,27.1649},
				{0.107308,0.0482255,6.00016,26.5809},
				{126.732,40.9353,4,15.7632},
				{-1.29807,-2.08763,5.6061,21.9216}},
			       {{0.0863814,0.0662239,78.8143,22.7513},
				{0.0819542,0.0353392,4.05495,26.0594},
				{123.313,39.0737,4,17.3371},
				{-1.23367,-2.10145,4.37322,21.3527}},
			       {{0.136301,0.0624967,4.06521,39.5042},
				{0.0573779,0.0258018,4.99407,24.1826},
				{119.322,74.2621,7.26117,10},
				{-4.36463,-2.59454,7.47078,36.4093}},
			       {{0.0273149,-0.009257,20.2085,24.9454},
				{0.0689756,0.0280117,4.01384,30.0473},
				{128.494,53.5087,4.35928,17.9939},
				{-1.59758,-1.5448,4.0001,29.9412}},
			       {{0.114001,0.0792582,37.7721,36.9068},
				{0.075125,0.0282635,4.00108,29.5526},
				{138.358,59.4015,4,10},
				{-0.996918,-1.03077,4.00001,27.4224}},
			       {{0.108582,0.0801866,79.5092,28.1964},
				{0.068627,0.0341774,16.8015,27.416},
				{130.348,48.6141,4,17.9232},
				{-0.926971,-0.686175,4,20.191}}};

double params_Momentum_CD[3][2]={{-0.107267,0.00193308},
				 {0.0267032,0.000545869},
				 {-2.31836,0.0177897}};

double Function_Erf(double x, double A, double B, double C, double D){
  return A - B*(1+erf(((-x+D)/C))); 
}

double Function_Trig(double x, double A, double B, double C, double D){
  return A + B*sin((x*2*M_PI/C)+D); 
}

double Function_TrigFixedPeriod(double x, double A, double B, double C){
  return A + B*sin((x*2*M_PI/180)+C); 
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
      params[i] = Function_Erf(theta,params_Momentum_FD[sector-1][i][0],params_Momentum_FD[sector-1][i][1],params_Momentum_FD[sector-1][i][2],params_Momentum_FD[sector-1][i][3]);
    }    
    mom+=Function_Trig(phi,params[0],params[1],params[2],params[3]);
  }
  else if(p->getRegion()==CD){
    double params[3];
    for(int i = 0; i < 3; i++){
      params[i] = Function_Algebraic1(theta,params_Momentum_CD[i][0],params_Momentum_CD[i][1]);
    }    
    mom+=Function_TrigFixedPeriod(phi,params[0],params[1],params[2]);
  }
  else{
    cout<<"Problem\n\n\n\n\n\n\n\n\n";
    exit(-2);
  }
  v3.SetMag(mom);
  p4.SetXYZM(v3.X(),v3.Y(),v3.Z(),M);  
}
*/

/*
double GetThetaCorrected(clas12::region_particle* p){
  TVector3 v(p->par()->getPx(),p->par()->getPy(),p->par()->getPz());
  double theta = v.Theta()*180/M_PI;
  double phi = v.Phi()*180/M_PI;
  if(p->getRegion()==FD){
    int sector = p->getSector();
    double shift = 0;
    shift += (sector==1)?0:(sector==2)?60:(sector==3)?120:(sector==4 && phi>0)?180:(sector==4 && phi<0)?-180:(sector==5)?-120:(sector==6)?-60:0;
    phi-=shift;

    double params[4];
    for(int i = 0; i < 4; i++){
      params[i] = Function_Erf(theta,params_Theta_FD[sector-1][i][0],params_Theta_FD[sector-1][i][1],params_Theta_FD[sector-1][i][2],params_Theta_FD[sector-1][i][3]);
    }
    return Function_Trig(phi,params[0],params[1],params[2],params[3]);
  }
  else if(p->getRegion()==CD){
    double params[3];
      for(int i = 0; i < 3; i++){
	params[i] = Function_Algebraic2(theta,params_Theta_CD[i][0],params_Theta_CD[i][1]);
      }
    return Function_TrigFixedPeriod(phi,params[0],params[1],params[2]);    
  }
  cout<<"Problem\n\n\n\n\n\n\n\n\n";
  return -1000.00;
}


void SetLorentzVectorCorrectedTheta(TLorentzVector &p4, clas12::region_particle* p){
  TVector3 vInit(p->par()->getPx(),p->par()->getPy(),p->par()->getPz());
  double newtheta = vInit.Theta() + (GetThetaCorrected(p)*M_PI/180);
  TVector3 vCorr;
  vCorr.SetMagThetaPhi(vInit.Mag(),newtheta,vInit.Phi());
  p4.SetXYZM(vCorr.X(),vCorr.Y(),vCorr.Z(),p4.M());
}


double FDCorr(double mom, double theta){
  double A[2]= {-0.000695124,-0.000355869};
  double B[2]= {0.00182181,7.77933e-05};
  double C[3]= {0.000266541,0.424055,49.068};  
  double a = Function_Algebraic1(theta,A[0],A[1]);
  double b = Function_Algebraic1(theta,B[0],B[1]);
  double c = Function_Algebraic4(theta,C[0],C[1],C[2]);
  double correction = Function_Algebraic3(mom,a,b,c);
  return correction;
}

double CDCorr(double mom, double theta){
  double A[2]= {-0.00555986,-6.06201e-05};
  double B[2]= {0.00695634,8.24535e-05};
  double C[2]= {0.00155103,1.74283e-05};
  double a = Function_Algebraic1(theta,A[0],A[1]);
  double b = Function_Algebraic1(theta,B[0],B[1]);
  double c = Function_Algebraic1(theta,C[0],C[1]);
  double correction = Function_Algebraic3(mom,a,b,c);
  return correction;
}

void GetMomEnergyLossProtonsCorrected(TLorentzVector &p4, clas12::region_particle* p){
  TVector3 v3 = p4.Vect();
  double mom = v3.Mag();
  double theta = v3.Theta()*180/M_PI;
  double M = p4.M();        
  if(p->getRegion()==FD){
    mom+=FDCorr(mom,theta);
  }
  else if(p->getRegion()==CD){
    mom+=CDCorr(mom,theta);
  }
  else{
    cout<<"Problem\n\n\n\n\n\n\n\n\n";
    exit(-2);
  }
  v3.SetMag(mom);
  p4.SetXYZM(v3.X(),v3.Y(),v3.Z(),M);  
}

*/
