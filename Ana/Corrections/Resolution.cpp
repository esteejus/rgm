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
#include "TRandom3.h"

using namespace std;
using namespace clas12;

const double c = 29.9792458;
/*
TRandom3 * thisRand = new TRandom3(0);;
double SmearFD[6][6]={{1.27164,-0.103322,0.00374144,0.22998,0.00155623,0.000888604},
		      {0.716006,-0.0386559,0.00203629,0.557439,-0.0323154,0.00163883},
		      {1.10877,-0.0853149,0.00335627,-0.465281,0.0735502,-0.000816378},
		      {0.916793,-0.0645254,0.0026928,0.267525,-0.00692347,0.00114717},
		      {1.18955,-0.0993153,0.00359098,0.262575,-0.00198241,0.000950544},
		      {0.486033,-0.0162104,0.0015784,-0.143655,0.0362722,0.000199701}};

double SmearCD[8]={10.049,3.07365,0.461349,1.47116,-0.832429,1.45479,1.67907,3.90999};

double Quad(double x, double A, double B, double C){
  return A + B*x + C*x*x; 
}

double DoubQuad(double x, double A, double B, double C, double D, double E, double F){
  double X = Quad(x,A,B,C)*Quad(x,A,B,C) - Quad(x,D,E,F)*Quad(x,D,E,F);
  if(X>0.01){
    return sqrt(X); 
  }
  return 0.1;
}

double Trig3(double x, double A, double B, double C, double D, double E, double F, double G){
  return A + B*sin((x*1*M_PI/180)+C) + D*sin((x*2*M_PI/180)+E) + F*sin((x*3*M_PI/180)+G); 
}

double DoubTrig3(double x, double A, double B, double C, double D, double E, double F, double G, double H){
  double X = Trig3(x,A,B,C,D,E,F,G)*Trig3(x,A,B,C,D,E,F,G) - H*H;
  if(X>0.01){
    return sqrt(X); 
  }
  return 0.1;
}

void SetLorentzVector_MomentumSimulationSmear(TLorentzVector &p4, clas12::region_particle* p){
  TVector3 v3 = p4.Vect();
  double mom = v3.Mag();
  double theta = v3.Theta()*180/M_PI;
  double phi = v3.Phi()*180/M_PI;
  double M = p4.M();        
  if(p->getRegion()==FD){
    int sector = p->getSector();
    double smear = 0.01*DoubQuad(theta,SmearFD[sector-1][0],SmearFD[sector-1][1],SmearFD[sector-1][2],SmearFD[sector-1][3],SmearFD[sector-1][4],SmearFD[sector-1][5]);
    mom*=thisRand->Gaus(1.0,smear);
  }
  else if(p->getRegion()==CD){
    double smear = 0.01*DoubTrig3(phi,SmearCD[0],SmearCD[1],SmearCD[2],SmearCD[3],SmearCD[4],SmearCD[5],SmearCD[6],SmearCD[7]);
    mom*=thisRand->Gaus(1.0,smear);
  }
  else{
    cout<<"Problem\n\n\n\n\n\n\n\n\n";
    exit(-2);
  }
  v3.SetMag(mom);
  p4.SetXYZM(v3.X(),v3.Y(),v3.Z(),M);  
}
*/

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
vector<double> bE_ThetaCD = {35,38,41,44,47,50,53,56,59,62,65,68,71,74,77,80,83,86,89};
vector<double> bE_PhiCD = {-180,-160,-140,-120,-100,-80,-60,-40,-20,0,20,40,60,80,100,120,140,160,180};
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

  TH1D * h_Q2FD_ep = new TH1D("Q2FD_ep","Q^{2} (e,e'p);Q^{2};Counts",100,0,5);
  hist_list.push_back(h_Q2FD_ep);
  TH1D * h_thetaFD_ep = new TH1D("thetaFD_ep","#theta_{q} (e,e'p);#theta_{q};Counts",100,5,55);
  hist_list.push_back(h_thetaFD_ep);
  TH1D * h_DEpFD_ep = new TH1D("DEpFD_ep","#Delta E' (e,e'p);#Delta #E';Counts",100,-0.5,0.5);
  hist_list.push_back(h_DEpFD_ep);
  TH1D * h_phidiffFD_ep = new TH1D("phidiffFD_ep","#Delta #phi (e,e'p);#Delta #phi;Counts",100,-15,15);
  hist_list.push_back(h_phidiffFD_ep);
  TH1D * h_angleFD_ep = new TH1D("angleFD_ep","#theta (e,e'p);#theta;Counts",100,0,50);
  hist_list.push_back(h_angleFD_ep);

  TH1D * h_e_corr_binSector_binTheta[6][14];
  for(int j=1; j<=6; j++){
    for(int i=0; i<14; i++){
      int min = bE_Theta[i];
      int max = bE_Theta[i+1];
      sprintf(temp_name,"ecorr_sector_%d_theta_%d",j,i);
      sprintf(temp_title,"#Delta p Sector %d (%d< #theta < %d);#Delta p [GeV];Counts",j,min,max);
      h_e_corr_binSector_binTheta[j-1][i] = new TH1D(temp_name,temp_title,50,-0.05,0.05);
      hist_list.push_back(h_e_corr_binSector_binTheta[j-1][i]);
    }
  }

  TH1D * h_p_corr_binSector_binTheta[6][14];
  for(int j=1; j<=6; j++){
    for(int i=0; i<14; i++){
      int min = bE_Theta[i];
      int max = bE_Theta[i+1];
      sprintf(temp_name,"pcorr_sector_%d_theta_%d",j,i);
      sprintf(temp_title,"#Delta p Sector %d (%d< #theta < %d);#Delta p [GeV];Counts",j,min,max);
      h_p_corr_binSector_binTheta[j-1][i] = new TH1D(temp_name,temp_title,50,-0.05,0.05);
      hist_list.push_back(h_p_corr_binSector_binTheta[j-1][i]);
    }
  }
  
  TH1D * h_corr_binThetaCD[18];
  for(int i=0; i<18; i++){
    int min = bE_ThetaCD[i];
    int max = bE_ThetaCD[i+1];
    sprintf(temp_name,"corr_theta_%d",i);
    sprintf(temp_title,"#Delta p (%d< #theta < %d);#Delta p;Counts",min,max);
    h_corr_binThetaCD[i] = new TH1D(temp_name,temp_title,50,-0.15,0.15);
    hist_list.push_back(h_corr_binThetaCD[i]);
  }

  TH1D * h_corr_binPhiCD[18];
  for(int i=0; i<18; i++){
    int min = bE_PhiCD[i];
    int max = bE_PhiCD[i+1];
    sprintf(temp_name,"corr_phi_%d",i);
    sprintf(temp_title,"#Delta p (%d< #phi < %d);#Delta p;Counts",min,max);
    h_corr_binPhiCD[i] = new TH1D(temp_name,temp_title,50,-0.35,0.35);
    hist_list.push_back(h_corr_binPhiCD[i]);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TH1D * h_e_corr_Smear_binSector_binTheta[6][14];
  for(int j=1; j<=6; j++){
    for(int i=0; i<14; i++){
      int min = bE_Theta[i];
      int max = bE_Theta[i+1];
      sprintf(temp_name,"ecorr_Smear_sector_%d_theta_%d",j,i);
      sprintf(temp_title,"#Delta p Sector %d (%d< #theta < %d);#Delta p [GeV];Counts",j,min,max);
      h_e_corr_Smear_binSector_binTheta[j-1][i] = new TH1D(temp_name,temp_title,50,-0.05,0.05);
      hist_list.push_back(h_e_corr_Smear_binSector_binTheta[j-1][i]);
    }
  }

  TH1D * h_p_corr_Smear_binSector_binTheta[6][14];
  for(int j=1; j<=6; j++){
    for(int i=0; i<14; i++){
      int min = bE_Theta[i];
      int max = bE_Theta[i+1];
      sprintf(temp_name,"pcorr_Smear_sector_%d_theta_%d",j,i);
      sprintf(temp_title,"#Delta p Sector %d (%d< #theta < %d);#Delta p [GeV];Counts",j,min,max);
      h_p_corr_Smear_binSector_binTheta[j-1][i] = new TH1D(temp_name,temp_title,50,-0.05,0.05);
      hist_list.push_back(h_p_corr_Smear_binSector_binTheta[j-1][i]);
    }
  }
  
  TH1D * h_corr_Smear_binThetaCD[18];
  for(int i=0; i<18; i++){
    int min = bE_ThetaCD[i];
    int max = bE_ThetaCD[i+1];
    sprintf(temp_name,"corr_Smear_theta_%d",i);
    sprintf(temp_title,"#Delta p (%d< #theta < %d);#Delta p;Counts",min,max);
    h_corr_Smear_binThetaCD[i] = new TH1D(temp_name,temp_title,50,-0.15,0.15);
    hist_list.push_back(h_corr_Smear_binThetaCD[i]);
  }

  TH1D * h_corr_Smear_binPhiCD[18];
  for(int i=0; i<18; i++){
    int min = bE_PhiCD[i];
    int max = bE_PhiCD[i+1];
    sprintf(temp_name,"corr_Smear_phi_%d",i);
    sprintf(temp_title,"#Delta p (%d< #phi < %d);#Delta p;Counts",min,max);
    h_corr_Smear_binPhiCD[i] = new TH1D(temp_name,temp_title,50,-0.35,0.35);
    hist_list.push_back(h_corr_Smear_binPhiCD[i]);
  }
  
  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Sumw2();
    hist_list[i]->GetXaxis()->CenterTitle();
    hist_list[i]->GetYaxis()->CenterTitle();
  }

  int counter = 0;
  //while(chain.Next())
  while(chain.Next() && (counter<10000000000))
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
	  TLorentzVector el_smear = el;
	  if(isMC==1){
	    SetLorentzVector_MomentumSimulationSmear(el_smear,electrons[0]);
	  }
	  
	  double Eprime = ( mN*beam_E )/( beam_E*(1-cos(el.Theta())) + mN );
	  TVector3 vel_fromAngle;
	  vel_fromAngle.SetMagThetaPhi(Eprime,el.Theta(),el.Phi());
	  TVector3 vp_fromAngle = vbeam-vel_fromAngle;
	  double Delta_Eprime = ( mN*beam_E )/( beam_E*(1-cos(el.Theta())) + mN ) - el.P();
	  double eres = Delta_Eprime/Eprime;

	  double Delta_Eprime_smear = ( mN*beam_E )/( beam_E*(1-cos(el.Theta())) + mN ) - el_smear.P();
	  double eres_smear = Delta_Eprime_smear/Eprime;

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

	  if(protons.size() <= 0){continue;}
	  //if(protons[0]->getRegion()!=FD){continue;}

	  GetLorentzVector_ReconVector(proton_ptr,protons[0]);
	  if(isMC==0){
	    SetLorentzVector_ThetaCorrection(proton_ptr,protons[0]);
	  }
	  SetLorentzVector_EnergyLossCorrection(proton_ptr,protons[0]);
	  if(isMC==0){
	    SetLorentzVector_MomentumCorrection(proton_ptr,protons[0]);
	  }
	  TLorentzVector proton_ptr_smear = proton_ptr;
	  if(isMC==1){
	    SetLorentzVector_MomentumSimulationSmear(proton_ptr_smear,protons[0]);
	  }
	  
	  double mom_q = q.P();
	  double phi_q = q.Phi()*180/M_PI;
	  double mom_p = proton_ptr.P();
	  double Delta_Mom_FromAngle = vp_fromAngle.Mag()-mom_p;
	  double pres = Delta_Mom_FromAngle/vp_fromAngle.Mag();

	  double mom_p_smear = proton_ptr_smear.P();
	  double Delta_Mom_FromAngle_smear = vp_fromAngle.Mag()-mom_p_smear;
	  double pres_smear = Delta_Mom_FromAngle_smear/vp_fromAngle.Mag();

	  
	  double theta_p = proton_ptr.Theta()*180/M_PI;
	  double phi_p = proton_ptr.Phi()*180/M_PI;
	  double Delta_mom = (mom_p-mom_q)/mom_q;
	  double Delta_theta = theta_p-theta_q;
	  //double Delta_phi = getPhiDiff(q.Vect(),proton_ptr.Vect());;
	  double Delta_phi = (phi_p-phi_q);
	  int sector_p = protons[0]->getSector();

	  if(Delta_phi<-180){Delta_phi+=360;}
	  else if(Delta_phi>180){Delta_phi-=360;}
	  if(fabs(Delta_Eprime)>0.15){continue;}

	  h_DEpFD_ep->Fill(Delta_Eprime,wep);
	  h_phidiffFD_ep->Fill(Delta_phi,wep);
	  h_angleFD_ep->Fill(q.Vect().Angle(proton_ptr.Vect())*180/M_PI,wep);
	  if(fabs(Delta_phi)>3){continue;}
	  //if((q.Vect().Angle(proton_ptr.Vect())*180/M_PI)>5){continue;}
	  h_Q2FD_ep->Fill(Q2,wep);	  
	  h_thetaFD_ep->Fill(theta_q,wep);

	  if(binX(bE_Theta,theta_e)!=-1){
	    h_e_corr_binSector_binTheta[sector_e-1][binX(bE_Theta,theta_e)]->Fill(eres,wep);
	    h_e_corr_Smear_binSector_binTheta[sector_e-1][binX(bE_Theta,theta_e)]->Fill(eres_smear,wep);
	  }
	    
	  if(protons[0]->getRegion()==FD){
	    if(binX(bE_Theta,theta_p)!=-1){
	      h_p_corr_binSector_binTheta[sector_p-1][binX(bE_Theta,theta_p)]->Fill(pres,wep);
	      h_p_corr_Smear_binSector_binTheta[sector_p-1][binX(bE_Theta,theta_p)]->Fill(pres_smear,wep);
	    }
	  }
	  if(protons[0]->getRegion()==CD){
	    if(binX(bE_ThetaCD,theta_p)!=-1){
	      h_corr_binThetaCD[binX(bE_ThetaCD,theta_p)]->Fill(pres,wep);
	      h_corr_Smear_binThetaCD[binX(bE_ThetaCD,theta_p)]->Fill(pres_smear,wep);
	    }
	    if(theta_p<62){
	    if(binX(bE_PhiCD,phi_p)!=-1){
	      h_corr_binPhiCD[binX(bE_PhiCD,phi_p)]->Fill(pres,wep);
	      h_corr_Smear_binPhiCD[binX(bE_PhiCD,phi_p)]->Fill(pres_smear,wep);
	    }
	    }
	  }
	  

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
  h_thetaFD_ep->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  TGraph * g_ep_sigma[6];
  TGraph * g_e_sigma[6];
  for(int j = 0; j < 6 ; j++){
    g_e_sigma[j] = new TGraph;
    myCanvas->Divide(3,4);
    for(int i = 0; i < 11; i++){
      myCanvas->cd(i+1);
      h_e_corr_binSector_binTheta[j][i]->Draw();      
      TF1 * f_thetabin = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.5,0.5,3);
      f_thetabin->SetParameter(0,h_e_corr_binSector_binTheta[j][i]->GetMaximum());
      f_thetabin->SetParameter(1,0);
      f_thetabin->SetParLimits(1,-0.075,0.075);
      f_thetabin->SetParameter(2,0.05);
      f_thetabin->SetParLimits(2,0.01,0.15);
      TFitResultPtr point = h_e_corr_binSector_binTheta[j][i]->Fit(f_thetabin,"SrBeqn","",-0.2,0.06);
      f_thetabin->Draw("SAME");
      double x = (bE_Theta[i]+bE_Theta[i+1])/2;
      if(point!=-1){
	g_e_sigma[j]->SetPoint(g_e_sigma[j]->GetN(),x,1000*point->Parameter(2));
	//g_ep_sigma[j]->SetPoint(g_ep_sigma[j]->GetN(),x,1000*point->Parameter(2));
      }
    }
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();
  }

  
  TGraph * g_p_sigma[6];
  for(int j = 0; j < 6 ; j++){
    g_p_sigma[j] = new TGraph;
    myCanvas->Divide(3,3);
    for(int i = 7; i < 14; i++){
      myCanvas->cd(i-6);
      h_p_corr_binSector_binTheta[j][i]->Draw();
      TF1 * f_thetabin = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.5,0.5,3);
      f_thetabin->SetParameter(0,h_p_corr_binSector_binTheta[j][i]->GetMaximum());
      f_thetabin->SetParameter(1,0);
      f_thetabin->SetParLimits(1,-0.15,0.15);
      f_thetabin->SetParameter(2,0.1);
      f_thetabin->SetParLimits(2,0.01,0.25);
      TFitResultPtr point = h_p_corr_binSector_binTheta[j][i]->Fit(f_thetabin,"SrBeqn","",-0.5,0.5);
      f_thetabin->Draw("SAME");
      double x = (bE_Theta[i]+bE_Theta[i+1])/2;
      if(point!=-1){
	g_p_sigma[j]->SetPoint(g_p_sigma[j]->GetN(),x,1000*point->Parameter(2));
	//g_ep_sigma[j]->SetPoint(g_ep_sigma[j]->GetN(),x,1000*point->Parameter(2));
      }
    }
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();
  }
  
  TGraph * g_sigma = new TGraph;
  myCanvas->Divide(3,3);
  for(int i = 0; i < 9; i++){
    myCanvas->cd(i+1);
    h_corr_binThetaCD[i]->Draw();
    TF1 * f_thetabin = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.5,0.5,3);
    f_thetabin->SetParameter(0,h_corr_binThetaCD[i]->GetMaximum());
    f_thetabin->SetParameter(1,0);
    f_thetabin->SetParameter(2,0.2);
    TFitResultPtr point = h_corr_binThetaCD[i]->Fit(f_thetabin,"SrBeqn","",-0.5,0.5);
    f_thetabin->Draw("SAME");
    double x = (bE_ThetaCD[i]+bE_ThetaCD[i+1])/2;
    if(point!=-1){
      g_sigma->SetPoint(g_sigma->GetN(),x,1000*point->Parameter(2));
    }
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,3);
  for(int i = 9; i < 18; i++){
    myCanvas->cd(i-8);
    h_corr_binThetaCD[i]->Draw();
    TF1 * f_thetabin = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.5,0.5,3);
    f_thetabin->SetParameter(0,h_corr_binThetaCD[i]->GetMaximum());
    f_thetabin->SetParameter(1,0);
    f_thetabin->SetParameter(2,0.2);
    TFitResultPtr point = h_corr_binThetaCD[i]->Fit(f_thetabin,"SrBeqn","",-0.5,0.5);
    f_thetabin->Draw("SAME");
    double x = (bE_ThetaCD[i]+bE_ThetaCD[i+1])/2;
    if(point!=-1){
      g_sigma->SetPoint(g_sigma->GetN(),x,1000*point->Parameter(2));
    }
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  double x_ab1[2] = {7,45};
  double y_ab1[2] = {0,200};  
  TGraph * r_ab1 = new TGraph(2,x_ab1,y_ab1);
  r_ab1->SetLineColor(0);
  r_ab1->SetTitle("#sigma_{#Delta p} vs. #theta_{e} (e,e'p);#theta_{e}^{#circ};#sigma_{#Delta p} [MeV]");
  myCanvas->Divide(2,3);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    r_ab1->Draw();
    g_e_sigma[i]->Draw("SAME");
    g_p_sigma[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  /*
  myCanvas->Divide(2,3);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    r_ab1->Draw();
    g_ep_sigma[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  */
  myCanvas->Divide(1,1);
  g_sigma->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(3,3);
  for(int i = 9; i < 18; i++){
    myCanvas->cd(i-8);
    h_corr_binThetaCD[i]->Draw();
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  f->Close();


  return 0;
}
