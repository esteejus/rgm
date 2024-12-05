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
int tBinCD(double t)
{
  int b = (t - 35) / 5;
  if(b < 0){return -1;}
  if(b > 9){return -1;}
  return b;
}

int tBinFD(double t)
{
  int b = (t - 20) / 5;
  if(b < 0){return -1;}
  if(b > 4){return -1;}
  return b;
}

int tBinFDe(double t)
{
  int b = (t - 10) / 3;
  if(b < 0){return -1;}
  if(b > 8){return -1;}
  return b;
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

//double beam_E = 5.98636;
double beam_E = 5.984792;
double beam_E_sigma = 0.00299;
const double me = 0.000511;
const double mU = 0.9314941024;
const double m_4He = 4.00260325415 * mU - 2*me;


double G(double x, double N, double mu, double sigma){
  return (N/(sigma*sqrt(2*M_PI))) * exp(-0.5 * sq((x-mu)/sigma)) ; 
}

double Quad(double x, double A, double B, double C){
  return A + B*x + C*x*x; 
}


double cotan(double x){ return cos(x)/sin(x);}

double E(double x, double N, double tau){
  return N * exp( x / tau) ; 
}

double ThetaE(double ThetaP){
  double ThetaP_rad = ThetaP*M_PI/180;
  return (180/M_PI) * 2*atan( (1/((beam_E/mass_p)+1)) * cotan(ThetaP_rad) );
}


double ThetaP(double ThetaE){
  double ThetaE_rad = ThetaE*M_PI/180;
  return (180/M_PI) * atan( (1/((beam_E/mass_p)+1)) * cotan(ThetaE_rad/2) );
}

double SQ(double x){ return x*x;}

double ClosestPoint(double x, double ThetaEp, double ThetaPp){
  return sqrt(SQ(x-ThetaEp) + SQ(ThetaP(x)-ThetaPp));
}


void getMax(TH1D * hist, TF1 * func, double & mean_fit, double & sigma_fit){
  int N = hist->GetEntries();
  if(N<20000){hist->Rebin(2);}
  else if(N<2000){hist->Rebin(4);}
  double max = hist->GetMaximum();
  double mode = hist->GetBinCenter(hist->GetMaximumBin());
  double diff = fabs(mode - hist->GetMean());
  if(mode<0){mode=0;}
  double stddev = hist->GetStdDev();

  func->SetParameter(0,max);
  func->SetParameter(1,mode);
  func->SetParLimits(1,mode-diff,mode+diff);
  func->SetParameter(2,stddev);
  func->SetParLimits(2,0.0005,stddev);
  TFitResultPtr point = hist->Fit(func,"SrBeqn","",mode-stddev,mode+0.5*stddev);
  mean_fit = point->Parameter(1);
  sigma_fit = point->Parameter(2);
}

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

void getFunction(TGraphErrors * g_mygraph, TF1 * f_myfunc, TFitResultPtr & p_mypoint){  
  f_myfunc->SetLineColor(3);
  f_myfunc->SetLineWidth(1);
  f_myfunc->SetParameter(0,0);
  f_myfunc->SetParLimits(0,-1,1);
  f_myfunc->SetParameter(1,0);
  f_myfunc->SetParLimits(1,-1,1);
  f_myfunc->SetParameter(2,0);
  f_myfunc->SetParLimits(2,-1,1);
  f_myfunc->SetParameter(3,0);
  f_myfunc->SetParLimits(3,-1,1);
  p_mypoint = g_mygraph->Fit(f_myfunc,"SrBeqn","",-40,40);
      
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
  vector<double> x_rad = {11.5,14.5,17.5,20.5,23.5,26.5,29.5,32.5,35.5};  
  vector<double> y_rad = {12.24,14.7628,14.7636,15.2256,17.0205,22.3517,23.9472,29.5068};  
  TGraph * g_rad = new TGraph;
  g_rad->SetName("g_rad");
  g_rad->SetTitle("#Delta E' Shift Due to Radiation;#theta;#Delta E'");
  for(int i = 0; i < 8; i++){
    g_rad->SetPoint(g_rad->GetN(),x_rad[i],y_rad[i]);
  }
  TF1 * f_rad = new TF1("f_rad",[&](double *x, double *p){ return Quad(x[0],p[0],p[1],p[2]); },10,40,3);
  f_rad->SetParameter(0,5.0);
  f_rad->SetParameter(1,1.0);
  f_rad->SetParameter(2,1.0);
  TFitResultPtr p_rad = g_rad->Fit(f_rad,"SrBeqn","",10,40);

  //////////////////////////////////
  char temp_name[100];
  char temp_title[100];

  vector<TH1*> hist_list;
  
  TH1D * h_E_Res[14];
  for(int i=0; i<14; i++){
    int min = bE_Theta[i];
    int max = bE_Theta[i+1];
    sprintf(temp_name,"h_E_Res_%d",i);
    sprintf(temp_title,"Counts vs. #Delta E' (%d< #theta < %d);#Delta E';Counts",min,max);
    h_E_Res[i] = new TH1D(temp_name,temp_title,100,-0.15,0.15);
  }
    
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
	  
	  //SetLorentzVectorCorrectedTheta(el_corrected,electrons[0]);
	  //SetLorentzVectorCorrectedMomentum(el_corrected,electrons[0]);


	  double Eprime = ( mN*beam_E )/( beam_E*(1-cos(el.Theta())) + mN );
	  TVector3 vel_fromAngle;
	  vel_fromAngle.SetMagThetaPhi(Eprime,el.Theta(),el.Phi());
	  TVector3 vp_fromAngle = vbeam-vel_fromAngle;
	  double Delta_Eprime = ( mN*beam_E )/( beam_E*(1-cos(el.Theta())) + mN ) - el.P();
	  ///
	  double Eprime_corrected = ( mN*beam_E )/( beam_E*(1-cos(el_corrected.Theta())) + mN );
	  TVector3 vel_fromAngle_corrected;
	  vel_fromAngle_corrected.SetMagThetaPhi(Eprime_corrected,el_corrected.Theta(),el_corrected.Phi());
	  TVector3 vp_fromAngle_corrected = vbeam-vel_fromAngle_corrected;
	  double Delta_Eprime_corrected = ( mN*beam_E )/( beam_E*(1-cos(el_corrected.Theta())) + mN ) - el_corrected.P();
	  ///
	  
	  int sector_e = electrons[0]->getSector();
	  TLorentzVector q = beam - el;
          double Q2        = -q.M2();
	  double omega = q.E();
          double xB       = Q2/(2 * mass_p * (beam.E() - el.E()) );
	  double omega_diff = omega - (Q2/(2*mN));

	  double vtz_e = electrons[0]->par()->getVz();
	  double WSq = (mN*mN) - Q2 + (2*omega*mN);
	  double W = sqrt(WSq);
	  double phi_e = el.Phi() * 180/M_PI;
	  double theta_e = el.Theta()*180/M_PI;

	  double theta_e_corrected = el_corrected.Theta()*180/M_PI;
	  double theta_q = q.Theta()*180/M_PI;
	  double shift_e = 0;
	  shift_e += (sector_e==1)?0:(sector_e==2)?60:(sector_e==3)?120:(sector_e==4 && phi_e>0)?180:(sector_e==4 && phi_e<0)?-180:(sector_e==5)?-120:(sector_e==6)?-60:0;

	  if(protons.size() <= 0){continue;}
	  GetLorentzVector_ReconVector(proton_ptr,protons[0]);
	  SetLorentzVector_EnergyLossCorrection(proton_ptr,protons[0]);

	  double mom_q = q.P();
	  double phi_q = q.Phi()*180/M_PI;
	  
	  double mom_p = proton_ptr.P();
	  double Delta_Mom_FromAngle = vp_fromAngle.Mag()-mom_p;
	  double mom_p_corrected = proton_ptr_corrected.P();
	  double Delta_Mom_FromAngle_corrected = vp_fromAngle_corrected.Mag()-mom_p_corrected;


	  double theta_p = proton_ptr.Theta()*180/M_PI;
	  double phi_p = proton_ptr.Phi()*180/M_PI;
	  
	  //double mom_t = beta * gamma * mN;
	  
	  double Delta_mom = (mom_p-mom_q)/mom_q;
	  //cout<<Delta_mom<<endl;
	  double Delta_theta = theta_p-theta_q;
	  double Delta_phi = (phi_p-phi_q);
	  if(Delta_phi<-180){Delta_phi+=360;}
	  else if(Delta_phi>180){Delta_phi-=360;}

	  if(fabs(Delta_Eprime)>0.15){continue;}
	  if(fabs(Delta_phi)>3){continue;}
	  
	    
	  if(binX(bE_Theta,theta_e)!=-1){
	    h_E_Res[binX(bE_Theta,theta_e)]->Fill(Delta_Eprime,wep);
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

  
  TGraph * g_mode = new TGraph;
  TGraph * g_sigma = new TGraph;
  TGraph * g_mode_corr = new TGraph;
  TGraph * g_sigma_corr = new TGraph;
  myCanvas->Divide(3,4);
  for(int i = 0; i < 11; i++){
    double mode, sigma, mode_corr, sigma_corr;
    
    TF1 * f_thetabin = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.15,0.15,3);
    TF1 * f_thetabin_corr = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.15,0.15,3);
    getMax(h_E_Res[i],f_thetabin,mode,sigma);

    myCanvas->cd(i+1);
    h_E_Res[i]->SetLineColor(1);
    h_E_Res[i]->Draw();      
    f_thetabin->SetLineColor(1);
    f_thetabin->Draw("SAME");
    
    double x = (bE_Theta[i]+bE_Theta[i+1])/2;
    g_mode->SetPoint(g_mode->GetN(),x,1000*mode);
    g_sigma->SetPoint(g_sigma->GetN(),x,1000*sigma);

    
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  cout<<"vector<double> x_rad = {";
  for(int i = 0; i < 11; i++){
    double x = (bE_Theta[i]+bE_Theta[i+1])/2;    
    cout<<x<<",";
  }
  cout<<"};\n";
  
  cout<<"vector<double> y_rad = {";
  for(int i = 0; i < 11; i++){
    cout<<g_mode->GetY()[i]<<",";
  }
  cout<<"};\n";

  cout<<"vector<double> e_rad = {";
  for(int i = 0; i < 11; i++){
    cout<<g_sigma->GetY()[i]<<",";
  }
  cout<<"};\n";

  
  TF1 * f_radnew = new TF1("f_radnew",[&](double *x, double *p){ return Quad(x[0],p[0],p[1],p[2]); },8,40,3);
  f_radnew->SetParameter(0,5.0);
  f_radnew->SetParameter(1,1.0);
  f_radnew->SetParameter(2,1.0);
  TFitResultPtr p_radnew = g_mode->Fit(f_radnew,"SrBeqn","",8,32);

  /*
  double x_ab1[2] = {10,35};
  double y_ab1[2] = {10,40};  
  TGraph * r_ab1 = new TGraph(2,x_ab1,y_ab1);
  r_ab1->SetLineColor(0);
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  r_ab1->Draw();
  g_mode->SetLineColor(1);
  g_mode->Draw();
  f_rad->Draw("SAME");
  f_radnew->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  double x_ab2[2] = {10,35};
  double y_ab2[2] = {0,60};  
  TGraph * r_ab2 = new TGraph(2,x_ab2,y_ab2);
  r_ab2->SetLineColor(0);
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  r_ab2->Draw();
  g_sigma->SetLineColor(1);
  g_sigma->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  */
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  f->Close();


  return 0;
}
