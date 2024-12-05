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

auto db=TDatabasePDG::Instance();
double mass_p = db->GetParticle(2212)->Mass();
double mass_n = db->GetParticle(2112)->Mass();
double mD = 1.8756;
double beam_E = 5.984792;
double beam_E_sigma = 0.00299;
const double me = 0.000511;
const double mU = 0.9314941024;
const double m_4He = 4.00260325415 * mU - 2*me;

const double c = 29.9792458;

void SetMom(TLorentzVector &p4,double mom){
  TVector3 v3 = p4.Vect();
  double theta = v3.Theta()*180/M_PI;
  double phi = v3.Phi()*180/M_PI;
  double M = p4.M();        
  v3.SetMag(mom);
  p4.SetXYZM(v3.X(),v3.Y(),v3.Z(),M);  
}

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
vector<double> bE_MomFD = {0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0};
vector<double> bE_MomCD = {0.5,1.0,1.3,1.6,2.0,2.5,3.0};
vector<double> bE_ThetaCD = {38,41,44,47,50,53,56,59,62,65,68,71,74,77,80,83,86,89,92};
vector<double> bE_Theta = {8,10,12,14,16,18,20,23,26,29,32,35,45};
vector<double> bE_Phi = {-35,-15,-5,0,5,10,15,25,35};
int binX(vector<double> XS, double X){
  for(int i = 0; i <= XS.size(); i++){if(X<XS[i]){return i-1;}}
  return -1;
}

double G(double x, double N, double mu, double sigma){
  return (N/(sigma*sqrt(2*M_PI))) * exp(-0.5 * sq((x-mu)/sigma)) ; 
}

double cotan(double x){ return cos(x)/sin(x);}

double E(double x, double N, double tau){
  return N * exp( x / tau) ; 
}


double getExp(TLorentzVector balance_ptr, TLorentzVector par){
  double theta_bpar = balance_ptr.Vect().Angle(par.Vect());
  double Eb = balance_ptr.E();
  double Pb = balance_ptr.P();
  double K = ((mass_n*mass_n) - balance_ptr.M2() - par.M2()) / 2;
  double a = Pb*Pb*cos(theta_bpar)*cos(theta_bpar) - Eb*Eb;
  double b = -2 * K * Pb * cos(theta_bpar);
  double c = K*K - Eb*Eb*par.M2();
  double x_min = (-b - sqrt(b*b - 4*a*c))/(2*a);
  double x_max = (-b + sqrt(b*b - 4*a*c))/(2*a);
  return x_min;
}

double getExpProton(double eMom, double eTheta, double ePhi, double pTheta, double pPhi){
  TVector3 v3e;
  v3e.SetMagThetaPhi(eMom,eTheta,ePhi);
  TLorentzVector vLe;
  vLe.SetXYZM(v3e.X(),v3e.Y(),v3e.Z(),me);

  TLorentzVector beam(0,0,beam_E,beam_E);
  TLorentzVector deut_ptr(0,0,0,mD);
  TLorentzVector balance_ptr = beam + deut_ptr - vLe;

  TVector3 vp;
  vp.SetMagThetaPhi(1,pTheta,pPhi);
  
  double theta_bpar = balance_ptr.Vect().Angle(vp);
  double Eb = balance_ptr.E();
  double Pb = balance_ptr.P();
  double K = ((mass_n*mass_n) - balance_ptr.M2() - mass_p*mass_p) / 2;
  double a = Pb*Pb*cos(theta_bpar)*cos(theta_bpar) - Eb*Eb;
  double b = -2 * K * Pb * cos(theta_bpar);
  double c = K*K - Eb*Eb*mass_p*mass_p;
  double x_min = (-b - sqrt(b*b - 4*a*c))/(2*a);
  double x_max = (-b + sqrt(b*b - 4*a*c))/(2*a);
  return x_min;
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

double ClosestPointMom(double x, double eMomP, double eThetaP, double ePhiP, double pMomP, double pThetaP, double pPhiP){
  return sqrt(SQ(x-eMomP) + SQ(getExpProton(x,eThetaP,ePhiP,pThetaP,pPhiP)-pMomP));
}

void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),rp->par()->getPz(),p4.M());
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
    
    TFitResultPtr gPoint = proj->Fit(gFit,"SrBeqn","",-0.75,0.75);
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
  TLorentzVector nucleus_ptr(0,0,0,m_4He);
  TLorentzVector deut_ptr(0,0,0,mD);
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector proton_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector ntr(0,0,0,db->GetParticle(2112)->Mass());
  reweighter newWeight(beam_E,2,2,kelly,"AV18");

  //////////////////////////////////
  char temp_name[100];
  char temp_title[100];

  
  vector<TH1*> hist_list;


  TH1D * h_mmissdiffFD = new TH1D("mmissdiffFD","m_{miss}-m_{n} (e,e'p_{FD});m_{miss}-m_{n} [GeV];Counts",100,-0.6,0.6);
  hist_list.push_back(h_mmissdiffFD);
  TH1D * h_mmissdiffCD = new TH1D("mmissdiffCD","m_{miss}-m_{n} (e,e'p_{CD});m_{miss}-m_{n} [GeV];Counts",100,-0.6,0.6);
  hist_list.push_back(h_mmissdiffCD);
  TH1D * h_mmissdiffcorrFD = new TH1D("mmissdiffcorrFD","m_{miss}-m_{n} (e,e'p_{FD});m_{miss}-m_{n} [GeV];Counts",100,-0.6,0.6);
  hist_list.push_back(h_mmissdiffcorrFD);
  TH1D * h_mmissdiffcorrCD = new TH1D("mmissdiffcorrCD","m_{miss}-m_{n} (e,e'p_{CD});m_{miss}-m_{n} [GeV];Counts",100,-0.6,0.6);
  hist_list.push_back(h_mmissdiffcorrCD);

  TH1D * h_mmissdiffCD_binThetaCD[18];
  TH1D * h_mmissdiffcorrCD_binThetaCD[18];
  for(int i=0; i<18; i++){
    int min = bE_ThetaCD[i];
    int max = bE_ThetaCD[i+1];
    sprintf(temp_name,"missdiffCD_%d",i);
    sprintf(temp_title,"m_{miss}-m_{n} (e,e'p_{CD}) (%d< #theta < %d);m_{miss}-m_{n} [GeV];Counts",min,max);
    h_mmissdiffCD_binThetaCD[i] = new TH1D(temp_name,temp_title,100,-0.6,0.6);
    hist_list.push_back(h_mmissdiffCD_binThetaCD[i]);
    
    sprintf(temp_name,"missdiffcorrCD_%d",i);
    h_mmissdiffcorrCD_binThetaCD[i] = new TH1D(temp_name,temp_title,100,-0.6,0.6);
    hist_list.push_back(h_mmissdiffcorrCD_binThetaCD[i]);
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
	  SetLorentzVector_ThetaCorrection(el,electrons[0]);
	  TLorentzVector el_corrected = el;
	  SetLorentzVector_MomentumCorrection(el_corrected,electrons[0]);

	  double Delta_Eprime = ( mN*beam_E )/( beam_E*(1-cos(el.Theta())) + mN ) - el.P();

	  	  	  
	  if(protons.size() <= 0){continue;}
	  GetLorentzVector_ReconVector(proton_ptr,protons[0]);
	  SetLorentzVector_ThetaCorrection(proton_ptr,protons[0]);
	  SetLorentzVector_EnergyLossCorrection(proton_ptr,protons[0]);
	  TLorentzVector proton_ptr_corrected = proton_ptr;
	  SetLorentzVector_MomentumCorrection(proton_ptr_corrected,protons[0]);

	  TLorentzVector miss = beam + deut_ptr - el - proton_ptr;
	  TLorentzVector miss_corrected = beam + deut_ptr - el_corrected - proton_ptr_corrected;
	
	  TLorentzVector q = beam - el;
	  double phi_q = q.Phi()*180/M_PI;
	  double phi_p = proton_ptr.Phi()*180/M_PI;
	  double theta_p = proton_ptr.Theta()*180/M_PI;
	  double Delta_phi = (phi_p-phi_q);

	  if(fabs(Delta_Eprime)>0.15){continue;}
	  if(fabs(Delta_phi)>7){continue;}
	  if(protons[0]->getRegion()==FD){
	    h_mmissdiffFD->Fill(miss.M()-mN,wep);
	    h_mmissdiffcorrFD->Fill(miss_corrected.M()-mN,wep);
	  }
	  else if(protons[0]->getRegion()==CD){
	    h_mmissdiffCD->Fill(miss.M()-mN,wep);
	    h_mmissdiffcorrCD->Fill(miss_corrected.M()-mN,wep);
	    if(binX(bE_ThetaCD,theta_p)!=-1){
	      h_mmissdiffCD_binThetaCD[binX(bE_ThetaCD,theta_p)]->Fill(miss.M()-mN,wep);
	      h_mmissdiffcorrCD_binThetaCD[binX(bE_ThetaCD,theta_p)]->Fill(miss_corrected.M()-mN,wep);	    }

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

  //////////////////////////////////////////////////

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmissdiffcorrFD->SetLineColor(2);
  h_mmissdiffcorrFD->Draw();
  h_mmissdiffFD->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmissdiffcorrCD->SetLineColor(2);
  h_mmissdiffcorrCD->Draw();
  h_mmissdiffCD->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,3);
  for(int i = 0; i < 9; i++){
    myCanvas->cd(i+1);
    h_mmissdiffcorrCD_binThetaCD[i]->SetLineColor(2);
    h_mmissdiffcorrCD_binThetaCD[i]->Draw();
    h_mmissdiffCD_binThetaCD[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,3);
  for(int i = 9; i < 18; i++){
    myCanvas->cd(i-8);
    h_mmissdiffcorrCD_binThetaCD[i]->SetLineColor(2);
    h_mmissdiffcorrCD_binThetaCD[i]->Draw();
    h_mmissdiffCD_binThetaCD[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  myCanvas->Divide(3,3);

  
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  f->Close();


  return 0;
}
