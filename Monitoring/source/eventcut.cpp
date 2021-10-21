#include "eventcut.h"

eventcut::eventcut(double E, char * filename)
{
  //Set beam energy
  Ebeam = E;
  vbeam.SetXYZ(0,0,Ebeam);

  //Set all cuts to false
  cutInfo cutStruct;
  cutStruct.docut = false;
  cutStruct.min = -1000000;
  cutStruct.max =  1000000;
  cutStruct.count =  0;
  cutStruct.label =  "";

  for(int i = e_nphe; i != fake; i++){
    cutName thisCut = static_cast<cutName>(i);
    cutmap[thisCut] = cutStruct;    
  }

  //Set the cuts from the text file
  set_cuts(filename);
  std::cout<<"\n\n\n\n\n\n\n\n\n Event selection class created from file: "<< filename <<"\n\n\n";
}

eventcut::eventcut(double E)
{
  Ebeam = E;
  vbeam.SetXYZ(0,0,Ebeam);
}

eventcut::~eventcut()
{
}

bool eventcut::getDoCut(cutName thisCut)
{
  return cutmap[thisCut].docut;
}

double eventcut::getCutMin(cutName thisCut)
{
  return cutmap[thisCut].min;
}

double eventcut::getCutMax(cutName thisCut)
{
  return cutmap[thisCut].max;
}

int eventcut::getCutCount(cutName thisCut)
{
  return cutmap[thisCut].count;
}

std::string eventcut::getCutLabel(cutName thisCut)
{
  return cutmap[thisCut].label;
}


void eventcut::set_cuts(char * filename)
{

  //Check if the file is open
  std::ifstream filestream(filename);
  if(!filestream.is_open()){
    std::cerr<< filename <<" failed to open. Aborting...\n";
    exit(-2);
  }
  
  //Loop over all lines in the file
  std::string line;
  while(getline(filestream,line))
    {
      //Parse the line into a cutName and a string with the cut values
      cutInfo cutStruct;
      cutStruct.docut=true;
      cutStruct.min = -1000000;
      cutStruct.max =  1000000;
      cutStruct.count =  0;
      cutStruct.label =  "";

      cutName thisCut = hashit(line.substr(0,line.find(": ")));
      std::string cut_values =line.erase(0, line.find(": ") + 2);	 

      //Most cuts are min and max
      //PID cuts are set as an integer
      //The TOF cut will be set with a string
      switch(thisCut)
	{
	case e_nphe:
	case e_calv:
	case e_calw:
	case e_SF:
	case e_mom:
	case l_theta:
	case l_thetalq:
	case l_chipid:
	case lsrc_Q2:
	case lsrc_xB:
	case lsrc_pmiss:
	case lsrc_mmiss:
	case lsrc_loq:
	case rsrc_mom:
	  std::string::size_type sz;
	  cutStruct.min = std::stod(cut_values,&sz);
	  cutStruct.max = std::stod(cut_values.substr(sz));
	  break;
	case l_pid:
	case rsrc_pid:
	  if(cut_values == "Proton"){
	    cutStruct.count = proton_number;
	  }
	  else if(cut_values == "Neutron"){
	    cutStruct.count = neutron_number;
	  }
	  else{
	    cutStruct.count = std::stoi(cut_values);	    
	  } 
	  break;
	case l_scint:
	  cutStruct.label = cut_values;
	  break;
	default:
	  std::cerr<<"This is an invalid cut:\n"
		   <<line.substr(0,line.find(": "))<<std::endl
		   <<"Aborting...\n";
	  exit(-2);
	  break;
	}
      cutmap[thisCut] = cutStruct;
    }
  filestream.close();
  
}

void eventcut::print_cuts()
{

  //Electron Cuts
  if(cutmap[e_nphe].docut){
    std::cout<<"#Photo-electrons: min="<<cutmap[e_nphe].min<<" max="<<cutmap[e_nphe].max<<std::endl;
  }
  if(cutmap[e_calv].docut){
    std::cout<<"Calorimeter hit position V: min="<<cutmap[e_calv].min<<" max="<<cutmap[e_calv].max<<std::endl;
  }
  if(cutmap[e_calw].docut){
    std::cout<<"Calorimeter hit position W: min="<<cutmap[e_calw].min<<" max="<<cutmap[e_calw].max<<std::endl;
  }
  if(cutmap[e_SF].docut){
    std::cout<<"Sampling Fraction: min="<<cutmap[e_SF].min<<" max="<<cutmap[e_SF].max<<std::endl;
  }
  if(cutmap[e_mom].docut){
    std::cout<<"Electron Momentum: min="<<cutmap[e_mom].min<<" max="<<cutmap[e_mom].max<<std::endl;
  }

  //Lead Nucleon Cuts
  if(cutmap[l_pid].docut){
    std::cout<<"Lead PID: pid#="<<cutmap[l_pid].count<<std::endl;
  }
  if(cutmap[l_scint].docut){
    std::cout<<"Lead particle scintillator: label="<<cutmap[l_scint].label<<std::endl;
  }
  if(cutmap[l_theta].docut){
    std::cout<<"Lead theta degrees: min="<<cutmap[l_theta].min<<" max="<<cutmap[l_theta].max<<std::endl;
  }
  if(cutmap[l_thetalq].docut){
    std::cout<<"Angle between q and lead momentum Degrees: min="<<cutmap[l_thetalq].min<<" max="<<cutmap[l_thetalq].max<<std::endl;
  }
  if(cutmap[l_chipid].docut){
    std::cout<<"Lead chi2 PID: min="<<cutmap[l_chipid].min<<" max="<<cutmap[l_chipid].max<<std::endl;
  }

  //SRC (e,e'N) Cuts
  if(cutmap[lsrc_Q2].docut){
    std::cout<<"Q2: min="<<cutmap[lsrc_Q2].min<<" max="<<cutmap[lsrc_Q2].max<<std::endl;
  }
  if(cutmap[lsrc_xB].docut){
    std::cout<<"xB: min="<<cutmap[lsrc_xB].min<<" max="<<cutmap[lsrc_xB].max<<std::endl;
  }
  if(cutmap[lsrc_pmiss].docut){
    std::cout<<"p_miss: min="<<cutmap[lsrc_pmiss].min<<" max="<<cutmap[lsrc_pmiss].max<<std::endl;
  }
  if(cutmap[lsrc_mmiss].docut){
    std::cout<<"m_miss: min="<<cutmap[lsrc_mmiss].min<<" max="<<cutmap[lsrc_mmiss].max<<std::endl;
  }
  if(cutmap[lsrc_loq].docut){
    std::cout<<"P/q: min="<<cutmap[lsrc_loq].min<<" max="<<cutmap[lsrc_loq].max<<std::endl;
  }


  //SRC (e,e'NN) Cuts
  if(cutmap[rsrc_pid].docut){
    std::cout<<"Recoil PID: pid#="<<cutmap[rsrc_pid].count<<std::endl;
  }
  if(cutmap[rsrc_mom].docut){
    std::cout<<"Recoil Nucleon Momentum: min="<<cutmap[rsrc_mom].min<<" max="<<cutmap[rsrc_mom].max<<std::endl;
  }

}

cutName eventcut::hashit(std::string cut_name)
{
  //Electron Cuts
  if(cut_name == "e_nphe"){ return e_nphe; }
  if(cut_name == "nphe"){ return e_nphe; }

  if(cut_name == "e_calv"){ return e_calv; }
  if(cut_name == "calv"){ return e_calv; }
  if(cut_name == "vcal"){ return e_calv; }

  if(cut_name == "e_calw"){ return e_calw; }
  if(cut_name == "calw"){ return e_calw; }
  if(cut_name == "wcal"){ return e_calw; }

  if(cut_name == "e_SF"){ return e_SF; }
  if(cut_name == "SF"){ return e_SF; }

  if(cut_name == "e_mom"){ return e_mom; }
  if(cut_name == "e_E"){ return e_mom; }
  if(cut_name == "E_e"){ return e_mom; }

  //Lead Nucleon Cuts
  if(cut_name == "l_pid"){ return l_pid; }
  if(cut_name == "lead_pid"){ return l_pid; }

  if(cut_name == "l_scint"){ return l_scint; }
  if(cut_name == "lead_scintillator"){ return l_scint; }

  if(cut_name == "l_theta"){ return l_theta; }
  if(cut_name == "lead_theta"){ return l_theta; }
  if(cut_name == "lead_angle"){ return l_theta; }

  if(cut_name == "l_thetalq"){ return l_thetalq; }
  if(cut_name == "theta_lq"){ return l_thetalq; }

  if(cut_name == "l_chipid"){ return l_chipid; }
  if(cut_name == "lead_chipid"){ return l_chipid; }

  //SRC (e,e'N) Cuts
  if(cut_name == "lsrc_Q2"){ return lsrc_Q2; }
  if(cut_name == "Q2"){ return lsrc_Q2; }

  if(cut_name == "lsrc_xB"){ return lsrc_xB; }
  if(cut_name == "xB"){ return lsrc_xB; }

  if(cut_name == "lsrc_pmiss"){ return lsrc_pmiss; }
  if(cut_name == "pmiss"){ return lsrc_pmiss; }

  if(cut_name == "lsrc_mmiss"){ return lsrc_mmiss; }
  if(cut_name == "mmiss"){ return lsrc_mmiss; }

  if(cut_name == "lsrc_loq"){ return lsrc_loq; }
  if(cut_name == "loq"){ return lsrc_loq; }

  //SRC (e,e'NN) Cuts
  if(cut_name == "rsrc_pid"){ return rsrc_pid; }
  if(cut_name == "recoil_pid"){ return rsrc_pid; }

  if(cut_name == "rsrc_mom"){ return rsrc_mom; }
  if(cut_name == "recoil_momentum"){ return rsrc_mom; }

  std::cerr<<"This is an invalid cut:\n"
	   <<cut_name<<std::endl
	   <<"Aborting...\n";
  exit(-2);
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
    if(!l_scintcut(c12,i)){ continue; }
    if(!l_thetacut(c12,i)){ continue; }
    if(!l_thetalqcut(c12,i)){ continue; }
    if(!l_chipidcut(c12,i)){ continue; }
    num_L++;
    index_L = i;
  }
  if(num_L != 1){return -1;}
  return index_L;
}

bool eventcut::leadSRCnucleoncut(const std::unique_ptr<clas12::clas12reader>& c12, int index_L)
{
  if(!lsrc_Q2cut(c12)){ return false; }
  if(!lsrc_xBcut(c12)){ return false; }
  if(!lsrc_pmisscut(c12,index_L)){ return false; }
  if(!lsrc_mmisscut(c12,index_L)){ return false; }
  if(!lsrc_loqcut(c12,index_L)){ return false; }
  return true;
}

int eventcut::recoilSRCnucleoncut(const std::unique_ptr<clas12::clas12reader>& c12, int index_L)
{
  int pid_L = cutmap[l_pid].count;
  int pid_R = cutmap[rsrc_pid].count;
  auto nucleons=c12->getByID(pid_R);
  int num_R = 0;
  int index_R = -1;
  for(int j = 0; j < nucleons.size(); j++){
    if((index_L==j) && (pid_L==pid_R)){ continue; }
    if(!rsrc_momcut(c12,j)){ continue; }
    num_R++;
    index_R = j;
  }
  if(num_R != 1){return -1;}
  return index_R;
  
}



//Electron Cuts
bool eventcut::e_nphecut(const std::unique_ptr<clas12::clas12reader>& c12)
{
  if(!cutmap[e_nphe].docut){ return true; }

  auto electrons=c12->getByID(11);
  if(electrons[0]->che(clas12::HTCC)->getNphe() <= cutmap[e_nphe].min){ return false; }  
  if(electrons[0]->che(clas12::HTCC)->getNphe() >= cutmap[e_nphe].max){ return false; }
  return true;
}
bool eventcut::e_calvcut(const std::unique_ptr<clas12::clas12reader>& c12)
{
  if(!cutmap[e_calv].docut){ return true; }

  auto electrons=c12->getByID(11);  
  if(electrons[0]->cal(clas12::PCAL)->getLv() <= cutmap[e_calv].min){ return false; }  
  if(electrons[0]->cal(clas12::PCAL)->getLv() >= cutmap[e_calv].max){ return false; }
  return true;
}
bool eventcut::e_calwcut(const std::unique_ptr<clas12::clas12reader>& c12)
{
  if(!cutmap[e_calw].docut){ return true; }

  auto electrons=c12->getByID(11);
  if(electrons[0]->cal(clas12::PCAL)->getLw() <= cutmap[e_calw].min){ return false; }  
  if(electrons[0]->cal(clas12::PCAL)->getLw() >= cutmap[e_calw].max){ return false; }
  return true;
}
bool eventcut::e_SFcut(const std::unique_ptr<clas12::clas12reader>& c12)
{

  if(!cutmap[e_SF].docut){ return true; }

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
  if(!cutmap[e_mom].docut){ return true; }
  auto electrons=c12->getByID(11);

  if(electrons[0]->getP() <= cutmap[e_mom].min){ return false; }  
  if(electrons[0]->getP() >= cutmap[e_mom].max){ return false; }
  return true;
}



//Lead Nucleon Cuts
bool eventcut::l_scintcut(const std::unique_ptr<clas12::clas12reader>& c12, int i)
{
  if(!cutmap[l_scint].docut){ return true; }

  auto leadnucleons=c12->getByID(cutmap[l_pid].count);

  bool FTOF1A = false;
  if(leadnucleons[i]->sci(clas12::FTOF1A)->getDetector() == 12){FTOF1A = true;}

  bool FTOF1B = false;
  if(leadnucleons[i]->sci(clas12::FTOF1B)->getDetector() == 12){FTOF1B = true;}

  bool FTOF2 = false;
  if(leadnucleons[i]->sci(clas12::FTOF2)->getDetector() == 12){FTOF2 = true;}

  bool CTOF = false;
  if(leadnucleons[i]->sci(clas12::CTOF)->getDetector() == 4){CTOF = true;}

  std::string ct = cutmap[l_scint].label;
  bool nameCorrect = false;
  if(ct.compare("FTOF1A")==0){ nameCorrect = true; }
  if(ct.compare("FTOF1B")==0){ nameCorrect = true; }
  if(ct.compare("FTOF1") ==0){ nameCorrect = true; }
  if(ct.compare("FTOF2") ==0){ nameCorrect = true; }
  if(ct.compare("FTOF")  ==0){ nameCorrect = true; }
  if(ct.compare("CTOF")  ==0){ nameCorrect = true; }
  if(ct.compare("TOF")   ==0){ nameCorrect = true; }

  if(!nameCorrect){
    std::cerr<<"Incorrect lead nucleon scintillator choice -"<< ct<<"-\n Aborting...";
    exit(-2);
  }
  if((ct.compare("FTOF1A")==0) && (FTOF1A)){return true;}
  if((ct.compare("FTOF1B")==0) && (FTOF1B)){return true;}
  if((ct.compare("FTOF1")==0)  && (FTOF1A || FTOF1B)){return true;}
  if((ct.compare("FTOF2")==0)  && (FTOF2)){return true;}
  if((ct.compare("FTOF")==0)   && (FTOF1A || FTOF1B || FTOF2)){return true;}
  if((ct.compare("CTOF")==0)   && (CTOF)){return true;}
  if((ct.compare("TOF")==0)    && (FTOF1A || FTOF1B || FTOF2 || CTOF)){return true;}
  
  return false;
}

bool eventcut::l_thetacut(const std::unique_ptr<clas12::clas12reader>& c12, int i)
{
  if(!cutmap[l_theta].docut){ return true; }

  auto leadnucleons=c12->getByID(cutmap[l_pid].count);
  double theta = leadnucleons[i]->getTheta() * 180 / M_PI;
  if(theta < cutmap[l_theta].min){ return false; }
  if(theta > cutmap[l_theta].max){ return false; }
  return true;
}
bool eventcut::l_thetalqcut(const std::unique_ptr<clas12::clas12reader>& c12, int i)
{
  if(!cutmap[l_thetalq].docut){ return true; }

  auto electrons=c12->getByID(11);
  TVector3 ve;
  ve.SetMagThetaPhi(electrons[0]->getP(),electrons[0]->getTheta(),electrons[0]->getPhi());

  auto leadnucleons=c12->getByID(cutmap[l_pid].count);
  TVector3 vL;
  vL.SetMagThetaPhi(leadnucleons[i]->getP(),leadnucleons[i]->getTheta(),leadnucleons[i]->getPhi());
  
  TVector3 vq = vbeam - ve;
  
  double thetalq = vq.Angle(vL) * 180 / M_PI;
  if(thetalq < cutmap[l_thetalq].min){ return false; }
  if(thetalq > cutmap[l_thetalq].max){ return false; }
  return true;
}

bool eventcut::l_chipidcut(const std::unique_ptr<clas12::clas12reader>& c12, int i)
{
  if(!cutmap[l_chipid].docut){ return true; }

  auto leadnucleons=c12->getByID(cutmap[l_pid].count);
  double Chi2Pid_L = leadnucleons[i]->par()->getChi2Pid();
  if(Chi2Pid_L < cutmap[l_chipid].min){ return false; }
  if(Chi2Pid_L > cutmap[l_chipid].max){ return false; }
  return true;
}



//SRC (e,e'N) Cuts
bool eventcut::lsrc_Q2cut(const std::unique_ptr<clas12::clas12reader>& c12)
{
  if(!cutmap[lsrc_Q2].docut){ return true; }

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
  if(!cutmap[lsrc_xB].docut){ return true; }

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

bool eventcut::lsrc_pmisscut(const std::unique_ptr<clas12::clas12reader>& c12, int i)
{
  if(!cutmap[lsrc_pmiss].docut){ return true; }

  auto electrons=c12->getByID(11);
  TVector3 ve;
  ve.SetMagThetaPhi(electrons[0]->getP(),electrons[0]->getTheta(),electrons[0]->getPhi());

  auto leadnucleons=c12->getByID(cutmap[l_pid].count);
  TVector3 vL;
  vL.SetMagThetaPhi(leadnucleons[i]->getP(),leadnucleons[i]->getTheta(),leadnucleons[i]->getPhi());
  
  TVector3 vmiss = vbeam - ve - vL;

  if(vmiss.Mag() < cutmap[lsrc_pmiss].min){ return false; }
  if(vmiss.Mag() > cutmap[lsrc_pmiss].max){ return false; }
  return true;
}

bool eventcut::lsrc_mmisscut(const std::unique_ptr<clas12::clas12reader>& c12, int i)
{

  if(!cutmap[lsrc_mmiss].docut){ return true; }

  auto electrons=c12->getByID(11);
  TVector3 ve;
  ve.SetMagThetaPhi(electrons[0]->getP(),electrons[0]->getTheta(),electrons[0]->getPhi());

  auto leadnucleons=c12->getByID(cutmap[l_pid].count);
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

bool eventcut::lsrc_loqcut(const std::unique_ptr<clas12::clas12reader>& c12, int i)
{

  if(!cutmap[lsrc_loq].docut){ return true; }

  auto electrons=c12->getByID(11);
  TVector3 ve;
  ve.SetMagThetaPhi(electrons[0]->getP(),electrons[0]->getTheta(),electrons[0]->getPhi());

  auto leadnucleons=c12->getByID(cutmap[l_pid].count);
  TVector3 vL;
  vL.SetMagThetaPhi(leadnucleons[i]->getP(),leadnucleons[i]->getTheta(),leadnucleons[i]->getPhi());
  
  TVector3 vq = vbeam - ve;
  double Loq = vL.Mag()/vq.Mag();
  
  if(Loq < cutmap[lsrc_loq].min){ return false; }
  if(Loq > cutmap[lsrc_loq].max){ return false; }
  return true;

}



//SRC (e,e'NN) Cuts
bool eventcut::rsrc_momcut(const std::unique_ptr<clas12::clas12reader>& c12, int j)
{

  if(!cutmap[rsrc_mom].docut){ return true; }

  auto recoilnucleons=c12->getByID(cutmap[rsrc_pid].count);
  if(recoilnucleons[j]->getP() < cutmap[rsrc_mom].min){ return false; }
  if(recoilnucleons[j]->getP() > cutmap[rsrc_mom].max){ return false; }
  return true;

}
