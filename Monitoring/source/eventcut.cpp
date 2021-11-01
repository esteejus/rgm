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

  for(int i = e_cuts; i != fake; i++){
    cutName thisCut = static_cast<cutName>(i);
    cutmap[thisCut] = cutStruct;    
  }

  //Set the cuts from the text file
  set_cuts(filename);
  std::cout<<"\n\n\n\n Event selection class created from file: "<< filename <<"\n\n\n";
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

std::string eventcut::getCutName(cutName thisCut)
{

  std::string myCutName;
  switch(thisCut)
    {
    case e_cuts:
      myCutName = "(e,e') Cuts";
      break;
    case e_nphe:
      myCutName = "#e_{HTCC}^{-}";
      break;
    case e_calv:
      myCutName = "V_{e,cal}";
      break;
    case e_calw:
      myCutName = "W_{e,cal}";
      break;
    case e_SF:
      myCutName = "SF_{e,cal}";
      break;
    case e_mom:
      myCutName = "p_{e}";
      break;
    case e_vtze:
      myCutName = "Vertex Z_{e}";
      break;
    case l_cuts:
      myCutName = "(e,e'N_{Lead}) Cuts";
      break;
    case l_pid:
      myCutName = "Lead PID";
      break;
    case l_scint:
      myCutName = "Lead Scintillator";
      break;
    case l_theta:
      myCutName = "#theta_{Lead}";
      break;
    case l_thetalq:
      myCutName = "#theta_{Lead,q}";
      break;
    case l_chipid:
      myCutName = "Lead #chi^{2}_{PID}";
      break;
    case l_vtzdiff:
      myCutName = "Vertex Z_{e} - Z_{Lead}";
      break;
    case l_phidiff:
      myCutName = "|#phi_{e} - #phi_{Lead}|";
      break;
    case lsrc_cuts:
      myCutName = "(e,e'N_{Lead,SRC}) Cuts";
      break;
    case lsrc_Q2:
      myCutName = "Q^{2}";
      break;
    case lsrc_xB:
      myCutName = "x_{B}";
      break;
    case lsrc_pmiss:
      myCutName = "p_{miss}";
      break;
    case lsrc_mmiss:
      myCutName = "m_{miss}";
      break;
    case lsrc_loq:
      myCutName = "p/q";
      break;
    case rsrc_cuts:
      myCutName = "(e,e'N_{Lead,SRC},N_{Recoil,SRC}) Cuts";
      break;
    case rsrc_pid:
      myCutName = "Recoil PID";
      break;
    case rsrc_mom:
      myCutName = "p_{Recoil}";
      break;
    case rsrc_chipid:
      myCutName = "Recoil #chi^{2}_{PID}";
      break;
    default:
      myCutName = "Unknown Cut";
      break;
    }
  return myCutName;
}

std::string eventcut::getCutInformation(cutName thisCut)
{

  std::string min = std::to_string(cutmap[thisCut].min);
  min.erase(min.find_last_not_of('0')+1,std::string::npos );
  std::string max = std::to_string(cutmap[thisCut].max);
  max.erase(max.find_last_not_of('0')+1,std::string::npos );
  std::string cutonoff = (cutmap[thisCut].docut) ? "ON" : "OFF";

  std::string myCutInformation;
  switch(thisCut)
    {
    case e_nphe:
    case e_calv:
    case e_calw:
    case e_SF:
    case e_mom:
    case e_vtze:
    case l_theta:
    case l_thetalq:
    case l_chipid:
    case l_vtzdiff:
    case l_phidiff:
    case lsrc_Q2:
    case lsrc_xB:
    case lsrc_pmiss:
    case lsrc_mmiss:
    case lsrc_loq:
    case rsrc_mom:
    case rsrc_chipid:	  
      myCutInformation = getCutName(thisCut)+": min="+min+" max="+max;
      break;
    case l_pid:
    case rsrc_pid:
      myCutInformation = getCutName(thisCut)+": "+std::to_string(cutmap[thisCut].count);
      break;
    case l_scint:
      myCutInformation = getCutName(thisCut)+": "+cutmap[thisCut].label;
      break;
    case e_cuts:
    case l_cuts:
    case lsrc_cuts:
    case rsrc_cuts:
      myCutInformation = getCutName(thisCut)+": "+cutonoff;
      break;
    default:
      break;
    }
  return myCutInformation;
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
	case e_vtze:
	case l_theta:
	case l_thetalq:
	case l_chipid:
	case l_vtzdiff:
	case l_phidiff:
	case lsrc_Q2:
	case lsrc_xB:
	case lsrc_pmiss:
	case lsrc_mmiss:
	case lsrc_loq:
	case rsrc_mom:
	case rsrc_chipid:
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
	case e_cuts:
	case l_cuts:
	case lsrc_cuts:
	case rsrc_cuts:
	  if(cut_values == "ON"){
	    cutStruct.docut = true;
	  }
	  else if(cut_values == "OFF"){
	    cutStruct.docut = false;
	  }
	  else{
	    cutStruct.docut = false;
	  } 
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
  if(cutmap[e_cuts].docut){
    print_cut_loop(e_cuts,l_cuts);
  }
  else{
    std::cout<<getCutInformation(e_cuts)<<std::endl;
  }
  std::cout<<std::endl;

  if(cutmap[l_cuts].docut){
    print_cut_loop(l_cuts,lsrc_cuts);
  }
  else{
    std::cout<<getCutInformation(l_cuts)<<std::endl;
  }
  std::cout<<std::endl;

  if(cutmap[lsrc_cuts].docut){
    print_cut_loop(lsrc_cuts,rsrc_cuts);
  }
  else{
    std::cout<<getCutInformation(lsrc_cuts)<<std::endl;
  }
  std::cout<<std::endl;
  
  if(cutmap[rsrc_cuts].docut){
    print_cut_loop(rsrc_cuts,fake);
  }
  else{
    std::cout<<getCutInformation(rsrc_cuts)<<std::endl;
  }
  std::cout<<std::endl;
  
  std::cout<<std::endl<<std::endl;
}

void eventcut::print_cut_loop(cutName startCut, cutName endCut)
{
  for(int i = startCut; i != endCut; i++){
    cutName thisCut = static_cast<cutName>(i);
    if(cutmap[thisCut].docut){
      std::cout<<getCutInformation(thisCut)<<std::endl;
    }
  }  
}

void eventcut::print_cut_onPDF(TLatex &myText, cutName thisCut, double &line)
{
  if(cutmap[thisCut].docut){
    myText.DrawLatex(0.2,line,getCutInformation(thisCut).c_str());
    line-=0.1;
  }
}

cutName eventcut::hashit(std::string cut_name)
{
  //Electron Cuts
  if(cut_name == "e_cuts"){ return e_cuts; }
  if(cut_name == "(e,e')"){ return e_cuts; }
  if(cut_name == "(e,e') Cuts"){ return e_cuts; }
  if(cut_name == "Electron Cuts"){ return e_cuts; }

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

  if(cut_name == "e_vtze"){ return e_vtze; }
  if(cut_name == "e_vtz"){ return e_vtze; }
  if(cut_name == "evtz"){ return e_vtze; }

  //Lead Nucleon Cuts
  if(cut_name == "l_cuts"){ return l_cuts; }
  if(cut_name == "(e,e'N_{Lead})"){ return l_cuts; }
  if(cut_name == "(e,e'N_{Lead}) Cuts"){ return l_cuts; }
  if(cut_name == "Lead Cuts"){ return l_cuts; }

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

  if(cut_name == "l_vtzdiff"){ return l_vtzdiff; }

  if(cut_name == "l_phidiff"){ return l_phidiff; }

  //SRC (e,e'N) Cuts
  if(cut_name == "lsrc_cuts"){ return lsrc_cuts; }
  if(cut_name == "(e,e'N_{Lead,SRC})"){ return lsrc_cuts; }
  if(cut_name == "(e,e'N_{Lead,SRC}) Cuts"){ return lsrc_cuts; }
  if(cut_name == "Lead SRC Cuts"){ return lsrc_cuts; }

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
  if(cut_name == "rsrc_cuts"){ return rsrc_cuts; }
  if(cut_name == "(e,e'N_{Lead,SRC}N_{Recoil,SRC})"){ return rsrc_cuts; }
  if(cut_name == "(e,e'N_{Lead,SRC}N_{Recoil,SRC}) Cuts"){ return rsrc_cuts; }
  if(cut_name == "Recoil SRC Cuts"){ return rsrc_cuts; }

  if(cut_name == "rsrc_pid"){ return rsrc_pid; }
  if(cut_name == "recoil_pid"){ return rsrc_pid; }

  if(cut_name == "rsrc_mom"){ return rsrc_mom; }
  if(cut_name == "recoil_momentum"){ return rsrc_mom; }

  if(cut_name == "rsrc_chipid"){ return rsrc_chipid; }

  std::cerr<<"This is an invalid cut:\n"
	   <<cut_name<<std::endl
	   <<"Aborting...\n";
  exit(-2);
}


bool eventcut::electroncut(const std::unique_ptr<clas12::clas12reader>& c12)
{
  if(!cutmap[e_cuts].docut){ return true; }
  auto electrons=c12->getByID(11);
  if(electrons.size()!=1){ return false;}
  if(!e_nphecut(c12)){ return false; }
  if(!e_calvcut(c12)){ return false; }
  if(!e_calwcut(c12)){ return false; }
  if(!e_SFcut(c12)){ return false; }
  if(!e_momcut(c12)){ return false; }
  if(!e_vtzecut(c12)){ return false; }
  return true;
}

int eventcut::leadnucleoncut(const std::unique_ptr<clas12::clas12reader>& c12)
{
  if(!cutmap[l_cuts].docut){ return 0; }
  int pid = cutmap[l_pid].count;
  auto nucleons=c12->getByID(pid);
  int num_L = 0;
  int index_L = -1;
  for(int i = 0; i < nucleons.size(); i++){
    if(!l_scintcut(c12,i)){ continue; }
    if(!l_thetacut(c12,i)){ continue; }
    if(!l_thetalqcut(c12,i)){ continue; }
    if(!l_chipidcut(c12,i)){ continue; }
    if(!l_vtzdiffcut(c12,i)){ continue; }
    if(!l_phidiffcut(c12,i)){ continue; }
    num_L++;
    index_L = i;
  }
  if(num_L != 1){return -1;}
  return index_L;
}

bool eventcut::leadSRCnucleoncut(const std::unique_ptr<clas12::clas12reader>& c12, int index_L)
{
  if(!cutmap[lsrc_cuts].docut){ return true; }
  if(!lsrc_Q2cut(c12)){ return false; }
  if(!lsrc_xBcut(c12)){ return false; }
  if(!lsrc_pmisscut(c12,index_L)){ return false; }
  if(!lsrc_mmisscut(c12,index_L)){ return false; }
  if(!lsrc_loqcut(c12,index_L)){ return false; }
  return true;
}

int eventcut::recoilSRCnucleoncut(const std::unique_ptr<clas12::clas12reader>& c12, int index_L)
{
  if(!cutmap[rsrc_cuts].docut){ return 0; }
  int pid_L = cutmap[l_pid].count;
  int pid_R = cutmap[rsrc_pid].count;
  auto nucleons=c12->getByID(pid_R);
  int num_R = 0;
  int index_R = -1;
  for(int j = 0; j < nucleons.size(); j++){
    if((index_L==j) && (pid_L==pid_R)){ continue; }
    if(!rsrc_momcut(c12,j)){ continue; }
    if(!rsrc_chipidcut(c12,j)){ continue; }
    num_R++;
    index_R = j;
  }
  if(num_R != 1){return -1;}
  return index_R;
  
}



//Electron Cuts
bool eventcut::e_nphecut(const std::unique_ptr<clas12::clas12reader>& c12)
{
  auto electrons=c12->getByID(11);
  return inRange(electrons[0]->che(clas12::HTCC)->getNphe(),e_nphe);  
}
bool eventcut::e_calvcut(const std::unique_ptr<clas12::clas12reader>& c12)
{
  auto electrons=c12->getByID(11);  
  return inRange(electrons[0]->cal(clas12::PCAL)->getLv(),e_calv);  
}
bool eventcut::e_calwcut(const std::unique_ptr<clas12::clas12reader>& c12)
{
  auto electrons=c12->getByID(11);
  return inRange(electrons[0]->cal(clas12::PCAL)->getLw(),e_calw);  
}
bool eventcut::e_SFcut(const std::unique_ptr<clas12::clas12reader>& c12)
{
  auto electrons=c12->getByID(11);

  TVector3 ve;
  ve.SetMagThetaPhi(electrons[0]->getP(),electrons[0]->getTheta(),electrons[0]->getPhi());
  double EoP_e =  (electrons[0]->cal(clas12::PCAL)->getEnergy() +  electrons[0]->cal(clas12::ECIN)->getEnergy() +  electrons[0]->cal(clas12::ECOUT)->getEnergy()) / ve.Mag();  
  return inRange(EoP_e,e_SF);  
}
bool eventcut::e_momcut(const std::unique_ptr<clas12::clas12reader>& c12)
{
  auto electrons=c12->getByID(11);
  return inRange(electrons[0]->getP(),e_mom);  
}
bool eventcut::e_vtzecut(const std::unique_ptr<clas12::clas12reader>& c12)
{
  auto electrons=c12->getByID(11);
  return inRange(electrons[0]->par()->getVz(),e_vtze);  
}



//Lead Nucleon Cuts
bool eventcut::l_scintcut(const std::unique_ptr<clas12::clas12reader>& c12, int i)
{
  if(!cutmap[l_scint].docut){ return true; }

  auto leadnucleons=c12->getByID(cutmap[l_pid].count);

  bool FTOF1A = (leadnucleons[i]->sci(clas12::FTOF1A)->getDetector() == 12);
  bool FTOF1B = (leadnucleons[i]->sci(clas12::FTOF1B)->getDetector() == 12);
  bool FTOF2 = (leadnucleons[i]->sci(clas12::FTOF2)->getDetector() == 12);
  bool CTOF = (leadnucleons[i]->sci(clas12::CTOF)->getDetector() == 4);

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
  auto leadnucleons=c12->getByID(cutmap[l_pid].count);
  double theta = leadnucleons[i]->getTheta() * 180 / M_PI;
  return inRange(theta,l_theta);  
}
bool eventcut::l_thetalqcut(const std::unique_ptr<clas12::clas12reader>& c12, int i)
{
  auto electrons=c12->getByID(11);
  TVector3 ve;
  ve.SetMagThetaPhi(electrons[0]->getP(),electrons[0]->getTheta(),electrons[0]->getPhi());

  auto leadnucleons=c12->getByID(cutmap[l_pid].count);
  TVector3 vL;
  vL.SetMagThetaPhi(leadnucleons[i]->getP(),leadnucleons[i]->getTheta(),leadnucleons[i]->getPhi());
  
  TVector3 vq = vbeam - ve;
  
  double thetalq = vq.Angle(vL) * 180 / M_PI;
  return inRange(thetalq,l_thetalq);  
}
bool eventcut::l_chipidcut(const std::unique_ptr<clas12::clas12reader>& c12, int i)
{
  auto leadnucleons=c12->getByID(cutmap[l_pid].count);
  return inRange(leadnucleons[i]->par()->getChi2Pid(),l_chipid);  
}
bool eventcut::l_vtzdiffcut(const std::unique_ptr<clas12::clas12reader>& c12, int i)
{
  auto electrons=c12->getByID(11);
  double vtze = electrons[0]->par()->getVz();  
  auto leadnucleons=c12->getByID(cutmap[l_pid].count);
  double vtzl = leadnucleons[i]->par()->getVz();  
  return inRange(vtze-vtzl,l_vtzdiff);  
}
bool eventcut::l_phidiffcut(const std::unique_ptr<clas12::clas12reader>& c12, int i)
{
  auto electrons=c12->getByID(11);
  double e_phi = electrons[0]->getPhi() * 180 / M_PI;
  auto leadnucleons=c12->getByID(cutmap[l_pid].count);
  double p_phi = leadnucleons[i]->getPhi() * 180 / M_PI;
  double phidiff;

  if(e_phi>p_phi){
    if((e_phi-p_phi)<=180){
      phidiff = (e_phi-p_phi);
    }
    else{
      phidiff = 360 - (e_phi-p_phi);
    }
  }
  else{
    if((p_phi-e_phi)<=180){
      phidiff = (p_phi-e_phi);
    }
    else{
      phidiff = 360 - (p_phi-e_phi);
    }
  }

  return inRange(phidiff,l_phidiff);  
}


//SRC (e,e'N) Cuts
bool eventcut::lsrc_Q2cut(const std::unique_ptr<clas12::clas12reader>& c12)
{
  auto electrons=c12->getByID(11);
  TVector3 ve;
  ve.SetMagThetaPhi(electrons[0]->getP(),electrons[0]->getTheta(),electrons[0]->getPhi());
  TVector3 vq = vbeam - ve;
  double nu = Ebeam - ve.Mag();
  double Q2 = vq.Mag2() - (nu*nu);

  return inRange(Q2,lsrc_Q2);  
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
  
  return inRange(xB,lsrc_xB);  
}

bool eventcut::lsrc_pmisscut(const std::unique_ptr<clas12::clas12reader>& c12, int i)
{
  auto electrons=c12->getByID(11);
  TVector3 ve;
  ve.SetMagThetaPhi(electrons[0]->getP(),electrons[0]->getTheta(),electrons[0]->getPhi());

  auto leadnucleons=c12->getByID(cutmap[l_pid].count);
  TVector3 vL;
  vL.SetMagThetaPhi(leadnucleons[i]->getP(),leadnucleons[i]->getTheta(),leadnucleons[i]->getPhi());
  
  TVector3 vmiss = vbeam - ve - vL;

  return inRange(vmiss.Mag(),lsrc_pmiss);  
}

bool eventcut::lsrc_mmisscut(const std::unique_ptr<clas12::clas12reader>& c12, int i)
{
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

  return inRange(mmiss,lsrc_mmiss);  
}

bool eventcut::lsrc_loqcut(const std::unique_ptr<clas12::clas12reader>& c12, int i)
{
  auto electrons=c12->getByID(11);
  TVector3 ve;
  ve.SetMagThetaPhi(electrons[0]->getP(),electrons[0]->getTheta(),electrons[0]->getPhi());

  auto leadnucleons=c12->getByID(cutmap[l_pid].count);
  TVector3 vL;
  vL.SetMagThetaPhi(leadnucleons[i]->getP(),leadnucleons[i]->getTheta(),leadnucleons[i]->getPhi());
  
  TVector3 vq = vbeam - ve;
  double Loq = vL.Mag()/vq.Mag();
  
  return inRange(Loq,lsrc_loq);
}

//SRC (e,e'NN) Cuts
bool eventcut::rsrc_momcut(const std::unique_ptr<clas12::clas12reader>& c12, int j)
{
  auto recoilnucleons=c12->getByID(cutmap[rsrc_pid].count);
  return inRange(recoilnucleons[j]->getP(),rsrc_mom);
}
bool eventcut::rsrc_chipidcut(const std::unique_ptr<clas12::clas12reader>& c12, int j)
{
  auto recoilnucleons=c12->getByID(cutmap[rsrc_pid].count);
  return inRange(recoilnucleons[j]->par()->getChi2Pid(),rsrc_chipid);
}



//General Cut
bool eventcut::inRange(double x, cutName thisCut)
{
  if(!cutmap[thisCut].docut){ return true; }

  if(x < cutmap[thisCut].min){ return false; }
  if(x > cutmap[thisCut].max){ return false; }
  return true;

}
