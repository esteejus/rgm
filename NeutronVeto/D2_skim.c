void D2_skim(TString inName, TString outName, bool keep_good)
{

  // charge is 0 for neutrons and 1 for protons
  TString hipo_name = "2gev_hipo/" + outName + ".hipo"; 
  clas12root::HipoChainWriter chain(hipo_name.Data()); // output hipo
  chain.Add(inName.Data()); // input hipo file
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();

  auto config_c12=chain.GetC12Reader();


  // create output root file and tree
  TString root_name = "2gev_root/" + outName + ".root";
  TFile * f = new TFile(root_name.Data(),"RECREATE");
  TTree * ntree = new TTree("T","NeutronTree");

  // create output txt file
  TString txt_name = "2gev_txt/" + outName + ".txt";
  ofstream outtxt(txt_name.Data());


  Int_t nhits;
  double px, py, pz, energy, time, path;
  //Float_t energy[100] = {-1};
  Int_t sec[100] = {-1};
  Int_t lay[100] = {-1};
  int event;
  //int sec, lay, event;
  ntree->Branch("px",&px,"momentum x/D");
  ntree->Branch("py",&py,"momentum y/D");
  ntree->Branch("pz",&pz,"momentum z/D");
  ntree->Branch("nhits",&nhits,"number of hits/I");
  ntree->Branch("sec",sec,"sec[10]/I");
  ntree->Branch("lay",lay,"lay[10]/I");
  ntree->Branch("energy",&energy,"energy/D");
  ntree->Branch("event",&event,"event/I");



  // REC::Scintillator
  auto rec_scint = config_c12->addBank("REC::Scintillator");
  auto scint_detector = config_c12->getBankOrder(rec_scint,"detector");
  auto scint_sector = config_c12->getBankOrder(rec_scint,"sector");
  auto scint_layer = config_c12->getBankOrder(rec_scint,"layer");
  auto scint_energy = config_c12->getBankOrder(rec_scint,"energy");
  
  
  // CND hits
  auto cnd_hits = config_c12->addBank("CND::hits");
  auto cnd_id = config_c12->getBankOrder(cnd_hits,"id");
  auto cnd_status = config_c12->getBankOrder(cnd_hits,"status");
  auto cnd_trkID = config_c12->getBankOrder(cnd_hits,"trkID");
  auto cnd_sector = config_c12->getBankOrder(cnd_hits,"sector");
  auto cnd_layer = config_c12->getBankOrder(cnd_hits,"layer");
  auto cnd_component = config_c12->getBankOrder(cnd_hits,"component");
  auto cnd_energy = config_c12->getBankOrder(cnd_hits,"energy");
  auto cnd_time = config_c12->getBankOrder(cnd_hits,"time");
  auto cnd_energy_unc = config_c12->getBankOrder(cnd_hits,"energy_unc");
  auto cnd_time_unc = config_c12->getBankOrder(cnd_hits,"time_unc");
  auto cnd_x = config_c12->getBankOrder(cnd_hits,"x");
  auto cnd_y = config_c12->getBankOrder(cnd_hits,"y");
  auto cnd_z = config_c12->getBankOrder(cnd_hits,"z");
  auto cnd_x_unc = config_c12->getBankOrder(cnd_hits,"x_unc");
  auto cnd_y_unc = config_c12->getBankOrder(cnd_hits,"y_unc");
  auto cnd_z_unc = config_c12->getBankOrder(cnd_hits,"z_unc");
  auto cnd_tx = config_c12->getBankOrder(cnd_hits,"tx");
  auto cnd_ty = config_c12->getBankOrder(cnd_hits,"ty");
  auto cnd_tz = config_c12->getBankOrder(cnd_hits,"tz");
  auto cnd_tlength = config_c12->getBankOrder(cnd_hits,"tlength");
  auto cnd_pathlength = config_c12->getBankOrder(cnd_hits,"pathlength");
  auto cnd_indexLadc = config_c12->getBankOrder(cnd_hits,"indexLadc");
  auto cnd_indexRadc = config_c12->getBankOrder(cnd_hits,"indexRadc");
  auto cnd_indexLtdc = config_c12->getBankOrder(cnd_hits,"indexLtdc");
  auto cnd_indexRtdc = config_c12->getBankOrder(cnd_hits,"indexRtdc");

  // CND adc/tdc
  auto cnd_adc = config_c12->addBank("CND::adc");
  auto cnd_tdc = config_c12->addBank("CND::tdc");

  // CND clusters
  auto cnd_clusters = config_c12->addBank("CND::clusters");
  auto clust_id = config_c12->getBankOrder(cnd_clusters,"id");
  auto clust_sector = config_c12->getBankOrder(cnd_clusters,"sector");
  auto clust_layer = config_c12->getBankOrder(cnd_clusters,"layer");
  auto clust_component = config_c12->getBankOrder(cnd_clusters,"component");
  auto clust_nhits = config_c12->getBankOrder(cnd_clusters,"nhits");
  auto clust_energy = config_c12->getBankOrder(cnd_clusters,"energy");
  auto clust_x = config_c12->getBankOrder(cnd_clusters,"x");
  auto clust_y = config_c12->getBankOrder(cnd_clusters,"y");
  auto clust_z = config_c12->getBankOrder(cnd_clusters,"z");
  auto clust_time = config_c12->getBankOrder(cnd_clusters,"time");
  auto clust_status = config_c12->getBankOrder(cnd_clusters,"status");
  auto clust_size = config_c12->getBankOrder(cnd_clusters,"size");

  // CTOF hits
  auto ctof_hits = config_c12->addBank("CTOF::hits");
  auto ctof_id = config_c12->getBankOrder(ctof_hits,"id");
  auto ctof_layer = config_c12->getBankOrder(ctof_hits,"layer");
  auto ctof_sector = config_c12->getBankOrder(ctof_hits,"sector");
  auto ctof_component = config_c12->getBankOrder(ctof_hits,"component");
  auto ctof_energy = config_c12->getBankOrder(ctof_hits,"energy");
  auto ctof_x = config_c12->getBankOrder(ctof_hits,"x");
  auto ctof_y = config_c12->getBankOrder(ctof_hits,"y");
  auto ctof_z = config_c12->getBankOrder(ctof_hits,"z");

  // CTOF adc/tdc
  auto ctof_adc = config_c12->addBank("CTOF::adc");
  auto ctof_tdc = config_c12->addBank("CTOF::tdc");

  // CTOF clusters
  auto ctof_clusters = config_c12->addBank("CTOF::clusters");
  auto ctof_clus_size = config_c12->getBankOrder(ctof_clusters,"size");
  auto ctof_clus_sector = config_c12->getBankOrder(ctof_clusters,"sector");
  auto ctof_clus_layer = config_c12->getBankOrder(ctof_clusters,"layer");
  auto ctof_clus_component = config_c12->getBankOrder(ctof_clusters,"component");
  auto ctof_clus_energy = config_c12->getBankOrder(ctof_clusters,"energy");
  auto ctof_clus_time = config_c12->getBankOrder(ctof_clusters,"time"); // try this one!!!
  auto ctof_clus_status = config_c12->getBankOrder(ctof_clusters,"status");


  // ScintExtras
  auto scintextras = config_c12->addBank("RECHB::ScintExtras");
  auto scint_dedx = config_c12->getBankOrder(scintextras,"dedx");
  auto scint_size = config_c12->getBankOrder(scintextras,"size");
  auto scint_layermult = config_c12->getBankOrder(scintextras,"layermult");
  

  // other banks
  auto rec_part = config_c12->addBank("REC::Particle");

 
  int counter = 0;
  auto& c12=chain.C12ref();



  // histos

  // generated momentum
  TH2D * h_px = new TH2D("px","px;generated px;reconstructed px",100,0,1.5,100,0,1.5);
    h_px->SetOption("colz");
  TH2D * h_py = new TH2D("py","py;generated py;reconstructed py",100,0,1.5,100,0,1.5);
    h_py->SetOption("colz");
  TH2D * h_pz = new TH2D("pz","pz;generated pz;reconstructed pz",100,0,1.5,100,0,1.5);
    h_pz->SetOption("colz");
  TH2D * h_p = new TH2D("p","p;generated p;reconstructed p",100,0,1.5,100,0,1.5);
    h_p->SetOption("colz");



  // proton stuff
  TH1D * h_psize = new TH1D("psize","Number of Protons in Event",10,0,10);
  TH2D * h_dbeta_p = new TH2D("dbeta_p","#Delta #beta vs proton momentum",50,0,3,50,-0.2,0.2);
    h_dbeta_p->SetOption("colz");

  // neutron stuff
  TH1D * h_nsize = new TH1D("nsize","Number of Neutrons in Event",10,0,10);


  // reconstructed momentum
  TH2D * h_pangles = new TH2D("pangles","Proton Angles;phi;theta",48,-180,180,45,0,180);
    h_pangles->SetOption("colz");
  TH2D * h_nangles = new TH2D("nangles","Neutron Angles;phi;theta",48,-180,180,45,0,180);
    h_nangles->SetOption("colz");
  TH1D * h_pxminuspx = new TH1D("pxminuspx","px_{n}-px_{miss};Counts",100,-0.5,0.5);
  TH1D * h_pyminuspy = new TH1D("pyminuspy","py_{n}-py_{miss};Counts",100,-0.5,0.5);
  TH1D * h_pzminuspz = new TH1D("pzminuspz","pz_{n}-pz_{miss};Counts",100,-0.5,0.5);
  TH1D * h_pminusp = new TH1D("pminusp","p_{n}-p_{gen};Counts",100,-0.5,0.5);
  TH2D * h_pvsp = new TH2D("pvsp","Momentum Resolution;p_{miss} (GeV/c);p_{measured} (GeV/c)",100,0,1,100,0,1);
    h_pvsp->SetOption("colz");
  TH1D * h_cos0 = new TH1D("cos0","Cosine of angle between generated and reconstructed p",50,-1.1,1.1);
  TH1D * h_hitsec = new TH1D("hitsec","CND hit sector",25,0,25);
  TH1D * h_hitlay = new TH1D("hitlay","CND hit layer",10,0,4);
  TH1D * h_cos1 = new TH1D("cos1","Cosine of angle between generated p and cluster hit",50,-1.1,1.1); 
  TH2D * h_dpp = new TH2D("dpp","Momentum Resolution;p_{generated} (GeV/c);#Delta p/p",100,0,1,100,-0.4,0.4);
    h_dpp->SetOption("colz");
  TH1D * h_energy = new TH1D("energy","Neutron Energy Deposition;Energy (MeV);Counts",100,0,50);
  TH2D * h_sec_phi = new TH2D("sec_phi","Sector vs Phi of CND hits;phi (deg);Sector",90,0,360,25,0,25);
    h_sec_phi->SetOption("colz");
  TH1D * h_mmiss = new TH1D("mmiss","Missing Mass",50,0.5,1.5);
  TH2D * h_mmiss_xb = new TH2D("mmiss_xb","Missing Mass vs x_{B}",50,0,3,50,0.5,1.5);
    h_mmiss_xb->SetOption("colz");
  TH2D * h_theta_beta = new TH2D("theta_beta","Neutron theta vs beta;#beta;#theta",50,-0.1,1.1,55,35,145);
    h_theta_beta->SetOption("colz");
  TH2D * h_p_theta = new TH2D("p_theta","Neutron Momentum vs Theta;#theta;p (GeV/c)",55,35,145,50,0,1.2);
    h_p_theta->SetOption("colz");
  TH2D * h_pmiss_thetamiss = new TH2D("pmiss_thetamiss","pmiss vs #theta_{pmiss};#theta_{pmiss};pmiss",90,0,180,50,0,1.2);
    h_pmiss_thetamiss->SetOption("colz");


  // good n / bad n set
  TH2D * h_nangles2 = new TH2D("nangles2","Neutron Angles;phi;theta",48,-180,180,45,0,180);
    h_nangles2->SetOption("colz");
  TH1D * h_pxminuspx2 = new TH1D("pxminuspx2","px_{n}-px_{miss};Counts",100,-0.5,0.5);
  TH1D * h_pyminuspy2 = new TH1D("pyminuspy2","py_{n}-py_{miss};Counts",100,-0.5,0.5);
  TH1D * h_pzminuspz2 = new TH1D("pzminuspz2","pz_{n}-pz_{miss};Counts",100,-0.5,0.5);
  TH1D * h_pminusp2 = new TH1D("pminusp2","p_{n}-p_{gen};Counts",100,-0.5,0.5);
  TH2D * h_pvsp2 = new TH2D("pvsp2","Momentum Resolution;p_{miss} (GeV/c);p_{measured} (GeV/c)",100,0,1,100,0,1);
    h_pvsp2->SetOption("colz");
  TH1D * h_cos02 = new TH1D("cos02","Cosine of angle between generated and reconstructed p",50,-1.1,1.1);
  TH2D * h_dpp2 = new TH2D("dpp2","Momentum Resolution;p_{generated} (GeV/c);#Delta p/p",100,0,1,100,-0.4,0.4);
    h_dpp2->SetOption("colz");
  TH1D * h_mmiss2 = new TH1D("mmiss2","Missing Mass",50,0.5,1.5);
  TH2D * h_mmiss_xb2 = new TH2D("mmiss_xb2","Missing Mass vs x_{B}",50,0,3,50,0.5,1.5);
    h_mmiss_xb2->SetOption("colz");
  TH1D * h_energy2 = new TH1D("energy2","Neutron Energy Deposition;Energy (MeV);Counts",100,0,50);
  TH2D * h_theta_beta2 = new TH2D("theta_beta2","Neutron theta vs beta;#beta;#theta",50,-0.1,1.1,55,35,145);
    h_theta_beta2->SetOption("colz");
  TH2D * h_p_theta2 = new TH2D("p_theta2","Neutron Momentum vs Theta;#theta;p (GeV/c)",55,35,145,50,0,1.2);
    h_p_theta2->SetOption("colz");
  TH2D * h_pmiss_thetamiss2 = new TH2D("pmiss_thetamiss2","pmiss vs #theta_{pmiss};#theta_{pmiss};pmiss",90,0,180,50,0,1.2);
    h_pmiss_thetamiss2->SetOption("colz");


const double mP = 0.93828;
const double mN = 0.939;
const double mD = 1.8756;



int numevent = 0;
  //while(chain.Next() && numevent<200)
  while(chain.Next())
  {

    // define particles
    TVector3 p(0.,0.,0.);
    TVector3 pe(0.,0.,0.);


    // identify particles from REC::Particle
    //if (!myCut.electroncut(c12)) {continue;}
    auto elec=c12->getByID(11);
    auto prot = c12->getByID(2212);
    auto neut = c12->getByID(2112);
    auto allParticles=c12->getDetParticles();
    if (elec.size()!=1) {continue;}
    if (prot.size()!=1) {continue;}
    event = c12->runconfig()->getEvent() << '\n';

    // reject particles with the wrong PID
    bool trash = 0;
    for (int i=0; i<allParticles.size(); i++)
    {
      int pid = allParticles[i]->par()->getPid();
      if (pid!=2112 && pid!=11 && pid!=2212) {trash=1;}
    }
    if (trash==1) {continue;}


    numevent = numevent + 1;

    // electron momentum
    double Ebeam = 2.07052;
    double pe_x = elec[0]->par()->getPx();
    double pe_y = elec[0]->par()->getPy();
    double pe_z = elec[0]->par()->getPz();
    pe.SetXYZ(pe_x,pe_y,pe_z);
    double vze = elec[0]->par()->getVz();
    TVector3 pb(0,0,Ebeam);
    TVector3 pq = pb - pe;
    double nu = Ebeam - pe.Mag();
    double QSq = pq.Mag() - (nu*nu);
    double xB = QSq / (2*mN*nu);



  double starttime = c12->event()->getStartTime();

  
  // FIND PROTON
  // using skimmed data, so I don't need to apply theta cuts on the protons
  h_psize->Fill(prot.size());
  int p_index = -1;
  TVector3 pp;
  for (int i=0; i<prot.size(); i++)
  {
    pp.SetMagThetaPhi(prot[i]->getP(),prot[i]->getTheta(),prot[i]->getPhi());
    double dbeta = prot[i]->par()->getBeta() - pp.Mag()/sqrt(pp.Mag2()+mP*mP);
    h_dbeta_p->Fill(pp.Mag(),dbeta);
    double vzp = prot[i]->par()->getVz();
    double chipid = prot[i]->par()->getChi2Pid();
    if ((vzp-vze)<-4. || (vzp-vze)>4.) {continue;}
    if (chipid<-3. || chipid>3.) {continue;}
    //bool is_FD = (prot[i]->getRegion()==FD);
    //if (!is_FD) {continue;}
    //if (pp.Theta()*180./M_PI>35) {continue;}
    //if (pp.Theta()*180./M_PI<40 || pp.Theta()*180./M_PI>140) {continue;}

    if (dbeta<-0.05 || dbeta>0.05) {continue;}
    if (pp.Mag()<0.2 || pp.Mag()>1.0) {continue;}
    p_index = i;
  }
  // NOT YET OPTIMIZED - what do I do if there are two protons?

  if (p_index<0) {continue;}
  pp.SetMagThetaPhi(prot[p_index]->getP(),prot[p_index]->getTheta(),prot[p_index]->getPhi());



  double p_theta = pp.Theta()*180./M_PI;
  h_pangles->Fill(pp.Phi()*180./M_PI,p_theta);




  // missing momentum
  TVector3 pmiss = pq-pp;


  // PRINT BANK INFO //
  // LOOP OVER NEUTRONS
  h_nsize->Fill(neut.size());
  for (int i=0; i<neut.size(); i++)
  {

  // get neutron momentum
  double pn_x = neut[i]->par()->getPx();
  double pn_y = neut[i]->par()->getPy();
  double pn_z = neut[i]->par()->getPz();
  TVector3 pn;
  pn.SetXYZ(pn_x,pn_y,pn_z);

  // figure out what layer the hit is in
  bool is_CND1 = (neut[i]->sci(CND1)->getLayer()==1);
  bool is_CND2 = (neut[i]->sci(CND2)->getLayer()==2);
  bool is_CND3 = (neut[i]->sci(CND3)->getLayer()==3);

  int num_hits_inc = 0;
  if (is_CND1) {num_hits_inc = num_hits_inc + 1;}
  if (is_CND2) {num_hits_inc = num_hits_inc + 1;}
  if (is_CND3) {num_hits_inc = num_hits_inc + 1;}

  // put REC::Scintillator information
  int pindex, sector, layer, component, status, size, clusterid;
  char layermult;
  double time, energy, path, x, y, z, dedx;
  double beta = neut[i]->par()->getBeta();

  if (is_CND1)
  {
    pindex = neut[i]->sci(CND1)->getPindex();
    sector = neut[i]->sci(CND1)->getSector();
    layer =  neut[i]->sci(CND1)->getLayer();
    component =  neut[i]->sci(CND1)->getComponent();
    time =   neut[i]->sci(CND1)->getTime() - starttime;
    energy = neut[i]->sci(CND1)->getEnergy();
    path =   neut[i]->sci(CND1)->getPath();
    status = neut[i]->sci(CND1)->getStatus();
    x =      neut[i]->sci(CND1)->getX();
    y =      neut[i]->sci(CND1)->getY();
    z =      neut[i]->sci(CND1)->getZ();
    dedx =   neut[i]->sci(CND1)->getDedx();
    size =   neut[i]->sci(CND1)->getSize();
    //layermult = neut[i]->sci(CND1)->getLayerMultiplicity();
  }

  if (is_CND3)
  {
    pindex = neut[i]->sci(CND3)->getPindex();
    sector = neut[i]->sci(CND3)->getSector();
    layer =  neut[i]->sci(CND3)->getLayer();
    component =  neut[i]->sci(CND3)->getComponent();
    time =   neut[i]->sci(CND3)->getTime() - starttime;
    energy = neut[i]->sci(CND3)->getEnergy();
    path =   neut[i]->sci(CND3)->getPath();
    status = neut[i]->sci(CND3)->getStatus();
    x =      neut[i]->sci(CND3)->getX();
    y =      neut[i]->sci(CND3)->getY();
    z =      neut[i]->sci(CND3)->getZ();
    dedx =   neut[i]->sci(CND3)->getDedx();
    size =   neut[i]->sci(CND3)->getSize();
    //layermult = neut[i]->sci(CND3)->getLayerMultiplicity();
  }

  if (is_CND2)
  {
    pindex = neut[i]->sci(CND2)->getPindex();
    sector = neut[i]->sci(CND2)->getSector();
    layer =  neut[i]->sci(CND2)->getLayer();
    component =  neut[i]->sci(CND2)->getComponent();
    time =   neut[i]->sci(CND2)->getTime() - starttime;
    energy = neut[i]->sci(CND2)->getEnergy();
    path =   neut[i]->sci(CND2)->getPath();
    status = neut[i]->sci(CND2)->getStatus();
    x =      neut[i]->sci(CND2)->getX();
    y =      neut[i]->sci(CND2)->getY();
    z =      neut[i]->sci(CND2)->getZ();
    dedx =   neut[i]->sci(CND2)->getDedx();
    size =   neut[i]->sci(CND2)->getSize();
    //layermult = neut[i]->sci(CND2)->getLayerMultiplicity();
  }
  // PROBLEM: this gives preference to 2nd-layer hits
  if (!is_CND1 && !is_CND2 && !is_CND3)
  {
    sector = (neut[i]->sci(CTOF)->getComponent())/2; // rounded down, ctof component mapped onto cnd sector
    layer = 0;
    component = 1; // value doesn't matter - not used
    time =   neut[i]->sci(CTOF)->getTime() - starttime;
    energy = neut[i]->sci(CTOF)->getEnergy();
    path =   neut[i]->sci(CTOF)->getPath();
    //status = neut[i]->sci(CTOF)->getStatus();
    x =      neut[i]->sci(CTOF)->getX();
    y =      neut[i]->sci(CTOF)->getY();
    z =      neut[i]->sci(CTOF)->getZ();
    dedx =   energy/3; // getDedx() is 0 (true for sim - is it also true for data?)
    //size = neut[i]->sci(CTOF)->getSize();
    //layermult = neut[i]->sci(CTOF)->getLayerMultiplicity();
  }

  double Ep = sqrt(mN*mN + pp.Mag2());
  double Emiss = Ebeam + mD - pe.Mag() - Ep;
  double mmiss = sqrt((Emiss*Emiss) - pmiss.Mag2());


  // BASIC NEUTRONS CUTS
  double n_theta = pn.Theta()*180./M_PI;
  if (n_theta<40 || n_theta>140) {continue;}
  if (pn_x==0 || pn_y==0 || pn_z==0) {continue;}
  if (pn.Mag()<0.2) {continue;}
  if (energy<3) {continue;}
  //if (pmiss.Theta()*180./M_PI<40 || pmiss.Theta()*180./M_PI>140) {continue;}
  if (pmiss.Mag()<0.2 || pmiss.Mag()>1.2) {continue;}

  //if (xB<0.6) {continue;}
  double cos0 = pmiss.Dot(pn) / (pmiss.Mag()*pn.Mag());


  // fill histos
  h_nangles->Fill(pn.Phi()*180./M_PI,n_theta);
  h_cos0->Fill(pmiss.Dot(pn) / (pmiss.Mag()*pn.Mag()));
  h_pxminuspx->Fill(pn_x-pmiss.X());
  h_pyminuspy->Fill(pn_y-pmiss.Y());
  h_pzminuspz->Fill(pn_z-pmiss.Z());
  h_pminusp->Fill(pn.Mag()-pmiss.Mag());
  h_pvsp->Fill(pmiss.Mag(),pn.Mag());
  h_energy->Fill(energy);
  h_dpp->Fill(pmiss.Mag(),(pmiss.Mag()-pn.Mag())/pmiss.Mag());
  h_mmiss->Fill(mmiss);
  h_mmiss_xb->Fill(xB,mmiss);
  h_theta_beta->Fill(beta,n_theta);
  h_p_theta->Fill(n_theta,pn.Mag());
  h_pmiss_thetamiss->Fill(pmiss.Theta()*180./M_PI,pmiss.Mag());
  





  // CTOF BANK INFO
  /*double ctof_event_energy = 0;
  bool ctof_hit_match = 0;
  int ctof_hit_component = -1;

  for (int j=0; j<c12->getBank(ctof_hits)->getRows(); j++)
  {
    if (c12->getBank(ctof_hits)->getRows()==0) {std::cout << "skipping\n";continue;}
    ctof_hit_component = c12->getBank(ctof_hits)->getInt(ctof_component,j);
    ctof_event_energy = ctof_event_energy = c12->getBank(ctof_hits)->getFloat(ctof_energy,j);
    //check nonzero, just in expected range
    if (abs(sector*2-ctof_component)<2) {ctof_hit_match = 1;}
    // check if there is a ctof hit in a sector corresponding to the cnd hit
  }*/



/*  // CTOF CLUSTERS BANK INFO
  double ctofclus_event_energy = 0;
  int ctofclus_component = -1;
  int ctofclus_nearby = 0;
  double ctof_time = 0;
  for (int j=0; j<c12->getBank(ctof_clusters)->getRows(); j++)
  {
    if (c12->getBank(ctof_clusters)->getRows()==0) {continue;}
    ctofclus_component = c12->getBank(ctof_clusters)->getInt(ctof_clus_component,j);
    if (abs(sector*2-ctofclus_component)<3)
    {
      ctofclus_event_energy = ctofclus_event_energy + c12->getBank(ctof_clusters)->getFloat(ctof_clus_energy,j);
      ctofclus_nearby = ctofclus_nearby + 1;
      ctof_time = ctof_time + c12->getBank(ctof_clusters)->getFloat(ctof_clus_time,j);
    }
    else if (abs(sector*2-ctofclus_component)>44)
    {
      ctofclus_event_energy = ctofclus_event_energy + c12->getBank(ctof_clusters)->getFloat(ctof_clus_energy,j);
      ctofclus_nearby = ctofclus_nearby + 1;
      ctof_time = ctof_time + c12->getBank(ctof_clusters)->getFloat(ctof_clus_time,j);
    }
  }
  if (ctofclus_nearby>0) {ctof_time = ctof_time/ctofclus_nearby;}*/



  /*// CND CLUSTERS BANK INFO
  int hits_nearby7 = 0;
  double cluster_energy7 = 0;
  //int layer_width = 0;
  //double cluster_dedx7 = 0;
  //int sector_depth = 0;
  for (int j=0; j<c12->getBank(cnd_clusters)->getRows(); j++)
  {
    int cluster_sector = c12->getBank(cnd_clusters)->getInt(clust_sector,j);
    int cluster_layer = c12->getBank(cnd_clusters)->getInt(clust_layer,j);
    // count nearby hit and nearby energy
    if (abs(cluster_sector-sector)<4)
    {
      hits_nearby7 = hits_nearby7 + 1;
      cluster_energy7 = cluster_energy7 + c12->getBank(cnd_clusters)->getFloat(clust_energy,j);
    }
    else if (abs(cluster_sector-sector)>20)
    {
      hits_nearby7 = hits_nearby7 + 1;
      cluster_energy7 = cluster_energy7 + c12->getBank(cnd_clusters)->getFloat(clust_energy,j);
    }

    // measure width of response in neutron layer
    if (abs(cluster_sector-sector)<3 && cluster_layer==layer) { layer_width = layer_width + 1; }
    else if (abs(cluster_sector-sector)==22 && cluster_layer==layer) { layer_width = layer_width + 1; }
    else if (abs(cluster_sector-sector)==23 && cluster_layer==layer) { layer_width = layer_width + 1; }

    // count how many layers in sector have response
    if (cluster_sector==sector) { sector_depth = sector_depth + 1; }
  }
  
  // calculate average energy per hit
  //if (hits_nearby7>0) {cluster_dedx7 = cluster_energy7/hits_nearby7; }*/



  int hits_nearby7 = 0; int ctof_nearby7 = 0;
  double cluster_energy7 = 0; double ctof_energy7 = 0;
  //std::cout << "NEUTRON: " << '\t' << (is_CND1||is_CND2||is_CND3) << '\t' << sector << '\t' << layer << '\n';
  for (int j=0; j<c12->getBank(rec_scint)->getRows(); j++)
  {
    int rec_detector = c12->getBank(rec_scint)->getInt(scint_detector,j);
    if (rec_detector!=3 && rec_detector!=4) {continue;}
    int rec_sector = c12->getBank(rec_scint)->getInt(scint_sector,j);
    int rec_layer = c12->getBank(rec_scint)->getInt(scint_layer,j);
    double rec_energy = c12->getBank(rec_scint)->getFloat(scint_energy,j);
    //std::cout << rec_detector << '\t' << rec_sector << '\t' << rec_layer << '\t' << rec_energy << '\n';
    if (rec_detector==3 && (abs(rec_sector-sector)<4)) // hits in CND
    {
      hits_nearby7 = hits_nearby7 + 1;
      cluster_energy7 = cluster_energy7 + rec_energy;
    }
    else if (rec_detector==3 && (abs(rec_sector-sector)>20)) // hits in CND, boundary
    {
      hits_nearby7 = hits_nearby7 + 1;
      cluster_energy7 = cluster_energy7 + rec_energy;
    }
    else if (rec_detector==4 && (abs(rec_sector-2*sector)<2)) // hits in CTOF
    {
      ctof_nearby7 = hits_nearby7 + 1;
      ctof_energy7 = cluster_energy7 + rec_energy;
    }
    else if (rec_detector==4 && abs(rec_sector-2*sector)>44) // hits in CTOF, boundary
    {
      ctof_nearby7 = hits_nearby7 + 1;
      ctof_energy7 = cluster_energy7 + rec_energy;
    }
    //std::cout << "CND hits nearby = " << hits_nearby7 << '\t' << cluster_energy7 << '\n';
    //std::cout << "CTOF hits nearby = " << ctof_nearby7 << '\t' << ctof_energy7 << '\n';

  }
  //std::cout << '\n';


  // option 1: loop over rec::scintillator banks (probably more efficient)
  // option 2: loop over sectors near neutron, loop over all particles, count how many have correct sector/layer







  // Determine whether to write to "good neutron" or "bad neutron" file


  //bool good_N = (cos0<-0.9 && p.Mag()>0.1 && abs(pn_x-pmiss.X())/pmiss.X()<0.2 && abs(pn_y-pmiss.Y())/pmiss.Y()<0.2 && abs(pn_z-pmiss.Z())/pmiss.Z()<0.2 && abs(pn.Mag()-pmiss.Mag())/pmiss.Mag()<0.1);

  //bool good_N = (cos0<-0.8 && p.Mag()>0.1 && abs(pn_x-pmiss.X())/pmiss.X()<0.3 && abs(pn_y-pmiss.Y())/pmiss.Y()<0.3 && abs(pn_z-pmiss.Z())/pmiss.Z()<0.3 && abs(pn.Mag()-pmiss.Mag())/pmiss.Mag()<0.3);
  bool good_N = (cos0>0.9 && abs(pmiss.Mag()-pn.Mag())<0.1 && mmiss< 1.05);

  //bool bad_N =  (cos0>-0.8 || abs(pn_x-pmiss.X())>0.2 || abs(pn_y-pmiss.Y())>0.2 || abs(pn_z-pmiss.Z())>0.2 || abs(pn.Mag()-pmiss.Mag())>0.2);
  bool bad_N = (cos0<0.8 && abs(pmiss.Mag()-pn.Mag())>0.2 && mmiss>1.15);
 // cos0 histogram filled earlier!
 // shoudl cos0 be 1 or -1 with this definition?

  bool keep_this_one = keep_good ? good_N : bad_N;

  if (keep_this_one)
  {
    // all neutrons - print features
    outtxt << layer << ' ';
    outtxt << energy << ' ';
    //outtxt << dedx << ' ';
    outtxt << z << ' ';
    outtxt << beta << ' ';
    outtxt << hits_nearby7 << ' ';
    outtxt << cluster_energy7 << ' ';
    outtxt << pn.Theta()*180./M_PI << ' ';
    outtxt << ctof_energy7 << ' ';
    outtxt << ctof_nearby7 << ' ';
    //outtxt << ctof_time << ' ';
    outtxt << '\n';

  h_nangles2->Fill(pn.Phi()*180./M_PI,n_theta);
  h_cos02->Fill(pmiss.Dot(pn) / (pmiss.Mag()*pn.Mag()));
  h_pxminuspx2->Fill(pn_x-pmiss.X());
  h_pyminuspy2->Fill(pn_y-pmiss.Y());
  h_pzminuspz2->Fill(pn_z-pmiss.Z());
  h_pminusp2->Fill(pn.Mag()-pmiss.Mag());
  h_pvsp2->Fill(pmiss.Mag(),pn.Mag());
  h_dpp2->Fill(pmiss.Mag(),(pmiss.Mag()-pn.Mag())/pmiss.Mag());
  h_mmiss2->Fill(mmiss);
  h_mmiss_xb2->Fill(xB,mmiss);
  h_energy2->Fill(energy);
  h_theta_beta2->Fill(beta,n_theta);
  h_p_theta2->Fill(n_theta,pn.Mag());
  h_pmiss_thetamiss2->Fill(pmiss.Theta()*180./M_PI,pmiss.Mag());
  }

  }



    chain.WriteEvent();
    ntree->Fill();

    counter++;

  }

  std::cout << '\n' <<counter << " events counted!\n\n";



  // write histograms
  h_psize->Write();
  h_p->Write();
  h_pxminuspx->Write();
  h_pyminuspy->Write();
  h_pzminuspz->Write();
  h_pminusp->Write();
  h_pvsp->Write();
  h_pangles->Write();
  h_nangles->Write();
  h_px->Write();
  h_py->Write();
  h_pz->Write();
  h_cos0->Write();
  h_cos1->Write();
  h_energy->Write();
  h_dpp->Write();
  h_mmiss->Write();
  h_mmiss_xb->Write();
  h_dbeta_p->Write();
  h_theta_beta->Write();
  h_p_theta->Write();
  h_pmiss_thetamiss->Write();



  h_nangles2->Write();
  h_cos02->Write();
  h_pxminuspx2->Write();
  h_pyminuspy2->Write();
  h_pzminuspz2->Write();
  h_pminusp2->Write();
  h_pvsp2->Write();
  h_dpp2->Write();
  h_mmiss2->Write();
  h_mmiss_xb2->Write();
  h_energy2->Write();
  h_theta_beta2->Write();
  h_p_theta2->Write();
  h_pmiss_thetamiss2->Write();


  outtxt.close();
  ntree->Write();
  f->Close();


} // end macro
