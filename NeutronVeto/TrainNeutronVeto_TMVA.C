/// \file
/// \ingroup tutorial_tmva
/// \notebook -nodraw
/// This macro provides examples for the training and testing of the
/// TMVA classifiers.
///
/// As input data is used a toy-MC sample consisting of four Gaussian-distributed
/// and linearly correlated input variables.
/// The methods to be used can be switched on and off by means of booleans, or
/// via the prompt command, for example:
///
///     root -l ./TrainNeutronVeto_TMVA.C\(\"Fisher,Likelihood\"\)
///
/// (note that the backslashes are mandatory)
/// If no method given, a default set of classifiers is used.
/// The output file "TMVA.root" can be analysed with the use of dedicated
/// macros (simply say: root -l <macro.C>), which can be conveniently
/// invoked through a GUI that will appear at the end of the run of this macro.
/// Launch the GUI via the command:
///
///     root -l ./TMVAGui.C
///
/// You can also compile and run the example with the following commands
///
///     make
///     ./TrainNeutronVeto_TMVA <Methods>
///
/// where: `<Methods> = "method1 method2"` are the TMVA classifier names
/// example:
///
///     ./TrainNeutronVeto_TMVA Fisher LikelihoodPCA BDT
///
/// If no method given, a default set is of classifiers is used
///
/// - Project   : TMVA - a ROOT-integrated toolkit for multivariate data analysis
/// - Package   : TMVA
/// - Root Macro: TrainNeutronVeto_TMVA
///
/// \macro_output
/// \macro_code
/// \author Andreas Hoecker


#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"


#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"



int TrainNeutronVeto_TMVA( TString myMethodList = "" )
{
   // Methods to be processed can be given as an argument; use format:
   //
   //     mylinux~> root -l TrainNeutronVeto_TMVA.C\(\"myMethod1,myMethod2,myMethod3\"\)

   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // Cut optimisation
   Use["Cuts"]            = 1;
   //
   // 1-dimensional likelihood ("naive Bayes estimator")
   Use["Likelihood"]      = 1;
   //
   // Linear Discriminant Analysis
   Use["LD"]              = 1; // Linear Discriminant identical to Fisher
   Use["Fisher"]          = 0;
   //
   // Neural Networks (all are feed-forward Multilayer Perceptrons)
   Use["MLP"]             = 0; // Recommended ANN
   //
   // Boosted Decision Trees
   Use["BDT"]             = 1; // uses Adaptive Boost
   Use["BDTG"]            = 1; // uses Gradient Boost
   //
   // Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"]         = 1;

   // Categorized methods
   Use["MLPcat"]          = 1; // what does 0/1 mean???


   // other models
   Use["PDERS"]           = 1;
   Use["SVM"]             = 1;
   Use["KNN"]             = 1;
   Use["HMatrix"]         = 1;


   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TrainNeutronVeto_TMVA" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return 1;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // Here the preparation phase begins


/*   // simulation
   TChain * chS = new TChain("T");
   TChain * chB = new TChain("T");
   TString inputDirectory = "/lustre19/expphy/volatile/clas12/users/erins/neutron-veto/sim_rgm_bknd/ana_root/";
   for (int i=1; i<100; i++)
   {
     TString filenameS = inputDirectory + "ana2_neutrons_" + i + ".root";
     chS->Add(filenameS.Data());
     TString filenameB = inputDirectory + "ana2_protons_" + i + ".root";
     chB->Add(filenameB.Data());
   }

   TTree *signalTree     = chS;
   TTree *background     = chB;
   // use 5000 events to train*/
   

/*   // D 2 GeV
   TChain * chS = new TChain("T");
   TChain * chB = new TChain("T");
   TString inputDirectory = "/lustre19/expphy/volatile/clas12/users/erins/neutron-veto/d_2gev/2gev_root/";
   vector<int> vec = {5567,5568,5569,5570,5572,5573,5574,5575,5576,5777,5578,5579,5580,5581,5583,5586,5587,5588,5589,5590,5591,5592,5593,5594,5595,5598,5599,5600,5601,5602,5603,5604,5606,5607,5608,5609,5610,5611,5612,5613,5614,5615,5616,5617,5618,5619,5620,5622,5623,5624,5625,5626,5627};
   //vector<int> vec = {5568};
   for (int i=0; i<vec.size(); i++)
   {
     TString filenameS = inputDirectory + "goodn_pCD_01" + vec[i] + ".root";
     chS->Add(filenameS.Data());
     TString filenameB = inputDirectory + "ppipn_pCD_01" + vec[i] + ".root";
     chB->Add(filenameB.Data());
   }

   // read trees (from 2 GeV or simulation)

   TTree *signalTree     = chS;
   TTree *background     = chB;*/


   // D 6 GeV
   TChain * chS = new TChain("T");
   TChain * chB = new TChain("T");
   TString inputDirectory = "/lustre19/expphy/volatile/clas12/users/erins/neutron-veto/d_6gev/6gev_root/";
   vector<int> vec = {5045,5046,5047,5049,5050,5051,5052,5053,5054,5055,5056,5057,5058,5059,5060,5061,5062,5065,5066,5067,5072,5073,5074,5075,5077,5078,5079,5081,5082,5093,5094,5095,5096,5097,5098,5099,5100,5101,5102,5103,5104,5105,5106,5434,5435,5436,5437,5439,5441,5442,5443,5444,5445,5447,5448,5449,5450,5451,5452,5454,5455,5456};
   //vector<int> vec = {5045,5046,5047,5049,5050,5051,5052,5053,5054,5055,5056};
   for (int i=0; i<vec.size(); i++)
   {
std::cout << "NOW READING " << vec[i] << '\n';
     TString filenameS = inputDirectory + "goodn_e5_pCD_01" + vec[i] + ".root";
     chS->Add(filenameS.Data());
     TString filenameB = inputDirectory + "ppipn_e5_CD_01" + vec[i] + ".root";
     chB->Add(filenameB.Data());
   }
   TTree *signalTree     = chS;
   TTree *background     = chB;





   // read trees (from 2 GeV or simulation)


   //TTree *signalTree     = (TTree*)inputS->Get("T");
   //TTree *background     = (TTree*)inputB->Get("T");




   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName( "TMVA.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   // Create the factory object. Later you can choose the methods
   // whose performance you'd like to investigate. The factory is
   // the only TMVA object you have to interact with
   //
   // The first argument is the base of the name of all the
   // weightfiles in the directory weight/
   //
   // The second argument is the output file for the training results
   // All TMVA output can be suppressed by removing the "!" (not) in
   // front of the "Silent" argument in the option string
   TMVA::Factory *factory = new TMVA::Factory( "TrainNeutronVeto_TMVA", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

   TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset_sim_may_JUNK");
   // If you wish to modify default settings
   // (please check "src/Config.h" to see all available global options)
   //
   //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
   dataloader->AddVariable( "energy", 'F' );
   dataloader->AddVariable( "layermult", 'F' );
   dataloader->AddVariable( "size", 'F' );
   dataloader->AddVariable( "cnd_hits", 'F' );
   dataloader->AddVariable( "cnd_energy", 'F' );
   dataloader->AddVariable( "ctof_energy", 'F' );
   dataloader->AddVariable( "ctof_hits", 'F' );
   dataloader->AddVariable( "angle_diff", 'F' );

   // You can add so-called "Spectator variables", which are not used in the MVA training,
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
   // input variables, the response values of all trained MVAs, and the spectator variables

   dataloader->AddSpectator( "momentum",  "Momentum", "units", 'F' );


   // global event weights per tree (see below for setting event-wise weights)
   Double_t signalWeight     = 1.0;
   Double_t backgroundWeight = 1.0;

   // You can add an arbitrary number of signal or background trees
   dataloader->AddSignalTree    ( signalTree,     signalWeight );
   dataloader->AddBackgroundTree( background, backgroundWeight );


   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

   // Tell the dataloader how to use the training and testing events
   //
   // If no numbers of events are given, half of the events in the tree are used
   // for training, and the other half for testing:
   //
   //    dataloader->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
   //
   // To also specify the number of testing events, use:
   //
   //    dataloader->PrepareTrainingAndTestTree( mycut,
   //         "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );
   dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=5000:nTrain_Background=5000:SplitMode=Random:NormMode=NumEvents:!V" );

   // ### Book MVA methods
   //
   // Please lookup the various method configuration options in the corresponding cxx files, eg:
   // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
   // it is possible to preset ranges in the option string in which the cut optimisation should be done:
   // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

   // Cut optimisation
   if (Use["Cuts"])
      factory->BookMethod( dataloader, TMVA::Types::kCuts, "Cuts",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );

   // Test the multi-dimensional probability density estimator
   // here are the options strings for the MinMax and RMS methods, respectively:
   //
   // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
   if (Use["MLP"])
      factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=800:HiddenLayers=N+5:TestRate=5:!UseRegulator" );


   if (Use["BDT"])  // Adaptive Boost
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=4:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );


   // Decorrelated likelihood
   if (Use["LikelihoodD"])
     factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "LikelihoodD",
                            "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" );
 

   if (Use["PDERS"])
       factory->BookMethod( dataloader, TMVA::Types::kPDERS, "PDERS",
                            "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );

   if (Use["KNN"])
       factory->BookMethod( dataloader, TMVA::Types::kKNN, "KNN",
                            "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );


  // H-Matrix (chi2-squared) method
   if (Use["HMatrix"])
       factory->BookMethod( dataloader, TMVA::Types::kHMatrix, "HMatrix", "!H:!V:VarTransform=None" );


  // Fisher discriminant
   if (Use["Fisher"])
       factory->BookMethod( dataloader, TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );


  if (Use["BDTG"]) // Gradient Boost
       factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG",
                            "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );



  if (Use["RuleFit"])
       factory->BookMethod( dataloader, TMVA::Types::kRuleFit, "RuleFit",
                            "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );  //


   // --------------------------------------------------------------------------------------------------
   //  Now you can optimize the setting (configuration) of the MVAs using the set of training events
   // STILL EXPERIMENTAL and only implemented for BDT's !
   //
   //     factory->OptimizeAllMethods("SigEffAtBkg0.01","Scan");
   //     factory->OptimizeAllMethods("ROCIntegral","FitGA");
   //
   // --------------------------------------------------------------------------------------------------

   // Now you can tell the factory to train, test, and evaluate the MVAs
   //
   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TrainNeutronVeto_TMVA is done!" << std::endl;

   delete factory;
   delete dataloader;
   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );

   return 0;
}

int main( int argc, char** argv )
{
   // Select methods (don't look at this code - not of interest)
   TString methodList;
   for (int i=1; i<argc; i++) {
      TString regMethod(argv[i]);
      if(regMethod=="-b" || regMethod=="--batch") continue;
      if (!methodList.IsNull()) methodList += TString(",");
      methodList += regMethod;
   }
   return TrainNeutronVeto_TMVA(methodList);
}
