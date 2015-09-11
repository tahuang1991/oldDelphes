// -*- C++ -*-
//
//// Package:    DiHiggsWWAnalyzer
//// class MMC
//
// Description: [one line class summary]
// use to run MMC algorithm. 
//  Implementation:
//       [Notes on implementation]
//       */
//       //
//       // Original Author:  tao huang
//       //         Created:  Thurs, 04 06 2015
//       // $Id$
//       //
//       //
#ifndef DiHiggsWWAnalyzer_MMC_h
#define DiHiggsWWAnalyzer_MMC_h

//std lib
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "math.h"  


//root lib
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TDirectory.h"

typedef std::pair<float, float> EtaPhi;

class MMC{

    public:
    //constructor
    MMC(TLorentzVector mu1_lorentz, TLorentzVector mu2_lorentz, TLorentzVector bbbar_lorentz,TLorentzVector totjets_lorentz,
	TLorentzVector met_lorentz, TLorentzVector nu1_lorentz, TLorentzVector nu2_lorentz, TLorentzVector bbbar_genp_lorentz,
	TLorentzVector h2tohh_lorentz, int onshellMarker, bool simulation,// only for simulation 
	int ievent, bool weightfromonshellnupt_func, bool weightfromonshellnupt_hist, bool weightfromonoffshellWmass_hist, 
        int iterations, std::string RefPDFfile, bool useMET, int verbose_ =0
	);
    MMC();
    ~MMC();

    private:
      TTree *mmctree;
      TH1F MMC_h2Mass;
 
   //runMMC
   public: 
      bool runMMC();
      TH1F getMMCh2();
      TTree* getMMCTree();
      //TH1F* getMMCNeutrio_onshell1();
      //TH1F* getMMCNeutrio_onshell2();
      //TH1F* getMMCNeutrio_offshell1();
      //TH1F* getMMCNeutrio_offshell2();

   private:
      void initTree(TTree* mmctree);
        
     
      float genEtaGuass(float mean, float rms);
      float genPhiFlat();
      EtaPhi generatenu1_etaphi();
      float nu1pt_onshellW(EtaPhi nu1_etaphi, TLorentzVector* mu1lorentz, float wMass);
      bool  nulorentz_offshellW(TLorentzVector* jetlorentz, TLorentzVector* mu1lorentz, 
			       TLorentzVector* mu2lorentz, TLorentzVector* nu1lorentz, 
 			       TLorentzVector* nu2lorentz, int control, float hMass);
      bool  nulorentz_offshellW(TVector2* met, TLorentzVector* mu1lorentz, 
			       TLorentzVector* mu2lorentz, TLorentzVector* nu1lorentz, 
 			       TLorentzVector* nu2lorentz, int control, float hMass);
      bool checkSolution(TLorentzVector* jetslorentz,
                          TLorentzVector* mu1lorentz,
                          TLorentzVector* mu2lorentz,
                          TLorentzVector* nu1lorentz, int control, float hMass); 
      bool cutsCheck();
      void assignMuLorentzVec(int control);  
          
    private:
      float onshellWMassRandomWalk(float x0, float step, float random);
      float onshellWMassRandomWalk(float x0, float step, float random, TH1F* hist);
      float onshellWMassPDF(float wmass);
  
    private:
      //void readoutonshellnuptPDF(TFile*, TH1F *);
      TH1F* readoutonshellWMassPDF();
      TH1F* readoutoffshellWMassPDF();
      TH2F* readoutonoffshellWMassPDF();
      TH1F* readoutonshellnuptPDF();
 
    private:
      float weightfromhist(TH1F* pdf, float x); 
      float weightfromhist(TH2F* pdf, float x, float y, bool whole=true); 
      float weightfromonshellnupt(float nupt); 
   
    private:
      bool weightfromonshellnupt_func_;
      bool weightfromonshellnupt_hist_;
      bool weightfromonoffshellWmass_hist_;
      bool useMET_;   

    private:
      TLorentzVector calculateMET(); 
    public:
      void printTrueLorentz();
      void printMMCresult(); 

    private:
      int iev;
      int onshellMarker;
      bool simulation;
      int iterations_;
      int seed_;
      std::string RefPDFfile_;
      int verbose;

   

    private:
      TLorentzVector* mmc_mu1_lorentz;
      TLorentzVector* mmc_mu2_lorentz;
      TLorentzVector* mmc_bjets_lorentz;
      TLorentzVector* mmc_totjets_lorentz;
      TLorentzVector* nu1_lorentz_true;
      TLorentzVector* nu2_lorentz_true;
      TLorentzVector* onshellW_lorentz_true;
      TLorentzVector* offshellW_lorentz_true;
      TVector2* mmcmet_vec2;
      TLorentzVector* htoWW_lorentz_true;
      TLorentzVector* htoBB_lorentz_true;
      TLorentzVector* h2tohh_lorentz_true;
      
      TLorentzVector* mu_onshellW_lorentz;
      TLorentzVector* mu_offshellW_lorentz;
      TLorentzVector* jets_lorentz;
      TVector2* met_vec2;
      TLorentzVector* nu_onshellW_lorentz;
      TLorentzVector* nu_offshellW_lorentz;
      TLorentzVector* offshellW_lorentz;
      TLorentzVector* onshellW_lorentz;
      TLorentzVector* htoWW_lorentz;
      TLorentzVector* htoBB_lorentz;
      TLorentzVector* h2tohh_lorentz;
    private:
      //branches
      float eta_mean;
      float eta_rms;
      float eta_gen; 
      float phi_gen;
      float wmass_gen;
      float hmass_gen;
       
      
      int control;
      float weight;
      float weight1;//extra weight
      float weight2;//extra weight
      float weight3;//extra weight
 
      float mu_onshellW_Eta;
      float mu_onshellW_Phi;
      float mu_onshellW_Pt;
      float mu_onshellW_E;
      float mu_offshellW_Eta;
      float mu_offshellW_Phi;
      float mu_offshellW_Pt;
      float mu_offshellW_E;
      float nu_onshellW_Eta;
      float nu_onshellW_Phi;
      float nu_onshellW_Pt;
      float nu_onshellW_E;
      float nu_offshellW_Eta;
      float nu_offshellW_Phi;
      float nu_offshellW_Pt;
      float nu_offshellW_E;
      
      float onshellW_Eta;
      float onshellW_Phi;
      float onshellW_Pt;
      float onshellW_E;
      float onshellW_Mass;
      float offshellW_Eta;
      float offshellW_Phi;
      float offshellW_Pt;
      float offshellW_E;
      float offshellW_Mass;
     
      float htoBB_jets_Eta;
      float htoBB_jets_Phi;
      float htoBB_jets_Pt;
      float htoBB_jets_E;
      float htoBB_jets_Mass;
      float htoBB_Eta;
      float htoBB_Phi;
      float htoBB_Pt;
      float htoBB_E;
      float htoBB_Mass;
      float htoWW_Eta;
      float htoWW_Phi;
      float htoWW_Pt;
      float htoWW_E;
      float htoWW_Mass;

      float MMCmet_E;
      float MMCmet_Phi;
      float MMCmet_Px;
      float MMCmet_Py;

      float h2tohh_Eta;
      float h2tohh_Phi;
      float h2tohh_Pt;
      float h2tohh_E;
      float h2tohh_Mass;


      float met;
      float met_phi;
      float met_px;
      float met_py;

      float eta_nuoffshellW_true;
      float phi_nuoffshellW_true;
      float pt_nuoffshellW_true;
      float eta_nuonshellW_true;
      float phi_nuonshellW_true;
      float pt_nuonshellW_true;
      float mass_offshellW_true;
      float mass_onshellW_true;
      float mass_htoWW_true;
      float pt_h2tohh_true;
      float mass_h2tohh_true;

};
   






#endif
