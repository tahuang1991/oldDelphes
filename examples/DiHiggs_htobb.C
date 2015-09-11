
#include "TTree.h"
#include "TSystem.h"


class ExRootTreeReader;
class ExRootResult;
class GenParticle;
class Jet;
class Track;
class Tower;

void printGenParticle(GenParticle *genP)
{
    cout << " genP Id " << genP->PID <<" Pt " << genP->PT << " M1 "<< genP->M1<<" M2 "<< genP->M2;
    cout << " D1 "<< genP->D1 <<" D2 "<<genP->D2 << endl;
}

void printJet(Jet *jet)
{

  GenParticle *particle;
  Muon *muon;

  Track *track;
  Tower *tower;

  Jet *jet;
  TObject *object;
  TLorentzVector momentum;
      momentum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
      TRefArray constituentarray(jet->Constituents);
      TRefArray particlearray(jet->Particles);

      cout<<"Looping over jet constituents. Jet pt: "<<jet->PT<<", eta: "<<jet->Eta<<", phi: "<<jet->Phi<<endl;      

      // Loop over all jet's constituents
      for(Int_t j = 0; j < jet->Constituents.GetEntriesFast(); ++j)
      {
        object = jet->Constituents.At(j);
        // Check if the constituent is accessible
        if(object == 0) continue;

        if(object->IsA() == GenParticle::Class())
        {
          particle = (GenParticle*) object;
          cout << "    GenPart pt: " << particle->PT << ", eta: " << particle->Eta << ", phi: " << particle->Phi << endl;
          momentum += particle->P4();
        }
        else if(object->IsA() == Track::Class())
        {
          track = (Track*) object;
          cout << "    Track pt: " << track->PT << ", eta: " << track->Eta << ", phi: " << track->Phi << endl;
          momentum += track->P4();
        }
        else if(object->IsA() == Tower::Class())
        {
          tower = (Tower*) object;
          cout << "    Tower pt: " << tower->ET << ", eta: " << tower->Eta << ", phi: " << tower->Phi << endl;
          momentum += tower->P4();
        }
        else if(object->IsA() == Muon::Class())
        {
          muon = (Muon*) object;
          cout << "    Muon pt: " << muon->PT << ", eta: " << muon->Eta << ", phi: " << muon->Phi << endl;
          momentum += muon->P4();
        }
      }
      cout << " constituent sum pt:  " << momentum.Pt() <<" eta "<< momentum.Eta()  <<"  phi " << momentum.Phi() << std::endl;


      for (Int_t j =0; j<jet->Particles.GetEntries();  j++){
     		GenParticle *p_tmp((GenParticle) particlearray.At(j));
		printGenParticle(p_tmp);
	}
}
//get the finalsate of genp 
void getFinalState(GenParticle* genp, TClonesArray *branchParticle)
{
     cout << "before getFinalState "; printGenParticle(genp);
       while ((genp->D1>0 && ((GenParticle*)branchParticle->At(genp->D1))->PID == genp->PID)
	    || (genp->D2 >0 && ((GenParticle*)branchParticle->At(genp->D2))->PID == genp->PID)){
 	if ( genp->D1>0 && ((GenParticle*)branchParticle->At(genp->D1))->PID == genp->PID) 
 		genp = (GenParticle*)branchParticle->At(genp->D1);
        else   genp = (GenParticle*)branchParticle->At(genp->D2);
        cout << "during getFinalState "; printGenParticle(genp);
       }

}

//------------------------------------------------------------------------------


void AnalyseEvents(ExRootTreeReader *treeReader, TTree *evtree)
{

//create branches 
  Float_t b1_px =0;
  Float_t b1_py =0;
  Float_t b1_pz =0;
  Float_t b1_energy =0;
  Float_t b2_px =0;
  Float_t b2_py=0;
  Float_t b2_pz=0;
  Float_t b2_energy=0;
  Float_t htobb_px=0;
  Float_t htobb_py=0;
  Float_t htobb_pz=0;
  Float_t htobb_energy=0;
  Float_t htobb_mass =0;
  Bool_t htobb=false; 
   
  Float_t genbjet_px=0;
  Float_t genbjet_py=0;
  Float_t genbjet_pz=0;
  Float_t genbjet_energy=0;
  Float_t genbbarjet_px=0;
  Float_t genbbarjet_py=0;
  Float_t genbbarjet_pz=0;
  Float_t genbbarjet_energy=0;
  Bool_t hasgenbjet=false;
  Bool_t hasgenbbarjet=false;
  
  Float_t dR_bjet = 2.0;
  Float_t dR_bbarjet = 2.0;
  Float_t bjet_px=0;
  Float_t bjet_py=0;
  Float_t bjet_pz=0;
  Float_t bjet_energy=0;
  Float_t bbarjet_px=0;
  Float_t bbarjet_py=0;
  Float_t bbarjet_pz=0;
  Float_t bbarjet_energy=0;
  Bool_t hasbjet=false;
  Bool_t hasbbarjet=false;

  Float_t mu1_px =0;
  Float_t mu1_py =0;
  Float_t mu1_pz =0;
  Float_t mu1_energy =0;
  Float_t mu2_px =0;
  Float_t mu2_py =0;
  Float_t mu2_pz =0;
  Float_t mu2_energy =0;
  Float_t nu1_px =0;
  Float_t nu1_py =0;
  Float_t nu1_pz =0;
  Float_t nu1_energy =0;
  Float_t nu2_px =0;
  Float_t nu2_py =0;
  Float_t nu2_pz =0;
  Float_t nu2_energy =0;
  
  Bool_t Wtomu2nu2=false;
  Bool_t Wtomu1nu1=false;
  Bool_t htoWW=false;
  
  Float_t Muon1_px = 0;
  Float_t Muon1_py = 0;
  Float_t Muon1_pz = 0;
  Float_t Muon1_energy = 0;
  Float_t Muon2_px = 0;
  Float_t Muon2_py = 0;
  Float_t Muon2_pz = 0;
  Float_t Muon2_energy = 0;
  
  Float_t htoWW_px =0;
  Float_t htoWW_py =0;
  Float_t htoWW_pz =0;
  Float_t htoWW_energy =0;
  Float_t htoWW_mass = 0;
  Float_t met = 0;
  Float_t met_phi = 0;
  Float_t met_px = 0;
  Float_t met_py = 0;

  Float_t h2tohh_mass =0;
//additional cuts
  Bool_t hasMET = false;
  Bool_t hastwomuons = false;
  Bool_t hasMuon1 = false;
  Bool_t hasMuon2 = false;
  
  evtree->Branch("b1_px",&b1_px, "b1_px/F");
  evtree->Branch("b1_py",&b1_py, "b1_py/F");
  evtree->Branch("b1_pz",&b1_pz, "b1_pz/F");
  evtree->Branch("b1_energy",&b1_energy, "b1_energy/F");
  evtree->Branch("b2_px",&b2_px, "b2_px/F");
  evtree->Branch("b2_py",&b2_py, "b2_py/F");
  evtree->Branch("b2_pz",&b2_pz, "b2_pz/F");
  evtree->Branch("b2_energy",&b2_energy, "b2_energy/F");
  evtree->Branch("htobb_px",&htobb_px, "htobb_px/F");
  evtree->Branch("htobb_py",&htobb_py, "htobb_py/F");
  evtree->Branch("htobb_pz",&htobb_pz, "htobb_pz/F");
  evtree->Branch("htobb_energy",&htobb_energy, "htobb_energy/F");
  evtree->Branch("htobb_mass",&htobb_mass, "htobb_mass/F");
  evtree->Branch("htobb",&htobb, "htobb/B");
  
 
  evtree->Branch("bjet_px",&bjet_px, "bjet_px/F");
  evtree->Branch("bjet_py",&bjet_py, "bjet_py/F");
  evtree->Branch("bjet_pz",&bjet_pz, "bjet_pz/F");
  evtree->Branch("bjet_energy",&bjet_energy, "bjet_energy/F");
  evtree->Branch("bbarjet_px",&bbarjet_px, "bbarjet_px/F");
  evtree->Branch("bbarjet_py",&bbarjet_py, "bbarjet_py/F");
  evtree->Branch("bbarjet_pz",&bbarjet_pz, "bbarjet_pz/F");
  evtree->Branch("bbarjet_energy",&bbarjet_energy, "bbarjet_energy/F");
  evtree->Branch("dR_bjet", &dR_bjet,"dR_bjet/F");  
  evtree->Branch("dR_bbarjet", &dR_bbarjet,"dR_bbarjet/F");  

  evtree->Branch("hasbjet",&hasbjet, "hasbjet/B");
  evtree->Branch("hasbbarjet",&hasbbarjet, "hasbbarjet/B");
  evtree->Branch("hasMET",&hasMET, "hasMET/B");
  evtree->Branch("hasMuon1",&hasMuon1, "hasMuon1/B");
  evtree->Branch("hasMuon2",&hasMuon2, "hasMuon2/B");
  evtree->Branch("hastwomuons",&hastwomuons, "hastwomuons/B");
  
  evtree->Branch("mu1_px",&mu1_px, "mu1_px/F");
  evtree->Branch("mu1_py",&mu1_py, "mu1_py/F");
  evtree->Branch("mu1_pz",&mu1_pz, "mu1_pz/F");
  evtree->Branch("mu1_energy",&mu1_energy, "mu1_energy/F");
  evtree->Branch("nu1_px",&nu1_px, "nu1_px/F");
  evtree->Branch("nu1_py",&nu1_py, "nu1_py/F");
  evtree->Branch("nu1_pz",&nu1_pz, "nu1_pz/F");
  evtree->Branch("nu1_energy",&nu1_energy, "nu1_energy/F");
  evtree->Branch("mu2_px",&mu2_px, "mu2_px/F");
  evtree->Branch("mu2_py",&mu2_py, "mu2_py/F");
  evtree->Branch("mu2_pz",&mu2_pz, "mu2_pz/F");
  evtree->Branch("mu2_energy",&mu2_energy, "mu2_energy/F");
  evtree->Branch("nu2_px",&nu2_px, "nu2_px/F");
  evtree->Branch("nu2_py",&nu2_py, "nu2_py/F");
  evtree->Branch("nu2_pz",&nu2_pz, "nu2_pz/F");
  evtree->Branch("nu2_energy",&nu2_energy, "nu2_energy/F");
  evtree->Branch("htoWW_energy",&htoWW_energy);
  evtree->Branch("htoWW_px",&htoWW_px,"htoWW_px/F");
  evtree->Branch("htoWW_py",&htoWW_py,"htoWW_px/F");
  evtree->Branch("htoWW_pz",&htoWW_pz,"htoWW_pz/F");
  evtree->Branch("htoWW_mass",&htoWW_mass,"htoWW_mass/F");

  evtree->Branch("Muon1_px",&Muon1_px, "Muon1_px/F");
  evtree->Branch("Muon1_py",&Muon1_py, "Muon1_py/F");
  evtree->Branch("Muon1_pz",&Muon1_pz, "Muon1_pz/F");
  evtree->Branch("Muon1_energy",&Muon1_energy, "Muon1_energy/F");
  evtree->Branch("Muon2_px",&Muon2_px, "Muon2_px/F");
  evtree->Branch("Muon2_py",&Muon2_py, "Muon2_py/F");
  evtree->Branch("Muon2_pz",&Muon2_pz, "Muon2_pz/F");
  evtree->Branch("Muon2_energy",&Muon2_energy, "Muon2_energy/F");

  evtree->Branch("met",&met,"met/F");
  evtree->Branch("met_phi",&met_phi,"met_phi/F");
  evtree->Branch("met_px",&met_px,"met_px/F");
  evtree->Branch("met_py",&met_py,"met_py/F");
	 
  evtree->Branch("h2tohh_mass",&h2tohh_mass,"h2tohh_mass/F");
	//loop over events 1. find h->b bbar; 2 find two bjets after cuts; 3 fill Tree
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");

  TClonesArray *branchTrack = treeReader->UseBranch("Track");
  TClonesArray *branchTower = treeReader->UseBranch("Tower");

  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
  TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");

  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;

  Long64_t entry;
  Jet *jet, *b1jet, *b2jet;
  MissingET *Met;
  GenParticle *particle, *genP, *genhiggs1, *genhiggs2, *genhtobb, *genb1, *genb2;
  GenParticle *genhtoWW, *genW1, *genW2, *genmu1, *genmu2, *gennu1, *gennu2;
  Electron *electron;
  Photon *photon;
  Muon *muon, *muon1, *muon2;

  Track *track;
  Tower *tower;

  TLorentzVector momentum;
  

  Int_t i;
  Int_t j;
  // Loop over all events

  for(entry = 0; entry < allEntries; ++entry)
  //for(entry = 0; entry < 1000; ++entry)
  {
    htobb = false;
    hasbjet = false;
    hasbbarjet = false;
    Wtomu1nu1 = false;
    Wtomu2nu2 = false;
    htoWW = false;
    hastwomuons = false;
    hasMET = false; 
    hasMuon1 = false; 
    hasMuon2 = false; 
     
    dR_bjet = 2.0;
    dR_bbarjet = 2.0;
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    cout <<" event id " << entry << endl;
    /*for (i =0; i < branchParticle->GetEntries(); ++i){
	genP = (GenParticle*) branchParticle->At(i);
        //defualt M1 M2 D1 D2 is -1;
        cout << " genP Id " << genP->PID <<" Pt " << genP->PT << " M1 "<< genP->M1<<" M2 "<< genP->M2 << " D1 "<< genP->D1 <<" D2 "<<genP->D2<< endl;
        if ( genP->M1 >= 0 && genP->M1 <branchParticle->GetEntries()){
        	GenParticle *M1P= (GenParticle*) branchParticle->At(genP->M1);
        	cout <<" M1 Id " << M1P->PID <<" Pt " << M1P->PT << std::endl;
	}
        if ( genP->D1 >= 0 && genP->D1 <branchParticle->GetEntries()){
        	GenParticle *D1P= (GenParticle*) branchParticle->At(genP->D1);
        	cout <<" D1 Id " << D1P->PID <<" Pt " << D1P->PT << std::endl;
         }
     }*/
     genP =  (GenParticle*) branchParticle->At(0); printGenParticle(genP);
     h2tohh_mass = genP->Mass;
     genhiggs1 =  (GenParticle*) branchParticle->At(genP->D1); 
     //printGenParticle(genhiggs1); 
     genhiggs2 =  (GenParticle*) branchParticle->At(genP->D2); 
     //printGenParticle(genhiggs2);
     while ((genhiggs1->D1>0 && ((GenParticle*)branchParticle->At(genhiggs1->D1))->PID == genhiggs1->PID) 
	    || (genhiggs1->D2>0 && ((GenParticle*)branchParticle->At(genhiggs1->D2))->PID == genhiggs1->PID)){
	if (genhiggs1->D1>0 && ((GenParticle*)branchParticle->At(genhiggs1->D1))->PID == genhiggs1->PID) 
		genhiggs1 = (GenParticle*)branchParticle->At(genhiggs1->D1);
        else   genhiggs1 = (GenParticle*)branchParticle->At(genhiggs1->D2);
     }
     while ((genhiggs2->D1>0 && ((GenParticle*)branchParticle->At(genhiggs2->D1))->PID == genhiggs2->PID) 
	    || (genhiggs2->D2>0 && ((GenParticle*)branchParticle->At(genhiggs2->D2))->PID == genhiggs2->PID)){
	if (genhiggs2->D1>0 && ((GenParticle*)branchParticle->At(genhiggs2->D1))->PID == genhiggs2->PID) 
		genhiggs2 = (GenParticle*)branchParticle->At(genhiggs2->D1);
        else   genhiggs2 = (GenParticle*)branchParticle->At(genhiggs2->D2);
     }
    if (abs(((GenParticle*)branchParticle->At(genhiggs1->D1))->PID) ==5 ) {
		//printGenParticle(genhiggs1);
		//printGenParticle((GenParticle*)branchParticle->At(genhiggs1->D1));
		//printGenParticle((GenParticle*)branchParticle->At(genhiggs1->D2));
		htobb = true;
		genhtobb = genhiggs1;
    }
    else if (abs(((GenParticle*)branchParticle->At(genhiggs2->D1))->PID) ==5 ) {
		//printGenParticle(genhiggs2);
		//printGenParticle((GenParticle*)branchParticle->At(genhiggs2->D1));
		//printGenParticle((GenParticle*)branchParticle->At(genhiggs2->D2));
		htobb = true;
		genhtobb = genhiggs2;
     }
  
    if (abs(((GenParticle*)branchParticle->At(genhiggs2->D1))->PID) == 24){
		htoWW = true;
		genhtoWW = genhiggs2;
    }else if (abs(((GenParticle*)branchParticle->At(genhiggs1->D1))->PID) == 24){
		htoWW = true;
		genhtoWW = genhiggs1;
     }
    
    if (htobb){
		if (((GenParticle*)branchParticle->At(genhtobb->D1))->PID == 5){
 			genb1 = (GenParticle*)branchParticle->At(genhtobb->D1);	
 			genb2 = (GenParticle*)branchParticle->At(genhtobb->D2);	
		} 
  		else {
 			genb1 = (GenParticle*)branchParticle->At(genhtobb->D2);	
 			genb2 = (GenParticle*)branchParticle->At(genhtobb->D1);	
		}
    //move on to "final state" bquarks
     while ((genb1->D1>0 && ((GenParticle*)branchParticle->At(genb1->D1))->PID == genb1->PID 
	    )|| (genb1->D2 >0 && ((GenParticle*)branchParticle->At(genb1->D2))->PID == genb1->PID)){
	if ( genb1->D1>0 && ((GenParticle*)branchParticle->At(genb1->D1))->PID == genb1->PID) 
		genb1 = (GenParticle*)branchParticle->At(genb1->D1);
        else   genb1 = (GenParticle*)branchParticle->At(genb1->D2);
     }
     while ((genb2->D1>0 && ((GenParticle*)branchParticle->At(genb2->D1))->PID == genb2->PID 
	    )|| (genb2->D2 >0 && ((GenParticle*)branchParticle->At(genb2->D2))->PID == genb2->PID)){
	if ( genb2->D1>0 && ((GenParticle*)branchParticle->At(genb2->D1))->PID == genb2->PID) 
		genb2 = (GenParticle*)branchParticle->At(genb2->D1);
        else   genb2 = (GenParticle*)branchParticle->At(genb2->D2);
     }
           b1_px = genb1->Px;
	   b1_py = genb1->Py;
	   b1_pz = genb1->Pz;
	   b1_energy = genb1->E;
	   b2_px = genb2->Px;
	   b2_py = genb2->Py;
	   b2_pz = genb2->Pz;
	   b2_energy = genb2->E;
           htobb_px = genhtobb->Px; htobb_py= genhtobb->Py; htobb_pz = genhtobb->Pz; htobb_energy = genhtobb->E;
	   htobb_mass = genhtobb->Mass;
    }

    if (htoWW){
		if (((GenParticle*)branchParticle->At(genhtoWW->D1))->PID == -24){
 			genW1 = (GenParticle*)branchParticle->At(genhtoWW->D1);	//to muon(-13)
 			genW2 = (GenParticle*)branchParticle->At(genhtoWW->D2);	
		} 
  		else {
 			genW1 = (GenParticle*)branchParticle->At(genhtoWW->D2);	
 			genW2 = (GenParticle*)branchParticle->At(genhtoWW->D1);	
		}
       while ((genW1->D1>0 && ((GenParticle*)branchParticle->At(genW1->D1))->PID == genW1->PID)
	    || (genW1->D2 >0 && ((GenParticle*)branchParticle->At(genW1->D2))->PID == genW1->PID)){
 	if ( genW1->D1>0 && ((GenParticle*)branchParticle->At(genW1->D1))->PID == genW1->PID) 
 		genW1 = (GenParticle*)branchParticle->At(genW1->D1);
        else   genW1 = (GenParticle*)branchParticle->At(genW1->D2);
       }
       while ((genW2->D1>0 && ((GenParticle*)branchParticle->At(genW2->D1))->PID == genW2->PID)
	    || (genW2->D2 >0 && ((GenParticle*)branchParticle->At(genW2->D2))->PID == genW2->PID)){
	if ( genW2->D1>0 && ((GenParticle*)branchParticle->At(genW2->D1))->PID == genW2->PID) 
		genW2 = (GenParticle*)branchParticle->At(genW2->D1);
        else   genW2 = (GenParticle*)branchParticle->At(genW2->D2);
       }

       if (genW1->D1>0 && ((GenParticle*)branchParticle->At(genW1->D1))->PID == -13){
	 	genmu1 = (GenParticle*)branchParticle->At(genW1->D1);
		gennu1 = (GenParticle*)branchParticle->At(genW1->D2);
		Wtomu1nu1 = true;
       }else if (genW1->D2 >0 && ((GenParticle*)branchParticle->At(genW1->D2))->PID == -13){
	 	genmu1 = (GenParticle*)branchParticle->At(genW1->D2);
		gennu1 = (GenParticle*)branchParticle->At(genW1->D1);
		Wtomu1nu1 = true;
	}
       if (genW2->D1>0 && ((GenParticle*)branchParticle->At(genW2->D1))->PID == 13){
	 	genmu2 = (GenParticle*)branchParticle->At(genW2->D1);
		gennu2 = (GenParticle*)branchParticle->At(genW2->D2);
		Wtomu2nu2 = true;
	} else if (genW2->D1>0 && ((GenParticle*)branchParticle->At(genW2->D2))->PID == 13){
	 	genmu2 = (GenParticle*)branchParticle->At(genW2->D2);
		gennu2 = (GenParticle*)branchParticle->At(genW2->D1);
		Wtomu2nu2 = true;
	}
    }
   	
    if (Wtomu1nu1){
    	getFinalState(genmu1, branchParticle);	
    	getFinalState(gennu1, branchParticle);	
	mu1_px = genmu1->Px; mu1_py = genmu1->Py; mu1_pz = genmu1->Pz; mu1_energy = genmu1->E;
	nu1_px = gennu1->Px; nu1_py = gennu1->Py; nu1_pz = gennu1->Pz; nu1_energy = gennu1->E;
   	cout << "mu1 from W "; printGenParticle(genmu1);
    }
    if (Wtomu2nu2){
    	getFinalState(genmu2, branchParticle);	
    	getFinalState(gennu2, branchParticle);	
	mu2_px = genmu2->Px; mu2_py = genmu2->Py; mu2_pz = genmu2->Pz; mu2_energy = genmu2->E;
	nu2_px = gennu2->Px; nu2_py = gennu2->Py; nu2_pz = gennu2->Pz; nu2_energy = gennu2->E;
        cout << "mu2 from W "; printGenParticle(genmu2);
    }
 
   //loop all jets 
    for (i =0;  i<  branchJet->GetEntries(); i++)
    {
      jet = (Jet*) branchJet->At(i);
      if (jet->PT < 30 || abs(jet->Eta)>5.0) continue;
      TLorentzVector jet_p4 = jet->P4();
      if (htobb && jet_p4.DeltaR(genb1->P4()) < dR_bjet) {
      	b1jet = jet;
	dR_bjet = jet_p4.DeltaR(genb1->P4());
        bjet_px = jet_p4.Px(); bjet_py = jet_p4.Py(); bjet_pz=jet_p4.Pz(); bjet_energy = jet_p4.Energy();
	hasbjet = true;
        
       }
      if (htobb && jet_p4.DeltaR(genb2->P4()) < dR_bbarjet){
	b2jet = jet;
        dR_bbarjet = jet_p4.DeltaR(genb2->P4());
        bbarjet_px = jet_p4.Px(); bbarjet_py = jet_p4.Py(); bbarjet_pz=jet_p4.Pz(); bbarjet_energy = jet_p4.Energy();
	hasbbarjet = true;
       }
      //cout <<" jet Pt " << jet->PT <<" Refparticles size " << particlearray.GetEntries() << " Constituents size " << constituentarray.GetEntries() << endl;
      //printJet(jet); 
    }
      // bjet should be different from bbarjet
    if (hasbjet && hasbbarjet && b1jet == b2jet){
	hasbjet = false;
	hasbbarjet = false;
     }
    else if (hasbjet && hasbbarjet){
     //cout <<" bjet Pt " << b1jet->PT  <<"  bbarjet Pt " << b2jet->PT << endl; 
    }
    
    // Analyse missing ET
    if(branchMissingET->GetEntriesFast() > 0)
    {
      Met = (MissingET*) branchMissingET->At(0);
      met = Met->MET;
      met_phi = Met->Phi;
      met_px = Met->P4().Px();
      met_py = Met->P4().Py();
      if (met > 20) hasMET = true;
    }

    // Loop over all Muon in event
    if (Wtomu1nu1 && Wtomu2nu2){
    if (genmu1->PT > 10 && genmu2->PT > 10 && fabs(genmu1->Eta)<2.4 && fabs(genmu2->Eta)<2.4) hastwomuons =true;
    }
    for(i = 0; i < branchMuon->GetEntriesFast(); ++i)
    {
      muon = (Muon*) branchMuon->At(i);
      particle = (GenParticle*) muon->Particle.GetObject();
//      printGenParticle(particle);
      if (Wtomu1nu1 && particle == genmu1) {
	muon1 = muon;
	hasMuon1 = true;
      }
      else if (Wtomu2nu2 && particle == genmu2) {
	muon2 = muon;
        hasMuon2 = true;
      }
      cout <<" muon eta " << muon->Eta << " phi " << muon->Phi << " Pt "<< muon->PT << endl; 
    }
    if (hasMuon1){
	Muon1_px = muon1->P4().Px(); Muon1_py = muon1->P4().Py(); Muon1_pz = muon1->P4().Pz(); Muon1_energy = muon1->P4().E();
	}
    if (hasMuon2){
	Muon2_px = muon2->P4().Px(); Muon2_py = muon2->P4().Py(); Muon2_pz = muon2->P4().Pz(); Muon2_energy = muon2->P4().E();
	}
   
//fill branches
    if (htobb || (Wtomu1nu1 && Wtomu2nu2)) evtree->Fill();
    }
    //evtree->Fill();
}



//------------------------------------------------------------------------------

//void DiHiggs_htobb(const char *inputFile)
void DiHiggs_htobb()
{
  gSystem->Load("libDelphes");

  TChain *chain = new TChain("Delphes");
  //const char *inputFile("/home/taohuang/Herwig++/Delphes-3.2.0/delphes_output.root");
  const char *inputFile("/fdata/hepx/store/user/taohuang/Hhh/delphes320_B3_100k.root");
  chain->Add(inputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);


  TTree *evtree = new TTree("evtree","event tree");
  
  
  AnalyseEvents(treeReader, evtree);
  //evtree->Print();
  TFile *file = new TFile("/fdata/hepx/store/user/taohuang/Hhh/DiHiggs_0702_100k.root","recreate");
  evtree->Write();
  file->Close();

  cout << "** Exiting..." << endl;

  delete treeReader;
  delete chain;
}

