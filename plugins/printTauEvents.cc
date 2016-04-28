/*
 * \file printTauEvents.cc
 *
 * \author I. Ojalvo
 * Written for miniAOD
 */

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "L1Trigger/Stage3Ntuplizer/interface/printTauEvents.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/Math/interface/deltaR.h"
#include <iostream>
#include <fstream>

using namespace edm;
using std::cout;
using std::endl;
using std::vector;

printTauEvents::printTauEvents( const ParameterSet & cfg ) :
  vtxLabel_(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("vertices"))),
  discriminatorMu_(consumes<reco::PFTauDiscriminator>(cfg.getParameter<edm::InputTag>("recoTauDiscriminatorMu"))),
  discriminatorIso_(consumes<reco::PFTauDiscriminator>(cfg.getParameter<edm::InputTag>("recoTauDiscriminatorIso"))),
  tauSrc_(consumes<vector<pat::Tau> >(cfg.getParameter<edm::InputTag>("recoTau")))
  {

  file0.open ("events-DecayMode0.txt");
  file1.open ("events-DecayMode1.txt");
  file10.open ("events-DecayMode10.txt");

  recoPt_              = cfg.getParameter<double>("recoPtCut");
  maxRecoPt_           = cfg.getParameter<double>("maxRecoPt");

  efficiencyTree = tfs_->make<TTree>("EfficiencyTree", "Efficiency Tree");
  efficiencyTree->Branch("run",    &run,     "run/I");
  efficiencyTree->Branch("lumi",   &lumi,    "lumi/I");
  efficiencyTree->Branch("event",  &event,   "event/I");
  efficiencyTree->Branch("nvtx",          &nvtx,         "nvtx/I");
  
  efficiencyTree->Branch("recoPt",    &recoPt,   "recoPt/D");
  efficiencyTree->Branch("decayMode", &decayMode,   "decayMode/I");
  
  efficiencyTree->Branch("recoEta",       &recoEta,   "recoEta/D");
  efficiencyTree->Branch("recoPhi",       &recoPhi,   "recoPhi/D");
  
  }

void printTauEvents::beginJob( const EventSetup & es) {
}


void printTauEvents::analyze( const Event& evt, const EventSetup& es )
 {

  run = evt.id().run();
  lumi = evt.id().luminosityBlock();
  event = evt.id().event();
  std::cout<<run<<":"<<lumi<<":"<<evt.id().event()<<std::endl;
  edm::Handle<reco::VertexCollection> vertices;
  if(evt.getByToken(vtxLabel_, vertices)){
    nvtx = vertices->size();
  }
  
  Handle<L1CaloRegionCollection> regions;

  edm::Handle<reco::PFTauDiscriminator> discriminatorIso;
  edm::Handle<reco::PFTauDiscriminator> discriminatorMu;

  std::vector<reco::PFTauRef> goodTausRef;
  std::vector<pat::Tau> goodTaus;

  Handle<vector<pat::Tau> > taus;
  if(evt.getByToken(tauSrc_, taus)){//Begin Getting Reco Taus
    for ( unsigned iTau = 0; iTau < taus->size(); ++iTau ) {
      //reco::PFTauRef tauCandidate(taus, iTau);
      pat::Tau tau = taus->at(iTau);
      // std::cout<<"reco tau pt "<< tau.pt()<< " eta "<<tau.eta() << " phi " << tau.phi() << " ... three hits 10: " << tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") << " ... mu loose gt0: " << tau.tauID("againstMuonLoose3")<< " ... decaymodefinding: " <<tau.tauID("decayModeFinding") <<std::endl;
      // discriminate by decaymodefiniding, antiMu, Iso, 
      if(tau.tauID("decayModeFinding")>0.5&&tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits")<20&&tau.tauID("againstMuonLoose3")>0){
	if( tau.decayMode()>-1){
	  
	  std::cout<<"tauPT = "<< tau.pt() <<" eta= " <<tau.eta() <<" phi= " << tau.phi()<<" tau.decayMode() "<<tau.decayMode()<<std::endl;
	  std::cout<<run<<":"<<lumi<<":"<<evt.id().event()<<std::endl;
	  //get rid of the cruft for analysis to save disk space
	  if(tau.pt() > recoPt_ && tau.pt() < maxRecoPt_) {
	    goodTaus.push_back(tau);

	    if(tau.decayMode()==0){
	      file0 <<"tauPT = "<< tau.pt() <<" eta= " <<tau.eta() <<" phi= " << tau.phi()<<std::endl;
	      file0 <<run<<":"<<lumi<<":"<<evt.id().event()<<std::endl;
	    }
	    if(tau.decayMode()==1){
	      file1 <<"tauPT = "<< tau.pt() <<" eta= " <<tau.eta() <<" phi= " << tau.phi()<<std::endl;
	      file1 <<run<<":"<<lumi<<":"<<evt.id().event()<<std::endl;
	    }
	    if(tau.decayMode()==10){
	      file10 <<"tauPT = "<< tau.pt() <<" eta= " <<tau.eta() <<" phi= " << tau.phi()<<std::endl;
	      file10 <<run<<":"<<lumi<<":"<<evt.id().event()<<std::endl;
	    }
	  }
	}
      }
    }
  }//End Getting Reco Taus
  else
    std::cout<<"Error getting reco taus"<<std::endl;

  //If there isn't at least 1 good reco tau don't bother filling
  if(goodTaus.size()>0){
    
    ////Make efficiencies
    for(unsigned int i = 0; i < goodTaus.size(); i++){
      
      pat::Tau recoTau = goodTaus.at(i);
      decayMode = -999;
      
      ////Fill Reco Objects
      recoPt  = recoTau.pt();
      recoEta = recoTau.eta();
      recoPhi = recoTau.phi();
      decayMode = recoTau.decayMode();

      efficiencyTree->Fill();
    }

  }

 }


void printTauEvents::endJob() {
}

printTauEvents::~printTauEvents(){
  file0.close();  file1.close();  file10.close();
}

DEFINE_FWK_MODULE(printTauEvents);
