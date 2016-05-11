/*
 * \file L1TRatesAndEffJets.cc
 *
 * \author I. Ojalvo
 * Written for miniAOD
 */

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "L1Trigger/Stage3Ntuplizer/interface/L1TRatesAndEffJets.h"
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

bool compareByPt (l1extra::L1JetParticle i,l1extra::L1JetParticle j) { return(i.pt()>j.pt()); };

L1TRatesAndEffJets::L1TRatesAndEffJets( const ParameterSet & cfg ) :
  rctSource_(cfg.getParameter<edm::InputTag>("rctSource")),
  pfCandsToken_(consumes<vector<pat::PackedCandidate> >(cfg.getParameter<edm::InputTag>("pfCands"))),
  ecalSrc_(consumes<EcalTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("ecalDigis"))),
  hcalSrc_(consumes<HcalTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("hcalDigis"))),
  vtxLabel_(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("vertices"))),
  jetSrc_(consumes<vector<pat::Jet> >(cfg.getParameter<edm::InputTag>("recoJets"))),
  l1ExtraJetSource_(consumes<vector <l1extra::L1JetParticle> >(cfg.getParameter<edm::InputTag>("l1ExtraJetSource")))
  {

    folderName_          = cfg.getUntrackedParameter<std::string>("folderName");
    recoPt_              = cfg.getParameter<double>("recoPtCut");
    folder               = tfs_->mkdir(folderName_);
    efficiencyTree = folder.make<TTree>("EfficiencyTree", "Efficiency Tree");
    efficiencyTree->Branch("run",    &run,     "run/I");
    efficiencyTree->Branch("lumi",   &lumi,    "lumi/I");
    efficiencyTree->Branch("event",  &event,   "event/I");
    efficiencyTree->Branch("nvtx",          &nvtx,         "nvtx/D");
    
    efficiencyTree->Branch("recoPt",    &recoPt,   "recoPt/D");
    
    efficiencyTree->Branch("jetPt",      &jetPt, "jetPt/D");

    efficiencyTree->Branch("recoEta",       &recoEta,   "recoEta/D");
    efficiencyTree->Branch("jetEta",     &jetEta, "jetEta/D");
    
    efficiencyTree->Branch("recoPhi",       &recoPhi,   "recoPhi/D");
    efficiencyTree->Branch("jetPhi",     &jetPhi, "jetPhi/D");

    efficiencyTree->Branch("l1Matched",  &l1Matched, "l1Matched/I");

    nEvents     = folder.make<TH1F>( "nEvents"  , "nEvents", 2,  0., 1. );
    
    jet_pt      = folder.make<TH1F>( "jet_pt"  , "p_{t}", 100,  0., 100. );
    jet_eta      = folder.make<TH1F>( "jet_eta"  , "#eta", 100,  -3, 3 );
    jet_phi   = folder.make<TH1F>( "jet_phi"  , "#phi", 100,  -4, 4 );

    jet_pt_jet1         = folder.make<TH1F>( "jet_pt_jet1"         , "p_{t}", 150,  0., 150. ); 
    jet_pt_jet2         = folder.make<TH1F>( "jet_pt_jet2"         , "p_{t}", 150,  0., 150. ); 
    jet_pt_jet3         = folder.make<TH1F>( "jet_pt_jet3"         , "p_{t}", 150,  0., 150. ); 
    jet_pt_jet4         = folder.make<TH1F>( "jet_pt_jet4"         , "p_{t}", 150,  0., 150. ); 

    jet_pt_jet1_eta2p4  = folder.make<TH1F>( "jet_pt_jet1_eta2p4"  , "p_{t}", 150,  0., 150. ); 
    jet_pt_jet2_eta2p4  = folder.make<TH1F>( "jet_pt_jet2_eta2p4"  , "p_{t}", 150,  0., 150. ); 
    jet_pt_jet3_eta2p4  = folder.make<TH1F>( "jet_pt_jet3_eta2p4"  , "p_{t}", 150,  0., 150. ); 
    jet_pt_jet4_eta2p4  = folder.make<TH1F>( "jet_pt_jet4_eta2p4"  , "p_{t}", 150,  0., 150. ); 

    recoJet_pt   = folder.make<TH1F>( "recoJet_pt"  , "p_{t}", 100,  0., 100. );
    recoJet_eta  = folder.make<TH1F>( "recoJet_eta"  , "eta", 100,  -3, 3. );
    recoJet_phi  = folder.make<TH1F>( "recoJet_phi"  , "phi", 100,  -4, 4. );

  }

void L1TRatesAndEffJets::beginJob( const EventSetup & es) {
   std::cout<<"begin job..."<<std::endl;
}

void L1TRatesAndEffJets::analyze( const Event& evt, const EventSetup& es )
 {
   std::cout<<"Analyzing..."<<std::endl;
   nEvents->Fill(1);
   
   run = evt.id().run();
   lumi = evt.id().luminosityBlock();
   event = evt.id().event();
   edm::Handle<reco::VertexCollection> vertices;   
   //edm::Handle < BXVector<l1t::Tau> > stage1IsoTaus;

  edm::Handle < vector<l1extra::L1JetParticle> > l1ExtraJets;
  
  std::vector<pat::Tau> goodJets;
  
  edm::Handle<vector<pat::PackedCandidate> >pfCands;
  edm::Handle<EcalTrigPrimDigiCollection> ecalTPGs;
  edm::Handle<HcalTrigPrimDigiCollection> hcalTPGs;
  

  ESHandle<L1CaloHcalScale> hcalScale;
  es.get<L1CaloHcalScaleRcd>().get(hcalScale);

  //get pfcandidate particles
  if(!evt.getByToken(pfCandsToken_, pfCands)){
    std::cout<<"Error Getting PFCandidates"<<std::endl;
  }

  //Make Rates
  // loop over taus
  Handle<vector<pat::Jet> > jets;
  if(evt.getByToken(jetSrc_, jets)){//Begin Getting Reco Taus
    for (const pat::Jet &jet : *jets) {
      recoJet_pt->Fill( jet.pt() );
      recoJet_eta->Fill( jet.eta() );
      recoJet_phi->Fill( jet.phi() );
      //get rid of the cruft for analysis to save disk space
      if(jet.pt() > recoPt_ ) {
	goodJets.push_back(jet);

      }
    }
  }
  else
    std::cout<<"Error getting reco taus"<<std::endl;

  ////
  vector<l1extra::L1JetParticle> l1JetSorted;
  vector<l1extra::L1JetParticle> l1JetSortedEtaRestricted2p1;
  vector<l1extra::L1JetParticle> l1JetSortedEtaRestricted2p4;
  if(evt.getByToken(l1ExtraJetSource_, l1ExtraJets)){
    std::cout<<"found rlx stage 3 jets"<<std::endl;
    for( vector<l1extra::L1JetParticle>::const_iterator l1Jet = l1ExtraJets->begin(); l1Jet != l1ExtraJets->end(); l1Jet++ ){
      l1JetSorted.push_back(*l1Jet);
      if(abs(l1Jet->eta()) < 2.4) l1JetSortedEtaRestricted2p4.push_back(*l1Jet);
    }
  }
  
  std::sort(l1JetSorted.begin(),l1JetSorted.end(),compareByPt);
  std::sort(l1JetSortedEtaRestricted2p4.begin(),l1JetSortedEtaRestricted2p4.end(),compareByPt);
  std::cout<<"making jets"<<std::endl;

  //Begin Making Rate Plots
  for( auto l1Jet : l1JetSorted ) {
    if( l1Jet.pt() > 150 )
      jet_pt->Fill( 150 );
    else jet_pt->Fill( l1Jet.pt() );
    jet_eta->Fill( l1Jet.eta() );
    jet_phi->Fill( l1Jet.phi() );
    //std::cout<<"rlx jet pt "<< l1Jet->pt() << " rlx jet eta "<< l1Jet->eta()<<" rlx jet phi "<< l1Jet->phi()<<std::endl;
  }

  if(l1JetSorted.size()>0)
    jet_pt_jet1->Fill(l1JetSorted.at(0).pt());
  if(l1JetSorted.size()>1)
    jet_pt_jet2->Fill(l1JetSorted.at(1).pt());
  if(l1JetSorted.size()>2)
    jet_pt_jet3->Fill(l1JetSorted.at(2).pt());
  if(l1JetSorted.size()>3)
    jet_pt_jet4->Fill(l1JetSorted.at(3).pt());

  if(l1JetSortedEtaRestricted2p4.size()>0)
    jet_pt_jet1_eta2p4->Fill(l1JetSortedEtaRestricted2p4.at(0).pt());
  if(l1JetSortedEtaRestricted2p4.size()>1)
    jet_pt_jet2_eta2p4->Fill(l1JetSortedEtaRestricted2p4.at(1).pt());
  if(l1JetSortedEtaRestricted2p4.size()>2)
    jet_pt_jet3_eta2p4->Fill(l1JetSortedEtaRestricted2p4.at(2).pt());
  if(l1JetSortedEtaRestricted2p4.size()>3)
    jet_pt_jet4_eta2p4->Fill(l1JetSortedEtaRestricted2p4.at(3).pt());


  //End Making Rate Plots
  std::cout<<"looking through the jets"<<std::endl;

  //If there isn't at least 1 good reco jet don't bother doing all the work
  if(goodJets.size()>0){
    double deltaR_ = 0.5;
    
    ////Make efficiencies
    for(unsigned int i = 0; i < goodJets.size(); i++){
      nvtx = 0;
      if(evt.getByToken(vtxLabel_, vertices)){
	nvtx = (double)vertices->size();
	//std::cout<<"nvtx "<<nvtx<<std::endl;
      }

      pat::Jet recoJet = goodJets.at(i);
      
      ////Fill Reco Objects
      recoPt  = recoJet.pt();
      recoEta = recoJet.eta();
      recoPhi = recoJet.phi();
      
      //Fill L1 Objects
      l1Matched = -1; 
      jetPt = 0; jetEta = -99; jetPhi = -99; 
      
      for(uint32_t k = 0; k<jetSorted.size(); k++){
	double dR = deltaR( recoJet.p4(), isoJetSorted.at(k).p4());
	if( dR < deltaR_){
	  jetPt  = l1JetSorted.at(k).pt();
	  jetEta = l1JetSorted.at(k).eta();
	  jetPhi = l1JetSorted.at(k).phi();
	  l1Matched = 1;
	  break;
	}
      }
      if(l1Matched==0)
	std::cout<<"Not Matched!"<<std::endl;
      std::cout<<"nvtx "<<nvtx<<std::endl;
      efficiencyTree->Fill();
    }
  }

}


void L1TRatesAndEffJets::endJob() {
}

L1TRatesAndEffJets::~L1TRatesAndEffJets(){

}

DEFINE_FWK_MODULE(L1TRatesAndEffJets);
