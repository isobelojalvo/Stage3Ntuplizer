/*
 * \file L1TRatesAndEff.cc
 *
 * \author I. Ojalvo
 * Written for miniAOD
 */

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "L1Trigger/Stage3Ntuplizer/interface/L1TRatesAndEff.h"
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

L1TRatesAndEff::L1TRatesAndEff( const ParameterSet & cfg ) :
  rctSource_(cfg.getParameter<edm::InputTag>("rctSource")),
  pfCandsToken_(consumes<vector<pat::PackedCandidate> >(cfg.getParameter<edm::InputTag>("pfCands"))),
  ecalSrc_(consumes<EcalTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("ecalDigis"))),
  hcalSrc_(consumes<HcalTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("hcalDigis"))),
  vtxLabel_(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("vertices"))),
  discriminatorMu_(consumes<reco::PFTauDiscriminator>(cfg.getParameter<edm::InputTag>("recoTauDiscriminatorMu"))),
  discriminatorIso_(consumes<reco::PFTauDiscriminator>(cfg.getParameter<edm::InputTag>("recoTauDiscriminatorIso"))),
  tauSrc_(consumes<vector<pat::Tau> >(cfg.getParameter<edm::InputTag>("recoTau"))),
  l1ExtraIsoTauSource_(consumes<vector <l1extra::L1JetParticle> >(cfg.getParameter<edm::InputTag>("l1ExtraIsoTauSource"))),
  l1ExtraTauSource_(consumes<vector <l1extra::L1JetParticle> >(cfg.getParameter<edm::InputTag>("l1ExtraTauSource"))),
  l1Stage2TauSource_(consumes<BXVector<l1t::Tau> > (cfg.getParameter<edm::InputTag>("stage2TauSource"))),
  l1Stage1TauSource_(consumes<BXVector<l1t::Tau> > (cfg.getParameter<edm::InputTag>("stage1TauSource"))),
  l1Stage1IsoTauSource_(consumes<BXVector<l1t::Tau> > (cfg.getParameter<edm::InputTag>("stage1IsoTauSource"))),
  l1GctTauSource_(consumes<vector<L1GctJetCand> > (cfg.getParameter<edm::InputTag>("gctTauJets"))),
  //  vector<L1GctJetCand>                  "gctDigis"                  "tauJets"         "RAW2DIGI"
  regionSource_(consumes<vector <L1CaloRegion> >(cfg.getParameter<edm::InputTag>("UCTRegion")))
  {

    folderName_          = cfg.getUntrackedParameter<std::string>("folderName");
    recoPt_              = cfg.getParameter<double>("recoPtCut");
    folder               = tfs_->mkdir(folderName_);
    //folder->cd();
    //efficiencyTree = tfs_->make<TTree>("EfficiencyTree", "Efficiency Tree");
    efficiencyTree = folder.make<TTree>("EfficiencyTree", "Efficiency Tree");
    efficiencyTree->Branch("run",    &run,     "run/I");
    efficiencyTree->Branch("lumi",   &lumi,    "lumi/I");
    efficiencyTree->Branch("event",  &event,   "event/I");
    efficiencyTree->Branch("nvtx",          &nvtx,         "nvtx/D");
    
    efficiencyTree->Branch("recoPt",    &recoPt,   "recoPt/D");
    efficiencyTree->Branch("decayMode", &decayMode,   "decayMode/I");
    
    efficiencyTree->Branch("tauEtaEcalEnt", &tauEtaEcalEnt,"tauEtaEcalEnt/D");
    efficiencyTree->Branch("tauPhiEcalEnt", &tauPhiEcalEnt,"tauPhiEcalEnt/D");
    efficiencyTree->Branch("rawEcal",       &rawEcal,      "rawEcal/D");
    efficiencyTree->Branch("rawHcal",       &rawHcal,      "rawHcal/D");
    efficiencyTree->Branch("ecal",          &ecal,         "ecal/D");
    efficiencyTree->Branch("hcal",          &hcal,         "hcal/D");
    efficiencyTree->Branch("jetEt",         &jetEt,        "jetEt/D");
    efficiencyTree->Branch("jetEta",        &jetEta,       "jetEta/D");
    efficiencyTree->Branch("jetPhi",        &jetPhi,       "jetPhi/D");
    efficiencyTree->Branch("pfCandsEt",     &pfCandsEt,    "pfCandsEt/D");
    efficiencyTree->Branch("signalCandsEt", &signalCandsEt,"signalCandsEt/D");
    efficiencyTree->Branch("isoCandsEt",    &isoCandsEt,   "isoCandsEt/D");
    
    efficiencyTree->Branch("max3ProngDeltaR",     &max3ProngDeltaR,    "max3ProngDeltaR/D");
    efficiencyTree->Branch("minProngPt",          &minProngPt,         "minProngPt/D");
    efficiencyTree->Branch("midProngPt",          &midProngPt,         "midProngPt/D");
    efficiencyTree->Branch("maxProngPt",          &maxProngPt,         "maxProngPt/D");
    efficiencyTree->Branch("n3ProngCands",        &n3ProngCands,       "maxProngPt/I");
    
    efficiencyTree->Branch("hcal",          &hcal,         "hcal/D");
    
    efficiencyTree->Branch("ecalTPG2x2",    &TPGE2x2,      "ecalTPG2x2/D");
    efficiencyTree->Branch("hcalTPG2x2",    &TPGH2x2,      "hcalTPG2x2/D");
    efficiencyTree->Branch("TPG2x2",        &TPG2x2,       "TPG2x2/D");
    
    efficiencyTree->Branch("ecalTPG5x5",    &TPGE5x5,      "ecalTPG5x5/D");
    efficiencyTree->Branch("hcalTPG5x5",    &TPGH5x5,      "hcalTPG5x5/D");
    efficiencyTree->Branch("TPG5x5",        &TPG5x5,       "TPG5x5/D");

    efficiencyTree->Branch("ecalTPG6x6",    &TPGE6x6,      "ecalTPG6x6/D");
    efficiencyTree->Branch("hcalTPG6x6",    &TPGH6x6,      "hcalTPG6x6/D");
    efficiencyTree->Branch("TPG6x6",        &TPG6x6,       "TPG6x6/D");
    
    efficiencyTree->Branch("ecalTPG7x7",    &TPGE7x7,      "ecalTPG7x7/D");
    efficiencyTree->Branch("hcalTPG7x7",    &TPGH7x7,      "hcalTPG7x7/D");
    efficiencyTree->Branch("TPG7x7",        &TPG7x7,       "TPG7x7/D");
    

    efficiencyTree->Branch("isoTauPt",      &isoTauPt, "isoTauPt/D");
    efficiencyTree->Branch("rlxTauPt",      &rlxTauPt, "rlxTauPt/D");
    
    efficiencyTree->Branch("recoEta",       &recoEta,   "recoEta/D");
    efficiencyTree->Branch("isoTauEta",     &isoTauEta, "isoTauEta/D");
    efficiencyTree->Branch("rlxTauEta",     &rlxTauEta, "rlxTauEta/D");
    
    efficiencyTree->Branch("recoPhi",       &recoPhi,   "recoPhi/D");
    efficiencyTree->Branch("isoTauPhi",     &isoTauPhi, "isoTauPhi/D");
    efficiencyTree->Branch("rlxTauPhi",     &rlxTauPhi, "rlxTauPhi/D");

    efficiencyTree->Branch("l1IsoMatched",  &l1IsoMatched, "l1IsoMatched/I");
    efficiencyTree->Branch("l1RlxMatched",  &l1RlxMatched, "l1RlxMatched/I");

    nEvents     = folder.make<TH1F>( "nEvents"  , "nEvents", 2,  0., 1. );
    
    isoTau_pt   = folder.make<TH1F>( "isoTau_pt"  , "p_{T}", 100,  0., 100. );
    tau_pt      = folder.make<TH1F>( "tau_pt"  , "p_{t}", 100,  0., 100. );

    isoTau_eta   = folder.make<TH1F>( "isoTau_eta"  , "#eta", 100,  -3, 3 );
    tau_eta      = folder.make<TH1F>( "tau_eta"  , "#eta", 100,  -3, 3 );

    isoTau_phi   = folder.make<TH1F>( "isoTau_phi"  , "#phi", 100,  -4, 4 );
    tau_phi      = folder.make<TH1F>( "tau_phi"  , "#phi", 100,  -4, 4 );

    tau_pt_diTau         = folder.make<TH1F>( "tau_pt_diTau"         , "p_{t}", 150,  0., 150. );    
    tau_pt_diTau_eta2p4  = folder.make<TH1F>( "tau_pt_diTau_eta2p4"  , "p_{t}", 150,  0., 150. );    
    tau_pt_diTau_eta2p1  = folder.make<TH1F>( "tau_pt_diTau_eta2p1"  , "p_{t}", 150,  0., 150. );    

    isoTau_pt_diTau         = folder.make<TH1F>( "isoTau_pt_diTau"         , "p_{t}", 150,  0., 150. );    
    isoTau_pt_diTau_eta2p4  = folder.make<TH1F>( "isoTau_pt_diTau_eta2p4"  , "p_{t}", 150,  0., 150. );    
    isoTau_pt_diTau_eta2p1  = folder.make<TH1F>( "isoTau_pt_diTau_eta2p1"  , "p_{t}", 150,  0., 150. );    

    recoTau_pt   = folder.make<TH1F>( "recoTau_pt"  , "p_{t}", 100,  0., 100. );
    recoTau_eta  = folder.make<TH1F>( "recoTau_eta"  , "eta", 100,  -3, 3. );
    recoTau_phi  = folder.make<TH1F>( "recoTau_phi"  , "phi", 100,  -4, 4. );

    regionHitEta   = folder.make<TH1F>( "regionHit_eta"  , "eta", 16, 1, 16. );
    regionHitPhi   = folder.make<TH1F>( "regionHit_phi"  , "phi", 16, 1, 16. );
    regionTotal   = folder.make<TH1F>( "regionHit_total"  , "fullmap", 16, 1, 16. );

    regionEta   = folder.make<TH1F>( "region_eta"  , "eta", 22, 1, 22. );
    regionPhi   = folder.make<TH1F>( "region_phi"  , "phi", 72, 1, 72. );
    regionPt   = folder.make<TH1F>( "region_pt"  , "pt", 100, 0, 100. );

    regionEtaFine   = folder.make<TH1F>( "region_eta_Fine"  , "eta", 88, 1, 88. );
    regionPhiFine   = folder.make<TH1F>( "region_phi_Fine"  , "phi", 72, 1, 72. );

  }

void L1TRatesAndEff::beginJob( const EventSetup & es) {
   std::cout<<"begin job..."<<std::endl;
}

void L1TRatesAndEff::analyze( const Event& evt, const EventSetup& es )
 {
   std::cout<<"Analyzing..."<<std::endl;
   nEvents->Fill(1);
   
   run = evt.id().run();
   lumi = evt.id().luminosityBlock();
   event = evt.id().event();
   edm::Handle<reco::VertexCollection> vertices;   
   Handle<L1CaloRegionCollection> regions;
   
   edm::Handle<reco::PFTauDiscriminator> discriminatorIso;
   edm::Handle<reco::PFTauDiscriminator> discriminatorMu;
   
   edm::Handle < L1GctJetCandCollection > l1IsoTauJets;
   edm::Handle < L1GctJetCandCollection > l1TauJets;
   edm::Handle < BXVector<l1t::Tau> > stage2Taus;
   edm::Handle < BXVector<l1t::Tau> > stage1Taus;
   edm::Handle < BXVector<l1t::Tau> > stage1IsoTaus;
   
   edm::Handle < vector<l1extra::L1JetParticle> > l1ExtraTaus;
   edm::Handle < vector<l1extra::L1JetParticle> > l1ExtraIsoTaus;
   
   edm::Handle < vector <L1GctJetCand> > l1GctTaus;
   
   std::vector<reco::PFTauRef> goodTausRef;
   std::vector<pat::Tau> goodTaus;
   
   edm::Handle<vector<pat::PackedCandidate> >pfCands;
   edm::Handle<EcalTrigPrimDigiCollection> ecalTPGs;
   edm::Handle<HcalTrigPrimDigiCollection> hcalTPGs;
   
   
  if(!evt.getByToken(ecalSrc_, ecalTPGs))
    std::cout<<"ERROR GETTING THE ECAL TPGS"<<std::endl;
  if(!evt.getByToken(hcalSrc_, hcalTPGs))
    std::cout<<"ERROR GETTING THE HCAL TPGS"<<std::endl;

  ESHandle<L1CaloHcalScale> hcalScale;
  es.get<L1CaloHcalScaleRcd>().get(hcalScale);

  //get pfcandidate particles
  if(!evt.getByToken(pfCandsToken_, pfCands)){
    std::cout<<"Error Getting PFCandidates"<<std::endl;
  }

  //Make Rates
  // loop over taus
  Handle<vector<pat::Tau> > taus;
  if(evt.getByToken(tauSrc_, taus)){//Begin Getting Reco Taus
    for ( unsigned iTau = 0; iTau < taus->size(); ++iTau ) {
      pat::Tau tau = taus->at(iTau);
      // discriminate by decaymodefiniding, antiMu, Iso, 
      if(tau.tauID("decayModeFinding")>0.5&&tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits")<20&&tau.tauID("againstMuonLoose3")>0){
	if( tau.decayMode()>-1){
	  
	  recoTau_pt->Fill( tau.pt() );
	  recoTau_eta->Fill( tau.eta() );
	  recoTau_phi->Fill( tau.phi() );

	  //get rid of the cruft for analysis to save disk space
	  if(tau.pt() > recoPt_ ) {
	    goodTaus.push_back(tau);

	  }
	}
      }
    }
  }//End Getting Reco Taus
  else
    std::cout<<"Error getting reco taus"<<std::endl;

  ////
  vector<l1extra::L1JetParticle> rlxTauSorted;
  vector<l1extra::L1JetParticle> rlxTauSortedEtaRestricted2p1;
  vector<l1extra::L1JetParticle> rlxTauSortedEtaRestricted2p4;
  if(evt.getByToken(l1ExtraTauSource_, l1ExtraTaus)){
    std::cout<<"found rlx stage 3 taus"<<std::endl;
    for( vector<l1extra::L1JetParticle>::const_iterator rlxTau = l1ExtraTaus->begin(); rlxTau != l1ExtraTaus->end(); rlxTau++ ){
      rlxTauSorted.push_back(*rlxTau);
      if(abs(rlxTau->eta()) < 2.4) rlxTauSortedEtaRestricted2p4.push_back(*rlxTau);
      if(abs(rlxTau->eta()) < 2.1) rlxTauSortedEtaRestricted2p1.push_back(*rlxTau);
    }
  }
  else if(evt.getByToken(l1Stage2TauSource_, stage2Taus)){
    std::cout<<"found rlx stage 2 taus size: "<< stage2Taus->size() <<std::endl;
    for(BXVector<l1t::Tau>::const_iterator rlxTau = stage2Taus->begin(); rlxTau != stage2Taus->end(); rlxTau++ ) {
      //make this int l1extra::L1JetParticle
      if(rlxTau->hwIso() > 0){
	std::cout<<"rlxTau Pt "<<rlxTau->pt()<<" eta: "<<rlxTau->eta()<<" phi: "<<rlxTau->phi()<<std::endl;
	l1extra::L1JetParticle tempJet(rlxTau->p4());
	rlxTauSorted.push_back(tempJet);
	if(abs(rlxTau->eta()) < 2.4) rlxTauSortedEtaRestricted2p4.push_back(tempJet);
	if(abs(rlxTau->eta()) < 2.1) rlxTauSortedEtaRestricted2p1.push_back(tempJet);
      }
    }
  }
  else if(evt.getByToken(l1Stage1TauSource_, stage1Taus)){
    std::cout<<"found rlx stage 1 taus size: "<< stage1Taus->size() <<std::endl;
    //std::cout<<"1. pt: "<< stage1Taus->at(0,0).pt()<<std::endl;
    for(BXVector<l1t::Tau>::const_iterator rlxTau = stage1Taus->begin(); rlxTau != stage1Taus->end(); rlxTau++ ) {
      //make this int l1extra::L1JetParticle
      //std::cout<<"rlxTau Pt "<<rlxTau->pt()<<" eta: "<<rlxTau->eta()<<" phi: "<<rlxTau->phi()<<std::endl;
      l1extra::L1JetParticle tempJet(rlxTau->p4());
      rlxTauSorted.push_back(tempJet);
      if(abs(rlxTau->eta()) < 2.4) rlxTauSortedEtaRestricted2p4.push_back(tempJet);
      if(abs(rlxTau->eta()) < 2.1) rlxTauSortedEtaRestricted2p1.push_back(tempJet);
    }
  }
  else if(evt.getByToken(l1GctTauSource_, l1GctTaus)){
    std::cout<<"found gct taus size: "<< l1GctTaus->size() <<std::endl;
    for(vector<L1GctJetCand>::const_iterator gctTau = l1GctTaus->begin(); gctTau != l1GctTaus->end(); gctTau++) {
      //L1CaloRegion region = L1CaloRegion( gctTau->rank(), false, gctTau->regionId().rctCrate(), gctTau->regionId().rctRegion());
      
      //UCTRegionProcess uctRegion(region);
      float pt = gctTau->rank();
      float eta = convertRCTEta(gctTau->regionId().ieta());
      float phi = convertRCTPhi(gctTau->regionId().iphi());
      //std::cout<<"rank "<<gctTau->rank()<<" gctEta "<<gctTau->regionId().ieta()<<" gctPhi "<<gctTau->regionId().iphi()<<std::endl;

      if(gctTau->regionId().ieta()>17 || gctTau->regionId().ieta()<4) continue;
      if(gctTau->regionId().iphi()>19){std::cout<<"region phi is out of bounds!!"<<std::endl; continue;}
      
      if(gctTau->rank()>0){
 	//std::cout<<"rank "<<gctTau->rank()<<" gctEta "<<gctTau->regionId().ieta()<<" gctPhi "<<gctTau->regionId().iphi()<<std::endl;
        //std::cout<<" pt "<<pt<< " eta "<< eta<< " phi "<<phi<<std::endl;
	
	reco::LeafCandidate::PolarLorentzVector tempLorentz;
	tempLorentz.SetPt(pt);
	tempLorentz.SetEta(eta);
	tempLorentz.SetPhi(phi);

	//std::cout<<"temp lorentz pt  "<<pt<<std::endl;
	//std::cout<<"temp lorentz eta "<<eta<<std::endl;
	//std::cout<<"temp lorentz phi "<<phi<<std::endl;
	
	l1extra::L1JetParticle tempJet(tempLorentz);
	//std::cout<<"temp jet pt  "<<tempJet.pt()<<std::endl;
	//std::cout<<"temp jet eta "<<tempJet.eta()<<std::endl;
	//std::cout<<"temp jet phi "<<tempJet.phi()<<std::endl;

	rlxTauSorted.push_back(tempJet);
	if(abs(eta) < 2.4) rlxTauSortedEtaRestricted2p4.push_back(tempJet);
	if(abs(eta) < 2.1) rlxTauSortedEtaRestricted2p1.push_back(tempJet);
      }
    }
  }
  //else
  //std::cout<<"did not get relaxed taus"<<std::endl;
  
  std::sort(rlxTauSorted.begin(),rlxTauSorted.end(),compareByPt);
  std::sort(rlxTauSortedEtaRestricted2p4.begin(),rlxTauSortedEtaRestricted2p4.end(),compareByPt);
  std::sort(rlxTauSortedEtaRestricted2p1.begin(),rlxTauSortedEtaRestricted2p1.end(),compareByPt);
  std::cout<<"rlxTauSorted size "<<rlxTauSorted.size() <<std::endl;


  //Begin Making Rate Plots
  for( auto rlxTau : rlxTauSorted ) {
    if( rlxTau.pt() > 150 )
      tau_pt->Fill( 150 );
    else tau_pt->Fill( rlxTau.pt() );
    tau_eta->Fill( rlxTau.eta() );
    tau_phi->Fill( rlxTau.phi() );
    std::cout<<"rlx tau pt "<< rlxTau.pt() << " rlx tau eta "<< rlxTau.eta()<<" rlx tau phi "<< rlxTau.phi()<<std::endl;
  }

  if(rlxTauSorted.size()>1)
    tau_pt_diTau->Fill(rlxTauSorted.at(1).pt());
  if(rlxTauSortedEtaRestricted2p4.size()>1)
    tau_pt_diTau_eta2p4->Fill(rlxTauSortedEtaRestricted2p4.at(1).pt());
  if(rlxTauSortedEtaRestricted2p1.size()>1)
    tau_pt_diTau_eta2p1->Fill(rlxTauSortedEtaRestricted2p1.at(1).pt());

  //Now for the iso taus
  vector<l1extra::L1JetParticle> isoTauSorted;
  vector<l1extra::L1JetParticle> isoTauSortedEtaRestricted2p1;
  vector<l1extra::L1JetParticle> isoTauSortedEtaRestricted2p4;
  if(evt.getByToken(l1ExtraIsoTauSource_, l1ExtraIsoTaus)){
    for( vector<l1extra::L1JetParticle>::const_iterator isoTau = l1ExtraIsoTaus->begin(); isoTau != l1ExtraIsoTaus->end(); isoTau++ ) {
      isoTauSorted.push_back(*isoTau);
      if(abs(isoTau->eta()) < 2.4) isoTauSortedEtaRestricted2p4.push_back(*isoTau);
      if(abs(isoTau->eta()) < 2.1) isoTauSortedEtaRestricted2p1.push_back(*isoTau);
    }
  }//in case we're looking at stage 2
  else if(evt.getByToken(l1Stage2TauSource_, stage2Taus)){
    for(BXVector<l1t::Tau>::const_iterator isoTau = stage2Taus->begin(); isoTau != stage2Taus->end(); isoTau++ ) {
      //make this int l1extra::L1JetParticle
      if(isoTau->hwIso() < 1){
	//std::cout<<"isoTau Pt "<<isoTau->pt()<<" eta: "<<isoTau->eta()<<" phi: "<<isoTau->phi()<<std::endl;
	l1extra::L1JetParticle tempJet(isoTau->p4());
	isoTauSorted.push_back(tempJet);
	if(abs(isoTau->eta()) < 2.4) isoTauSortedEtaRestricted2p4.push_back(tempJet);
	if(abs(isoTau->eta()) < 2.1) isoTauSortedEtaRestricted2p1.push_back(tempJet);
      }
    }
  }
  else if(evt.getByToken(l1Stage1IsoTauSource_, stage1IsoTaus)){
    std::cout<<"found stage 1 taus size: "<< stage1IsoTaus->size() <<std::endl;
    for(BXVector<l1t::Tau>::const_iterator isoTau = stage1IsoTaus->begin(); isoTau != stage1IsoTaus->end(); isoTau++ ) {
      //make this int l1extra::L1JetParticle
      std::cout<<"isoTau Pt "<<isoTau->pt()<<" eta: "<<isoTau->eta()<<" phi: "<<isoTau->phi()<<std::endl;
      l1extra::L1JetParticle tempJet(isoTau->p4());
      isoTauSorted.push_back(tempJet);
      if(abs(isoTau->eta()) < 2.4) isoTauSortedEtaRestricted2p4.push_back(tempJet);
      if(abs(isoTau->eta()) < 2.1) isoTauSortedEtaRestricted2p1.push_back(tempJet);
    }
  }

  std::sort(isoTauSorted.begin(),isoTauSorted.end(),compareByPt);
  std::sort(isoTauSortedEtaRestricted2p4.begin(),isoTauSortedEtaRestricted2p4.end(),compareByPt);
  std::sort(isoTauSortedEtaRestricted2p1.begin(),isoTauSortedEtaRestricted2p1.end(),compareByPt);

  std::cout<<"making taus"<<std::endl;
  for( auto isoTau : isoTauSorted ) {
    if( isoTau.pt() > 150 )
      isoTau_pt->Fill( 150 );
    else
      isoTau_pt->Fill( isoTau.pt() );
    isoTau_eta->Fill( isoTau.eta() );
    isoTau_phi->Fill( isoTau.phi() );
  }
  if(isoTauSorted.size()>1)
    isoTau_pt_diTau->Fill(isoTauSorted.at(1).pt());
  if(isoTauSortedEtaRestricted2p4.size()>1)
    isoTau_pt_diTau_eta2p4->Fill(isoTauSortedEtaRestricted2p4.at(1).pt());
  if(isoTauSortedEtaRestricted2p1.size()>1)
    isoTau_pt_diTau_eta2p1->Fill(isoTauSortedEtaRestricted2p1.at(1).pt());

  //End Making Rate Plots
  std::cout<<"making regions"<<std::endl;
  //************* Get Regions and make region plots
  if(!evt.getByToken(regionSource_,regions)){
    std::cout<<"ERROR GETTING THE REGIONS!!!"<<std::endl;}
  else{
    for(vector<L1CaloRegion>::const_iterator region = regions->begin(); region != regions->end(); ++region){
      UCTRegionProcess uctRegion(*region);
      if(region->et()>0){
 	float pt = (region->et());
	
	float eta = uctRegion.getFineRecoEta();
	
	float phi = uctRegion.getFineRecoPhi();
	regionEta->Fill(eta);
	regionPhi->Fill(phi);
	regionPt->Fill(pt);

      }
    }
  }

  bool testMode = false;
  std::cout<<"looking through the taus"<<std::endl;

  //If there isn't at least 1 good reco tau don't bother doing all the work
  if(goodTaus.size()>0){
    double deltaR_ = 1;
    
    ////Make TPG Maps
    //ecal
    //vector<vector<double> > eTowerETMap(75, vector<double>(65));
    double eTowerETMap[73][57]={{0}};
    //hcal
    //vector<vector<double> > hTowerETMap(75, vector<double>(65));
    double hTowerETMap[73][57]={{0}};

    initializeECALTPGMap( ecalTPGs, eTowerETMap, testMode );
    initializeHCALTPGMap( hcalTPGs, hcalScale, hTowerETMap, testMode);
    std::cout<<"Finished initializing the tpg maps"<<std::endl;
    ////Make efficiencies
    for(unsigned int i = 0; i < goodTaus.size(); i++){
      //nvtx = 0;
      if(evt.getByToken(vtxLabel_, vertices)){
	nvtx = (double)vertices->size();
	//std::cout<<"nvtx "<<nvtx<<std::endl;
      }

      pat::Tau recoTau = goodTaus.at(i);
      tauEtaEcalEnt =-999, tauPhiEcalEnt =-999;
      decayMode = -999, jetEt = -999, jetEta = -999, jetPhi = -999, rawEcal = 0, rawHcal = 0,  ecal = 0, hcal = 0;
      
      ////Fill Reco Objects
      recoPt  = recoTau.pt();
      recoEta = recoTau.eta();
      recoPhi = recoTau.phi();
      
      decayMode = recoTau.decayMode();

      //Fill L1 Objects
      l1IsoMatched = -1; l1RlxMatched = -1;
      isoTauPt = 0; isoTauEta = -99; isoTauPhi = -99; 
      rlxTauPt = 0; rlxTauEta = -99; rlxTauPhi = -99;
      
      
      for(uint32_t k = 0; k<isoTauSorted.size(); k++){
	double dR = deltaR( recoTau.p4(), isoTauSorted.at(k).p4());
	if( dR < deltaR_){
	  isoTauPt  = isoTauSorted.at(k).pt();
	  isoTauEta = isoTauSorted.at(k).eta();
	  isoTauPhi = isoTauSorted.at(k).phi();
	  l1IsoMatched = 1;
	  break;
	}
      }
      std::cout<<"rlxTauSorted size "<<rlxTauSorted.size() <<std::endl;

      for( auto rlxTau : rlxTauSorted ) {
	std::cout<<"rlx tau pt "<< rlxTau.pt() << " rlx tau eta "<< rlxTau.eta()<<" rlx tau phi "<< rlxTau.phi()<<std::endl;
	double dR = deltaR( recoTau.p4(), rlxTau.p4());
	if(dR < deltaR_){
	  rlxTauPt  = rlxTau.pt();
	  rlxTauEta = rlxTau.eta();
	  rlxTauPhi = rlxTau.phi();
	  std::cout<< "rlxTauPt" << rlxTauPt <<" rlxTauEta "<< rlxTauEta<<"rlxTauPhi"<<rlxTauPhi<< " DR "<< dR<<std::endl;
	  l1RlxMatched = 1;
	  break;
	}
      }
      
      /*
      if(l1RlxMatched<1){
	std::cout<<"Not Matched!"<<std::endl;
	std::cout<<"recoPt "<<recoPt<< " recoEta "<< recoEta<< " recoPhi "<<recoPhi<<std::endl;
      }
      else{
	std::cout<<"Matched"<<std::endl;
	std::cout<<"recoPt "<<recoPt<< " recoEta "<< recoEta<< " recoPhi "<<recoPhi<<std::endl;
	std::cout<< "rlxTauPt" << rlxTauPt <<" rlxTauEta "<< rlxTauEta<<"rlxTauPhi"<<rlxTauPhi<<std::endl;
	}*/
      efficiencyTree->Fill();
    }
    
  }

}

/*
 * Return DeltaR of three Prong
 * Min prong Pt
 * nSignal Cands
 */
void 
L1TRatesAndEff::getThreeProngInfo(const pat::Tau & tau, double &maxDeltaR, double &minProngPt, double &midProngPt, double &maxProngPt, int &nCands){
  maxDeltaR = 0;
  minProngPt = tau.signalPFCands().at(0)->pt();
  midProngPt = tau.signalPFCands().at(0)->pt();
  maxProngPt = tau.signalPFCands().at(0)->pt();
  int i = 0;
  for(i = 0; i < (int) tau.signalPFCands().size(); i++){
    double dR = deltaR( tau.p4(), tau.signalPFCands().at(i)->p4());
    double ptProng = tau.signalPFCands().at(i)->pt();
    if(dR>maxDeltaR)
      maxDeltaR = dR;

    if(minProngPt>ptProng)
      minProngPt = ptProng;

    if(maxProngPt<ptProng){
      midProngPt = maxProngPt;
      maxProngPt = ptProng;
    }
    
  }
  nCands = i;
}

/*
void
L1TRatesAndEff::getRawEcalHcalEnergy(const pat::PackedCandidate pfCand, double &rawEcal, double &rawHcal, double &ecal, double &hcal){
  //edm::Ptr<reco::PFCandidate> PFCand = tau.leadPFCand();
  //std::cout<<"Pt "<<pfCand->pt()<<std::endl;
  if(pfCand.isNonnull()){
    rawEcal = pfCand->rawEcalEnergy();
    rawHcal = pfCand->rawHcalEnergy();
    ecal = pfCand->ecalEnergy();
    hcal = pfCand->hcalEnergy();
  }
  }*/

//four vector addition vs linear sum
double
L1TRatesAndEff::getPFCandsEt(const std::vector<pat::PackedCandidate> pfCands){
  double etTotal = 0;
  for (uint32_t i = 0; i < pfCands.size(); i++ ) {
    etTotal +=pfCands.at(i).et();
  }
  return etTotal;
}

double
L1TRatesAndEff::getPFCandsEtEtaPhi(edm::Handle<std::vector<pat::PackedCandidate> >& pfCands, const pat::Tau & tau, double dR){
  double etTotal = 0;
  for (uint32_t i = 0; i < pfCands->size(); i++ ) {
    if(reco::deltaR(tau, pfCands->at(i)) < dR){
      etTotal +=pfCands->at(i).et();
    }
  }
  return etTotal;
}


int
L1TRatesAndEff::get2x2TPGs(const int maxTPGPt_eta,const int maxTPGPt_phi, const double eTowerETMap[73][57], const double hTowerETMap[73][57], double &TPGe5x5_, double &TPGh5x5_ ){
  
  //TPG5x5_gcteta_ = twrEta2RegionEta(maxTPGPt_eta);
  //TPG5x5_tpgeta_ = maxTPGPt_eta;
  for (int j = -1; j < 2; ++j) {//phi
    for (int k = -1; k < 2; ++k) { //eta
      int tpgsquarephi= maxTPGPt_phi+j;
      int tpgsquareeta= maxTPGPt_eta+k;
      if (tpgsquarephi==-1) {tpgsquarephi=71;}
      if (tpgsquarephi==-2) {tpgsquarephi=70;}
      if (tpgsquarephi==-3) {tpgsquarephi=69;}
      if (tpgsquarephi==-4) {tpgsquarephi=68;}
      if (tpgsquarephi==-5) {tpgsquarephi=67;}
      if (tpgsquarephi==72) {tpgsquarephi=0;}
      if (tpgsquarephi==73) {tpgsquarephi=1;}
      if (tpgsquarephi==74) {tpgsquarephi=2;}
      if (tpgsquarephi==75) {tpgsquarephi=3;}
      if (tpgsquarephi==76) {tpgsquarephi=4;}
      if (tpgsquareeta>55 || tpgsquareeta<0) {continue;}//No Eta values beyond
      TPGh5x5_ += hTowerETMap[tpgsquarephi][tpgsquareeta];
      TPGe5x5_ += eTowerETMap[tpgsquarephi][tpgsquareeta];
    }
  }
  //std::cout<<"TPGe5x5_ "<<TPGe5x5_<<" TPGh5x5_ "<<TPGh5x5_<<std::endl;
  return (TPGe5x5_ + TPGh5x5_);
}


int
L1TRatesAndEff::get5x5TPGs(const int maxTPGPt_eta,const int maxTPGPt_phi, const double eTowerETMap[73][57], const double hTowerETMap[73][57], double &TPGe5x5_, double &TPGh5x5_ ){
  
  //TPG5x5_gcteta_ = twrEta2RegionEta(maxTPGPt_eta);
  //TPG5x5_tpgeta_ = maxTPGPt_eta;
  for (int j = -5; j < 6; ++j) {//phi
    for (int k = -5; k < 6; ++k) { //eta
      int tpgsquarephi= maxTPGPt_phi+j;
      int tpgsquareeta= maxTPGPt_eta+k;
      if (tpgsquarephi==-1) {tpgsquarephi=71;}
      if (tpgsquarephi==-2) {tpgsquarephi=70;}
      if (tpgsquarephi==-3) {tpgsquarephi=69;}
      if (tpgsquarephi==-4) {tpgsquarephi=68;}
      if (tpgsquarephi==-5) {tpgsquarephi=67;}
      if (tpgsquarephi==72) {tpgsquarephi=0;}
      if (tpgsquarephi==73) {tpgsquarephi=1;}
      if (tpgsquarephi==74) {tpgsquarephi=2;}
      if (tpgsquarephi==75) {tpgsquarephi=3;}
      if (tpgsquarephi==76) {tpgsquarephi=4;}
      if (tpgsquareeta>55 || tpgsquareeta<0) {continue;}//No Eta values beyond
      TPGh5x5_ += hTowerETMap[tpgsquarephi][tpgsquareeta];
      TPGe5x5_ += eTowerETMap[tpgsquarephi][tpgsquareeta];
    }
  }
  //std::cout<<"TPGe5x5_ "<<TPGe5x5_<<" TPGh5x5_ "<<TPGh5x5_<<std::endl;
  return (TPGe5x5_ + TPGh5x5_);
}

int
L1TRatesAndEff::get6x6TPGs(const int maxTPGPt_eta,const int maxTPGPt_phi, const double eTowerETMap[73][57], const double hTowerETMap[73][57], double &TPGe6x6_, double &TPGh6x6_ ){

  for (int j = -6; j < 7; ++j) {//phi
    for (int k = -6; k < 7; ++k) { //eta
      int tpgsquarephi= maxTPGPt_phi+j;
      int tpgsquareeta= maxTPGPt_eta+k;
      if (tpgsquarephi==-1) {tpgsquarephi=71;}
      if (tpgsquarephi==-2) {tpgsquarephi=70;}
      if (tpgsquarephi==-3) {tpgsquarephi=69;}
      if (tpgsquarephi==-4) {tpgsquarephi=68;}
      if (tpgsquarephi==-5) {tpgsquarephi=67;}
      if (tpgsquarephi==72) {tpgsquarephi=0;}
      if (tpgsquarephi==73) {tpgsquarephi=1;}
      if (tpgsquarephi==74) {tpgsquarephi=2;}
      if (tpgsquarephi==75) {tpgsquarephi=3;}
      if (tpgsquarephi==76) {tpgsquarephi=4;}
      if (tpgsquareeta>55 || tpgsquareeta<0) {continue;}//No Eta values beyond
      TPGh6x6_ += hTowerETMap[tpgsquarephi][tpgsquareeta];
      TPGe6x6_ += eTowerETMap[tpgsquarephi][tpgsquareeta];
    }
  }
  //std::cout<<"TPGe6x6_ "<<TPGe6x6_<<" TPGh6x6_ "<<TPGh6x6_<<std::endl;
  return (TPGe6x6_ + TPGh6x6_);
}


int
L1TRatesAndEff::get7x7TPGs(const int maxTPGPt_eta,const int maxTPGPt_phi, const double eTowerETMap[73][57], const double hTowerETMap[73][57], double &TPGe7x7_, double &TPGh7x7_ ){
  for (int j = -7; j < 8; ++j) {//phi
    for (int k = -7; k < 8; ++k) { //eta
      int tpgsquarephi= maxTPGPt_phi+j;
      int tpgsquareeta= maxTPGPt_eta+k;
      if (tpgsquarephi==-1) {tpgsquarephi=71;}
      if (tpgsquarephi==-2) {tpgsquarephi=70;}
      if (tpgsquarephi==-3) {tpgsquarephi=69;}
      if (tpgsquarephi==-4) {tpgsquarephi=68;}
      if (tpgsquarephi==-5) {tpgsquarephi=67;}
      if (tpgsquarephi==72) {tpgsquarephi=0;}
      if (tpgsquarephi==73) {tpgsquarephi=1;}
      if (tpgsquarephi==74) {tpgsquarephi=2;}
      if (tpgsquarephi==75) {tpgsquarephi=3;}
      if (tpgsquarephi==76) {tpgsquarephi=4;}
      if (tpgsquareeta>55 || tpgsquareeta<0) {continue;}//No Eta values beyond
      TPGh7x7_ += hTowerETMap[tpgsquarephi][tpgsquareeta];
      TPGe7x7_ += eTowerETMap[tpgsquarephi][tpgsquareeta];
    }
  }
  //std::cout<<"TPGe7x7_ "<<TPGe7x7_<<" TPGh7x7_ "<<TPGh7x7_<<std::endl;
  return (TPGe7x7_ + TPGh7x7_);
}

  
/*
 * Get the ECAL TPGS create a TPG map for the event
 *
 */
  
void L1TRatesAndEff::initializeECALTPGMap(Handle<EcalTrigPrimDigiCollection> ecal, double eTowerETMap[73][57], bool testMode){
  
  //std::cout << "ECAL TPGS" << std::endl;
  for (size_t i = 0; i < ecal->size(); ++i) {
    int cal_ieta = (*ecal)[i].id().ieta();
    int cal_iphi = (*ecal)[i].id().iphi();
    int iphi = cal_iphi-1;
    int ieta = TPGEtaRange(cal_ieta);
    // TPG iPhi starts at 1 and goes to 72.  Let's index starting at zero.
    // TPG ieta ideal goes from 0-55.
    double LSB = 0.5;
    double et= (*ecal)[i].compressedEt()*LSB;
    //if(et>0)std::cout<<"et "<< et<<std::endl;
    if(testMode && iphi == 34 && ieta == 11){
      et = 40;
    }

    //if(et>0)
    //cout<<"Before filling eTower"
    //<<"ECAL ieta:"<<ieta<<" cal_ieta:"<< cal_ieta<<" iphi:"<<iphi<<" et:"<<et<<endl;

    if (iphi >= 0 && iphi <= 72 &&
	ieta >= 0 && ieta <= 55) {
      eTowerETMap[iphi][ieta] = et; 
    }

  }

}

 void L1TRatesAndEff::initializeHCALTPGMap(const Handle<HcalTrigPrimDigiCollection> hcal, 
					 const ESHandle<L1CaloHcalScale>  hcalScale, 
					 double hTowerETMap[73][57], bool testMode){
  for (size_t i = 0; i < hcal->size(); ++i) {
    HcalTriggerPrimitiveDigi tpg = (*hcal)[i];
    int cal_ieta = tpg.id().ieta();
    int cal_iphi = tpg.id().iphi();
    int iphi = cal_iphi-1;
    int ieta = TPGEtaRange(cal_ieta);
    short absieta = std::abs(tpg.id().ieta());
    short zside = tpg.id().zside();
    double energy = hcalScale->et(tpg.SOI_compressedEt(), absieta, zside); 
    //if(energy>0)std::cout<<"energy "<< energy<<std::endl;

    if(testMode && iphi == 34 && ieta == 12){
      energy = 40;
    }

    if (iphi >= 0 && iphi <= 71 &&
	ieta >= 0 && ieta <= 55) {

      //(*hcal)[i].SOI_compressedEt(), absieta, zside)*LSB; //*LSB
      //if(energy>0)
      //std::cout<<"hcal iphi "<<iphi<<" ieta "<<ieta<<" energy "<<energy<<std::endl;
      hTowerETMap[iphi][ieta] = energy;
      //TPGSum_ +=energy;
      //TPGH_ += energy;
      //double alpha_h = TPGSFp_[cal_ieta]; //v3
      //hCorrTowerETMap[cal_iphi][cal_ieta] = alpha_h*energy;
      //cTPGH_ += alpha_h*energy;
      //if (energy > 0) {
      //std::cout << "hcal eta/phi=" << ieta << "/" << iphi
      //<< " = (" << getEtaTPG(ieta) << "/" << getPhiTPG(iphi) << ") "
      //<< " et=" << (*hcal)[i].SOI_compressedEt()
      //<< " energy=" << energy
      //<< " rctEta="<< twrEta2RegionEta(cal_ieta) << " rctPhi=" << twrPhi2RegionPhi(cal_iphi)
      //<< " fg=" << (*hcal)[i].SOI_fineGrain() << std::endl;
      //}
      //if (energy>maxTPGHPt){
      //maxTPGHPt=energy;
      //maxTPGHPt_phi = cal_iphi; //this one starts at 0-72
      //maxTPGHPt_eta = cal_ieta; //this one is 0-54
      //} 
    }
    //else
      //std::cout<<"HCAL failed checks iphi "<<iphi<<" ieta "<<ieta<<std::endl;
  }//end HCAL TPG
}

void L1TRatesAndEff::endJob() {
}

L1TRatesAndEff::~L1TRatesAndEff(){

}

DEFINE_FWK_MODULE(L1TRatesAndEff);
