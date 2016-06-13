/*
 * \file triggerTestArea.cc
 *
 * \author I. Ojalvo
 * Written for miniAOD
 */

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "L1Trigger/Stage3Ntuplizer/interface/triggerTestArea.h"
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
//vector<L1CaloRegion>                  "l1tCaloLayer1Digis"

bool compareByPtTrig(l1extra::L1JetParticle i,l1extra::L1JetParticle j) { return(i.pt()>j.pt()); };
bool compareByPtTower(tower i,tower j) { return(i.EpH>j.EpH); };

triggerTestArea::triggerTestArea( const ParameterSet & cfg ) :
  rctSource_(cfg.getParameter<edm::InputTag>("rctSource")),
  pfCandsToken_(consumes<vector<pat::PackedCandidate> >(cfg.getParameter<edm::InputTag>("pfCands"))),
  ecalSrc_(consumes<EcalTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("ecalDigis"))),
  hcalSrc_(consumes<HcalTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("hcalDigis"))),
  vtxLabel_(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("vertices"))),
  discriminatorMu_(consumes<reco::PFTauDiscriminator>(cfg.getParameter<edm::InputTag>("recoTauDiscriminatorMu"))),
  discriminatorIso_(consumes<reco::PFTauDiscriminator>(cfg.getParameter<edm::InputTag>("recoTauDiscriminatorIso"))),
  tauSrc_(consumes<vector<pat::Tau> >(cfg.getParameter<edm::InputTag>("recoTau"))),
  jetSrc_(consumes<vector<pat::Jet> >(cfg.getParameter<edm::InputTag>("recoJets"))),
  l1ExtraIsoTauSource_(consumes<vector <l1extra::L1JetParticle> >(cfg.getParameter<edm::InputTag>("l1ExtraIsoTauSource"))),
  l1ExtraTauSource_(consumes<vector <l1extra::L1JetParticle> >(cfg.getParameter<edm::InputTag>("l1ExtraTauSource"))),
  l1Stage2TauSource_(consumes<BXVector<l1t::Tau> > (cfg.getParameter<edm::InputTag>("stage2TauSource"))),
  regionSource_(consumes<vector <L1CaloRegion> >(cfg.getParameter<edm::InputTag>("UCTRegion")))
  {

    folderName_          = cfg.getUntrackedParameter<std::string>("folderName");
    recoPt_              = cfg.getParameter<double>("recoPtCut");
    folder               = tfs_->mkdir("recoTaus");
    efficiencyTreeTaus = folder.make<TTree>("EfficiencyTreeTaus", "Efficiency Tree Taus");
    efficiencyTreeTaus->Branch("run",    &run,     "run/I");
    efficiencyTreeTaus->Branch("lumi",   &lumi,    "lumi/I");
    efficiencyTreeTaus->Branch("event",  &event,   "event/I");
    efficiencyTreeTaus->Branch("nvtx",   &nvtx,         "nvtx/I");
    
    efficiencyTreeTaus->Branch("recoPt",    &recoPt,   "recoPt/D");
    efficiencyTreeTaus->Branch("recoEta",   &recoEta,   "recoEta/D");
    efficiencyTreeTaus->Branch("recoPhi",   &recoPhi,   "recoPhi/D");
    efficiencyTreeTaus->Branch("decayMode", &decayMode,   "decayMode/I");
    //efficiencyTreeTaus->Branch("jetEt",         &jetEt,        "jetEt/D");
    //efficiencyTreeTaus->Branch("jetEta",        &jetEta,       "jetEta/D");
    //efficiencyTreeTaus->Branch("jetPhi",        &jetPhi,       "jetPhi/D");
    efficiencyTreeTaus->Branch("rlxTauPt",    &rlxTauPt,    "rlxTauPt/D");
    efficiencyTreeTaus->Branch("rlxTauEta",   &rlxTauEta,   "rlxTauEta/D");
    efficiencyTreeTaus->Branch("rlxTauPhi",   &rlxTauPhi,   "rlxTauPhi/D");

    efficiencyTreeTaus->Branch("isoTauPt",    &isoTauPt,    "isoTauPt/D");
    efficiencyTreeTaus->Branch("isoTauEta",   &isoTauEta,   "isoTauEta/D");
    efficiencyTreeTaus->Branch("isoTauPhi",   &isoTauPhi,   "isoTauPhi/D");

    //to calculate
    efficiencyTreeTaus->Branch("ecalTPG3x3",        &TPGE3x3,          "ecalTPG3x3/D");
    efficiencyTreeTaus->Branch("hcalTPG3x3",        &TPGH3x3,          "hcalTPG3x3/D");
    efficiencyTreeTaus->Branch("TPG3x3",            &TPG3x3,           "TPG3x3/D");
    efficiencyTreeTaus->Branch("EoHTPG3x3",         &EoHTPG3x3,        "EoHTPG3x3/D");
    efficiencyTreeTaus->Branch("HoETPG3x3",         &HoETPG3x3,        "HoETPG3x3/D");
    efficiencyTreeTaus->Branch("HoEpHTPG3x3",       &HoEpHTPG3x3,      "HoEpHTPG3x3/D");
    efficiencyTreeTaus->Branch("EoEpHTPG3x3",       &EoEpHTPG3x3,      "EoEpHTPG3x3/D");

    efficiencyTreeTaus->Branch("EoH_highPtTower",   &EoH_highPtTower,  "EoH_highPtTower/D");
    efficiencyTreeTaus->Branch("HoE_highPtTower",   &HoE_highPtTower,  "HoE_highPtTower/D");
    efficiencyTreeTaus->Branch("EoEpH_highPtTower", &EoEpH_highPtTower, "EoEpH_highPtTower/D");
    efficiencyTreeTaus->Branch("HoEpH_highPtTower", &HoEpH_highPtTower, "HoEpH_highPtTower/D");

    efficiencyTreeTaus->Branch("ecalTPG_nearTower", &TPGE_nearTower,      "ecalTPG_nearTower/D");
    efficiencyTreeTaus->Branch("hcalTPG_nearTower", &TPGH_nearTower,      "hcalTPG_nearTower/D");
    efficiencyTreeTaus->Branch("TPG_nearTower",     &TPG_nearTower,       "TPG_nearTower/D");

    efficiencyTreeTaus->Branch("EoH_nearTower",     &EoH_nearTower,    "EoH_nearTower/D");
    efficiencyTreeTaus->Branch("HoE_nearTower",     &HoE_nearTower,    "HoE_nearTower/D");
    efficiencyTreeTaus->Branch("EoEpH_nearTower",   &EoEpH_nearTower,  "EoEpH_nearTower/D");
    efficiencyTreeTaus->Branch("HoEpH_nearTower",   &HoEpH_nearTower,  "HoEpH_nearTower/D");
						    
    efficiencyTreeTaus->Branch("EoH_eAndhTower",    &EoH_eAndhTower,   "EoH_eAndhTower/D");
    efficiencyTreeTaus->Branch("HoE_eAndhTower",    &HoE_eAndhTower,   "HoE_eAndhTower/D");
    efficiencyTreeTaus->Branch("EoEpH_eAndhTower",  &EoEpH_eAndhTower, "EoEpH_eAndhTower/D");
    efficiencyTreeTaus->Branch("HoEpH_eAndhTower",  &HoEpH_eAndhTower, "HoEpH_eAndhTower/D");

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
    regionPt    = folder.make<TH1F>( "region_pt"  , "pt", 100, 0, 100. );

    regionEtaFine   = folder.make<TH1F>( "region_eta_Fine"  , "eta", 88, 1, 88. );
    regionPhiFine   = folder.make<TH1F>( "region_phi_Fine"  , "phi", 72, 1, 72. );


    ////////////////////////////////Now for the jets
    folderJets               = tfs_->mkdir("recoJets");
    //cout<<"afterfolder jets path "<<(gDirectory->GetPath())<<endl;
    //folderJets.cd();
    efficiencyTreeJets = folderJets.make<TTree>("EfficiencyTreeJets","Efficiency Tree Jets");
    efficiencyTreeJets->Branch("run",    &runJ,     "run/I");
    efficiencyTreeJets->Branch("lumi",   &lumiJ,    "lumi/I");
    efficiencyTreeJets->Branch("event",  &eventJ,   "event/I");
    efficiencyTreeJets->Branch("nvtx",   &nvtxJ,         "nvtx/I");
    
    efficiencyTreeJets->Branch("recoPt",    &recoPtJ,   "recoPt/D");
    efficiencyTreeJets->Branch("recoEta",   &recoEtaJ,   "recoEta/D");
    efficiencyTreeJets->Branch("recoPhi",   &recoPhiJ,   "recoPhi/D");
    efficiencyTreeJets->Branch("decayMode", &decayModeJ,   "decayMode/I");
    //efficiencyTreeJets->Branch("jetEt",         &jetEt,        "jetEt/D");
    //efficiencyTreeJets->Branch("jetEta",        &jetEta,       "jetEta/D");
    //efficiencyTreeJets->Branch("jetPhi",        &jetPhi,       "jetPhi/D");

    efficiencyTreeJets->Branch("rlxTauPt",    &rlxTauPtJ,    "rlxTauPt/D");
    efficiencyTreeJets->Branch("rlxTauEta",   &rlxTauEtaJ,   "rlxTauEta/D");
    efficiencyTreeJets->Branch("rlxTauPhi",   &rlxTauPhiJ,   "rlxTauPhi/D");

    efficiencyTreeJets->Branch("isoTauPt",    &isoTauPtJ,    "isoTauPt/D");
    efficiencyTreeJets->Branch("isoTauEta",   &isoTauEtaJ,   "isoTauEta/D");
    efficiencyTreeJets->Branch("isoTauPhi",   &isoTauPhiJ,   "isoTauPhi/D");

    //to calculate
    efficiencyTreeJets->Branch("ecalTPG3x3",    &TPGE3x3J,      "ecalTPG3x3/D");
    efficiencyTreeJets->Branch("hcalTPG3x3",    &TPGH3x3J,      "hcalTPG3x3/D");
    efficiencyTreeJets->Branch("TPG3x3",        &TPG3x3J,       "TPG3x3/D");
    efficiencyTreeJets->Branch("EoHTPG3x3",     &EoHTPG3x3J,    "EoHTPG3x3/D");
    efficiencyTreeJets->Branch("HoETPG3x3",     &HoETPG3x3J,    "HoETPG3x3/D");
    efficiencyTreeJets->Branch("HoEpHTPG3x3",   &HoEpHTPG3x3J,  "HoEpHTPG3x3/D");
    efficiencyTreeJets->Branch("EoEpHTPG3x3",   &EoEpHTPG3x3J,  "EoEpHTPG3x3/D");

    efficiencyTreeJets->Branch("EoH_highPtTower",  &EoH_highPtTowerJ, "EoH_highPtTower/D");
    efficiencyTreeJets->Branch("HoE_highPtTower",  &HoE_highPtTowerJ, "HoE_highPtTower/D");
    efficiencyTreeJets->Branch("EoEpH_highPtTower",&EoEpH_highPtTowerJ, "EoEpH_highPtTower/D");
    efficiencyTreeJets->Branch("HoEpH_highPtTower",&HoEpH_highPtTowerJ, "HoEpH_highPtTower/D");

    efficiencyTreeJets->Branch("ecalTPG_nearTower",    &TPGE_nearTowerJ,      "ecalTPG_nearTower/D");
    efficiencyTreeJets->Branch("hcalTPG_nearTower",    &TPGH_nearTowerJ,      "hcalTPG_nearTower/D");
    efficiencyTreeJets->Branch("TPG_nearTower",        &TPG_nearTowerJ,       "TPG_nearTower/D");

    efficiencyTreeJets->Branch("EoH_nearTower",  &EoH_nearTowerJ, "EoH_nearTower/D");
    efficiencyTreeJets->Branch("HoE_nearTower",  &HoE_nearTowerJ, "HoE_nearTower/D");
    efficiencyTreeJets->Branch("EoEpH_nearTower",&EoEpH_nearTowerJ, "EoEpH_nearTower/D");
    efficiencyTreeJets->Branch("HoEpH_nearTower",&HoEpH_nearTowerJ, "HoEpH_nearTower/D");

    efficiencyTreeJets->Branch("EoH_eAndhTower",  &EoH_eAndhTowerJ, "EoH_eAndhTower/D");
    efficiencyTreeJets->Branch("HoE_eAndhTower",  &HoE_eAndhTowerJ, "HoE_eAndhTower/D");
    efficiencyTreeJets->Branch("EoEpH_eAndhTower",&EoEpH_eAndhTowerJ, "EoEpH_eAndhTower/D");
    efficiencyTreeJets->Branch("HoEpH_eAndhTower",&HoEpH_eAndhTowerJ, "HoEpH_eAndhTower/D");


  }

void triggerTestArea::beginJob( const EventSetup & es) {
   std::cout<<"begin job..."<<std::endl;
}

void triggerTestArea::analyze( const Event& evt, const EventSetup& es )
 {
   std::cout<<"Analyzing..."<<std::endl;
   nEvents->Fill(1);
   
   run = evt.id().run();
   lumi = evt.id().luminosityBlock();
   event = evt.id().event();
   edm::Handle<reco::VertexCollection> vertices;
   if(evt.getByToken(vtxLabel_, vertices)){
     nvtx = vertices->size();
   }
   
   Handle<L1CaloRegionCollection> regions;
   
   edm::Handle<reco::PFTauDiscriminator> discriminatorIso;
   edm::Handle<reco::PFTauDiscriminator> discriminatorMu;
   
   edm::Handle < L1GctJetCandCollection > l1IsoTauJets;
   edm::Handle < L1GctJetCandCollection > l1TauJets;
   edm::Handle < BXVector<l1t::Tau> > stage2Taus;

  edm::Handle < vector<l1extra::L1JetParticle> > l1ExtraTaus;
  edm::Handle < vector<l1extra::L1JetParticle> > l1ExtraIsoTaus;
  
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
      if(tau.tauID("decayModeFinding")>0.5&&tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits")<20&&tau.tauID("againstMuonLoose3")>0.5&&tau.tauID("againstElectronTightMVA5")>0.5){
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

  /// get recojets
  Handle<vector<pat::Jet> > jets;
  std::vector<pat::Jet> recoJets;
  if(evt.getByToken(jetSrc_, jets)){//Begin Getting Reco Jets
    ///////////
    for (const pat::Jet &jet : *jets) {
      for(auto goodTau : goodTaus)
	//take all jets but cross clean them
	if(reco::deltaR(goodTau.p4(), jet.p4()) > 0.5)
	  recoJets.push_back(jet);
    }
  }  


  ////
  vector<l1extra::L1JetParticle> rlxTauSorted;
  vector<l1extra::L1JetParticle> rlxTauSortedEtaRestricted2p1;
  vector<l1extra::L1JetParticle> rlxTauSortedEtaRestricted2p4;
  if(evt.getByToken(l1ExtraTauSource_, l1ExtraTaus)){
    std::cout<<"found stage 3 taus"<<std::endl;
    for( vector<l1extra::L1JetParticle>::const_iterator rlxTau = l1ExtraTaus->begin(); rlxTau != l1ExtraTaus->end(); rlxTau++ ){
      rlxTauSorted.push_back(*rlxTau);
      if(abs(rlxTau->eta()) < 2.4) rlxTauSortedEtaRestricted2p4.push_back(*rlxTau);
      if(abs(rlxTau->eta()) < 2.1) rlxTauSortedEtaRestricted2p1.push_back(*rlxTau);
    }
  }
  
  std::sort(rlxTauSorted.begin(),rlxTauSorted.end(),compareByPtTrig);
  std::sort(rlxTauSortedEtaRestricted2p4.begin(),rlxTauSortedEtaRestricted2p4.end(),compareByPtTrig);
  std::sort(rlxTauSortedEtaRestricted2p1.begin(),rlxTauSortedEtaRestricted2p1.end(),compareByPtTrig);

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
  }

  std::sort(isoTauSorted.begin(),isoTauSorted.end(),compareByPtTrig);
  std::sort(isoTauSortedEtaRestricted2p4.begin(),isoTauSortedEtaRestricted2p4.end(),compareByPtTrig);
  std::sort(isoTauSortedEtaRestricted2p1.begin(),isoTauSortedEtaRestricted2p1.end(),compareByPtTrig);

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
  //std::cout<<"looking through the taus"<<std::endl;

  //folder.cd();
  //If there isn't at least 1 good reco tau or jet don't bother doing all the work
  if(goodTaus.size()>0||recoJets.size()>0){
    std::cout<<"goodTaus.size() "<<goodTaus.size()<<std::endl;
    std::cout<<"recoJets.size() "<<recoJets.size()<<std::endl;
    double deltaR_ = 0.5;
    ////Make TPG Maps
    //ecal
    double eTowerETMap[73][57]={{0}};
    //hcal
    double hTowerETMap[73][57]={{0}};

    initializeECALTPGMap( ecalTPGs, eTowerETMap, testMode );
    initializeHCALTPGMap( hcalTPGs, hcalScale, hTowerETMap, testMode);
    //std::cout<<"Finished initializing the tpg maps"<<std::endl;

    ////Make efficiencies
    for(unsigned int i = 0; i < goodTaus.size(); i++){
      //std::cout<<"start looking through taus"<<std::endl;
      pat::Tau recoTau = goodTaus.at(i);
      float towerEta = -10;
      float towerPhi = -10;

      ////Fill Reco Objects
      recoPt  = recoTau.pt();
      recoEta = recoTau.eta();
      recoPhi = recoTau.phi();
      TPGE3x3=0;TPGH3x3=0;TPG3x3=0;EoHTPG3x3=0;HoETPG3x3=0;HoEpHTPG3x3=0;EoEpHTPG3x3=0;
      EoH_highPtTower=0;HoE_highPtTower=0;EoEpH_highPtTower=0;HoEpH_highPtTower=0;
      TPGE_nearTower=0;TPGH_nearTower=0;TPG_nearTower=0;EoH_nearTower=0;HoE_nearTower=0;EoEpH_nearTower=0;HoEpH_nearTower=0;
      EoH_eAndhTower=0;HoE_eAndhTower=0;EoEpH_eAndhTower=0; HoEpH_eAndhTower=0;      
      
      //std::cout<<"recoEta "<<recoEta<<" recoPhi "<<recoPhi<<std::endl;
      if(recoEta <= -2.4 || recoEta >= 2.4 || recoPhi <= -3.2 || recoPhi > 3.2 ){
	std::cout<<"recoEta "<<recoEta<<" recoPhi "<<recoPhi<<std::endl;
	continue;
      }
      //std::cout<<"tau passed reco cuts"<<std::endl;
      towerEta = convertGenEta(recoEta);
      towerPhi = convertGenPhi(recoPhi);
      std::cout<<"towerEta "<<towerEta<<" towerPhi "<<towerPhi<<std::endl;

      decayMode = recoTau.decayMode();
      std::vector<tower> towers_nearTower;
      towers_nearTower.clear();
      TPG_nearTower = (float) get1x1TPGs(towerEta, towerPhi, eTowerETMap, hTowerETMap, TPGE_nearTower, TPGH_nearTower , towers_nearTower);

      if(TPGH_nearTower>0)
	EoH_nearTower =  TPGE_nearTower/TPGH_nearTower;
      if(TPGE_nearTower>0)
	HoE_nearTower =  TPGH_nearTower/TPGE_nearTower;
      if(TPGE_nearTower+TPGH_nearTower>0)
	HoEpH_nearTower = TPGH_nearTower/(TPGE_nearTower+TPGH_nearTower);
      if(TPGE_nearTower+TPGH_nearTower>0)
	EoEpH_nearTower = TPGE_nearTower/(TPGE_nearTower+TPGH_nearTower);

      std::vector<tower> towers;
      towers.clear();
      TPG3x3 = get3x3TPGs(towerEta,towerPhi, eTowerETMap, hTowerETMap, TPGE3x3, TPGH3x3 , towers);
      //std::cout<<"TPG3x3 "<<TPG3x3<<std::endl;
      TPG3x3 = TPGE3x3 + TPGH3x3;
      //std::cout<<"TPG3x3 "<<TPG3x3<<std::endl;
      if(TPGH3x3>0)
	EoHTPG3x3 =  TPGE3x3/TPGH3x3;
      if(TPGE3x3>0)
	HoETPG3x3 =  TPGH3x3/TPGE3x3;
      if((TPGE3x3+TPGH3x3)>0)
	HoEpHTPG3x3 = TPGH3x3/(TPGE3x3+TPGH3x3);
      if((TPGE3x3+TPGH3x3)>0)
	EoEpHTPG3x3 = TPGE3x3/(TPGE3x3+TPGH3x3);
      //sort by pt
      std::sort(towers.begin(),towers.end(),compareByPtTower);
      EoH_highPtTower = towers.at(0).EoH;
      HoE_highPtTower = towers.at(0).HoE;

      EoEpH_highPtTower = towers.at(0).EoEpH;
      HoEpH_highPtTower = towers.at(0).HoEpH;

      //sort by pt look for e and h > 0
      for(auto iTower : towers){
	if(iTower.ecalEt>0&&iTower.hcalEt>0){
	  EoH_eAndhTower   = iTower.EoH;
	  HoE_eAndhTower   = iTower.HoEpH;
          EoEpH_eAndhTower = iTower.EoEpH;
          HoEpH_eAndhTower = iTower.HoEpH;
	  break;
	}
      }

      //Fill L1 Objects
      isoTauPt = 0; isoTauEta = -99; isoTauPhi = -99; 
      rlxTauPt = 0; rlxTauEta = -99; rlxTauPhi = -99;
            
      for(uint32_t k = 0; k<isoTauSorted.size(); k++){
	double dR = deltaR( recoTau.p4(), isoTauSorted.at(k).p4());
	if( dR < deltaR_){
	  isoTauPt  = isoTauSorted.at(k).pt();
	  isoTauEta = isoTauSorted.at(k).eta();
	  isoTauPhi = isoTauSorted.at(k).phi();
	  break;
	}
      }
      for(uint32_t k = 0; k<rlxTauSorted.size(); k++){
	double dR = deltaR( recoTau.p4(), rlxTauSorted.at(k).p4());
	if(dR < deltaR_){
	  rlxTauPt  = rlxTauSorted.at(k).pt();
	  rlxTauEta = rlxTauSorted.at(k).eta();
	  rlxTauPhi = rlxTauSorted.at(k).phi();
	  break;
	}
      }


      //std::cout<<"filling tau tree"<<std::endl;
      efficiencyTreeTaus->Fill();
    }

        ////Make efficiencies
    for(unsigned int i = 0; i < recoJets.size(); i++){
      //std::cout<<"start looking through jets"<<std::endl;
      pat::Jet recoJet = recoJets.at(i);
      float towerEta = -10;
      float towerPhi = -10;

      ////Fill Reco Objects
      recoPtJ  = recoJet.pt();
      recoEtaJ = recoJet.eta();
      recoPhiJ = recoJet.phi();      
      TPGE3x3J=0;TPGH3x3J=0;TPG3x3J=0;EoHTPG3x3J=0;HoETPG3x3J=0;HoEpHTPG3x3J=0;EoEpHTPG3x3J=0;
      EoH_highPtTowerJ=0;HoE_highPtTowerJ=0;EoEpH_highPtTowerJ=0;HoEpH_highPtTowerJ=0;
      TPGE_nearTowerJ=0;TPGH_nearTowerJ=0;TPG_nearTowerJ=0;EoH_nearTowerJ=0;HoE_nearTowerJ=0;EoEpH_nearTowerJ=0;HoEpH_nearTowerJ=0;
      EoH_eAndhTowerJ=0;HoE_eAndhTowerJ=0;EoEpH_eAndhTowerJ=0; HoEpH_eAndhTowerJ=0;      
      if(recoEtaJ <= -2.4 || recoEtaJ >= 2.4 || recoPhiJ <= -3.2 || recoPhiJ > 3.2 )
	continue;
      towerEta = convertGenEta(recoEtaJ);
      towerPhi = convertGenPhi(recoPhiJ);
      //std::cout<<"putting all the towers in"<<std::endl;
      //if(towerEta <= -2.4 || towerEta >= 2.4 || towerPhi <= -3.2 || towerPhi > 3.2 )
      //continue;
      decayModeJ = -99;
      //now clear everything
      TPGE_nearTowerJ = 0; TPGH_nearTowerJ = 0; 
      TPG3x3J=0; TPGH3x3J=0; TPGE3x3J=0;
      EoHTPG3x3J=0;  HoETPG3x3J=0;  HoEpHTPG3x3J=0;  EoEpHTPG3x3J=0;
      
      TPG_highPtTowerJ=0; TPGH_highPtTowerJ=0; TPGE_highPtTowerJ=0;
      EoH_highPtTowerJ=0; HoE_highPtTowerJ=0; EoEpH_highPtTowerJ=0; HoEpH_highPtTowerJ=0;
      
      TPG_nearTowerJ=0; TPGH_nearTowerJ=0; TPGE_nearTowerJ=0;
      EoH_nearTowerJ=0; HoE_nearTowerJ=0; EoEpH_nearTowerJ=0; HoEpH_nearTowerJ=0;
      
      EoH_eAndhTowerJ=0; HoE_eAndhTowerJ=0; EoEpH_eAndhTowerJ=0; HoEpH_eAndhTowerJ=0;
      std::vector<tower> towers_nearTower;

      TPG_nearTowerJ = (float)get1x1TPGs(towerEta, towerPhi, eTowerETMap, hTowerETMap, TPGE_nearTowerJ, TPGH_nearTowerJ , towers_nearTower);

      if(TPGH_nearTowerJ>0)
	EoH_nearTowerJ =  ((float)TPGE_nearTowerJ)/((float)TPGH_nearTowerJ);
      if(TPGE_nearTowerJ>0)
	HoE_nearTowerJ =  ((float)TPGH_nearTowerJ)/((float)TPGE_nearTowerJ);
      if(TPGE_nearTowerJ+TPGH_nearTowerJ>0)
	HoEpH_nearTowerJ = ((float)TPGH_nearTowerJ)/((float)TPGE_nearTowerJ+(float)TPGH_nearTowerJ);
      if(TPGE_nearTowerJ+TPGH_nearTowerJ>0)
	EoEpH_nearTowerJ = ((float)TPGE_nearTowerJ)/((float)TPGE_nearTowerJ+(float)TPGH_nearTowerJ);

      std::vector<tower> towers;

      TPG3x3J = get3x3TPGs(towerEta,towerPhi, eTowerETMap, hTowerETMap, TPGE3x3J, TPGH3x3J, towers);
      TPG3x3J = TPGE3x3J + TPGH3x3J;
      if(TPGH3x3J>0)
	EoHTPG3x3J =  TPGE3x3J/TPGH3x3J;
      if(TPGE3x3J>0)
	HoETPG3x3J =  TPGH3x3J/TPGE3x3J;
      if((TPGE3x3J+TPGH3x3J)>0)
	HoEpHTPG3x3J = TPGH3x3J/(TPGE3x3J+TPGH3x3J);
      if((TPGE3x3J+TPGH3x3J)>0)
	EoEpHTPG3x3J = TPGE3x3J/(TPGE3x3J+TPGH3x3J);
      //sort by pt
      std::sort(towers.begin(),towers.end(),compareByPtTower);
      EoH_highPtTowerJ = (float)towers.at(0).EoH;
      HoE_highPtTowerJ = (float)towers.at(0).HoE;
      EoEpH_highPtTowerJ = (float)towers.at(0).EoEpH;
      HoEpH_highPtTowerJ = (float)towers.at(0).HoEpH;

      //sort by pt look for e and h > 0
      for(auto iTower : towers){
	if(iTower.ecalEt>0&&iTower.hcalEt>0){
	  EoH_eAndhTowerJ   = (float)iTower.EoH;
	  HoE_eAndhTowerJ   = (float)iTower.HoEpH;
          EoEpH_eAndhTowerJ = (float)iTower.EoEpH;
          HoEpH_eAndhTowerJ= (float)iTower.HoEpH;
	  break;
	}
      }
      
      //Fill L1 Objects
      isoTauPtJ = 0; isoTauEtaJ = -99; isoTauPhiJ = -99; 
      rlxTauPtJ = 0; rlxTauEtaJ = -99; rlxTauPhiJ = -99;
      
      for(uint32_t k = 0; k<isoTauSorted.size(); k++){
	double dR = deltaR( recoJet.p4(), isoTauSorted.at(k).p4());
	if( dR < deltaR_){
	  isoTauPtJ  = isoTauSorted.at(k).pt();
	  isoTauEtaJ = isoTauSorted.at(k).eta();
	  isoTauPhiJ = isoTauSorted.at(k).phi();
	  break;
	}
      }
      for(uint32_t k = 0; k<rlxTauSorted.size(); k++){
	double dR = deltaR( recoJet.p4(), rlxTauSorted.at(k).p4());
	if(dR < deltaR_){
	  rlxTauPtJ  = rlxTauSorted.at(k).pt();
	  rlxTauEtaJ = rlxTauSorted.at(k).eta();
	  rlxTauPhiJ = rlxTauSorted.at(k).phi();
	  break;
	}
      }

      efficiencyTreeJets->Fill();
    }
    
  }
}



/*
 * Return DeltaR of three Prong
 * Min prong Pt
 * nSignal Cands
 */
void 
triggerTestArea::getThreeProngInfo(const pat::Tau & tau, double &maxDeltaR, double &minProngPt, double &midProngPt, double &maxProngPt, int &nCands){
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


//four vector addition vs linear sum
double
triggerTestArea::getPFCandsEt(const std::vector<pat::PackedCandidate> pfCands){
  double etTotal = 0;
  for (uint32_t i = 0; i < pfCands.size(); i++ ) {
    etTotal +=pfCands.at(i).et();
  }
  return etTotal;
}

double
triggerTestArea::getPFCandsEtEtaPhi(edm::Handle<std::vector<pat::PackedCandidate> >& pfCands, const pat::Tau & tau, double dR){
  double etTotal = 0;
  for (uint32_t i = 0; i < pfCands->size(); i++ ) {
    if(reco::deltaR(tau, pfCands->at(i)) < dR){
      etTotal +=pfCands->at(i).et();
    }
  }
  return etTotal;
}


int
triggerTestArea::get3x3TPGs(const int maxTPGPt_eta,const int maxTPGPt_phi, const double eTowerETMap[73][57], const double hTowerETMap[73][57], double &TPGe5x5_, double &TPGh5x5_ , std::vector<tower> &towers){
  TPGh5x5_=0;
  TPGe5x5_=0;
  //TPG5x5_gcteta_ = twrEta2RegionEta(maxTPGPt_eta);
  //TPG5x5_tpgeta_ = maxTPGPt_eta;
  for (int j = -2; j < 1; ++j) {//phi
    //std::cout<<"iphi "<<maxTPGPt_phi+j<<std::endl;
    for (int k = -1; k < 2; ++k) { //eta
      int tpgsquarephi= maxTPGPt_phi+j;
      int tpgsquareeta= maxTPGPt_eta+k;
      //std::cout<<"eta "<<tpgsquareeta<<" phi "<<tpgsquarephi<<std::endl;
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
      //std::cout<<"hcal tpg "<<hTowerETMap[tpgsquarephi][tpgsquareeta]<<" ecal tpg "<<eTowerETMap[tpgsquarephi][tpgsquareeta]<<std::endl;
      tower temp;
      temp.ecalEt = eTowerETMap[tpgsquarephi][tpgsquareeta];
      temp.hcalEt = hTowerETMap[tpgsquarephi][tpgsquareeta];
      temp.EpH = (float)temp.ecalEt + (float)temp.hcalEt;
      if(temp.hcalEt>0)
	temp.EoH = ((float)temp.ecalEt)/(float)temp.hcalEt;
      else temp.EoH = 0;
      if(temp.ecalEt>0)
	temp.HoE = ((float)temp.hcalEt)/(float)temp.ecalEt;      
      else temp.HoE = 0;
      if(temp.ecalEt+temp.hcalEt > 0)
	temp.HoEpH = ((float)temp.hcalEt)/((float)temp.ecalEt+(float)temp.hcalEt);
      else temp.HoEpH = 0;
      if(temp.ecalEt+temp.hcalEt > 0)
	temp.EoEpH = ((float)temp.ecalEt)/((float)temp.ecalEt+(float)temp.hcalEt);
      else temp.EoEpH = 0;
      towers.push_back(temp);
    }
  }

  return (TPGe5x5_ + TPGh5x5_);
}


int
triggerTestArea::get1x1TPGs(const int maxTPGPt_eta,const int maxTPGPt_phi, const double eTowerETMap[73][57], const double hTowerETMap[73][57], double &TPGe5x5_, double &TPGh5x5_ , std::vector<tower> &towers){
  TPGh5x5_=0;
  TPGe5x5_=0;
  //TPG5x5_gcteta_ = twrEta2RegionEta(maxTPGPt_eta);
  //TPG5x5_tpgeta_ = maxTPGPt_eta;
  for (int j = 0; j < 1; ++j) {//phi
    for (int k = 0; k < 1; ++k) { //eta
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
      tower temp;
      temp.ecalEt = eTowerETMap[tpgsquarephi][tpgsquareeta];
      temp.hcalEt = hTowerETMap[tpgsquarephi][tpgsquareeta];
      temp.EpH = temp.ecalEt + temp.hcalEt;
      if(temp.hcalEt>0)
	temp.EoH = temp.ecalEt/temp.hcalEt;
      else temp.EoH = 0;
      if(temp.ecalEt>0)
	temp.HoE = temp.hcalEt/temp.ecalEt;      
      else temp.HoE = 0;
      if(temp.ecalEt+temp.hcalEt > 0)
	temp.HoEpH = temp.hcalEt/(temp.ecalEt+temp.hcalEt);
      else temp.HoEpH = 0;
      if(temp.ecalEt+temp.hcalEt > 0)
	temp.EoEpH = temp.ecalEt/(temp.ecalEt+temp.hcalEt);
      else temp.EoEpH = 0;

      towers.push_back(temp);
    }
  }

  return (TPGe5x5_ + TPGh5x5_);
}


int
triggerTestArea::get5x5TPGs(const int maxTPGPt_eta,const int maxTPGPt_phi, const double eTowerETMap[73][57], const double hTowerETMap[73][57], double &TPGe5x5_, double &TPGh5x5_ ){
  
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
triggerTestArea::get6x6TPGs(const int maxTPGPt_eta,const int maxTPGPt_phi, const double eTowerETMap[73][57], const double hTowerETMap[73][57], double &TPGe6x6_, double &TPGh6x6_ ){

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
triggerTestArea::get7x7TPGs(const int maxTPGPt_eta,const int maxTPGPt_phi, const double eTowerETMap[73][57], const double hTowerETMap[73][57], double &TPGe7x7_, double &TPGh7x7_ ){
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
  
void triggerTestArea::initializeECALTPGMap(Handle<EcalTrigPrimDigiCollection> ecal, double eTowerETMap[73][57], bool testMode){
  
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

 void triggerTestArea::initializeHCALTPGMap(const Handle<HcalTrigPrimDigiCollection> hcal, 
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

void triggerTestArea::endJob() {
}

triggerTestArea::~triggerTestArea(){

}

DEFINE_FWK_MODULE(triggerTestArea);
