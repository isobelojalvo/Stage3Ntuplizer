/*
 * \file L1TEventDisplayGenerator.cc
 *
 * \author I. Ojalvo
 * Written for miniAOD
 */

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "L1Trigger/Stage3Ntuplizer/interface/L1TEventDisplayGenerator.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "L1Trigger/L1TCaloLayer1/src/UCTGeometry.hh"
#include "L1Trigger/Stage3Ntuplizer/plugins/UCTRegionProcess.hh"
#include "L1Trigger/Stage3Ntuplizer/plugins/triggerGeometryTools.hh"

using namespace edm;
using std::cout;
using std::endl;
using std::vector;

L1TEventDisplayGenerator::L1TEventDisplayGenerator( const ParameterSet & cfg ) :
  rctSource_(cfg.getParameter<edm::InputTag>("rctSource")),
  packedPfCandsToken_(consumes<vector<pat::PackedCandidate> >(cfg.getParameter<edm::InputTag>("packedPfCands"))),
  pfCandsToken_(consumes<vector<reco::PFCandidate> >(cfg.getParameter<edm::InputTag>("pfCands"))),
  ecalSrc_(consumes<EcalTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("ecalDigis"))),
  hcalSrc_(consumes<HcalTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("hcalDigis"))),
  vtxLabel_(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("vertices"))),
  discriminatorMu_(consumes<reco::PFTauDiscriminator>(cfg.getParameter<edm::InputTag>("recoTauDiscriminatorMu"))),
  discriminatorIso_(consumes<reco::PFTauDiscriminator>(cfg.getParameter<edm::InputTag>("recoTauDiscriminatorIso"))),
  tauSrc_(consumes<vector<reco::PFTau> >(cfg.getParameter<edm::InputTag>("recoTaus"))),
  slimmedTauSrc_(consumes<vector<pat::Tau> >(cfg.getParameter<edm::InputTag>("slimmedTaus"))),
  jetSrc_(consumes<vector<pat::Jet> >(cfg.getParameter<edm::InputTag>("recoJets"))),
  l1ExtraIsoTauSource_(consumes<vector <l1extra::L1JetParticle> >(cfg.getParameter<edm::InputTag>("l1ExtraIsoTauSource"))),
  l1ExtraTauSource_(consumes<vector <l1extra::L1JetParticle> >(cfg.getParameter<edm::InputTag>("l1ExtraTauSource"))),
  l1ExtraJetSource_(consumes<vector <l1extra::L1JetParticle> >(cfg.getParameter<edm::InputTag>("l1ExtraJetSource"))),
  regionSource_(consumes<vector <L1CaloRegion> >(cfg.getParameter<edm::InputTag>("UCTRegion"))),
  ecalCaloSrc_(consumes<vector <reco::CaloCluster> >(cfg.getParameter<edm::InputTag>("ecalCaloClusters")))
  {
    folderName_          = cfg.getUntrackedParameter<std::string>("folderName");
    recoPt_              = cfg.getParameter<double>("recoPtCut");
    efficiencyTree = tfs_->make<TTree>("EfficiencyTree", "Efficiency Tree");

    efficiencyTree->Branch("hcalTpgs_Pt",  &hcalTpgs_Pt); 
    efficiencyTree->Branch("hcalTpgs_Eta", &hcalTpgs_Eta); 
    efficiencyTree->Branch("hcalTpgs_Phi", &hcalTpgs_Phi); 

    efficiencyTree->Branch("ecalTpgs_Pt",  &ecalTpgs_Pt); 
    efficiencyTree->Branch("ecalTpgs_Eta", &ecalTpgs_Eta); 
    efficiencyTree->Branch("ecalTpgs_Phi", &ecalTpgs_Phi); 

    efficiencyTree->Branch("sumTpgs_Pt",  &sumTpgs_Pt); 
    efficiencyTree->Branch("sumTpgs_Eta", &sumTpgs_Eta); 
    efficiencyTree->Branch("sumTpgs_Phi", &sumTpgs_Phi); 

    //putting bufsize at 32000 and changing split level to 0 so that the branch isn't split into multiple branches
    efficiencyTree->Branch("rlxTaus", "vector<TLorentzVector>", &rlxTaus, 32000, 0); 
    efficiencyTree->Branch("isoTaus", "vector<TLorentzVector>", &isoTaus, 32000, 0); 
    efficiencyTree->Branch("recoTaus", "vector<TLorentzVector>", &recoTaus, 32000, 0); 
    efficiencyTree->Branch("allRegions", "vector<TLorentzVector>", &allRegions, 32000, 0); 
    efficiencyTree->Branch("hcalTPGs", "vector<TLorentzVector>", &allHcalTPGs, 32000, 0); 
    efficiencyTree->Branch("ecalTPGs", "vector<TLorentzVector>", &allEcalTPGs, 32000, 0); 
    efficiencyTree->Branch("signalPFCands", "vector<TLorentzVector>", &signalPFCands, 32000, 0); 
    efficiencyTree->Branch("l1Jets", "vector<TLorentzVector>", &l1Jets, 32000, 0); 
    efficiencyTree->Branch("recoJets", "vector<TLorentzVector>", &recoJets, 32000, 0); 
    efficiencyTree->Branch("recoJetsDR", "vector<double>", &recoJetsDr, 32000, 0); 
    efficiencyTree->Branch("caloClusters", "vector<TLorentzVector>", &caloClusters, 32000, 0); 

    efficiencyTree->Branch("run",    &run,     "run/I");
    efficiencyTree->Branch("lumi",   &lumi,    "lumi/I");
    efficiencyTree->Branch("event",  &event,   "event/I");
    efficiencyTree->Branch("nvtx",   &nvtx,         "nvtx/I");

    efficiencyTree->Branch("decayMode", &decayMode,   "decayMode/I");
    
    efficiencyTree->Branch("tauEtaEcalEnt", &tauEtaEcalEnt,"tauEtaEcalEnt/D");
    efficiencyTree->Branch("tauPhiEcalEnt", &tauPhiEcalEnt,"tauPhiEcalEnt/D");

    efficiencyTree->Branch("recoPt",        &recoPt,   "recoPt/D");
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
    
  }

void L1TEventDisplayGenerator::beginJob( const EventSetup & es) {
}

//unsigned int const L1TEventDisplayGenerator::N_TOWER_PHI = 72;
//unsigned int const L1TEventDisplayGenerator::N_TOWER_ETA = 56;

void L1TEventDisplayGenerator::analyze( const Event& evt, const EventSetup& es )
 {

  run = evt.id().run();
  lumi = evt.id().luminosityBlock();
  event = evt.id().event();

  edm::Handle<reco::VertexCollection> vertices;
  if(evt.getByToken(vtxLabel_, vertices)){
    nvtx = (int) vertices->size();
    std::cout<<"nVertices "<<nvtx<<std::endl;
  }
  
  Handle<L1CaloRegionCollection> regions;

  edm::Handle<reco::PFTauDiscriminator> discriminatorIso;
  edm::Handle<reco::PFTauDiscriminator> discriminatorMu;

  edm::Handle < L1GctJetCandCollection > l1IsoTauJets;
  edm::Handle < L1GctJetCandCollection > l1TauJets;

  edm::Handle < vector<l1extra::L1JetParticle> > l1ExtraTaus;
  edm::Handle < vector<l1extra::L1JetParticle> > l1ExtraIsoTaus;
  edm::Handle < vector<reco::CaloCluster> > recoCaloClusters;

  edm::Handle < vector<l1extra::L1JetParticle> > l1ExtraJets;  

  std::vector<reco::PFTauRef> goodTausRef;
  std::vector<pat::Tau> goodTaus;
  
  edm::Handle<vector<pat::PackedCandidate> >pfCands;
  edm::Handle<EcalTrigPrimDigiCollection> ecalTPGs;
  edm::Handle<HcalTrigPrimDigiCollection> hcalTPGs;  

  /*Clear the vectors*/
  rlxTaus->clear(); 
  isoTaus->clear(); 
  recoTaus->clear(); 
  allRegions->clear(); 
  allEcalTPGs->clear(); 
  allHcalTPGs->clear(); 
  signalPFCands->clear(); 
  l1Jets->clear(); 
  recoJets->clear(); 
  caloClusters->clear(); 
  recoJetsDr->clear();
  


  if(!evt.getByToken(hcalSrc_, hcalTPGs))
    std::cout<<"ERROR GETTING THE HCAL TPGS"<<std::endl;

  ESHandle<L1CaloHcalScale> hcalScale;
  es.get<L1CaloHcalScaleRcd>().get(hcalScale);

  //get pfcandidate particles
  if(!evt.getByToken(packedPfCandsToken_, pfCands)){
    std::cout<<"Error Getting packed PFCandidates"<<std::endl;
  }

  //loop over jets
  Handle<vector<pat::Jet> > jets;
  if(evt.getByToken(jetSrc_, jets)){//Begin Getting Reco Jets
    ///////////
    for (const pat::Jet &jet : *jets) {
      //pat::Jet jet =  jets->at(iJet);
      TLorentzVector temp;
      temp.SetPtEtaPhiE(jet.pt(),jet.eta(),jet.phi(),jet.et());
      recoJets->push_back(temp);
      double jetDR=0;
      std::vector<reco::CandidatePtr> daus(jet.daughterPtrVector());
      std::sort(daus.begin(), daus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) 
		{ return p1->pt() > p2->pt(); }); // the joys of C++11

      for (unsigned int i2 = 0; i2 <  daus.size(); ++i2) {
	const pat::PackedCandidate &cand = dynamic_cast<const pat::PackedCandidate &>(*daus[i2]);
	float tempDR = reco::deltaR(jet.p4(), cand.p4());
	//std::cout<<"tempDR "<<tempDR<<std::endl;
	if(tempDR > jetDR){
	  jetDR = tempDR;
	}

	//printf("         constituent %3d: pt %6.2f, pdgId %+3d\n", i2,cand.pt(),cand.pdgId());
      }
      recoJetsDr->push_back(jetDR);
    }
  }  

  //std::cout<<"doing the taus"<<std::endl;
  // loop over taus
  Handle<vector<pat::Tau> > taus;
  if(evt.getByToken(slimmedTauSrc_, taus)){//Begin Getting Reco Taus
    for ( unsigned iTau = 0; iTau < taus->size(); ++iTau ) {
      pat::Tau tau =  taus->at(iTau);
      if(tau.tauID("decayModeFinding")>0.5&&tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits")<20&&tau.tauID("againstMuonLoose3")>0){
	TLorentzVector temp;
	temp.SetPtEtaPhiE(tau.pt(),tau.eta(),tau.phi(),tau.et());
	std::cout<<"decaymode "<<std::dec<<tau.decayMode()<<" pt: "<< tau.pt() << " eta: "<< tau.eta() << " phi: "<< tau.phi() <<std::endl;
	recoTaus->push_back(temp);
	
	//std::cout<<"# tau signal pf cands sigpf "<< tau.signalPFCands().size() <<std::endl;
	//std::cout<<"# tau signal pf cands sig "<< tau.signalCands().size() <<std::endl;
	//std::cout<<"# tau signal pf cands pi0 "<< tau.signalPiZeroCandidates().size() <<std::endl;
	//for(size_t i = 0; i < tau.signalCands().size(); i++){
	reco::CandidatePtrVector daus = tau.signalCands();
	std::sort(daus.begin(), daus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) 
		  { return p1->pt() > p2->pt(); }); // the joys of C++11
	
	for (unsigned int i2 = 0; i2 <  daus.size(); ++i2) {
	  const pat::PackedCandidate &cand = dynamic_cast<const pat::PackedCandidate &>(*daus[i2]);
	  TLorentzVector temp2;
	  temp2.SetPtEtaPhiE(cand.pt(),cand.eta(),cand.phi(),cand.et());
	  signalPFCands->push_back(temp2);
	}
      }
    }
  }//End Getting MiniAOD Taus
  else
    std::cout<<"Error getting reco taus"<<std::endl;
  //std::cout<<"end taus"<<std::endl;
  // loop over taus

  /*  
  Handle<vector<reco::PFTau> > tausReco;
  if(evt.getByToken(tauSrc_, tausReco)){//Begin Getting Reco Taus
    for ( unsigned iTau = 0; iTau < tausReco->size(); ++iTau ) {
      reco::PFTauRef tauCandidate(tausReco, iTau);
      reco::PFTau goodTau = tausReco->at(iTau);
      if(evt.getByToken(discriminatorMu_, discriminatorMu)){
	if(evt.getByToken(discriminatorIso_, discriminatorIso)){
	  //if( (*discriminatorIso)[tauCandidate] > 0.5 && tauCandidate->decayMode()>-1 && (*discriminatorMu)[tauCandidate] > 0.5 ){
	    TLorentzVector temp;
	    temp.SetPtEtaPhiE(tauCandidate->pt(),tauCandidate->eta(),tauCandidate->phi(),tauCandidate->et());
	    std::cout<<"decaymode "<<tauCandidate->decayMode()<<std::endl;
	    recoTaus->push_back(temp);
	    //std::cout<<"# tau signal pf cands sigpf "<< tau.signalPFCands().size() <<std::endl;
	    //std::cout<<"# tau signal pf cands sig "<< tau.signalCands().size() <<std::endl;
	    std::cout<<"# tau signal pf cands pi0 "<< goodTau.signalPFGammaCands().size() <<std::endl;
	    for(size_t i = 0; i < goodTau.signalPFGammaCands().size(); i++){
	      const reco::PFCandidate cand = goodTau.signalPFGammaCands().at(i);
	      std::cout<<"cand Pt "<<cand.pt()<< " eta "<<cand.eta()<<" phi "<<cand.phi()<<std::endl;
	      std::cout<<"ecal entrance  eta "<<cand.positionAtECALEntrance().eta()<<" phi "<<cand.positionAtECALEntrance().phi()<<std::endl;
	      TLorentzVector temp2;
	      temp2.SetPtEtaPhiE(cand.pt(),cand.eta(),cand.phi(),cand.et());
	      std::cout<<"signal pf cand "<<i<<" ... et: "<<cand.pt()<<std::endl;
	      signalPFCands->push_back(temp2);
	    }
	    //}else std::cout<<"tau failed cuts"<<std::endl;
	}else std::cout<<"isolation discriminator not found "<<std::endl;
      }else std::cout<<"muon discriminator not fond "<<std::endl;
    }
    //}
  }//End Getting RECO Taus
  else
    std::cout<<"Error getting reco taus"<<std::endl;
  */

  if(evt.getByToken(ecalCaloSrc_, recoCaloClusters)){
    for( vector<reco::CaloCluster>::const_iterator caloCluster = recoCaloClusters->begin(); 
	 caloCluster != recoCaloClusters->end(); 
	 caloCluster++ ) {
      //fill vector
      TLorentzVector temp ;
      temp.SetPtEtaPhiE(caloCluster->energy(),caloCluster->eta(),caloCluster->phi(),caloCluster->energy());
      caloClusters->push_back(temp);
    }
  }


  if(evt.getByToken(l1ExtraTauSource_, l1ExtraTaus)){
    for( vector<l1extra::L1JetParticle>::const_iterator rlxTau = l1ExtraTaus->begin(); 
	 rlxTau != l1ExtraTaus->end(); 
	 rlxTau++ ) {
      //fill vector
      TLorentzVector temp ;
      temp.SetPtEtaPhiE(rlxTau->pt(),rlxTau->eta(),rlxTau->phi(),rlxTau->et());
      // std::cout<< "rlxTau Pt: " <<rlxTau->pt()<<" Eta: " <<rlxTau->eta()<< " Phi: "<<rlxTau->phi()<<std::endl;
      rlxTaus->push_back(temp);
    }
  }

  if(evt.getByToken(l1ExtraIsoTauSource_ , l1ExtraIsoTaus)){
    for( vector<l1extra::L1JetParticle>::const_iterator isoTau = l1ExtraIsoTaus->begin();  
	 isoTau != l1ExtraIsoTaus->end();  
	 ++isoTau ) {
      //fill vector
      TLorentzVector temp ;
      temp.SetPtEtaPhiE(isoTau->pt(),isoTau->eta(),isoTau->phi(),isoTau->et());
      isoTaus->push_back(temp);
    }
  }

  if(evt.getByToken(l1ExtraJetSource_, l1ExtraJets)){
    for( vector<l1extra::L1JetParticle>::const_iterator l1Jet = l1ExtraJets->begin(); 
	 l1Jet != l1ExtraJets->end(); 
	 l1Jet++ ) {
      //fill vector
      TLorentzVector temp ;
      temp.SetPtEtaPhiE(l1Jet->pt(),l1Jet->eta(),l1Jet->phi(),l1Jet->et());
      l1Jets->push_back(temp);
    }
  }


  UCTGeometry g;
  if(!evt.getByToken(regionSource_,regions)){
    std::cout<<"ERROR GETTING THE REGIONS!!!"<<std::endl;}
  else{
    for(vector<L1CaloRegion>::const_iterator region = regions->begin(); region != regions->end(); ++region){
      UCTRegionProcess uctRegion(*region);
      if(region->et()>0){
	float pt = (region->et());
	
	float eta = uctRegion.getFineRecoEta();
	
	float phi = uctRegion.getFineRecoPhi();
	TLorentzVector temp ;
	temp.SetPtEtaPhiE(pt,eta,phi,pt);
	allRegions->push_back(temp);
      }
    }
  }

  if(!evt.getByToken(ecalSrc_, ecalTPGs))
    std::cout<<"ERROR GETTING THE ECAL TPGS"<<std::endl;
  else
    for (size_t i = 0; i < ecalTPGs->size(); ++i) {

      int cal_ieta = (*ecalTPGs)[i].id().ieta();
      int cal_iphi = (*ecalTPGs)[i].id().iphi();
      if(cal_iphi==0)
	std::cout<<"cal_phi is 0"<<std::endl;
      if(cal_ieta<-28)
	continue;
      if(cal_ieta>28)
	continue;
      int ieta = TPGEtaRange(cal_ieta);
      short zside = (*ecalTPGs)[i].id().zside();
      // TPG iPhi starts at 1 and goes to 72.  Let's index starting at zero.
      // TPG ieta ideal goes from 0-55.
      double LSB = 0.5;
      double et= (*ecalTPGs)[i].compressedEt()*LSB;
      if(ieta<0){
	std::cout<<"sorry, ieta less than 1 :("<<std::endl;
	std::cout<<"cal_ieta "<<cal_ieta<<" ieta "<<ieta<<std::endl;
      }
      float eta = getRecoEta(ieta, zside);
      float phi = getRecoPhi(cal_iphi);    
      //if(et>0)
      //std::cout<<"et "<<et<<std::endl;
      TLorentzVector temp ;
      temp.SetPtEtaPhiE(et,eta,phi,et);
      //if(et>5)
      //std::cout<<"Event Display tpg ecal pt() "<<temp.Pt()<< " eta " <<eta << " phi "<< phi <<std::endl;
      allEcalTPGs->push_back(temp);
    }

  if(!evt.getByToken(hcalSrc_, hcalTPGs))
    std::cout<<"ERROR GETTING THE HCAL TPGS"<<std::endl;
  else
    for (size_t i = 0; i < hcalTPGs->size(); ++i) {
      HcalTriggerPrimitiveDigi tpg = (*hcalTPGs)[i];
      int cal_ieta = tpg.id().ieta();
      int cal_iphi = tpg.id().iphi();
      if(cal_ieta>28)continue; 
      if(cal_ieta<-28)continue; 
      int ieta = TPGEtaRange(cal_ieta);
      short absieta = std::abs(tpg.id().ieta());
      short zside = tpg.id().zside();
      double et = hcalScale->et(tpg.SOI_compressedEt(), absieta, zside); 
      //if(et>0)
      //std::cout<<"HCAL ET "<<et<<std::endl;
      if(ieta<0){
	std::cout<<"sorry, ieta less than 1 :("<<std::endl;
	std::cout<<"cal_ieta "<<cal_ieta<<" ieta "<<ieta<<std::endl;
      }
      float eta = getRecoEta(ieta, zside);
      float phi = getRecoPhi(cal_iphi);    
      TLorentzVector temp ;
      temp.SetPtEtaPhiE(et,eta,phi,et);
      allHcalTPGs->push_back(temp);
    }

  efficiencyTree->Fill();
 }


//four vector addition vs linear sum
double
L1TEventDisplayGenerator::getPFCandsEt(const std::vector<pat::PackedCandidate> pfCands){
  double etTotal = 0;
  for (uint32_t i = 0; i < pfCands.size(); i++ ) {
    etTotal +=pfCands.at(i).et();
  }
  return etTotal;
}

double
L1TEventDisplayGenerator::getPFCandsEtEtaPhi(edm::Handle<std::vector<pat::PackedCandidate> >& pfCands, const pat::Tau & tau, double dR){
  double etTotal = 0;
  for (uint32_t i = 0; i < pfCands->size(); i++ ) {
    if(reco::deltaR(tau, pfCands->at(i)) < dR){
      etTotal +=pfCands->at(i).et();
    }
  }
  return etTotal;
}

int
L1TEventDisplayGenerator::get5x5TPGs(const int maxTPGPt_eta, 
				     const int maxTPGPt_phi, 
				     const double eTowerETMap[73][57], 
				     const double hTowerETMap[73][57], 
				     std::vector<double>* hcalTpgs_pt, 
				     std::vector<double>* hcalTpgs_eta, 
				     std::vector<double>* hcalTpgs_phi, 
				     std::vector<double>* ecalTpgs_pt, 
				     std::vector<double>* ecalTpgs_eta, 
				     std::vector<double>* ecalTpgs_phi,
				     std::vector<double>* sumTpgs_pt, 
				     std::vector<double>* sumTpgs_eta, 
				     std::vector<double>* sumTpgs_phi){
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
      hcalTpgs_pt->push_back(hTowerETMap[tpgsquarephi][tpgsquareeta]);
      hcalTpgs_eta->push_back(tpgEtaMap[k]);
      hcalTpgs_phi->push_back(tpgPhiMap[j]);

      ecalTpgs_pt->push_back(eTowerETMap[tpgsquarephi][tpgsquareeta]);
      ecalTpgs_eta->push_back(tpgEtaMap[tpgsquareeta]);
      ecalTpgs_phi->push_back(tpgPhiMap[tpgsquarephi]);

      sumTpgs_pt->push_back(eTowerETMap[tpgsquarephi][tpgsquareeta]+eTowerETMap[tpgsquarephi][tpgsquareeta]);
      sumTpgs_eta->push_back(tpgEtaMap[tpgsquareeta]);
      sumTpgs_phi->push_back(tpgPhiMap[tpgsquarephi]);
    }
  }
  //std::cout<<"TPGe5x5_ "<<TPGe5x5_<<" TPGh5x5_ "<<TPGh5x5_<<std::endl;
  //return (TPGe5x5_ + TPGh5x5_);
  return 1;
}

/*
 * Get the ECAL TPGS create a TPG map for the event
 *
 */
  
void L1TEventDisplayGenerator::initializeECALTPGMap(Handle<EcalTrigPrimDigiCollection> ecal, double eTowerETMap[73][57], bool testMode){
  
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

    if(testMode && iphi == 34 && ieta == 11){
      et = 40;
    }

    if (iphi >= 0 && iphi <= 72 &&
	ieta >= 0 && ieta <= 55) {
      eTowerETMap[iphi][ieta] = et; 
    }

  }

}

 void L1TEventDisplayGenerator::initializeHCALTPGMap(const Handle<HcalTrigPrimDigiCollection> hcal, 
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

void L1TEventDisplayGenerator::endJob() {
}

L1TEventDisplayGenerator::~L1TEventDisplayGenerator(){
}

DEFINE_FWK_MODULE(L1TEventDisplayGenerator);
