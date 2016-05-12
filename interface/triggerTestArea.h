#ifndef triggerTestArea_H
#define triggerTestArea_H

// system include files
#include <memory>
#include <unistd.h>
#include <iostream>
#include <fstream>

#include <iostream>
#include <fstream>
#include <vector>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

// GCT and RCT data formats
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctCollections.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/Tau.h"


#include <memory>
#include <math.h>
#include <vector>
#include <list>

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "CondFormats/L1TObjects/interface/L1CaloHcalScale.h"
#include "CondFormats/DataRecord/interface/L1CaloHcalScaleRcd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "L1Trigger/Stage3Ntuplizer/plugins/UCTRegionProcess.hh"

//
// class declaration
//
using std::vector;

struct tower{
  uint32_t ecalEt;
  uint32_t hcalEt;
  float EoH;
  float HoE;
  float HoEpH;
  float EoEpH;
  float EpH;
};

class triggerTestArea : public edm::EDAnalyzer {

 public:
  
  // Constructor
  triggerTestArea(const edm::ParameterSet& ps);
  
  // Destructor
  virtual ~triggerTestArea();

  edm::Service<TFileService> tfs_;

  std::ofstream file0, file1, file10;

  TH1F* nEvents;

  TH1F* isoTau_pt;
  TH1F* isoTau_eta;
  TH1F* isoTau_phi;

  TH1F* tau_pt;
  TH1F* tau_eta;
  TH1F* tau_phi;

  TH1F* tau_pt_diTau;
  TH1F* tau_pt_diTau_eta2p4;
  TH1F* tau_pt_diTau_eta2p1;

  TH1F* isoTau_pt_diTau;
  TH1F* isoTau_pt_diTau_eta2p4;
  TH1F* isoTau_pt_diTau_eta2p1;

  TH1F* recoTau_pt;
  TH1F* recoTau_eta;
  TH1F* recoTau_phi;

  TH1F* regionEta;
  TH1F* regionPhi;
  TH1F* regionPt;
  TH1F* regionEtaFine;
  TH1F* regionPhiFine;
  TH1F* regionTotal;

  TH1F* regionHitEta;
  TH1F* regionHitPhi;
  TTree* efficiencyTreeTaus;
  TTree* efficiencyTreeJets;
  TFileDirectory folder;
  TFileDirectory folderJets;

  int run, lumi, event, nvtx;
  int runJ, lumiJ, eventJ, nvtxJ;
  double isoTauPt, rlxTauPt, isoTauEta, rlxTauEta, isoTauPhi, rlxTauPhi;
  double recoPt, recoEta, recoPhi;
  double recoPtJ, recoEtaJ, recoPhiJ;
  int decayMode;
  int decayModeJ;
  //double tauEtaEcalEnt,tauPhiEcalEnt,rawEcal, rawHcal, ecal, hcal, jetEt, jetEta, jetPhi, nvtx;
  //double max3ProngDeltaR, minProngPt, maxProngPt, midProngPt; int n3ProngCands;
  //double pfCandsEt, signalCandsEt, isoCandsEt;
  double TPG3x3, TPGH3x3, TPGE3x3;
  double TPG3x3J, TPGH3x3J, TPGE3x3J;
  double EoHTPG3x3,  HoETPG3x3,  HoEpHTPG3x3,  EoEpHTPG3x3;
  double EoHTPG3x3J,  HoETPG3x3J,  HoEpHTPG3x3J,  EoEpHTPG3x3J;
  double TPG_highPtTower, TPGH_highPtTower, TPGE_highPtTower;
  double TPG_highPtTowerJ, TPGH_highPtTowerJ, TPGE_highPtTowerJ;
  double EoH_highPtTower, HoE_highPtTower, EoEpH_highPtTower, HoEpH_highPtTower;
  double EoH_highPtTowerJ, HoE_highPtTowerJ, EoEpH_highPtTowerJ, HoEpH_highPtTowerJ;

  double TPG_nearTower, TPGH_nearTower, TPGE_nearTower;
  double TPG_nearTowerJ, TPGH_nearTowerJ, TPGE_nearTowerJ;
  double EoH_nearTower, HoE_nearTower, EoEpH_nearTower, HoEpH_nearTower;
  double EoH_nearTowerJ, HoE_nearTowerJ, EoEpH_nearTowerJ, HoEpH_nearTowerJ;

  double EoH_eAndhTower, HoE_eAndhTower, EoEpH_eAndhTower, HoEpH_eAndhTower;
  double EoH_eAndhTowerJ, HoE_eAndhTowerJ, EoEpH_eAndhTowerJ, HoEpH_eAndhTowerJ;

  //double EoH_nearTower, HoE_nearTower, EoEpH_nearTower, HoEpH_nearTower;

  //double TPG5x5, TPGH5x5, TPGE5x5;
  //double TPG6x6, TPGH6x6, TPGE6x6;
  //double TPG7x7, TPGH7x7, TPGE7x7;

  void getThreeProngInfo(const pat::Tau & tau, double &maxDeltaR, double &minProngPt, double &midProngPt, double &maxProngPt, int &nCands);
  void getRawEcalHcalEnergy(const pat::PackedCandidate pfCand, double &rawEcal, double &rawHcal, double &ecal, double &hcal);
  double getPFCandsEt(const std::vector<pat::PackedCandidate> pfCands);
  double getPFCandsEtEtaPhi(edm::Handle<std::vector<pat::PackedCandidate> >& pfCands, const pat::Tau &tau, double dR);
  void initializeHCALTPGMap(const edm::Handle<HcalTrigPrimDigiCollection> hcal, const  edm::ESHandle<L1CaloHcalScale> hcalScale, double hTowerETMap[73][57], bool testMode = false);
  void initializeECALTPGMap(edm::Handle<EcalTrigPrimDigiCollection> ecal, double eTowerETMap[73][57], bool testMode = false);
  int get1x1TPGs(const int maxTPGPt_eta,const int maxTPGPt_phi, const double eTowerETMap[73][57], const double hTowerETMap[73][57], double &TPGe5x5_, double &TPGh5x5_ , std::vector<tower> &towers);
  int get3x3TPGs(const int maxTPGPt_eta,const int maxTPGPt_phi, const double eTowerETMap[73][57], const double hTowerETMap[73][57], double &TPGe2x2_, double &TPGh2x2_ , std::vector<tower> &towers);
  int get5x5TPGs(const int maxTPGPt_eta,const int maxTPGPt_phi, const double eTowerETMap[73][57], const double hTowerETMap[73][57], double &TPGe5x5_, double &TPGh5x5_ );
  int get6x6TPGs(const int maxTPGPt_eta,const int maxTPGPt_phi, const double eTowerETMap[73][57], const double hTowerETMap[73][57], double &TPGe6x6_, double &TPGh6x6_ );
  int get7x7TPGs(const int maxTPGPt_eta,const int maxTPGPt_phi, const double eTowerETMap[73][57], const double hTowerETMap[73][57], double &TPGe7x7_, double &TPGh7x7_ );

 protected:
  // Analyze
  void analyze(const edm::Event& evt, const edm::EventSetup& es);

  // BeginJob
  void beginJob(const edm::EventSetup &es);
  
  // EndJob
  void endJob(void);

  
 private:
  // ----------member data ---------------------------

  int nev_; // Number of events processed
  bool verbose_;
  std::ofstream logFile_;
  edm::InputTag rctSource_; 

  edm::EDGetTokenT<vector<pat::PackedCandidate> > pfCandsToken_;  
  edm::EDGetTokenT<L1CaloRegionCollection> L1RegionCollection;
  edm::EDGetTokenT<L1CaloEmCollection> L1EMCollection_;
  edm::EDGetTokenT<reco::VertexCollection> vertices_;
  edm::EDGetTokenT<EcalTrigPrimDigiCollection> ecalSrc_; 
  edm::EDGetTokenT<HcalTrigPrimDigiCollection> hcalSrc_;
  //edm::EDGetTokenT<double> recoPt_;
  //edm::EDGetTokenT<std::string> folderName_;
  edm::EDGetTokenT<reco::VertexCollection> vtxLabel_;
  edm::EDGetTokenT<reco::PFTauDiscriminator> discriminatorMu_;
  edm::EDGetTokenT<reco::PFTauDiscriminator> discriminatorIso_;
  edm::EDGetTokenT<vector<pat::Tau> > tauSrc_;
  edm::EDGetTokenT<vector<pat::Jet> > jetSrc_;
  //edm::EDGetTokenT<L1GctJetCandCollection> gctIsoTauJetsSource_;
  //edm::EDGetTokenT<L1GctJetCandCollection> gctTauJetsSource_;
  edm::EDGetTokenT<vector <l1extra::L1JetParticle> > l1ExtraIsoTauSource_;
  edm::EDGetTokenT<vector <l1extra::L1JetParticle> > l1ExtraTauSource_;
  edm::EDGetTokenT<BXVector <l1t::Tau> > l1Stage2TauSource_;
  edm::EDGetTokenT<vector <L1CaloRegion> > regionSource_;


  std::string folderName_;
  double recoPt_;



 int TPGEtaRange(int ieta){
   int iEta = 0;
   // So here, -28 becomes 0.  -1 be comes 27.  +1 becomes 28. +28 becomes 55.
   // And we have mapped [-28, -1], [1, 28] onto [0, 55]   
   if(ieta < 0)
     iEta = ieta + 28;
   else if(ieta > 0)
     iEta = ieta + 27;
   return iEta;
 }

  int convertGenEta(double inputEta) {
    const double tpgEtaValues[28] = {
      0.087,      
      0.174, // HB and inner HE bins are 0.348 wide
      0.261,
      0.348,
      0.435,
      0.522,
      0.609,
      0.696,
      0.783,
      0.870,
      0.957,
      1.044,
      1.131,
      1.218,
      1.305,
      1.392,
      1.479,
      1.566,
      1.653,
      1.74,
      1.827,
      1.914, // Last two HE bins are 0.432 and 0.828 wide
      2.001,
      2.217,
      2.391,
      2.565,
      2.739,
      //2.88,
      3
      //IGNORING HF
      //3.250, // HF bins are 0.5 wide
      //3.750,
      //4.250,
      //4.750
    };

    for (int n=1; n<=28; n++){
      //std::cout<<"inputEta "<<inputEta<< " n "<< n <<" tpgEtaValues[n-1] "<< tpgEtaValues[n-1] << " abs(inputEta)<tpgEtaValues[n-1]"<<std::endl;
      if (std::fabs(inputEta)<tpgEtaValues[n-1]) {
	//std::cout<<"found to be true"<<std::endl;
	//int tpgEta = n;
	//Positive eta is >28
	//negative eta is 0 to 27
	if(inputEta>0){
	  //std::cout<<"returning input eta >0 so + 28"<<std::endl;
	  return n + 26;}
	else{
	  //std::cout<<"returning input eta <0 so n"<<std::endl;
	  return (28 - n);}
	break;
      }
    }
    std::cout<<"OUT OF BOUNDS!!!!  inputeta: "<<inputEta<<std::endl;
    return -9;
  }

  //-pi < phi <= +pi,
  int convertGenPhi(double inputPhi){
    double posPhi[37];
    for(int n = 0; n < 37; n++)
      posPhi[n] = (0.087266) * n ;
    double negPhi[37];
    for(int n = 0; n < 37; n++){
      negPhi[n] = -3.14159 + 0.087266 * n + 0.087266;
      if(n == 36)
	negPhi[n] = 0;
    }

    //1 to 36 is 0 to pi
    if( 3.1416 > inputPhi && inputPhi >= 0){

      for(int n = 1; n < 37; n++){
	
	if(inputPhi <= posPhi[n-1]){
	  int tpgPhi = n - 1;
	  return tpgPhi;
	}
      }
    }

    //37 to 72 is -pi to 0
    else if(-3.1416 < inputPhi && inputPhi < 0){
      for(int n = 1; n < 37; n++){
	//std::cout<<"inputPhi "<<inputPhi<< " negPhi["<< n -1 <<"] "<< negPhi[n-1] << " n "<<n<<std::endl;
	if(inputPhi < negPhi[n-1]){
	  int tpgPhi = n + 36;
	  return tpgPhi;
	}}
    }
    std::cout<<"OUT OF BOUNDS!!!!  inputphi: "<<inputPhi<<std::endl;
    return -9;
  }


};

#endif
