#ifndef printTauEvents_H
#define printTauEvents_H

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
//#include "L1Trigger/L1TCaloLayer1/src/L1UCTCollections.h"

#include <memory>
#include <math.h>
#include <vector>
#include <list>

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "CondFormats/L1TObjects/interface/L1CaloHcalScale.h"
#include "CondFormats/DataRecord/interface/L1CaloHcalScaleRcd.h"
#include "FWCore/Framework/interface/ESHandle.h"

//
// class declaration
//
using std::vector;

class printTauEvents : public edm::EDAnalyzer {

 public:
  
  // Constructor
  printTauEvents(const edm::ParameterSet& ps);
  
  // Destructor
  virtual ~printTauEvents();

  edm::Service<TFileService> tfs_;

  std::ofstream file0, file1, file10;

  TTree* efficiencyTree;

  int run, lumi, event, nvtx;
  double recoPt, recoEta, recoPhi;
  int decayMode;
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
  edm::EDGetTokenT<reco::VertexCollection> vertices_;

  edm::EDGetTokenT<reco::VertexCollection> vtxLabel_;
  edm::EDGetTokenT<reco::PFTauDiscriminator> discriminatorMu_;
  edm::EDGetTokenT<reco::PFTauDiscriminator> discriminatorIso_;
  edm::EDGetTokenT<vector<pat::Tau> > tauSrc_;

  std::string folderName_;
  double recoPt_;
  double maxRecoPt_;

};

#endif
