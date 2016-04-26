#include <vector>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TSystem.h"
#include <TLorentzVector.h>
#include <TStyle.h>
#include "TLegend.h"
#include "TEllipse.h"
#include "TPaveText.h"
#include "TLine.h"
#include <sstream>
#include "Math/VectorUtil_Cint.h"

#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif

#ifdef __MAKECINT__
#pragma link C++ class vector<TLorentzVector>+;
#endif


void DrawRegionLines(){
  std::vector<TLine*> regionLines;
  float etaValues[17] = { -3, -2.088, -1.74, -1.392, -1.044, -0.696, -0.348, 0,
		0.348, 0.696, 1.044, 1.392, 1.74, 2.088, 3 };
  float phiValues[18] =
  {-2.965, -2.617, -2.268, -1.919, -1.570, -1.221, -0.872, -0.523, -0.174, 
      0.174, 0.523, 0.872, 1.221, 1.570, 1.919, 2.268, 2.617, 2.965};

  //eta lines
  for(int i = 0; i < 17; i++){
    TLine * line = new TLine(etaValues[i], -3.2, etaValues[i], 3.2); 
    line->SetLineColor(kBlue-7);
    line->SetLineStyle(1);
    regionLines.push_back(line);
  }

  //phi lines
  for(int i = 0; i < 18; i++){
    TLine * line = new TLine(-3, phiValues[i], 3, phiValues[i]); 
    line->SetLineColor(kBlue-7);
    line->SetLineStyle(1);
    regionLines.push_back(line);
  }

  for(size_t j = 0; j < regionLines.size(); j++){
    regionLines.at(j)->Draw();
  }
}

void DrawTowerLines(){
  std::vector<TLine*> TowerLines;
  float etaValues[59] = { -2.913, -2.739, -2.565, -2.391, -2.217, -2.088, -2.001, -1.914, -1.827, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.87, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.87, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.74, 1.827, 1.914, 2.001, 2.088, 2.217, 2.391, 2.565, 2.739, 2.913};

  float phiValues[73] =
    {-3.142, -3.054, -2.967, -2.880, -2.793, -2.705, -2.618, -2.531, -2.443, -2.356, -2.269, -2.182, -2.094, -2.007, -1.920, -1.833, -1.745, -1.658, -1.571, -1.484, -1.396, -1.309, -1.222, -1.134, -1.047, -0.960, -0.873, -0.785, -0.698, -0.611, -0.524, -0.436, -0.349, -0.262, -0.175, -0.087,
     0.000, 0.087, 0.175, 0.262, 0.349, 0.436, 0.524, 0.611, 0.698, 0.785, 0.873, 0.960, 1.047, 1.134, 1.222, 1.309, 1.396, 1.484, 1.571, 1.658, 1.745, 1.833, 1.920, 2.007, 2.094, 2.182, 2.269, 2.356, 2.443, 2.531, 2.618, 2.705, 2.793, 2.880, 2.967, 3.054, 3.142};
  
  //eta lines
  for(int i = 0; i < 59; i++){
    TLine * line = new TLine(etaValues[i], -3.2, etaValues[i], 3.2); 
    line->SetLineColor(kGray);
    line->SetLineStyle(3);
    TowerLines.push_back(line);
  }

  //phi lines
  for(int i = 0; i < 73; i++){
    TLine * line = new TLine(-3, phiValues[i], 3, phiValues[i]); 
    line->SetLineColor(kGray);
    line->SetLineStyle(3);
    TowerLines.push_back(line);
  }

  for(size_t j = 0; j < TowerLines.size(); j++){
    TowerLines.at(j)->Draw();
  }
}

void plotVector(int iEvent){
  gStyle->SetOptStat(0);
  //TFile *f = TFile::Open("l1TNtuple.root","READ");
  //TFile *f = TFile::Open("eventDisplay-2987803423.root","READ");
  //TFile *f = TFile::Open("eventDisplay-DecayMode1-Lumi1643.root","READ");
  //TFile *f = TFile::Open("eventDisplay-2988304974.root","READ");
  //TFile *f = TFile::Open("eventDisplay-ggH-546248.root","READ");
  //TFile *f = TFile::Open("highptTaus2.root","READ");
  TFile *f = TFile::Open("eventDisplay-HighL1LowReco-simTPGs.root","READ");
  //TFile *f = TFile::Open("eventDisplay-DecayMode0-Lumi1643-event2988093804.root","READ");
  
  if (!f) { return; }

  TTree *t = (TTree*) f->Get("l1NtupleProducer/EfficiencyTree");
  //f->GetObject("l1NtupleProducer/EfficiencyTree",t);

  std::vector<TLorentzVector> *vpx             = 0;
  std::vector<TLorentzVector> *vRecoTaus       = 0;
  std::vector<TLorentzVector> *vSignalTauCands = 0;
  std::vector<TLorentzVector> *vRecoJets       = 0;
  std::vector<double>         *vRecoJetsDr     = 0;
  std::vector<TLorentzVector> *vRlxTaus        = 0;
  std::vector<TLorentzVector> *vEcalTpgs       = 0;
  std::vector<TLorentzVector> *vHcalTpgs       = 0;
  std::vector<TLorentzVector> *vL1Jets         = 0;
  std::vector<TLorentzVector> *vCaloClusters   = 0;
  int event =0;

  // Create a new canvas.
  TCanvas *c1 = new TCanvas("c1","eta vs phi",200,10,700,700);
  c1->SetFillColor(0);
  c1->GetFrame()->SetFillColor(0);
  c1->GetFrame()->SetBorderSize(6);
  c1->GetFrame()->SetBorderMode(-1);

  const Int_t kUPDATE = 1000;

  TBranch *bvpx = 0;
  TBranch *bRecoTaus = 0;
  TBranch *bSignalTauCands = 0;
  TBranch *bRecoJets = 0;
  TBranch *bRecoJetsDr = 0;
  TBranch *bRlxTaus = 0;
  TBranch *bEcalTpgs = 0;
  TBranch *bHcalTpgs = 0;
  TBranch *bL1Jets = 0;
  TBranch *bCaloClusters = 0;
  TBranch *bEvent = 0;      

  t->SetBranchAddress("event",&event,&bEvent);
  t->SetBranchAddress("allRegions",&vpx,&bvpx);
  t->SetBranchAddress("recoTaus",&vRecoTaus,&bRecoTaus);
  t->SetBranchAddress("signalPFCands",&vSignalTauCands,&bSignalTauCands);
  t->SetBranchAddress("recoJets",&vRecoJets,&bRecoJets);
  t->SetBranchAddress("recoJetsDR",&vRecoJetsDr,&bRecoJetsDr);
  t->SetBranchAddress("rlxTaus",&vRlxTaus,&bRlxTaus);
  t->SetBranchAddress("ecalTPGs",&vEcalTpgs,&bEcalTpgs);
  t->SetBranchAddress("hcalTPGs",&vHcalTpgs,&bHcalTpgs);
  t->SetBranchAddress("l1Jets",&vL1Jets,&bL1Jets);
  t->SetBranchAddress("caloClusters",&vCaloClusters,&bCaloClusters);

  // Create one histograms
  TH1F   *h                = new TH1F("h","This is the eta distribution",100,-4,4);
  TH2F   *h2               = new TH2F("h2","Event 2988846758",68,-3,3,72,-3.142,3.142);
  TGraph *h2RecoTaus       = new TGraph();
  TGraph *h2SignalTauCands = new TGraph();
  TH2F   *h2RlxTaus        = new TH2F("h2RlxTaus","h2 title"  ,68,-3,3,72,-3.142,3.142);
  TH2F   *h2EcalTpgs       = new TH2F("h2EcalTpgs","h2 title",68,-3,3,72,-3.142,3.142);
  TH2F   *h2HcalTpgs       = new TH2F("h2HcalTpgs","h2 title",68,-3,3,72,-3.142,3.142);
  TH2F   *h2CaloClusters   = new TH2F("h2CaloClusters","h2 title",68,-3,3,72,-3.142,3.142);
  std::vector<TEllipse*>  jets; 
  std::vector<TPaveText*> l1JetText;

  std::vector<TEllipse*>  recoJets; 
  std::vector<TPaveText*> recoJetText;
  std::vector<TPaveText*> ecalTpgText;
  std::vector<TPaveText*> hcalTpgText;
  std::vector<TPaveText*> rlxTauText; 

  h->SetFillColor(48);

  //for (Int_t i = 0; i < 25000; i++) {
  int i = iEvent;
  Long64_t tentry = t->LoadTree(i);
  std::cout<<"i "<<i<< " tentry "<< tentry << std::endl;
  bEvent->GetEntry(tentry);
  bvpx->GetEntry(tentry);
  bRecoTaus->GetEntry(tentry);
  bSignalTauCands->GetEntry(tentry);
  bRecoJets->GetEntry(tentry);
  bRecoJetsDr->GetEntry(tentry);
  bRlxTaus->GetEntry(tentry);
  bEcalTpgs->GetEntry(tentry);
  bHcalTpgs->GetEntry(tentry);
  bL1Jets->GetEntry(tentry);
  bCaloClusters->GetEntry(tentry);

  //get the event number
  char* name = new char[30];
  sprintf(name,"Event %u",event);
  std::cout<<event<<std::endl;
  std::cout<<name<<std::endl;
  h2 = new TH2F("h2",name,68,-3,3,72,-3.142,3.142);
  delete name;
  int k = 0;
  for (UInt_t j = 0; j < vpx->size(); ++j) {
    if(vpx->at(j).Pt()>0)
      h->Fill(vpx->at(j).Eta());

    if(vpx->at(j).Pt()>1){
      h2->Fill(vpx->at(j).Eta(),vpx->at(j).Phi(),vpx->at(j).Pt());
      k++; 
      if(vpx->at(j).Pt()>10)
	std::cout<<"region eta "<<vpx->at(j).Eta()<<" region phi "<<vpx->at(j).Phi()<< " pt "<<vpx->at(j).Pt()/2 <<" GeV"<<std::endl;
    }
  }
  std::cout<<"n Regions above pt 0 :"<<k<<std::endl;

  for (UInt_t j = 0; j < vRecoTaus->size(); ++j) {
    //h2RecoTaus->Fill(vRecoTaus->at(j).Eta(),vRecoTaus->at(j).Phi(),vRecoTaus->at(j).Et());
    h2RecoTaus->SetPoint(j,vRecoTaus->at(j).Eta(),vRecoTaus->at(j).Phi());
    std::cout<<"recoTau "<< i<< " Pt: "<<vRecoTaus->at(j).Pt()<<std::endl;
  }    
  
  for (UInt_t j = 0; j < vSignalTauCands->size(); ++j) {
    //h2RecoTaus->Fill(vRecoTaus->at(j).Eta(),vRecoTaus->at(j).Phi(),vRecoTaus->at(j).Et());
    //if(vSignalTauCands->at(j).Pt()>20)
    h2SignalTauCands->SetPoint(j,vSignalTauCands->at(j).Eta(),vSignalTauCands->at(j).Phi());
    std::cout<<"vSignalTauCands->at(j).Pt() "<<vSignalTauCands->at(j).Pt()<<std::endl;
  }    
  
  for (UInt_t j = 0; j < vRlxTaus->size(); ++j) {
    double eta = vRlxTaus->at(j).Eta();
    double phi = vRlxTaus->at(j).Phi();
    double pt = vRlxTaus->at(j).Pt();

    h2RlxTaus->Fill(eta,phi,pt);
    std::cout<<"vRlxTaus->at(j).Pt() "<<pt<< " eta: "<<eta <<" phi "<< phi <<std::endl;

    if(pt>10){
      std::ostringstream strs;
      strs << pt;
      std::string text = strs.str();
      eta += 0.05;
      phi += 0.05;
      TPaveText *tempText = new TPaveText( eta, phi, eta+0.05, phi+0.05 );
      tempText->AddText(text.c_str());
      tempText->SetFillColor(0);
      tempText->SetLineColor(0);
      tempText->SetShadowColor(0);
      tempText->SetTextColor(kGreen);
      rlxTauText.push_back(tempText);
    }
  }

  for (UInt_t j = 0; j < vpx->size(); ++j) {
    if(ROOT::Math::VectorUtil::DeltaR(vRecoTaus->at(0), vpx->at(j))<0.2){
      std::cout<<"Matched Reco pt: "<< vRecoTaus->at(0).Pt() <<" region pt: "<< vpx->at(j).Pt()/2 <<" GeV... resolution: "<<(vRecoTaus->at(0).Pt()-(vpx->at(j).Pt()/2))/vRecoTaus->at(0).Pt()<<std::endl;
    }
  }

  

  for (UInt_t j = 0; j < vEcalTpgs->size(); ++j) {
    double eta = vEcalTpgs->at(j).Eta();
    double phi = vEcalTpgs->at(j).Phi();
    double pt  = vEcalTpgs->at(j).Pt();
    h2EcalTpgs->Fill(eta, phi, pt);

    if(pt>10){
      std::cout<<"vEcalTpgs->at(j).Pt() "<<vEcalTpgs->at(j).Pt()
	       <<" eta "<<vEcalTpgs->at(j).Eta()
	       <<" phi "<<vEcalTpgs->at(j).Phi()<<std::endl;
      
      
      std::ostringstream strs;
      strs << pt;
      std::string text = strs.str();
      eta += 0.01;
      phi += 0.01;
      TPaveText *tempText = new TPaveText( eta, phi, eta+0.1, phi+0.1 );
      tempText->AddText(text.c_str());
      tempText->SetFillColor(0);
      tempText->SetLineColor(0);
      tempText->SetShadowColor(0);
      tempText->SetTextColor(kBlue);
      ecalTpgText.push_back(tempText);
    }
  }
  
  for (UInt_t j = 0; j < vHcalTpgs->size(); ++j) {
    h2HcalTpgs->Fill(vHcalTpgs->at(j).Eta(), vHcalTpgs->at(j).Phi(), vHcalTpgs->at(j).Pt());

    if(vHcalTpgs->at(j).Pt()>10){
      std::cout<<"vHcalTpgs->at(j).Pt() "<<vHcalTpgs->at(j).Pt()
	       <<" eta "<<vHcalTpgs->at(j).Eta()
	       <<" phi "<<vHcalTpgs->at(j).Phi()<<std::endl;
    }
    //std::cout<<"vHcalTpgs->at(j).Pt() "<<vHcalTpgs->at(j).Pt()<<std::endl;
  }    
  
  for (UInt_t j = 0; j < vCaloClusters->size(); ++j) {

    h2CaloClusters->Fill(vCaloClusters->at(j).Eta(), vCaloClusters->at(j).Phi(), vCaloClusters->at(j).Pt());

    //std::cout<<"vHcalTpgs->at(j).Pt() "<<vHcalTpgs->at(j).Pt()<<std::endl;
  }    
  
  for (UInt_t j = 0; j < vL1Jets->size(); ++j) {
    
    //std::cout<<"jets eta "<<vL1Jets->at(j).Eta()<<" phi: "<<vL1Jets->at(j).Phi()<<std::endl;
    double eta = vL1Jets->at(j).Eta();
    double phi = vL1Jets->at(j).Phi();
    double pt  = vL1Jets->at(j).Pt();
    TEllipse *temp = new TEllipse(eta, phi, 0.5, 0.5);
    jets.push_back(temp);
    
    std::ostringstream strs;
    strs << pt << " GeV";
    std::string text = strs.str();
    
    eta += 0.5;
    phi += 0.0;
    TPaveText *tempText = new TPaveText( eta + 0.05, phi-0.1, eta+0.3, phi+0.2 );
    tempText->AddText(text.c_str());
    tempText->SetFillColor(0);
    tempText->SetLineColor(0);
    tempText->SetShadowColor(0);
    tempText->SetTextColor(kRed);
    l1JetText.push_back(tempText);
  }    

    for (UInt_t j = 0; j < vRecoJets->size(); ++j) {
      //std::cout<<"jets eta "<<vL1Jets->at(j).Eta()<<" phi: "<<vL1Jets->at(j).Phi()<<std::endl;
      double eta = vRecoJets->at(j).Eta();
      double phi = vRecoJets->at(j).Phi();
      double pt  = vRecoJets->at(j).Pt();
      double dr  = vRecoJetsDr->at(j);
      TEllipse *temp = new TEllipse(eta, phi, dr, dr);
      recoJets.push_back(temp);

      std::ostringstream strs;
      strs << pt << " GeV";
      std::string text = strs.str();

      eta += dr;
      phi += dr;
      TPaveText *tempText = new TPaveText( eta-0.1, phi, eta+0.3, phi+0.2 );
      tempText->AddText(text.c_str());
      tempText->SetFillColor(0);
      tempText->SetLineColor(0);
      tempText->SetShadowColor(0);
      tempText->SetTextColor(kBlue);
      recoJetText.push_back(tempText);
    }    



  h2->GetXaxis()->SetAxisColor(17);
  h2->GetYaxis()->SetAxisColor(17);

  h->Draw(); h2->Draw("BOX");
  DrawRegionLines();
  DrawTowerLines();
  h2RecoTaus->SetMarkerColor(kOrange);
  h2RecoTaus->SetFillColor(kOrange);
  h2RecoTaus->SetMarkerStyle(33);
  h2RecoTaus->Draw("SAME p0");

  h2SignalTauCands->SetMarkerColor(kBlack);
  h2SignalTauCands->SetFillColor(kBlack);
  h2SignalTauCands->SetMarkerStyle(3);
  h2SignalTauCands->Draw("SAME p0");


  h2HcalTpgs->SetFillColor(kMagenta);
  h2HcalTpgs->Draw("SAME BOX");
  //h2EcalTpgs->SetFillColor(kRed);
  h2EcalTpgs->SetFillColorAlpha(kRed, 0.75);
  h2EcalTpgs->Draw("SAME BOX");

  //h2CaloClusters->SetFillColor(kOrange);
  h2CaloClusters->SetFillStyle(3144);
  h2CaloClusters->SetFillColorAlpha(kOrange, 0.75);
  h2CaloClusters->Draw("SAME BOX");

  //h2RlxTaus->SetFillColor(kGreen);
  h2RlxTaus->SetFillColorAlpha(kGreen, 0.75);
  h2RlxTaus->Draw("SAME BOX");

  float xR=0.8;
  TLegend *l = new TLegend(xR,0.8,xR+0.2,1.0);
  l->AddEntry(h2,"Regions","F");
  l->AddEntry(h2RecoTaus,"Reco Tau","P");
  l->AddEntry(h2RlxTaus,"Rlx Taus","F");
  l->AddEntry(h2EcalTpgs,"ECAL TPGs","F");
  l->AddEntry(h2HcalTpgs,"HCAL TPGs","F");
  l->AddEntry(h2CaloClusters,"caloClusters","F");
  l->Draw();
  h2->GetXaxis()->SetTitle("eta");
  h2->GetYaxis()->SetTitle("phi");
  //c1->RedrawAxis();
  for (UInt_t j = 0; j < vL1Jets->size(); ++j) {
    //if(vL1Jets->at(j).Pt()<15)continue;
    jets.at(j)->SetLineColor(kRed);  
    jets.at(j)->SetLineWidth(1);  
    jets.at(j)->SetLineStyle(2);  
    jets.at(j)->SetFillColor(0);  
    jets.at(j)->SetFillColorAlpha(kWhite, 0.01);
    jets.at(j)->Draw();
    //l1JetText.at(j)->Draw();
    l1JetText.at(j)->SetFillColorAlpha(kWhite, 0.01);
    l1JetText.at(j)->SetLineColorAlpha(kWhite, 0.01);
  }

  for (UInt_t j = 0; j < vRecoJets->size(); ++j) {
    recoJets.at(j)->SetLineColor(kBlue);  
    recoJets.at(j)->SetLineWidth(2);  
    recoJets.at(j)->SetLineStyle(6);  
    recoJets.at(j)->SetFillColor(0);  
    recoJets.at(j)->SetFillColorAlpha(kWhite, 0.01);
    recoJets.at(j)->Draw();
    //recoJetText.at(j)->Draw();
    recoJetText.at(j)->SetFillColorAlpha(kWhite, 0.01);
    recoJetText.at(j)->SetLineColorAlpha(kWhite, 0.01);
  }

  
  for (UInt_t j = 0; j < ecalTpgText.size(); ++j) {
    //ecalTpgText.at(j)->Draw();
    ecalTpgText.at(j)->SetFillColorAlpha(kWhite, 0.01);
    ecalTpgText.at(j)->SetLineColorAlpha(kWhite, 0.001);
    }

  for (UInt_t j = 0; j < rlxTauText.size(); ++j) {
    rlxTauText.at(j)->Draw();
    rlxTauText.at(j)->SetFillColorAlpha(kWhite, 0.01);
    rlxTauText.at(j)->SetLineColorAlpha(kWhite, 0.001);
    rlxTauText.at(j)->SetBorderSize(0);
    }

  // this gives you the value taking into account the zoom
  Double_t xmin = c1->GetUxmin(); 
  Double_t xmax = c1->GetUxmax();
  Double_t ymin = c1->GetUymin(); 
  Double_t ymax = c1->GetUymax();

  // and this gives you the "real" value of the axis (no zoom)
  Double_t xminLimit = h2->GetXaxis()->GetXmin();  
  Double_t xmaxLimit = h2->GetXaxis()->GetXmax();
  Double_t yminLimit = h2->GetMinimum();
  Double_t ymaxLimit = h2->GetMaximum();
  float recoTauEta = vRecoTaus->at(0).Eta();
  float recoTauPhi = vRecoTaus->at(0).Phi();
  float recoTauPt = vRecoTaus->at(0).Pt();

  TPaveText *tempText = new TPaveText( recoTauEta+0.05, recoTauPhi+0.05, recoTauEta + 0.15, recoTauPhi + 0.15 );
  std::ostringstream strs;
  strs << recoTauPt;
  std::string text = strs.str();
  tempText->AddText(text.c_str());
  tempText->SetFillColor(0);
  tempText->SetLineColor(0);
  tempText->SetShadowColor(0);
  tempText->SetTextColor(kOrange);
  tempText->SetFillColorAlpha(kWhite, 0.01);
  tempText->SetLineColorAlpha(kWhite, 0.01);
  tempText->SetBorderSize(0);

  h2->GetXaxis()->SetRangeUser(recoTauEta - 0.55, recoTauEta + 0.55);
  h2->GetYaxis()->SetRangeUser(recoTauPhi - 0.55, recoTauPhi + 0.55);
  tempText->Draw();

  char* saveFile = new char[100];
  sprintf(saveFile,"/Users/isobelojalvo/Documents/work/Analysis/Level1/TopologyStudies/highPtTaus-Apr24/Event-%u.png",event);
  
  c1->SaveAs(saveFile);

  //time to clean up!
  /* 
  delete f;
  delete vpx      ; 
  delete vRecoTaus    ;
  delete vSignalTauCands;
  delete vRecoJets      ;
  delete vRecoJetsDr    ;
  delete vRlxTaus       ;
  delete vEcalTpgs      ;
  delete vHcalTpgs      ;
  delete vL1Jets        ;
  delete vCaloClusters  ;
  delete c1;
  jets.clear();       
  l1JetText.clear();       
  recoJets.clear();       
  recoJetText.clear();       
  ecalTpgText.clear();       
  hcalTpgText.clear();       
  rlxTauText.clear();       
  */
  //delete l;
  //delete saveFile;
}
