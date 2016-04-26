#ifndef UCTRegionProcess_hh
#define UCTRegionProcess_hh

#include "DataFormats/L1CaloTrigger/interface/L1CaloRegion.h"

#include <vector>
#include <iostream>

#include "L1Trigger/L1TCaloLayer1/src/UCTRegion.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTGeometry.hh"
//#include "L1Trigger/Stage3Ntuplizer/plugins/triggerGeometryTools.hh"
class UCTRegionProcess: L1CaloRegion {
public:

  //UCTRegionProcess();

  UCTRegionProcess(L1CaloRegion inRegion)
    : L1CaloRegion(inRegion){
    regionSummary = inRegion.raw();
    //rctEta = inRegion.id().ieta();
    //rctPhi = inRegion.id().iphi();
    
  };

  ~UCTRegionProcess();

  // Packed data access

  const uint32_t rawData() const {return regionSummary;}
  const uint32_t location() const {return ((regionSummary & LocationBits) >> LocationShift);}
 
 //fix me
  const int hitCaloEta() const {
    uint32_t highestTowerLocation = 0;
    return highestTowerLocation;
  }

  const int hitCaloPhi() const {
    uint32_t highestTowerLocation = 0;
    return highestTowerLocation;
  }

  const UCTTowerIndex hitTowerIndex() const {
    return UCTTowerIndex(hitCaloEta(), hitCaloPhi());
  }

  //input region->id().ieta(), region->id().iphi(), raw data

  uint32_t rctEta = this->id().ieta();
  uint32_t rctPhi = this->id().iphi();

  float getFineRecoEta(){
    int towereta = -999;
    float eta = -999;

    uint32_t fineLocationEtaBit = ((location() & 0xC)>>2);

    int rctEta = (this->id().ieta());

    towereta = (rctEta - 11) * 4 + fineLocationEtaBit ;
    if(rctEta<11){
      towereta = -((10 - rctEta ) * 4 + fineLocationEtaBit +1) ;
    }
    towereta = towereta + 28;
    //std::cout<<"rctEta "<< rctEta<< " raw " << this->raw() <<" fineLocationetabit "<<fineLocationEtaBit << " tower eta "<<towereta<< " towerEtaMap " << towerEtaMap[towereta]<<std::endl;
    if((towereta > -1 && towereta < 57))
      {
	eta = towerEtaMap[towereta];
	return eta;
      }
    else{
      eta = -999;
      std::cout<<"Error!!! towereta out of bounds in UCTREgionProcess.h "<<std::endl;
      std::cout<<"rctEta "<<rctEta<<" towereta "<<towereta<<" hit2BitEta "<<fineLocationEtaBit<<std::endl;
      std::cout<<"exiting..."<<std::endl;
      exit(0);
    }
    return eta;

  }


  //towerphi here is 0 to 71
  float getFineRecoPhi(){
    int towerphi = -999;
    float phi = -999;
    int rctPhi = (this->id().iphi());

    uint32_t fineLocationPhiBit = ((location() & 0x3));
    towerphi = rctPhi*4 + fineLocationPhiBit;

    //std::cout<<"rctPhi "<<rctPhi<<" towerphi "<<towerphi<<" hit2BitPhi "<<fineLocationPhiBit<<" "<< std::endl;    

    if(towerphi>-1&&towerphi<72)
      phi = towerPhiRCTMap[towerphi];
    else{
      phi = -999;
      std::cout<<"Error!!! towerphi out of bounds in UCTREgionProcess.h "<<std::endl;
      std::cout<<"rctPhi "<<rctPhi<<" towerphi "<<towerphi<<" hit2BitPhi "<<fineLocationPhiBit<<std::endl;
      exit(0);}

    return phi;
  }



//from tower eta -28 to 28
float towerEtaMap[57]=   { 
    -2.913, //-2.739, 
    -2.565, -2.391, 2.217, 	//switch to smaller trigger towers here    
    -2.0445, -1.9575, -1.8705, -1.7835, -1.6965, 
    -1.6095, -1.5225, -1.4355, -1.3485, -1.2615, 
    -1.1745, -1.0875, -1.0005, -0.9135, -0.8265, 
    -0.7395, -0.6525, -0.5655, -0.4785, -0.3915, 
    -0.3045, -0.2175, -0.1305, -0.0435, 0.0435, 
    0.1305, 0.2175, 0.3045, 0.3915, 0.4785, 
    0.5655, 0.6525, 0.7395, 0.8265, 0.9135, 
    1.0005, 1.0875, 1.1745, 1.2615, 1.3485, 
    1.4355, 1.5225, 1.6095, 1.6965, 1.7835, 
    1.8705, 1.9575, 2.0445, 2.217, 2.391, 
    2.565, //2.739,
    2.913 
  };

  //this is for mapping from RCT only
float towerPhiRCTMap[72]=                        
  {-0.131, -0.044, 0.044, 0.131, 0.218, 0.305, 0.393, 0.480, 0.567, 0.654, 0.742, 0.829, 0.916, 1.004, 1.091, 1.178, 1.265, 1.353, 1.440, 1.527, 1.614, 1.702, 1.789, 1.876, 1.963, 2.051, 2.138, 2.225, 2.313, 2.400, 2.487, 2.574, 2.662, 2.749, 2.836, 2.923, 3.011, 3.098,
    -3.098, -3.011, -2.923, -2.836, -2.749, -2.662, -2.574, -2.487, -2.400, -2.313, -2.225, -2.138, -2.051, -1.963, -1.876, -1.789, -1.702, -1.614, -1.527, -1.440, -1.353, -1.265, -1.178, -1.091, -1.004, -0.916, -0.829, -0.742, -0.654, -0.567, -0.480, -0.393, -0.305, -0.218};


  const uint32_t compressedData() const {return regionSummary;}

  // Access functions for convenience
  // Note that the bit fields are limited in hardware

  const uint32_t et() const {return (RegionETMask & regionSummary);}
  const bool isEGammaLike() const {return !((RegionEGVeto & regionSummary) == RegionEGVeto);}
  const bool isTauLike() const {return !((RegionTauVeto & regionSummary) == RegionTauVeto);}

  // More access functions
  /*
  const uint32_t getCrate() const {return crate;}
  const uint32_t getCard() const {return card;}
  const uint32_t getRegion() const {return region;}

  const bool isNegativeEta() const {return negativeEta;}

  const UCTTower* getTower(UCTTowerIndex t) const {
    return getTower(t.first, t.second);
  }

  friend std::ostream& operator<<(std::ostream&, const UCTRegion&);
  */
private:

  // No default constructor is needed

  //

  // No copy constructor is needed

  //UCTRegion(const UCTRegion&);

  // No equality operator is needed

  const UCTRegion& operator=(const UCTRegion&);

protected:

  // Helper functions

  //const UCTTower* getTower(uint32_t caloEta, uint32_t caloPhi) const;

  // Region location definition

  //uint32_t crate;
  //uint32_t card;
  //uint32_t region;
  //bool negativeEta;

  // Owned region level data 

  //std::vector<UCTTower*> towers;

  uint32_t regionSummary;

};

#endif
