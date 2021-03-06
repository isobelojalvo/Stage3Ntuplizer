#include "DataFormats/L1CaloTrigger/interface/L1CaloRegion.h"
#include "L1Trigger/Stage3Ntuplizer/plugins/UCTRegionProcess.hh"
#include <vector>
#include <iostream>

#include "L1Trigger/L1TCaloLayer1/src/UCTRegion.hh"

//constructor
//UCTRegionProcess::UCTRegionProcess(L1CaloRegion inRegion){
//  regionSummary = inRegion.raw();
//  rctEta = inRegion.id().ieta();
//  rctPhi = inRegion.id().iphi();
//  //std::cout<<"inRegion.raw() "<<inRegion.raw()<<std::endl;
//};

//destructor
UCTRegionProcess::~UCTRegionProcess(){};
