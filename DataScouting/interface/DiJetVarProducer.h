#ifndef DiJetVarProducer_h
#define DiJetVarProducer_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "TVector3.h"
#include "TLorentzVector.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

#include<vector>

class DiJetVarProducer : public edm::EDProducer { 
 public: 
  //explicit DiJetVarProducer(const edm::ParameterSet& iConfig, edm::Event& iEvent, const edm::EventSetup& iSetup) ;
  explicit DiJetVarProducer(const edm::ParameterSet& iConfig) ;
  ~DiJetVarProducer();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  
 private:


  edm::InputTag inputJetTag_; // input tag jets
  double wideJetDeltaR_; // Radius parameter for wide jets
  bool correctHLTJets_; // True (applies JEC to HLT jets)
  std::string recoJetCorrector;
};

#endif //DiJetVarProducer_h
