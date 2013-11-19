#ifndef DiJetVarProducer_h
#define DiJetVarProducer_h

#include "DQM/DataScouting/interface/ScoutingAnalyzerBase.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"

#include <vector>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TProfile.h"

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
  std::string recoJetCorrector;
  edm::InputTag metCollectionTag_;
  edm::InputTag metCleanCollectionTag_;
  double etawidejets_;
  double ptwidejets_;
  double detawidejets_;
  double dphiwidejets_;
  double maxEMfraction_;
  double maxHADfraction_;
  // trigger conditions
  triggerExpression::Evaluator * HLTpathMain_;
  triggerExpression::Evaluator * HLTpathMonitor_;
    
  // cache some data from the Event for faster access by the trigger conditions
  triggerExpression::Data triggerConfiguration_;

  TProfile *hPUCorr[30];
  TH1F *m_cutFlow;
  TH1F *m_MjjWide_finalSel;
  TH1F *m_MjjWide_finalSel_varbin;
  TH1F *m_MjjWide_finalSel_WithoutNoiseFilter;
  TH1F *m_MjjWide_finalSel_WithoutNoiseFilter_varbin;
  TH1F *m_MjjWide_deta_0p0_0p5;
  TH1F *m_MjjWide_deta_0p5_1p0;
  TH1F *m_MjjWide_deta_1p0_1p5;
  TH1F *m_MjjWide_deta_1p5_2p0;
  TH1F *m_MjjWide_deta_2p0_2p5;
  TH1F *m_MjjWide_deta_2p5_3p0;
  TH1F *m_MjjWide_deta_3p0_inf;
  TH1F *m_MjjWide_den_NOdeta;
  TH1F *m_MjjWide_num_NOdeta;
  TH1F *m_MjjWide_den_detaL4;
  TH1F *m_MjjWide_num_detaL4;
  TH1F *m_MjjWide_den_detaL3;
  TH1F *m_MjjWide_num_detaL3;
  TH1F *m_MjjWide_den_detaL2;
  TH1F *m_MjjWide_num_detaL2;
  TH1F *m_MjjWide_den;
  TH1F *m_MjjWide_num;
  TH1F *m_DetajjWide_finalSel;
  TH1F *m_DetajjWide;
  TH1F *m_DphijjWide_finalSel;
  TH1F *m_selJets_pt;
  TH1F *m_selJets_eta;
  TH1F *m_selJets_phi;
  TH1F *m_selJets_hadEnergyFraction;
  TH1F *m_selJets_emEnergyFraction;
  TH1F *m_selJets_towersArea;
  TH1F *m_metDiff;
  TH1F *m_metCases;
  TH1F *m_metCaseNoMetClean;
  TH1F *m_HT_inclusive;
  TH1F *m_HT_finalSel;

  TH2F *m_DetajjVsMjjWide;
  TH2F *m_DetajjVsMjjWide_rebin;
  TH2F *m_metVSmetclean;

};

#endif //DiJetVarProducer_h
