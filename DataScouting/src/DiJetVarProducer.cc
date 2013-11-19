#include "DQM/DataScouting/interface/DiJetVarProducer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"


#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

#include "DataFormats/Common/interface/AssociationMap.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"

#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TString.h"


#include <memory>
#include <vector>

//
// constructors and destructor
//
//DiJetVarProducer::DiJetVarProducer(const edm::ParameterSet& iConfig, edm::Event& iEvent, const edm::EventSetup& iSetup) :
DiJetVarProducer::DiJetVarProducer(const edm::ParameterSet& conf) :
  inputJetTag_(conf.getParameter<edm::InputTag>("inputJetTag")),
  wideJetDeltaR_(conf.getParameter<double>("wideJetDeltaR")),
  recoJetCorrector(conf.getParameter<std::string>("jetCorrections")),
  metCollectionTag_(conf.getUntrackedParameter<edm::InputTag>("metCollectionTag")),
  metCleanCollectionTag_(conf.getUntrackedParameter<edm::InputTag>("metCleanCollectionTag")),
  etawidejets_             (conf.getParameter<double>("etawidejets")),
  ptwidejets_              (conf.getParameter<double>("ptwidejets")),
  detawidejets_            (conf.getParameter<double>("detawidejets")),
  dphiwidejets_            (conf.getParameter<double>("dphiwidejets")),
  maxEMfraction_           (conf.getParameter<double>("maxEMfraction")),
  maxHADfraction_          (conf.getParameter<double>("maxHADfraction")),
  HLTpathMain_             (triggerExpression::parse( conf.getParameter<std::string>("HLTpathMain") )),
  HLTpathMonitor_          (triggerExpression::parse( conf.getParameter<std::string>("HLTpathMonitor") )),
  triggerConfiguration_    (conf.getParameterSet("triggerConfiguration"))

{
  //register your products
  //produces<std::vector<double> >("dijetvariables");
  edm::Service<TFileService> fs;
  TString s;
  for(int i=0; i<20; i++)
  {
    s.Form("%d",i+1);
    hPUCorr[i]  = fs->make<TProfile>("PUCorr_for_NPV:"+s,"Profile of PU Corrections versus eta for NPV "+s,20,-3.5,3.5,0,1);
  }
  
  /////////////////////////////////////////BOOK THE HISTOGRAMS/////////////////////////////////////////////
  // ==> TO BE UPDATED FOR sqrt(s)=8 TeV
  const int N_mass_bins=83;
  float massBins[N_mass_bins+1] = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7000};

  // 1D histograms
  m_cutFlow = fs->make<TH1F>("h1_cutFlow","Cut Flow",7, 0., 7.);
  m_cutFlow->Sumw2();
  m_cutFlow->GetXaxis()->SetBinLabel(1,"No cut");
  m_cutFlow->GetXaxis()->SetBinLabel(2,"N(WideJets)>=2");
  m_cutFlow->GetXaxis()->SetBinLabel(3,"|#eta|<2.5 , pT>30");
  m_cutFlow->GetXaxis()->SetBinLabel(4,"|#Delta#eta|<1.3");
  m_cutFlow->GetXaxis()->SetBinLabel(5,"JetID");
  m_cutFlow->GetXaxis()->SetBinLabel(6,"|#Delta#phi|>#pi/3");
  m_cutFlow->GetXaxis()->SetBinLabel(7,"|met-metClean|>0.1");

  m_MjjWide_finalSel = fs->make<TH1F>("h1_MjjWide_finalSel","M_{jj} WideJets (final selection)", 8000, 0., 8000.);
  m_MjjWide_finalSel->GetXaxis()->SetTitle("M_{jj} WideJets [GeV]");
  m_MjjWide_finalSel->GetYaxis()->SetTitle("Number of events");
  m_MjjWide_finalSel->Sumw2();

  m_MjjWide_finalSel_varbin = fs->make<TH1F>("h1_MjjWide_finalSel_varbin","M_{jj} WideJets (final selection)",N_mass_bins, massBins);
  m_MjjWide_finalSel_varbin->GetXaxis()->SetTitle("M_{jj} WideJets [GeV]");
  m_MjjWide_finalSel_varbin->GetYaxis()->SetTitle("Number of events");
  m_MjjWide_finalSel_varbin->Sumw2();

  m_MjjWide_finalSel_WithoutNoiseFilter = fs->make<TH1F>( "h1_MjjWide_finalSel_WithoutNoiseFilter","M_{jj} WideJets (final selection, without noise filters)",8000, 0., 8000.);
  m_MjjWide_finalSel_WithoutNoiseFilter->GetXaxis()->SetTitle("M_{jj} WideJets [GeV]");
  m_MjjWide_finalSel_WithoutNoiseFilter->GetYaxis()->SetTitle("Number of events");
  m_MjjWide_finalSel_WithoutNoiseFilter->Sumw2();


  m_MjjWide_finalSel_WithoutNoiseFilter_varbin = fs->make<TH1F>("h1_MjjWide_finalSel_WithoutNoiseFilter_varbin", "M_{jj} WideJets (final selection, without noise filters)", N_mass_bins, massBins);
  m_MjjWide_finalSel_WithoutNoiseFilter_varbin->GetXaxis()->SetTitle("M_{jj} WideJets [GeV]");
  m_MjjWide_finalSel_WithoutNoiseFilter_varbin->GetYaxis()->SetTitle("Number of events");
  m_MjjWide_finalSel_WithoutNoiseFilter_varbin->Sumw2();
                    
  m_MjjWide_deta_0p0_0p5 = fs->make<TH1F>("h1_MjjWide_deta_0p5_1p0", "M_{jj} WideJets (0.0<=#Delta#eta<0.5)", 8000, 0., 8000.);
  m_MjjWide_deta_0p0_0p5->GetXaxis()->SetTitle("M_{jj} WideJets [GeV]");
  m_MjjWide_deta_0p0_0p5->GetYaxis()->SetTitle("Number of events");
  m_MjjWide_deta_0p0_0p5->Sumw2();

  m_MjjWide_deta_0p5_1p0 = fs->make<TH1F>("h1_MjjWide_deta_0p5_1p0", "M_{jj} WideJets (0.5<=#Delta#eta<1.0)", 8000, 0., 8000.);
  m_MjjWide_deta_0p5_1p0->GetXaxis()->SetTitle("M_{jj} WideJets [GeV]");
  m_MjjWide_deta_0p5_1p0->GetYaxis()->SetTitle("Number of events");
  m_MjjWide_deta_0p5_1p0->Sumw2();

  m_MjjWide_deta_1p0_1p5 = fs->make<TH1F>("h1_MjjWide_deta_1p0_1p5", "M_{jj} WideJets (1.0<=#Delta#eta<1.5)", 8000, 0., 8000.);
  m_MjjWide_deta_1p0_1p5->GetXaxis()->SetTitle("M_{jj} WideJets [GeV]");
  m_MjjWide_deta_1p0_1p5->GetYaxis()->SetTitle("Number of events");
  m_MjjWide_deta_1p0_1p5->Sumw2();

  m_MjjWide_deta_1p5_2p0 = fs->make<TH1F>("h1_MjjWide_deta_1p5_2p0", "M_{jj} WideJets (1.5<=#Delta#eta<2.0)", 8000, 0., 8000.);
  m_MjjWide_deta_1p5_2p0->GetXaxis()->SetTitle("M_{jj} WideJets [GeV]");
  m_MjjWide_deta_1p5_2p0->GetYaxis()->SetTitle("Number of events");
  m_MjjWide_deta_1p5_2p0->Sumw2();

  m_MjjWide_deta_2p0_2p5 = fs->make<TH1F>("h1_MjjWide_deta_2p0_2p5", "M_{jj} WideJets (2.0<=#Delta#eta<2.5)", 8000, 0., 8000.);
  m_MjjWide_deta_2p0_2p5->GetXaxis()->SetTitle("M_{jj} WideJets [GeV]");
  m_MjjWide_deta_2p0_2p5->GetYaxis()->SetTitle("Number of events");
  m_MjjWide_deta_2p0_2p5->Sumw2();

  m_MjjWide_deta_2p5_3p0 = fs->make<TH1F>("h1_MjjWide_deta_2p5_3p0", "M_{jj} WideJets (2.5<=#Delta#eta<3.0)", 8000, 0., 8000.);
  m_MjjWide_deta_2p5_3p0->GetXaxis()->SetTitle("M_{jj} WideJets [GeV]");
  m_MjjWide_deta_2p5_3p0->GetYaxis()->SetTitle("Number of events");
  m_MjjWide_deta_2p5_3p0->Sumw2();

  m_MjjWide_deta_3p0_inf = fs->make<TH1F>("h1_MjjWide_deta_3p0_inf", "M_{jj} WideJets (3.0<=#Delta#eta<#infty)", 8000, 0., 8000.);
  m_MjjWide_deta_3p0_inf->GetXaxis()->SetTitle("M_{jj} WideJets [GeV]");
  m_MjjWide_deta_3p0_inf->GetYaxis()->SetTitle("Number of events");
  m_MjjWide_deta_3p0_inf->Sumw2();

  m_MjjWide_den_NOdeta = fs->make<TH1F>("h1_MjjWide_den_NOdeta","HLT Efficiency Studies (no deta cut)", 400, 0., 2000.);
  m_MjjWide_den_NOdeta->GetXaxis()->SetTitle("M_{jj} WideJets [GeV]");
  m_MjjWide_den_NOdeta->GetYaxis()->SetTitle("Number of events");
  m_MjjWide_den_NOdeta->Sumw2();
  
  m_MjjWide_num_NOdeta = fs->make<TH1F>("h1_MjjWide_num_NOdeta", "HLT Efficiency Studies (no deta cut)", 400, 0., 2000.);
  m_MjjWide_num_NOdeta->GetXaxis()->SetTitle("M_{jj} WideJets [GeV]");
  m_MjjWide_num_NOdeta->GetYaxis()->SetTitle("Number of events");
  m_MjjWide_num_NOdeta->Sumw2();

  m_MjjWide_den_detaL4 = fs->make<TH1F>("h1_MjjWide_den_detaL4","HLT Efficiency Studies (deta cut < 4.0)", 400, 0., 2000.);
  m_MjjWide_den_detaL4->GetXaxis()->SetTitle("M_{jj} WideJets [GeV]");
  m_MjjWide_den_detaL4->GetYaxis()->SetTitle("Number of events");
  m_MjjWide_den_detaL4->Sumw2();

  m_MjjWide_num_detaL4 = fs->make<TH1F>("h1_MjjWide_num_detaL4", "HLT Efficiency Studies (deta cut < 4.0)", 400, 0., 2000.);
  m_MjjWide_num_detaL4->GetXaxis()->SetTitle("M_{jj} WideJets [GeV]");
  m_MjjWide_num_detaL4->GetYaxis()->SetTitle("Number of events");
  m_MjjWide_num_detaL4->Sumw2();         

  m_MjjWide_den_detaL3 = fs->make<TH1F>( "h1_MjjWide_den_detaL3", "HLT Efficiency Studies (deta cut < 3.0)", 400, 0., 2000.);
  m_MjjWide_den_detaL3->GetXaxis()->SetTitle("M_{jj} WideJets [GeV]");
  m_MjjWide_den_detaL3->GetYaxis()->SetTitle("Number of events");
  m_MjjWide_den_detaL3->Sumw2();
  
  m_MjjWide_num_detaL3 = fs->make<TH1F>("h1_MjjWide_num_detaL3", "HLT Efficiency Studies (deta cut < 3.0)", 400, 0., 2000.);
  m_MjjWide_num_detaL3->GetXaxis()->SetTitle("M_{jj} WideJets [GeV]");
  m_MjjWide_num_detaL3->GetYaxis()->SetTitle("Number of events");
  m_MjjWide_num_detaL3->Sumw2();

  m_MjjWide_den_detaL2 = fs->make<TH1F>("h1_MjjWide_den_detaL2", "HLT Efficiency Studies (deta cut < 2.0)", 400, 0., 2000.);
  m_MjjWide_den_detaL2->GetXaxis()->SetTitle("M_{jj} WideJets [GeV]");
  m_MjjWide_den_detaL2->GetYaxis()->SetTitle("Number of events");
  m_MjjWide_den_detaL2->Sumw2();
  
  m_MjjWide_num_detaL2 = fs->make<TH1F>("h1_MjjWide_num_detaL2", "HLT Efficiency Studies (deta cut < 2.0)", 400, 0., 2000.);
  m_MjjWide_num_detaL2->GetXaxis()->SetTitle("M_{jj} WideJets [GeV]");
  m_MjjWide_num_detaL2->GetYaxis()->SetTitle("Number of events");
  m_MjjWide_num_detaL2->Sumw2(); 

  m_MjjWide_den = fs->make<TH1F>("h1_MjjWide_den", "HLT Efficiency Studies (default deta cut)", 400, 0., 2000.);
  m_MjjWide_den->GetXaxis()->SetTitle("M_{jj} WideJets [GeV]");
  m_MjjWide_den->GetYaxis()->SetTitle("Number of events");
  m_MjjWide_den->Sumw2();

  m_MjjWide_num = fs->make<TH1F>("h1_MjjWide_num", "HLT Efficiency Studies (default deta cut)", 400, 0., 2000.);
  m_MjjWide_num->GetXaxis()->SetTitle("M_{jj} WideJets [GeV]");
  m_MjjWide_num->GetYaxis()->SetTitle("Number of events");
  m_MjjWide_num->Sumw2();

  m_DetajjWide_finalSel = fs->make<TH1F>("h1_DetajjWide_finalSel", "#Delta#eta_{jj} WideJets (final selection)", 100, 0., 5.);
  m_DetajjWide_finalSel->GetXaxis()->SetTitle("#Delta#eta_{jj} WideJets");
  m_DetajjWide_finalSel->GetYaxis()->SetTitle("Number of events");
  m_MjjWide_den_detaL4->Sumw2();

  m_DetajjWide = fs->make<TH1F>("h1_DetajjWide", "#Delta#eta_{jj} WideJets (final selection except #Delta#eta cut)", 100, 0., 5.);
  m_DetajjWide->GetXaxis()->SetTitle("#Delta#eta_{jj} WideJets");
  m_DetajjWide->GetYaxis()->SetTitle("Number of events");
  m_DetajjWide->Sumw2();

  m_DphijjWide_finalSel = fs->make<TH1F>("h1_DphijjWide_finalSel", "#Delta#phi_{jj} WideJets (final selection)", 100, 0., TMath::Pi()+0.0001);
  m_DphijjWide_finalSel->GetXaxis()->SetTitle("#Delta#phi_{jj} WideJets [rad.]");
  m_DphijjWide_finalSel->GetYaxis()->SetTitle("Number of events");
  m_DphijjWide_finalSel->Sumw2();


  m_selJets_pt = fs->make<TH1F>("h1_selJets_pt", "Selected CaloJets", 500, 0., 5000.);
  m_selJets_pt->GetXaxis()->SetTitle("Jet Pt [GeV]");
  m_selJets_pt->GetYaxis()->SetTitle("Number of events");
  m_selJets_pt->Sumw2();

  m_selJets_eta = fs->make<TH1F>("h1_selJets_eta", "Selected CaloJets", 100, -5., 5.);
  m_selJets_eta->GetXaxis()->SetTitle("#eta");
  m_selJets_eta->GetYaxis()->SetTitle("Number of events");
  m_selJets_eta->Sumw2();

  m_selJets_phi = fs->make<TH1F>("h1_selJets_phi", "Selected CaloJets", 100, -TMath::Pi(), TMath::Pi());
  m_selJets_phi->GetXaxis()->SetTitle("#phi (rad.)");
  m_selJets_phi->GetYaxis()->SetTitle("Number of events");
  m_selJets_phi->Sumw2();

  m_selJets_hadEnergyFraction = fs->make<TH1F>("h1_selJets_hadEnergyFraction", "Selected CaloJets", 110, 0., 1.1);
  m_selJets_hadEnergyFraction->GetXaxis()->SetTitle("HAD Energy Fraction");
  m_selJets_hadEnergyFraction->GetYaxis()->SetTitle("Number of events");
  m_selJets_hadEnergyFraction->Sumw2();

  m_selJets_emEnergyFraction = fs->make<TH1F>("h1_selJets_emEnergyFraction", "Selected CaloJets", 110, 0., 1.1);
  m_selJets_emEnergyFraction->GetXaxis()->SetTitle("EM Energy Fraction");
  m_selJets_emEnergyFraction->GetYaxis()->SetTitle("Number of events");
  m_selJets_emEnergyFraction->Sumw2();

  m_selJets_towersArea = fs->make<TH1F>("h1_selJets_towersArea", "Selected CaloJets", 200, 0., 2.);
  m_selJets_towersArea->GetXaxis()->SetTitle("towers area");
  m_selJets_towersArea->GetYaxis()->SetTitle("Number of events");
  m_selJets_towersArea->Sumw2();

  m_metDiff = fs->make<TH1F>("h1_metDiff", "Met - MetCleaned", 500, -1000., 1000.);
  m_metDiff->GetXaxis()->SetTitle("met - metcleaned [GeV]");
  m_metDiff->GetYaxis()->SetTitle("Number of events");
  m_metDiff->Sumw2();

  m_metCases = fs->make<TH1F>( "h1_metCases", "Met cases", 3, 0., 3.);
  m_metCases->GetXaxis()->SetTitle("case");
  m_metCases->GetYaxis()->SetTitle("Number of events");
  m_metCases->Sumw2();

  m_metCases->GetXaxis()->SetBinLabel(1,"met , metclean");
  m_metCases->GetXaxis()->SetBinLabel(2,"met , !metclean");
  m_metCases->GetXaxis()->SetBinLabel(3,"!met , metclean");
  m_metCases->Sumw2();


  m_metCaseNoMetClean = fs->make<TH1F>( "h1_metCaseNoMetClean", "Met - MetCleaned", 1000, 0., 2000.);
  m_metCaseNoMetClean->GetXaxis()->SetTitle("MET [GeV]");
  m_metCaseNoMetClean->GetYaxis()->SetTitle("Number of events");
  m_metCaseNoMetClean->Sumw2();

  m_HT_inclusive = fs->make<TH1F>( "h1_HT_inclusive", "HT (inclusive)", 150, 0., 15000.);
  m_HT_inclusive->GetXaxis()->SetTitle("HT [GeV]");
  m_HT_inclusive->GetYaxis()->SetTitle("Number of events");
  m_HT_inclusive->Sumw2(); 

  m_HT_finalSel = fs->make<TH1F>( "h1_HT_finalSel", "HT (final selection)", 150, 0., 15000.);
  m_HT_finalSel->GetXaxis()->SetTitle("HT [GeV]");
  m_HT_finalSel->GetYaxis()->SetTitle("Number of events");
  m_HT_finalSel->Sumw2(); 


  // 2D histograms
  m_DetajjVsMjjWide = fs->make<TH2F>("h2_DetajjVsMjjWide", "#Delta#eta_{jj} vs M_{jj} WideJets", 8000, 0., 8000., 100, 0., 5.);
  m_DetajjVsMjjWide->GetXaxis()->SetTitle("M_{jj} WideJets [GeV]");
  m_DetajjVsMjjWide->GetYaxis()->SetTitle("#Delta#eta_{jj} WideJets");
  m_DetajjVsMjjWide->Sumw2();

  m_DetajjVsMjjWide_rebin = fs->make<TH2F>("h2_DetajjVsMjjWide_rebin", "#Delta#eta_{jj} vs M_{jj} WideJets", 400, 0., 8000., 50, 0., 5.);
  m_DetajjVsMjjWide_rebin->GetXaxis()->SetTitle("M_{jj} WideJets [GeV]");
  m_DetajjVsMjjWide_rebin->GetYaxis()->SetTitle("#Delta#eta_{jj} WideJets");
  m_DetajjVsMjjWide_rebin->Sumw2();

  m_metVSmetclean = fs->make<TH2F>("h2_metVSmetclean", "MET clean vs MET", 100, 0., 2000., 100, 0., 2000.);
  m_metVSmetclean->GetXaxis()->SetTitle("MET [GeV]");
  m_metVSmetclean->GetYaxis()->SetTitle("MET clean [GeV]");
  m_metVSmetclean->Sumw2();

/////////////////////////////////////////BOOK THE HISTOGRAMS/////////////////////////////////////////////

  produces<std::vector<math::PtEtaPhiMLorentzVector> >("widejets");
  
  LogDebug("") << "Input Jet Tag: " << inputJetTag_.encode() << " ";
  LogDebug("") << "Radius Parameter Wide Jet: " << wideJetDeltaR_ << ".";
}

DiJetVarProducer::~DiJetVarProducer()
{

}

// ------------ method called to produce the data  ------------ 
void
DiJetVarProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace std;
   using namespace edm;
   using namespace reco;
  
   // ## The output collections        
   //std::auto_ptr<std::vector<double> > dijetvariables(new std::vector<double>); 
   std::auto_ptr<std::vector<math::PtEtaPhiMLorentzVector> > widejets(new std::vector<math::PtEtaPhiMLorentzVector>); 

   // ## Get jet collection
   edm::Handle<reco::CaloJetCollection> calojets_handle;
   edm::Handle<double> dsRho;
   edm::Handle<reco::VertexCollection> recVtxs;

   iEvent.getByLabel(inputJetTag_,calojets_handle);
   iEvent.getByLabel("hltKT6CaloJets","rho",dsRho);
   iEvent.getByLabel("hltPixelVertices",recVtxs);

   TLorentzVector wj1_tmp;
   TLorentzVector wj2_tmp;
   TLorentzVector wj1;
   TLorentzVector wj2;
   TLorentzVector wdijet;

   double pileupcorr(0);
   double scaleL2L3(0);
   double DeltaR1(0);
   double DeltaR2(0);
   double thisHT(0);
   double MJJWide = -1;
   double DeltaEtaJJWide = -1;
   double DeltaPhiJJWide = -1;

	  // ## Event Selection
	  bool pass_nocut=false;
	  bool pass_twowidejets=false;
	  bool pass_etaptcuts=false;
	  bool pass_deta=false;
	  bool pass_JetIDtwojets=true;
	  bool pass_dphi=false;
	  bool pass_metFilter=true;
	  //--
	  bool pass_deta_L4=false;
	  bool pass_deta_L3=false;
	  bool pass_deta_L2=false;

	  bool pass_fullsel_NOdeta=false;
	  bool pass_fullsel_detaL4=false;
	  bool pass_fullsel_detaL3=false;
	  bool pass_fullsel_detaL2=false;
	  bool pass_fullsel=false;
	  // No cut
	  pass_nocut=true;

   // ## Get jet correction
   const JetCorrector* corrector = JetCorrector::getJetCorrector(recoJetCorrector, iSetup);       
   
   // ## Wide Jet Algorithm
   // At least two jets
   if( calojets_handle->size() >=2 )
     {

       TLorentzVector jet1, jet2;

       reco::CaloJetCollection::const_iterator j1 = calojets_handle->begin();
       reco::CaloJetCollection::const_iterator j2 = j1; ++j2; 
       

       jet1.SetPtEtaPhiM(j1->pt(),j1->eta(),j1->phi(),j1->mass());
       jet2.SetPtEtaPhiM(j2->pt(),j2->eta(),j2->phi(),j2->mass());  

       // Create wide jets (radiation recovery algorithm)


	   for(reco::CaloJetCollection::const_iterator it = calojets_handle->begin(); it != calojets_handle->end(); ++it)
		{
		   TLorentzVector currentJet;
	       CaloJet correctedJet = *it;
	       int i=0;
	        pileupcorr= 1-((*dsRho)*correctedJet.jetArea())/correctedJet.pt();
	        if(pileupcorr > 0. || pileupcorr < 1.) correctedJet.scaleEnergy(pileupcorr);

	        if(recVtxs->size()<=20)hPUCorr[recVtxs->size()-1]->Fill(correctedJet.eta(),pileupcorr,1);
	        scaleL2L3  = corrector->correction(correctedJet,iEvent,iSetup);
	        correctedJet.scaleEnergy(scaleL2L3);
	        
	        m_selJets_pt->Fill( it->pt() );
	        m_selJets_eta->Fill( it->eta() );
	        m_selJets_phi->Fill( it->phi() );
	        m_selJets_hadEnergyFraction->Fill( it->energyFractionHadronic() );
	        m_selJets_emEnergyFraction->Fill( it->emEnergyFraction() );
	        m_selJets_towersArea->Fill( it->towersArea() );
	        thisHT += it->pt();

	         // Jet id two leading jets
      		if( (i==0 || i==1) && (it->energyFractionHadronic()>maxHADfraction_ || it->emEnergyFraction()>maxEMfraction_) )
	     		  pass_JetIDtwojets=false;

       

	      	currentJet.SetPtEtaPhiM(correctedJet.pt(),correctedJet.eta(),correctedJet.phi(),correctedJet.mass());

		    DeltaR1 = currentJet.DeltaR(jet1);
		    DeltaR2 = currentJet.DeltaR(jet2);

		   if(DeltaR1 < DeltaR2 && DeltaR1 < wideJetDeltaR_){
		   	wj1_tmp += currentJet;
		   } 
		   else if(DeltaR2 < wideJetDeltaR_){
		     wj2_tmp += currentJet;
		    }
	 	} 
   
	   // HT
	   m_HT_inclusive->Fill(thisHT);        

	    // Re-order the wide jets in pT
	   if( wj1_tmp.Pt() > wj2_tmp.Pt()){
	    wj1 = wj1_tmp;
	    wj2 = wj2_tmp;
		}
	   else{
	    wj1 = wj2_tmp;
	    wj2 = wj1_tmp;
	  	}

	   // Create dijet system
	   wdijet = wj1 + wj2;
	   MJJWide = wdijet.M();
	   DeltaEtaJJWide = fabs(wj1.Eta()-wj2.Eta());
	   DeltaPhiJJWide = fabs(wj1.DeltaPhi(wj2));         

       // Put widejets in the container 
       

       math::PtEtaPhiMLorentzVector wj1math(wj1.Pt(), wj1.Eta(), wj1.Phi(), wj1.M());
       math::PtEtaPhiMLorentzVector wj2math(wj2.Pt(), wj2.Eta(), wj2.Phi(), wj2.M());
       widejets->push_back( wj1math );
       widejets->push_back( wj2math );
    }


	  // Two wide jets
	  	if( widejets->size() >= 2)
	    {
	      // Two wide jets
	      pass_twowidejets=true;

	      // Eta/pt cuts
	      if( fabs(wj1.Eta())<etawidejets_ && wj1.Pt()>ptwidejets_ && fabs(wj2.Eta())<etawidejets_ && wj2.Pt()>ptwidejets_){
	      	pass_etaptcuts=true;
	      }
	      
	      // Deta cut
	      if( DeltaEtaJJWide < detawidejets_ )
	      	pass_deta=true;

	      // Dphi cut
	      if( DeltaPhiJJWide > dphiwidejets_ )
	      	pass_dphi=true;

	      // Other Deta cuts
	      if( DeltaEtaJJWide < 4.0 )
		       pass_deta_L4=true;

	      if( DeltaEtaJJWide < 3.0 )
		       pass_deta_L3=true;

	      if( DeltaEtaJJWide < 2.0 )
		       pass_deta_L2=true;
	    }


  // met
  edm::Handle<reco::CaloMETCollection> calomet_handle;
  iEvent.getByLabel(metCollectionTag_,calomet_handle);

  // MET cleaned
  edm::Handle<reco::CaloMETCollection> calometClean_handle;
  iEvent.getByLabel(metCleanCollectionTag_,calometClean_handle);
 
  if( calomet_handle.isValid() && calometClean_handle.isValid() )
    {
      if( fabs ( (calomet_handle->front()).pt() - (calometClean_handle->front()).pt() ) > 0.1 )
	       pass_metFilter=false;

      m_metCases->Fill(0);
      m_metDiff->Fill( (calomet_handle->front()).pt() - (calometClean_handle->front()).pt() );
      m_metVSmetclean->Fill( (calomet_handle->front()).pt() , (calometClean_handle->front()).pt() );      
    }
  else if( calomet_handle.isValid() && !calometClean_handle.isValid() )
    {
      m_metCases->Fill(1);
      m_metCaseNoMetClean->Fill((calomet_handle->front()).pt());
    }
  else if( !calomet_handle.isValid() && calometClean_handle.isValid() )
    {
      m_metCases->Fill(2);
    }
     // Full selection (no deta cut)
  if( pass_nocut && pass_twowidejets && pass_etaptcuts && pass_JetIDtwojets && pass_dphi && pass_metFilter )
      pass_fullsel_NOdeta=true;

  // Full selection (various deta cuts)
  if( pass_nocut && pass_twowidejets && pass_etaptcuts && pass_JetIDtwojets && pass_dphi && pass_metFilter && pass_deta_L4 )
    pass_fullsel_detaL4=true;
  if( pass_nocut && pass_twowidejets && pass_etaptcuts && pass_JetIDtwojets && pass_dphi && pass_metFilter && pass_deta_L3 )
    pass_fullsel_detaL3=true;
  if( pass_nocut && pass_twowidejets && pass_etaptcuts && pass_JetIDtwojets && pass_dphi && pass_metFilter && pass_deta_L2 )
    pass_fullsel_detaL2=true;

  // Full selection (default deta cut)
  if( pass_nocut && pass_twowidejets && pass_etaptcuts && pass_deta && pass_JetIDtwojets && pass_dphi && pass_metFilter )
    pass_fullsel=true;


  // ## Fill Histograms 

  // Cut-flow plot
  if( pass_nocut )
    m_cutFlow->Fill(0);
  if( pass_nocut && pass_twowidejets )
    m_cutFlow->Fill(1);
  if( pass_nocut && pass_twowidejets && pass_etaptcuts )
    m_cutFlow->Fill(2);
  if( pass_nocut && pass_twowidejets && pass_etaptcuts && pass_deta )
    m_cutFlow->Fill(3);
  if( pass_nocut && pass_twowidejets && pass_etaptcuts && pass_deta && pass_JetIDtwojets )
    m_cutFlow->Fill(4);
  if( pass_nocut && pass_twowidejets && pass_etaptcuts && pass_deta && pass_JetIDtwojets && pass_dphi )
    m_cutFlow->Fill(5);
  if( pass_fullsel )
    m_cutFlow->Fill(6);

  // After full selection
  if( pass_fullsel ) 
    {
      // 1D histograms
      m_MjjWide_finalSel->Fill(MJJWide);
      m_MjjWide_finalSel_varbin->Fill(MJJWide);
      m_DetajjWide_finalSel->Fill(DeltaEtaJJWide);
      m_DphijjWide_finalSel->Fill(DeltaPhiJJWide);      
      m_HT_finalSel->Fill(thisHT);      
    }      

  // After full selection (without "noise" filters)
  if( pass_nocut && pass_twowidejets && pass_etaptcuts && pass_deta )
    {
      m_MjjWide_finalSel_WithoutNoiseFilter->Fill(MJJWide);
      m_MjjWide_finalSel_WithoutNoiseFilter_varbin->Fill(MJJWide);
    }

  // After full selection (except DeltaEta cut)
  if( pass_fullsel_NOdeta )
    {
      // 1D histograms
      m_DetajjWide->Fill(DeltaEtaJJWide);

      if( DeltaEtaJJWide >= 0.0 && DeltaEtaJJWide < 0.5 )
	       m_MjjWide_deta_0p0_0p5->Fill(MJJWide);
      if( DeltaEtaJJWide >= 0.5 && DeltaEtaJJWide < 1.0 )
	       m_MjjWide_deta_0p5_1p0->Fill(MJJWide);
      if( DeltaEtaJJWide >= 1.0 && DeltaEtaJJWide < 1.5 )
	       m_MjjWide_deta_1p0_1p5->Fill(MJJWide);
      if( DeltaEtaJJWide >= 1.5 && DeltaEtaJJWide < 2.0 )
	       m_MjjWide_deta_1p5_2p0->Fill(MJJWide);
      if( DeltaEtaJJWide >= 2.0 && DeltaEtaJJWide < 2.5 )
	       m_MjjWide_deta_2p0_2p5->Fill(MJJWide);
      if( DeltaEtaJJWide >= 2.5 && DeltaEtaJJWide < 3.0 )
	       m_MjjWide_deta_2p5_3p0->Fill(MJJWide);
      if( DeltaEtaJJWide >= 3.0 )
	       m_MjjWide_deta_3p0_inf->Fill(MJJWide);
      
      // 2D histograms
      m_DetajjVsMjjWide->Fill(MJJWide,DeltaEtaJJWide);            
      m_DetajjVsMjjWide_rebin->Fill(MJJWide,DeltaEtaJJWide);
    }
    // ## Get Trigger Info

  // HLT paths for DataScouting
  //  DST_HT250_v1
  //  DST_L1HTT_Or_L1MultiJet_v1
  //  DST_Mu5_HT250_v1
  //  DST_Ele8_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_HT250_v1

  int HLTpathMain_fired    = -1;
  int HLTpathMonitor_fired = -1;
  

  if (HLTpathMain_ and HLTpathMonitor_ and triggerConfiguration_.setEvent(iEvent, iSetup)) {
    // invalid HLT configuration, skip the processing
    
    // if the L1 or HLT configurations have changed, (re)initialize the filters (including during the first event)
       if (triggerConfiguration_.configurationUpdated()) {
	      HLTpathMain_->init(triggerConfiguration_);
	      HLTpathMonitor_->init(triggerConfiguration_);
	      
	      // log the expanded configuration
	      // std::cout << "HLT selector configurations updated" << std::endl;
	      // std::cout << "HLTpathMain:    " << *HLTpathMain_    << std::endl;
	      // std::cout << "HLTpathMonitor: " << *HLTpathMonitor_ << std::endl;
	    }
	    
	    HLTpathMain_fired    = (*HLTpathMain_)(triggerConfiguration_);
	    HLTpathMonitor_fired = (*HLTpathMonitor_)(triggerConfiguration_);
	    
	    // The OR of the two should always be "1"
	    // std::cout << *HLTpathMain_ << ": " << HLTpathMain_fired << " -- " << *HLTpathMonitor_ << ": " << HLTpathMonitor_fired << std::endl;
  	}
  
  // ## Trigger Efficiency Curves

  //denominator - full sel NO deta cut
	  if( pass_fullsel_NOdeta && HLTpathMonitor_fired == 1 )
	    {
	    	m_MjjWide_den_NOdeta->Fill(MJJWide);
            //numerator  
		    if( HLTpathMain_fired == 1){
			  m_MjjWide_num_NOdeta->Fill(MJJWide);
			}
	    }

	  //denominator - full sel deta < 4.0
	  if( pass_fullsel_detaL4 && HLTpathMonitor_fired == 1 )
	  {
		      m_MjjWide_den_detaL4->Fill(MJJWide);
			  //numerator  
		      if( HLTpathMain_fired == 1)
			{
			  m_MjjWide_num_detaL4->Fill(MJJWide);
			}
	   }

	  //denominator - full sel deta < 3.0
	  if( pass_fullsel_detaL3 && HLTpathMonitor_fired == 1 )
	  {
		      m_MjjWide_den_detaL3->Fill(MJJWide);
			  //numerator  
		      if( HLTpathMain_fired == 1)
			{
			  m_MjjWide_num_detaL3->Fill(MJJWide);
			}
	   }

	  //denominator - full sel deta < 2.0
	  if( pass_fullsel_detaL2 && HLTpathMonitor_fired == 1 )
	    {
		      m_MjjWide_den_detaL2->Fill(MJJWide);

		      //numerator  
		      if( HLTpathMain_fired == 1)
			{
			  m_MjjWide_num_detaL2->Fill(MJJWide);
			}
	    }
  
	  //denominator - full sel default deta cut (typically 1.3)
	  if( pass_fullsel && HLTpathMonitor_fired == 1 )
	    {
		      m_MjjWide_den->Fill(MJJWide);

		      //numerator  
		      if( HLTpathMain_fired == 1)
			{
			  m_MjjWide_num->Fill(MJJWide);
			}
	    }

   //
   // ## Put objects in the Event
   //iEvent.put(dijetvariables, "dijetvariables");
   iEvent.put(widejets, "widejets");
}



DEFINE_FWK_MODULE(DiJetVarProducer);
