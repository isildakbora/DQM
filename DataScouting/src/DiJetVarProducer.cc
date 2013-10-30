#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "DQM/DataScouting/interface/DiJetVarProducer.h"

#include "TVector3.h"
#include "TLorentzVector.h"

#include <memory>
#include <vector>

//
// constructors and destructor
//
//DiJetVarProducer::DiJetVarProducer(const edm::ParameterSet& iConfig, edm::Event& iEvent, const edm::EventSetup& iSetup) :
DiJetVarProducer::DiJetVarProducer(const edm::ParameterSet& iConfig) :
  inputJetTag_(iConfig.getParameter<edm::InputTag>("inputJetTag")),
  wideJetDeltaR_(iConfig.getParameter<double>("wideJetDeltaR")),
  correctHLTJets_(iConfig.getParameter<bool>("correctHLTJets")),
  recoJetCorrector(iConfig.getParameter<std::string>("jetCorrections"))

{
  //register your products
  //produces<std::vector<double> >("dijetvariables");
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
   iEvent.getByLabel(inputJetTag_,calojets_handle);
   iEvent.getByLabel("hltKT6CaloJets","rho",dsRho);

   //cout << "size: " << calojets_handle->size() << " rho: " << *dsRho << endl;
   
   // ## Get jet correction
   const JetCorrector* corrector = JetCorrector::getJetCorrector(recoJetCorrector, iSetup);       
   
   // ## Wide Jet Algorithm
   // At least two jets
   if( calojets_handle->size() >=2 )
     {
       TLorentzVector wj1_tmp;
       TLorentzVector wj2_tmp;
       TLorentzVector wj1;
       TLorentzVector wj2;
       TLorentzVector wdijet;
       
       //        // Loop over all the input jets
       //        for(reco::CaloJetCollection::const_iterator it = calojets_handle->begin(); it != calojets_handle->end(); ++it)
       //        	 {
       // 	   cout << "jet: " << it->pt() << " " << it->eta() << " " << it->phi() << endl;
       //        	 }
       
       // Find two leading jets
       TLorentzVector jet1, jet2;

       reco::CaloJetCollection::const_iterator j1 = calojets_handle->begin();
       reco::CaloJetCollection::const_iterator j2 = j1; ++j2; 
       

       jet1.SetPtEtaPhiM(j1->pt(),j1->eta(),j1->phi(),j1->mass());
       jet2.SetPtEtaPhiM(j2->pt(),j2->eta(),j2->phi(),j2->mass());  

       // Create wide jets (radiation recovery algorithm)
   double pileupcorr;
   double scaleL2L3;
   double DeltaR1;
   double DeltaR2;

   for(reco::CaloJetCollection::const_iterator it = calojets_handle->begin(); it != calojets_handle->end(); ++it)
	 {
	   TLorentzVector currentJet;
     CaloJet correctedJet = *it;
     if(correctHLTJets_)
       {
        std::cout << "jet pT=" << it->energy() << std::endl;
        pileupcorr= 1-((*dsRho)*correctedJet.jetArea())/correctedJet.pt();
        std::cout << "pileupcorr=" << pileupcorr << std::endl;
        if(pileupcorr > 0. || pileupcorr < 1.) correctedJet.scaleEnergy(pileupcorr);
        std::cout << "jet pT=" << correctedJet.energy() << std::endl;

        scaleL2L3  = corrector->correction(correctedJet,iEvent,iSetup);
        std::cout << "JEC Scale=" << scaleL2L3 <<std::endl;
        correctedJet.scaleEnergy(scaleL2L3);
        std::cout << "jet pT=" << correctedJet.energy() << std::endl;
       }

      currentJet.SetPtEtaPhiM(correctedJet.pt(),correctedJet.eta(),correctedJet.phi(),correctedJet.mass());

	    DeltaR1 = currentJet.DeltaR(jet1);
	    DeltaR2 = currentJet.DeltaR(jet2);

	   if(DeltaR1 < DeltaR2 && DeltaR1 < wideJetDeltaR_) {
	     wj1_tmp += currentJet;
	   } else if(DeltaR2 < wideJetDeltaR_) {
	     wj2_tmp += currentJet;
	   }
	 }       

       // Re-order the wide jets in pT
       if( wj1_tmp.Pt() > wj2_tmp.Pt() )	
	 { 
	   wj1 = wj1_tmp;
	   wj2 = wj2_tmp;
	 }
       else
	 {
	   wj1 = wj2_tmp;
	   wj2 = wj1_tmp;
	 }

       // Create dijet system
       wdijet = wj1 + wj2;

       //        cout << "j1 wide: " << wj1.Pt() << " " << wj1.Eta() << " " << wj1.Phi() << " " << wj1.M() << endl;
       //        cout << "j2 wide: " << wj2.Pt() << " " << wj2.Eta() << " " << wj2.Phi() << " " << wj2.M() << endl;
       //        cout << "MJJWide: " << wdijet.M() << endl;
       //        cout << "DeltaEtaJJWide: " << fabs(wj1.Eta()-wj2.Eta()) << endl;
       //        cout << "DeltaPhiJJWide: " << fabs(wj1.DeltaPhi(wj2)) << endl;
 
       //        // Put variables in the container
       //        dijetvariables->push_back( wdijet.M() );                 //0 = MJJWide
       //        dijetvariables->push_back( fabs(wj1.Eta()-wj2.Eta()) );  //1 = DeltaEtaJJWide       
       //        dijetvariables->push_back( fabs(wj1.DeltaPhi(wj2)) );    //2 = DeltaPhiJJWide             

       // Put widejets in the container 
       

       math::PtEtaPhiMLorentzVector wj1math(wj1.Pt(), wj1.Eta(), wj1.Phi(), wj1.M());
       math::PtEtaPhiMLorentzVector wj2math(wj2.Pt(), wj2.Eta(), wj2.Phi(), wj2.M());
       widejets->push_back( wj1math );
       widejets->push_back( wj2math );
     }
   //    else
   //      {
   //        // Put variables in the container
   //        dijetvariables->push_back( -1 );                //0 = MJJWide
   //        dijetvariables->push_back( -1 );                //1 = DeltaEtaJJWide       
   //        dijetvariables->push_back( -1 );                //2 = DeltaPhiJJWide       
   //      }

   // ## Put objects in the Event
   //iEvent.put(dijetvariables, "dijetvariables");
   iEvent.put(widejets, "widejets");

}



DEFINE_FWK_MODULE(DiJetVarProducer);