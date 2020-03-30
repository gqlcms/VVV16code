// system include files
#include <iostream>
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"  

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include <DataFormats/JetReco/interface/Jet.h>
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


#include "EDBRChannels.h"
#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include <TFormula.h>

#define Pi 3.141593
using namespace std;
//using namespace edm;

#include "VVV/VVVTreeMaker.h"
#include "VVV/VVVHLTInfo.h"
#include "VVV/VVV_GenInfo.h"



    //
    // constructors and destructor


bool
VVVTreeMaker::looseJetID( const pat::Jet& j ) {
    // refer to https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data
    double NHF = j.neutralHadronEnergyFraction();
    double NEMF = j.neutralEmEnergyFraction();
    double CHF = j.chargedHadronEnergyFraction();
    //double MUF = j.muonEnergyFraction();
    double CEMF = j.chargedEmEnergyFraction();
    int NumConst = j.chargedMultiplicity()+j.neutralMultiplicity();
    int NumNeutralParticle =j.neutralMultiplicity();
    int CHM = j.chargedMultiplicity();
    double eta = j.eta();
	return ( (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4) && abs(eta)<=3.0 ) || (NEMF<0.90 && NumNeutralParticle>10 && abs(eta)>3.0 )  ;
}

bool
VVVTreeMaker::tightJetID( const pat::Jet& j ) {
    // refer to https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data
    if(j.pt()>0.){
    double NHF = j.neutralHadronEnergyFraction();
    double NEMF = j.neutralEmEnergyFraction();
    double CHF = j.chargedHadronEnergyFraction();
    //double MUF = j.muonEnergyFraction();
    //double CEMF = j.chargedEmEnergyFraction();
    int NumConst = j.chargedMultiplicity()+j.neutralMultiplicity();
    int NumNeutralParticle =j.neutralMultiplicity();
    int CHM = j.chargedMultiplicity();
    double eta = j.eta();
    return ((  (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 ) || abs(eta)>2.4) && abs(eta)<=2.7 ) || (NHF<0.99 && NEMF>0.02 && NumNeutralParticle>2 && abs(eta)>2.7 && abs(eta)<=3.0 ) || (NEMF<0.90 && NHF>0.02 &&NumNeutralParticle>10 && abs(eta)>3.0) ) ;
}
else{
return (0);
    }
}

bool
VVVTreeMaker::tightJetIDpuppi( const pat::Jet& j ) {
    // refer to https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data
    if(j.pt()>0.){
    double NHF = j.neutralHadronEnergyFraction();
    double NEMF = j.neutralEmEnergyFraction();
    double CHF = j.chargedHadronEnergyFraction();
    //double MUF = j.muonEnergyFraction();
    int NumConst = j.chargedMultiplicity()+j.neutralMultiplicity();
    int NumNeutralParticle =j.neutralMultiplicity();
    int CHM = j.chargedMultiplicity();
    double eta = j.eta();
    return ((  (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 ) || (abs(eta)>2.4 && abs(eta)<=2.7) )) || (NHF<0.99 && abs(eta)>2.7 && abs(eta)<=3.0 ) || (NEMF<0.90 && NHF>0.02 &&NumNeutralParticle>2 && NumNeutralParticle<15 && abs(eta)>3.0) ) ;
}
else{
return (0);
    }
}

float
VVVTreeMaker::dEtaInSeed( const pat::Electron*  ele ){
    return ele->superCluster().isNonnull() && ele->superCluster()->seed().isNonnull() ? ele->deltaEtaSuperClusterTrackAtVtx() - ele->superCluster()->eta() + ele->superCluster()->seed()->eta() : std::numeric_limits<float>::max();

}


void VVVTreeMaker::initJetCorrFactors( void ){
    std::vector<JetCorrectorParameters> vPar;
    for ( std::vector<std::string>::const_iterator payloadBegin = jecAK8PayloadNames_.begin(), payloadEnd = jecAK8PayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
        JetCorrectorParameters pars(*ipayload);
        vPar.push_back(pars);
    }
    // Make the FactorizedJetCorrector
    jecAK8_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );

    vPar.clear();
    for ( std::vector<std::string>::const_iterator payloadBegin = jecAK8PayloadNamesGroomed_.begin(), payloadEnd = jecAK8PayloadNamesGroomed_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
        JetCorrectorParameters pars(*ipayload);
        vPar.push_back(pars);
    }
    // Make the FactorizedJetCorrector
    jecAK8Groomed_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
  
    vPar.clear();
    for ( std::vector<std::string>::const_iterator payloadBegin = jecAK8PayloadNamesGroomed_.begin(), payloadEnd = jecAK8PayloadNamesGroomed_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
        JetCorrectorParameters pars(*ipayload);
        vPar.push_back(pars);
    }
    jecAK8GroomedSD_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );

    vPar.clear();
    for ( std::vector<std::string>::const_iterator payloadBegin = jecAK8puppiPayloadNames_.begin(), payloadEnd = jecAK8puppiPayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
        JetCorrectorParameters pars(*ipayload);
        vPar.push_back(pars);
    }
    jecAK8puppi_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );

    vPar.clear();
    for ( std::vector<std::string>::const_iterator payloadBegin = jecAK8puppiPayloadNamesGroomed_.begin(), payloadEnd = jecAK8puppiPayloadNamesGroomed_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
        JetCorrectorParameters pars(*ipayload);
        vPar.push_back(pars);
    }
    jecAK8puppiGroomed_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );

    vPar.clear();
    for ( std::vector<std::string>::const_iterator payloadBegin = jecAK4PayloadNames_.begin(), payloadEnd = jecAK4PayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
        JetCorrectorParameters pars(*ipayload);
        vPar.push_back(pars);
    }
    // Make the FactorizedJetCorrector
    jecAK4_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );

    vPar.clear();
    for ( std::vector<std::string>::const_iterator payloadBegin = offsetCorrLabel_.begin(), payloadEnd = offsetCorrLabel_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
        JetCorrectorParameters pars(*ipayload);
        vPar.push_back(pars);
    }
    jecOffset_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
}


double VVVTreeMaker::getJEC( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ ){

    double jetCorrFactor = 1.;
    if ( fabs(rawJetP4.eta()) < jetCorrEtaMax ){
        jecAK4_->setJetEta( rawJetP4.eta() );
        jecAK4_->setJetPt ( rawJetP4.pt() );
        jecAK4_->setJetE  ( rawJetP4.energy() );
        jecAK4_->setJetPhi( rawJetP4.phi()    );
        jecAK4_->setJetA  ( jet.jetArea() );
        jecAK4_->setRho   ( *(rho_.product()) );
        jecAK4_->setNPV   ( nVtx );
        jetCorrFactor = jecAK4_->getCorrection();
    }
    reco::Candidate::LorentzVector corrJetP4 = rawJetP4;
    corrJetP4 *= jetCorrFactor;
    return jetCorrFactor;
}

double VVVTreeMaker::getJECOffset( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ ){

    double jetCorrFactor = 1.;
    if ( fabs(rawJetP4.eta()) < jetCorrEtaMax ){
        jecOffset_->setJetEta( rawJetP4.eta()     );
        jecOffset_->setJetPt ( rawJetP4.pt()      );
        jecOffset_->setJetE  ( rawJetP4.energy()  );
        jecOffset_->setJetPhi( rawJetP4.phi()     );
        jecOffset_->setJetA  ( jet.jetArea()      );
        jecOffset_->setRho   ( *(rho_.product())  );
        jecOffset_->setNPV   ( nVtx  );
        jetCorrFactor = jecOffset_->getCorrection();
    }

    reco::Candidate::LorentzVector corrJetP4 = rawJetP4;
    corrJetP4 *= jetCorrFactor;

    return jetCorrFactor;
}

//-------------------------------------------------------------------------------------------------------------------------------------//
//
// member functions
//
void VVVTreeMaker::addTypeICorr( edm::Event const & event ){
    TypeICorrMap_.clear();
    event.getByToken(jetToken_      , jets_    );
    event.getByToken(rhoToken_      , rho_     );
    //edm::Handle<double> rho_;
    //event.getByLabel("fixedGridRhoFastjetAll",rho_);
    //edm::Handle<reco::VertexCollection> vertices_;
    //event.getByLabel("offlineSlimmedPrimaryVertices", vertices_);
    //event.getByToken(vtxToken_, vertices_);
    edm::Handle<reco::VertexCollection> vertices_;
    event.getByToken(vtxToken_, vertices_);

    //event.getByToken(muonToken_     , muons_   );
    edm::Handle<edm::View<pat::Muon>> muons_;
    //event.getByLabel("slimmedMuons",muons_);
    event.getByToken(t1muSrc_,muons_);
    bool skipEM_                    = true;
    double skipEMfractionThreshold_ = 0.9;
    bool skipMuons_                 = true;
    
    std::string skipMuonSelection_string = "isGlobalMuon | isStandAloneMuon";
    StringCutObjectSelector<reco::Candidate>* skipMuonSelection_ = new StringCutObjectSelector<reco::Candidate>(skipMuonSelection_string,true);

    double jetCorrEtaMax_           = 9.9;
    double type1JetPtThreshold_     = 15.0; //10.0;

    double corrEx    = 0;
    double corrEy    = 0;
    double corrSumEt = 0;

    for (const pat::Jet &jet : *jets_) {

        double emEnergyFraction = jet.chargedEmEnergyFraction() + jet.neutralEmEnergyFraction();
        if ( skipEM_ && emEnergyFraction > skipEMfractionThreshold_ ) continue;

        reco::Candidate::LorentzVector rawJetP4 = jet.correctedP4(0);
        double corr = getJEC(rawJetP4, jet, jetCorrEtaMax_, jetCorrLabel_);

        if ( skipMuons_ ) {
            const std::vector<reco::CandidatePtr> & cands = jet.daughterPtrVector();
            for ( std::vector<reco::CandidatePtr>::const_iterator cand = cands.begin();cand != cands.end(); ++cand ) {
                const reco::PFCandidate *pfcand = dynamic_cast<const reco::PFCandidate *>(cand->get());
                const reco::Candidate *mu = (pfcand != 0 ? ( pfcand->muonRef().isNonnull() ? pfcand->muonRef().get() : 0) : cand->get());
                if ( mu != 0 && (*skipMuonSelection_)(*mu) ) {
                    reco::Candidate::LorentzVector muonP4 = (*cand)->p4();
                    rawJetP4 -= muonP4;
                }
            }
        }

        reco::Candidate::LorentzVector corrJetP4 = corr*rawJetP4;

        if ( corrJetP4.pt() > type1JetPtThreshold_ ) {
            reco::Candidate::LorentzVector tmpP4 = jet.correctedP4(0);
            corr = getJECOffset(tmpP4, jet, jetCorrEtaMax_, offsetCorrLabel_);
            reco::Candidate::LorentzVector rawJetP4offsetCorr = corr*rawJetP4;

            corrEx    -= (corrJetP4.px() - rawJetP4offsetCorr.px());
            corrEy    -= (corrJetP4.py() - rawJetP4offsetCorr.py());
            corrSumEt += (corrJetP4.Et() - rawJetP4offsetCorr.Et());
        }
    }
    TypeICorrMap_["corrEx"]    = corrEx;
    TypeICorrMap_["corrEy"]    = corrEy;
    TypeICorrMap_["corrSumEt"] = corrSumEt;
}
void VVVTreeMaker::addTypeICorr_user(edm::Event const& event) {
    TypeICorrMap_user_.clear();
    edm::Handle<pat::JetCollection> jets_;
    event.getByToken(t1jetSrc_userak4_, jets_);
    double corrEx_JEC         = 0;
    double corrEy_JEC         = 0;
    double corrSumEt_JEC      = 0;
    double corrEx_JEC_up      = 0;
    double corrEy_JEC_up      = 0;
    double corrSumEt_JEC_up   = 0;
    double corrEx_JEC_down    = 0;
    double corrEy_JEC_down    = 0;
    double corrSumEt_JEC_down = 0;
    
    double corrEx_JER         = 0;
    double corrEy_JER         = 0;
    double corrSumEt_JER      = 0;
    double corrEx_JER_up      = 0;
    double corrEy_JER_up      = 0;
    double corrSumEt_JER_up   = 0;
    double corrEx_JER_down    = 0;
    double corrEy_JER_down    = 0;
    double corrSumEt_JER_down = 0;
    for (const pat::Jet& jet : *jets_) {
        corrEx_JEC += jet.userFloat("corrEx_MET_JEC");
        corrEy_JEC += jet.userFloat("corrEy_MET_JEC");
        corrSumEt_JEC += jet.userFloat("corrSumEt_MET_JEC");
        corrEx_JEC_up += jet.userFloat("corrEx_MET_JEC_up");
        corrEy_JEC_up += jet.userFloat("corrEy_MET_JEC_up");
        corrSumEt_JEC_up += jet.userFloat("corrSumEt_MET_JEC_up");
        corrEx_JEC_down += jet.userFloat("corrEx_MET_JEC_down");
        corrEy_JEC_down += jet.userFloat("corrEy_MET_JEC_down");
        corrSumEt_JEC_down += jet.userFloat("corrSumEt_MET_JEC_down");
        corrEx_JER += jet.userFloat("corrEx_MET_JER");
        corrEy_JER += jet.userFloat("corrEy_MET_JER");
        corrSumEt_JER += jet.userFloat("corrSumEt_MET_JER");
        corrEx_JER_up += jet.userFloat("corrEx_MET_JER_up");
        corrEy_JER_up += jet.userFloat("corrEy_MET_JER_up");
        corrSumEt_JER_up += jet.userFloat("corrSumEt_MET_JER_up");
        corrEx_JER_down += jet.userFloat("corrEx_MET_JER_down");
        corrEy_JER_down += jet.userFloat("corrEy_MET_JER_down");
        corrSumEt_JER_down += jet.userFloat("corrSumEt_MET_JER_down");
    }
    TypeICorrMap_user_["corrEx_JEC"]         = corrEx_JEC;
    TypeICorrMap_user_["corrEy_JEC"]         = corrEy_JEC;
    TypeICorrMap_user_["corrSumEt_JEC"]      = corrSumEt_JEC;
    TypeICorrMap_user_["corrEx_JEC_up"]      = corrEx_JEC_up;
    TypeICorrMap_user_["corrEy_JEC_up"]      = corrEy_JEC_up;
    TypeICorrMap_user_["corrSumEt_JEC_up"]   = corrSumEt_JEC_up;
    TypeICorrMap_user_["corrEx_JEC_down"]    = corrEx_JEC_down;
    TypeICorrMap_user_["corrEy_JEC_down"]    = corrEy_JEC_down;
    TypeICorrMap_user_["corrSumEt_JEC_down"] = corrSumEt_JEC_down;
    
    TypeICorrMap_user_["corrEx_JER"]         = corrEx_JER;
    TypeICorrMap_user_["corrEy_JER"]         = corrEy_JER;
    TypeICorrMap_user_["corrSumEt_JER"]      = corrSumEt_JER;
    TypeICorrMap_user_["corrEx_JER_up"]      = corrEx_JER_up;
    TypeICorrMap_user_["corrEy_JER_up"]      = corrEy_JER_up;
    TypeICorrMap_user_["corrSumEt_JER_up"]   = corrSumEt_JER_up;
    TypeICorrMap_user_["corrEx_JER_down"]    = corrEx_JER_down;
    TypeICorrMap_user_["corrEy_JER_down"]    = corrEy_JER_down;
    TypeICorrMap_user_["corrSumEt_JER_down"] = corrSumEt_JER_down;
}


//-------------------------------------------------------------------------------------------------------------------------------------//
math::XYZTLorentzVector
VVVTreeMaker::getNeutrinoP4(double& MetPt, double& MetPhi, TLorentzVector& lep, int lepType){
    double leppt = lep.Pt();
    double lepphi = lep.Phi();
    double lepeta = lep.Eta();
    double lepenergy = lep.Energy();
    
    double metpt = MetPt;
    double metphi = MetPhi;
    
    double  px = metpt*cos(metphi);
    double  py = metpt*sin(metphi);
    double  pz = 0;
    double  pxl= leppt*cos(lepphi);
    double  pyl= leppt*sin(lepphi);
    double  pzl= leppt*sinh(lepeta);
    double  El = lepenergy;
    double  a = pow(MW_,2) + pow(px+pxl,2) + pow(py+pyl,2) - px*px - py*py - El*El + pzl*pzl;
    double  b = 2.*pzl;
    double  A = b*b -4.*El*El;
    double  B = 2.*a*b;
    double  C = a*a-4.*(px*px+py*py)*El*El;
    
    ///////////////////////////pz for fnal
    double M_mu =  0;
    
    //if(lepType==1)M_mu=0.105658367;//mu
    //if(lepType==0)M_mu=0.00051099891;//electron
    
    int type=2; // use the small abs real root
    
    a = MW_*MW_ - M_mu*M_mu + 2.0*pxl*px + 2.0*pyl*py;
    A = 4.0*(El*El - pzl*pzl);
    B = -4.0*a*pzl;
    C = 4.0*El*El*(px*px + py*py) - a*a;
    
    double tmproot = B*B - 4.0*A*C;
    
    if (tmproot<0) {
        //std::cout << "Complex root detected, taking real part..." << std::endl;
        pz = - B/(2*A); // take real part of complex roots
    }
    else {
        double tmpsol1 = (-B + sqrt(tmproot))/(2.0*A);
        double tmpsol2 = (-B - sqrt(tmproot))/(2.0*A);
        //std::cout << " Neutrino Solutions: " << tmpsol1 << ", " << tmpsol2 << std::endl;
        
        if (type == 0 ) {
            // two real roots, pick the one closest to pz of muon
            if (TMath::Abs(tmpsol2-pzl) < TMath::Abs(tmpsol1-pzl)) { pz = tmpsol2; }
            else { pz = tmpsol1; }
            // if pz is > 300 pick the most central root
            if ( abs(pz) > 300. ) {
                if (TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2) ) { pz = tmpsol1; }
                else { pz = tmpsol2; }
            }
        }
        if (type == 1 ) {
            // two real roots, pick the one closest to pz of muon
            if (TMath::Abs(tmpsol2-pzl) < TMath::Abs(tmpsol1-pzl)) { pz = tmpsol2; }
            else {pz = tmpsol1; }
        }
        if (type == 2 ) {
            // pick the most central root.
            if (TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2) ) { pz = tmpsol1; }
            else { pz = tmpsol2; }
        }
        /*if (type == 3 ) {
         // pick the largest value of the cosine
         TVector3 p3w, p3mu;
         p3w.SetXYZ(pxl+px, pyl+py, pzl+ tmpsol1);
         p3mu.SetXYZ(pxl, pyl, pzl );
         
         double sinthcm1 = 2.*(p3mu.Perp(p3w))/MW_;
         p3w.SetXYZ(pxl+px, pyl+py, pzl+ tmpsol2);
         double sinthcm2 = 2.*(p3mu.Perp(p3w))/MW_;
         
         double costhcm1 = sqrt(1. - sinthcm1*sinthcm1);
         double costhcm2 = sqrt(1. - sinthcm2*sinthcm2);
         
         if ( costhcm1 > costhcm2 ) { pz = tmpsol1; otherSol_ = tmpsol2; }
         else { pz = tmpsol2;otherSol_ = tmpsol1; }
         
         }*///end of type3
        
    }//endl of if real root
    
    //dont correct pt neutrino
    math::XYZTLorentzVector outP4(px,py,pz,sqrt(px*px+py*py+pz*pz));
    return outP4;
    
}//end neutrinoP4



// ------------ method called for each event  ------------
void
VVVTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   
    using namespace edm;
    setDummyValues(); //Initalize variables with dummy values

    nevent = iEvent.eventAuxiliary().event();
    run    = iEvent.eventAuxiliary().run();
    ls     = iEvent.eventAuxiliary().luminosityBlock();

    HLTStore(iEvent);
    GENStore(iEvent);


    edm::Handle<edm::View<pat::Jet> > hadronicVs;
    //iEvent.getByLabel(hadronicVSrc_.c_str(), hadronicVs);
    iEvent.getByToken(hadronicVSrc_, hadronicVs);
   
    edm::Handle<edm::View<pat::Jet> > hadronicVs_raw;
    //iEvent.getByLabel("slimmedJetsAK8", hadronicVs_raw);
    iEvent.getByToken(hadronicVSrc_raw_, hadronicVs_raw);

    edm::Handle<pat::JetCollection> hadronicVSoftDrop;
    iEvent.getByToken(hadronicVSoftDropSrc_, hadronicVSoftDrop);

    edm::Handle<edm::View<reco::Candidate> > leptonicVs;
    //iEvent.getByLabel(leptonicVSrc_.c_str(), leptonicVs);
    iEvent.getByToken(leptonicVSrc_, leptonicVs);

    edm::Handle<double> rho;
    //iEvent.getByLabel("fixedGridRhoFastjetAll",rho);

    iEvent.getByToken(rhoToken_      , rho     );
    double fastJetRho = *(rho.product());
    useless = fastJetRho;

    edm::Handle<edm::View<pat::Jet> > ak4jets;
    //iEvent.getByLabel(ak4jetsSrc_.c_str(), ak4jets);
    iEvent.getByToken(ak4jetsSrc_, ak4jets);
    
    edm::Handle<edm::View<reco::Candidate> > gravitons;
    //iEvent.getByLabel(gravitonSrc_.c_str(), gravitons);
    iEvent.getByToken(gravitonSrc_, gravitons);
    edm::Handle<edm::View<reco::Candidate> > metHandle;
    //iEvent.getByLabel(metSrc_.c_str(), metHandle);
    iEvent.getByToken(metSrc_, metHandle);
  
    edm::Handle<edm::View<pat::Muon>> loosemus;
    //iEvent.getByLabel(looseMuonSrc_.c_str(), loosemus);
    iEvent.getByToken(loosemuonToken_,loosemus);

    edm::Handle<edm::View<pat::Muon>> goodmus;
    //iEvent.getByLabel(goodMuSrc_.c_str(), goodmus);
    iEvent.getByToken(goodMuSrc_, goodmus);

    edm::Handle<edm::View<pat::Electron>> looseels;
    //iEvent.getByLabel(looseElectronsSrc_.c_str(), looseels);
    iEvent.getByToken(looseelectronToken_, looseels);

    edm::Handle<edm::View<reco::GenParticle> > genParticles;//define genParticle
    //iEvent.getByLabel(InputTag("prunedGenParticles"), genParticles);
    iEvent.getByToken(genSrc_, genParticles);

    edm::Handle<edm::View<pat::Muon>> mus;
    //iEvent.getByLabel("slimmedMuons",mus);
    iEvent.getByToken(MuSrc_, mus);
    edm::Handle<edm::View<pat::Electron>> eles;
    //iEvent.getByLabel("slimmedElectrons",eles);
    iEvent.getByToken(EleSrc_, eles);
    if (RunOnSig_||RunOnMC_){
        //  L1 prefiring
        edm::Handle< double > theprefweight;
        iEvent.getByToken(prefweight_token, theprefweight ) ;
        L1prefiring =(*theprefweight);
        edm::Handle< double > theprefweightup;
        iEvent.getByToken(prefweightup_token, theprefweightup ) ;
        L1prefiringup =(*theprefweightup);
        
        edm::Handle< double > theprefweightdown;
        iEvent.getByToken(prefweightdown_token, theprefweightdown ) ;
        L1prefiringdown =(*theprefweightdown);
        
        /*edm::Handle<LHEEventProduct> wgtsource;
        iEvent.getByToken(LheToken_, wgtsource);
        //std::cout<<"weight number "<<wgtsource->weights().size()<<std::endl;
        for ( int i=0; i<882; ++i) {
            pweight[i]= wgtsource->weights()[i].wgt/wgtsource->originalXWGTUP();
            //cout<<wgtsource->weights()[i].id<<"    "<<pweight[i]<<endl;
        }*/
	/*
        for ( int i=9; i<110; ++i) {
            pweight[i]= wgtsource->weights()[i+101].wgt/wgtsource->originalXWGTUP();
            //cout<<wgtsource->weights()[i].id<<"    "<<pweight[i]<<endl;
        }*/

        /*
        edm::Handle<LHERunInfoProduct> run;
        typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
        
        iEvent.getByToken(LhestrToken_,run);
        LHERunInfoProduct myLHERunInfoProduct = *(run.product());
        
        for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
            std::cout << iter->tag() << std::endl;
            std::vector<std::string> lines = iter->lines();
            for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
                std::cout << lines.at(iLine);
            }
        }
         */
        //iEvent.getByLabel("externalLHEProducer", wgtsource);
        //iEvent.getByLabel("source", wgtsource);

        edm::Handle<GenEventInfoProduct> genEvtInfo;
        //iEvent.getByLabel( "generator", genEvtInfo );
        iEvent.getByToken(GenToken_,genEvtInfo);

        //const std::vector<double>& evtWeights = genEvtInfo->weights();
        theWeight = genEvtInfo->weight();
        if(theWeight>0) nump = nump+1;
        if(theWeight<0) numm = numm+1;
	//cout<<theWeight<<endl;
        edm::Handle<std::vector<PileupSummaryInfo>>  PupInfo;
        //iEvent.getByLabel(edm::InputTag("slimmedAddPileupInfo"), PupInfo);
        iEvent.getByToken(PUToken_, PupInfo);
        std::vector<PileupSummaryInfo>::const_iterator PVI;
        for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
            nBX = PVI->getBunchCrossing();
            if(nBX == 0) { // "0" is the in-time crossing, negative values are the early crossings, positive are late
                npT = PVI->getTrueNumInteractions();
                npIT = PVI->getPU_NumInteractions();
            }
        }
    }
    //cout << "npT" << npT << " nBX" << nBX << endl;

    //filter
    iEvent.getByToken(noiseFilterToken_, noiseFilterBits_);
    const edm::TriggerNames &names = iEvent.triggerNames(*noiseFilterBits_);
    for (unsigned int i = 0, n = noiseFilterBits_->size(); i < n; ++i) {
        if (names.triggerName(i) == HBHENoiseFilter_Selector_)
            passFilter_HBHE_ = noiseFilterBits_->accept(i); // TO BE USED
        if (names.triggerName(i) == HBHENoiseIsoFilter_Selector_)
            passFilter_HBHEIso_ = noiseFilterBits_->accept(i); // TO BE USED
        if (names.triggerName(i) == GlobalHaloNoiseFilter_Selector_)
            passFilter_GlobalHalo_ = noiseFilterBits_->accept(i); // TO BE USED
        if (names.triggerName(i) == ECALDeadCellNoiseFilter_Selector_)
            passFilter_ECALDeadCell_ = noiseFilterBits_->accept(i); // under scrutiny
        if (names.triggerName(i) == GoodVtxNoiseFilter_Selector_)
            passFilter_GoodVtx_ = noiseFilterBits_->accept(i); // TO BE USED
        if (names.triggerName(i) == EEBadScNoiseFilter_Selector_)
            passFilter_EEBadSc_ = noiseFilterBits_->accept(i); // under scrutiny
    }

    edm::Handle<bool> badMuonResultHandle;
    edm::Handle<bool> badChargedHadronResultHandle;
    iEvent.getByToken(badMuon_Selector_, badMuonResultHandle);
    iEvent.getByToken(badChargedHadron_Selector_, badChargedHadronResultHandle);
    passFilter_badMuon_ = *badMuonResultHandle;
    passFilter_badChargedHadron_ = *badChargedHadronResultHandle;
    edm::Handle< bool > passecalBadCalibFilterUpdate ;
    iEvent.getByToken(ecalBadCalibFilterUpdate_token,passecalBadCalibFilterUpdate);
    passecalBadCalibFilterUpdate_ =  (*passecalBadCalibFilterUpdate );
 
    numCands = gravitons->size();
 





    //if(numCands != 0 ) {
    //const reco::Candidate& graviton  = gravitons->at(0);
    //cout<<hadronicVs->size()<<"  hadronicVs->size() "<<leptonicVs->size()<<endl;
    if((hadronicVs->size()!= 0 )  && (leptonicVs->size()!= 0) ){

       const reco::Candidate& leptonicV = leptonicVs->at(0);
       const reco::Candidate& hadronicV = hadronicVs->at(0);
       // const reco::Candidate& leptonicV = (*graviton.daughter("leptonicV"));
       const reco::Candidate& metCand = metHandle->at(0);
       const reco::Candidate& lepton = (*leptonicV.daughter(0));
       nLooseMu = loosemus->size();
       nLooseEle = looseels->size();

       edm::Handle<reco::VertexCollection> vertices;
       iEvent.getByToken(vtxToken_, vertices);
        // edm::Handle<reco::VertexCollection> vertices;
        // iEvent.getByLabel("offlineSlimmedPrimaryVertices", vertices);
        // iEvent.getByToken(vtxToken_, vertices);
        if (vertices->empty()) return; // skip the event if no PV found
        nVtx = vertices->size();
        reco::VertexCollection::const_iterator firstGoodVertex = vertices->end();
        for (reco::VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx) {
            // Replace isFake() for miniAOD because it requires tracks and miniAOD vertices don't have tracks:
            // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
            if (  /// !vtx->isFake() &&
                !(vtx->chi2()==0 && vtx->ndof()==0)
                &&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0
                && fabs(vtx->position().Z())<=24.0) {
                    firstGoodVertex = vtx;
                    break;
                }
        }
        if ( firstGoodVertex==vertices->end() ) return; // skip event if there are no good PVs
        // ***************************************************************** //
        // ************************* MET ********************** //
        iEvent.getByToken(metInputToken_ , METs_ );
        addTypeICorr(iEvent);
	if (RunOnMC_) addTypeICorr_user(iEvent);
        for (const pat::MET &met : *METs_) {            //         const float  rawPt    = met.shiftedPt(pat::MET::METUncertainty::NoShift, pat::MET::METUncertaintyLevel::Raw);
            //         const float  rawPhi   = met.shiftedPhi(pat::MET::METUncertainty::NoShift, pat::MET::METUncertaintyLevel::Raw);
            //         const float  rawSumEt = met.shiftedSumEt(pat::MET::METUncertainty::NoShift, pat::MET::METUncertaintyLevel::Raw);
            const float rawPt = met.uncorPt();
            const float rawPhi = met.uncorPhi();
            const float rawSumEt = met.uncorSumEt();
            TVector2 rawMET_;
            rawMET_.SetMagPhi (rawPt, rawPhi );
            Double_t rawPx = rawMET_.Px();
            Double_t rawPy = rawMET_.Py();
            Double_t rawEt = std::hypot(rawPx,rawPy);
            METraw_et = rawEt;
            METraw_phi = rawPhi;
            METraw_sumEt = rawSumEt;
            
            double pxcorr = rawPx+TypeICorrMap_["corrEx"];
            double pycorr = rawPy+TypeICorrMap_["corrEy"];
            double et     = std::hypot(pxcorr,pycorr);
            double sumEtcorr = rawSumEt+TypeICorrMap_["corrSumEt"];
            
            TLorentzVector corrmet;

            corrmet.SetPxPyPzE(pxcorr,pycorr,0.,et);
            MET_et = et;
            MET_phi = corrmet.Phi();
            MET_sumEt = sumEtcorr;
            useless = sumEtcorr;
            useless = rawEt;
            MET_corrPx = TypeICorrMap_["corrEx"];
            MET_corrPy = TypeICorrMap_["corrEy"];

            Double_t rawPtc       = met.corPt();
            Double_t rawPhic   = met.corPhi();
            Double_t rawSumEtc = met.corSumEt();
            TVector2 rawMET_c;
            rawMET_c.SetMagPhi (rawPtc, rawPhic );
            Double_t rawPxc = rawMET_c.Px();
            Double_t rawPyc = rawMET_c.Py();
            Double_t rawEtc = std::hypot(rawPxc,rawPyc);
            MET_et_m = rawEtc;
            MET_phi_m = rawPhic;
            MET_sumEt_m = rawSumEtc;
            MET_corrPx = TypeICorrMap_["corrEx"];
            MET_corrPy = TypeICorrMap_["corrEy"];

            if (RunOnMC_){ 
            double pxcorr_newo= rawPx+TypeICorrMap_user_["corrEx_JEC"];
            double pycorr_newo= rawPy+TypeICorrMap_user_["corrEy_JEC"];
            double et_newo     = std::hypot(pxcorr_newo,pycorr_newo);
	    MET_et_old=et_newo;
	    //cout<<MET_corrPx<<" MET_corrPx "<<MET_corrPy<<"  MET_corrPy  "<<TypeICorrMap_user_["corrEx_JEC"]<<"   "<<TypeICorrMap_user_["corrEy_JEC"]<<endl;
            // Marked for debug
            //------------------central value, correction from JetuserDataak4---------------------
            double pxcorr_new= rawPx+TypeICorrMap_user_["corrEx_JEC"]+TypeICorrMap_user_["corrEx_JER"];
            double pycorr_new= rawPy+TypeICorrMap_user_["corrEy_JEC"]+TypeICorrMap_user_["corrEy_JER"];
            double et_new     = std::hypot(pxcorr_new,pycorr_new);
            double sumEtcorr_new = rawSumEt+TypeICorrMap_user_["corrSumEt_JEC"]+TypeICorrMap_user_["corrSumEt_JER"];
            //----for JEC uncertainty study
            double pxcorr_JEC_up = rawPx+TypeICorrMap_user_["corrEx_JEC_up"]+TypeICorrMap_user_["corrEx_JER"];
            double pycorr_JEC_up = rawPy+TypeICorrMap_user_["corrEy_JEC_up"]+TypeICorrMap_user_["corrEy_JER"];
            double et_JEC_up     = std::hypot(pxcorr_JEC_up, pycorr_JEC_up);
            double sumEtcorr_JEC_up = rawSumEt+TypeICorrMap_user_["corrSumEt_JEC_up"]+TypeICorrMap_user_["corrSumEt_JER"];
            double pxcorr_JEC_down = rawPx+TypeICorrMap_user_["corrEx_JEC_down"]+TypeICorrMap_user_["corrEx_JER"];
            double pycorr_JEC_down = rawPy+TypeICorrMap_user_["corrEy_JEC_down"]+TypeICorrMap_user_["corrEy_JER"];
            double et_JEC_down     = std::hypot(pxcorr_JEC_down, pycorr_JEC_down);
            double sumEtcorr_JEC_down = rawSumEt+TypeICorrMap_user_["corrSumEt_JEC_down"]+TypeICorrMap_user_["corrSumEt_JER"];
            //----for JER uncertainty study
            double pxcorr_JER_up = rawPx+TypeICorrMap_user_["corrEx_JEC"]+TypeICorrMap_user_["corrEx_JER_up"];
            double pycorr_JER_up = rawPy+TypeICorrMap_user_["corrEy_JEC"]+TypeICorrMap_user_["corrEy_JER_up"];
            double et_JER_up     = std::hypot(pxcorr_JER_up, pycorr_JER_up);
            double sumEtcorr_JER_up = rawSumEt+TypeICorrMap_user_["corrSumEt_JEC"]+TypeICorrMap_user_["corrSumEt_JER_up"];
            double pxcorr_JER_down = rawPx+TypeICorrMap_user_["corrEx_JEC"]+TypeICorrMap_user_["corrEx_JER_down"];
            double pycorr_JER_down = rawPy+TypeICorrMap_user_["corrEy_JEC"]+TypeICorrMap_user_["corrEy_JER_down"];
            double et_JER_down     = std::hypot(pxcorr_JER_down,pycorr_JER_down);
            double sumEtcorr_JER_down = rawSumEt+TypeICorrMap_user_["corrSumEt_JEC"]+TypeICorrMap_user_["corrSumEt_JER_down"];
            //------------ 
            // Marked for debug
            MET_et_new= et_new;
            MET_et_JEC_up = et_JEC_up;
            MET_et_JEC_down = et_JEC_down;
            MET_et_JER_up = et_JER_up;
            MET_et_JER_down = et_JER_down;
            
            corrmet.SetPxPyPzE(pxcorr_new,pycorr_new,0.,et_new);
            MET_phi_new = corrmet.Phi();
            corrmet.SetPxPyPzE(pxcorr_JEC_up,pycorr_JEC_up,0.,et_JEC_up);
            MET_phi_JEC_up = corrmet.Phi();
            corrmet.SetPxPyPzE(pxcorr_JEC_down,pycorr_JEC_down,0.,et_JEC_down);
            MET_phi_JEC_down = corrmet.Phi();
            corrmet.SetPxPyPzE(pxcorr_JER_up,pycorr_JER_up,0.,et_JER_up);
            MET_phi_JER_up = corrmet.Phi();
            corrmet.SetPxPyPzE(pxcorr_JER_down,pycorr_JER_down,0.,et_JER_down);
            MET_phi_JER_down = corrmet.Phi();
            
            MET_sumEt_new = sumEtcorr_new;
            MET_sumEt_JEC_up = sumEtcorr_JEC_up;
            MET_sumEt_JEC_down = sumEtcorr_JEC_down;
            MET_sumEt_JER_up = sumEtcorr_JER_up;
            MET_sumEt_JER_down = sumEtcorr_JER_down;
            }// Marked for debug
            
        }
        // ***************************************************************** //
        
        /// For the time being, set these to 1
        triggerWeight=1.0;
        pileupWeight=1.0;
        
        double targetEvents = targetLumiInvPb_*crossSectionPb_;
        lumiWeight = targetEvents/originalNEvents_;

        ptlep1       = leptonicV.daughter(0)->pt();
        ptlep2       = leptonicV.daughter(1)->pt();
        etalep1      = leptonicV.daughter(0)->eta();
        etalep2      = leptonicV.daughter(1)->eta();
        philep1      = leptonicV.daughter(0)->phi();
        philep2      = leptonicV.daughter(1)->phi();
        lep          = std::max(abs(leptonicV.daughter(0)->pdgId()), abs(leptonicV.daughter(1)->pdgId()));
        double energylep1     = leptonicV.daughter(0)->energy();
        met          = metCand.pt();
        metPhi         = metCand.phi();
        //cout<<met<<" met candidate "<<metPhi<<endl;
        //candMass     = graviton.mass();
        ptVlep       = leptonicV.pt();
        yVlep        = leptonicV.eta();
        phiVlep      = leptonicV.phi();
        massVlep     = leptonicV.mass();
        mtVlep       = leptonicV.mt();
        TLorentzVector g_graviton, g_vhad, g_vlep;
        g_vlep.SetPtEtaPhiM(leptonicV.pt(),leptonicV.eta(),leptonicV.phi(),leptonicV.mass());
        g_vhad.SetPtEtaPhiM(hadronicV.pt(),hadronicV.eta(),hadronicV.phi(),hadronicV.mass());
        g_graviton = g_vlep + g_vhad ;
        candMass = g_graviton.Mag();

        ////////////////////////lep ID  ////////////////////////////////////
        if( leptonicV.daughter(0)->isMuon()||leptonicV.daughter(1)->isMuon()){

            const pat::Muon *mu1 = abs(leptonicV.daughter(0)->pdgId())==13 ?
                                                  (pat::Muon*)leptonicV.daughter(0):
                                                  (pat::Muon*)leptonicV.daughter(1);
            isHighPt = mu1->isHighPtMuon(vertices->at(0));
            trackIso = mu1->trackIso();
            muchaiso=mu1->pfIsolationR04().sumChargedHadronPt;
            muneuiso=mu1->pfIsolationR04().sumNeutralHadronEt;
            muphoiso=mu1->pfIsolationR04().sumPhotonEt;
            muPU=mu1->pfIsolationR04().sumPUPt;
            muisolation = (muchaiso+ std::max(0.0,muneuiso+muphoiso-0.5*muPU))/mu1->pt();

        }
        if( leptonicV.daughter(0)->isElectron()||leptonicV.daughter(1)->isElectron() ) {
            const pat::Electron *el1 = leptonicV.daughter(0)->isElectron() ?
                                                  (pat::Electron*)leptonicV.daughter(0):
                                                  (pat::Electron*)leptonicV.daughter(1);
            double etaSC1         = el1->superCluster()->eta();
            double d01            = (-1)*el1->gsfTrack()->dxy(firstGoodVertex->position());
            isHEEP = false;
            et = el1->energy()!=0. ? el1->et()/el1->energy()*el1->caloEnergy() : 0.;
            if( et > 35. ) {
                if( fabs(etaSC1) < 1.4442 ){
                    iso = el1->dr03EcalRecHitSumEt() + el1->dr03HcalDepth1TowerSumEt();
                    isoCut = 2 + 0.03*et + 0.28*fastJetRho;
                    if( el1->ecalDriven() == 1 && dEtaInSeed( el1 ) < 0.004 && el1->deltaPhiSuperClusterTrackAtVtx() < 0.06 &&
                        el1->hadronicOverEm() < (1./el1->superCluster()->energy()+0.05) &&
                        (el1->full5x5_e2x5Max()/el1->full5x5_e5x5() > 0.94 || el1->full5x5_e1x5()/el1->full5x5_e5x5() > 0.83) &&
                        el1->dr03TkSumPt() < 5. && el1->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 &&//numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 &&
                        iso < isoCut && fabs(d01) < 0.02 ) isHEEP = true;
                    }
                    if( fabs(etaSC1) > 1.566 && fabs(etaSC1) < 2.5 ){
                        iso = el1->dr03EcalRecHitSumEt() + el1->dr03HcalDepth1TowerSumEt();
                        if( et <= 50 )
                            isoCut = 2.5 + 0.28*fastJetRho;
                        else
                            isoCut = 2.5+0.03*(et-50.) + 0.28*fastJetRho;
                        if( el1->ecalDriven() == 1 && dEtaInSeed( el1 ) < 0.006 && el1->deltaPhiSuperClusterTrackAtVtx() < 0.06 &&
                            el1->hadronicOverEm() < (5./el1->superCluster()->energy()+0.05) && el1->full5x5_sigmaIetaIeta() < 0.03 &&
                            el1->dr03TkSumPt() < 5. && el1->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 &&//numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 &&
                            iso < isoCut && fabs(d01) < 0.05 ) isHEEP = true;
                    }
            }
        }

        /////////////////////////Leptonic Part///////////////////////////////
        TLorentzVector  glepton, gleptonicV,gleptonicV_new,gleptonicV_JEC_up,gleptonicV_JEC_down,gleptonicV_JER_up,gleptonicV_JER_down;
        glepton.SetPtEtaPhiE(ptlep1, etalep1, philep1, energylep1);
        math::XYZTLorentzVector neutrinoP4 = getNeutrinoP4(MET_et, MET_phi, glepton, 1);
        reco::CandidateBaseRef METBaseRef = metHandle->refAt(0);
        reco::ShallowCloneCandidate neutrino(METBaseRef, 0 , neutrinoP4);
        reco::CompositeCandidate WLeptonic;
        WLeptonic.addDaughter(lepton);
        WLeptonic.addDaughter(neutrino);
        AddFourMomenta addP4;
        addP4.set(WLeptonic);
        gleptonicV.SetPtEtaPhiM(WLeptonic.pt(),WLeptonic.eta(),WLeptonic.phi(),WLeptonic.mass());
        ptVlepJEC       = WLeptonic.pt();
        yVlepJEC        = WLeptonic.eta();
        phiVlepJEC      = WLeptonic.phi();
        massVlepJEC     = WLeptonic.mass();
        if (RunOnMC_){ 
        math::XYZTLorentzVector     neutrinoP4_new = getNeutrinoP4(MET_et_new, MET_phi_new, glepton, 1);
        math::XYZTLorentzVector     neutrinoP4_JEC_up = getNeutrinoP4(MET_et_JEC_up, MET_phi_JEC_up, glepton, 1);
        math::XYZTLorentzVector     neutrinoP4_JEC_down = getNeutrinoP4(MET_et_JEC_down, MET_phi_JEC_down, glepton, 1);
        math::XYZTLorentzVector     neutrinoP4_JER_up = getNeutrinoP4(MET_et_JER_up, MET_phi_JER_up, glepton, 1);
        math::XYZTLorentzVector     neutrinoP4_JER_down = getNeutrinoP4(MET_et_JER_down, MET_phi_JER_down, glepton, 1);
        reco::ShallowCloneCandidate neutrino_new(METBaseRef, 0, neutrinoP4_new);
        reco::ShallowCloneCandidate neutrino_JEC_up(METBaseRef, 0, neutrinoP4_JEC_up);
        reco::ShallowCloneCandidate neutrino_JEC_down(METBaseRef, 0, neutrinoP4_JEC_down);
        reco::ShallowCloneCandidate neutrino_JER_up(METBaseRef, 0, neutrinoP4_JER_up);
        reco::ShallowCloneCandidate neutrino_JER_down(METBaseRef, 0, neutrinoP4_JER_down);
        reco::CompositeCandidate    WLeptonic_new;
        reco::CompositeCandidate    WLeptonic_JEC_up;
        reco::CompositeCandidate    WLeptonic_JEC_down;
        reco::CompositeCandidate    WLeptonic_JER_up;
        reco::CompositeCandidate    WLeptonic_JER_down;
        WLeptonic_new.addDaughter(lepton);
        WLeptonic_new.addDaughter(neutrino_new);
        WLeptonic_JEC_up.addDaughter(lepton);
        WLeptonic_JEC_up.addDaughter(neutrino_JEC_up);
        WLeptonic_JEC_down.addDaughter(lepton);
        WLeptonic_JEC_down.addDaughter(neutrino_JEC_down);
        WLeptonic_JER_up.addDaughter(lepton);
        WLeptonic_JER_up.addDaughter(neutrino_JER_up);
        WLeptonic_JER_down.addDaughter(lepton);
        WLeptonic_JER_down.addDaughter(neutrino_JER_down);
        AddFourMomenta addP4_new;
        addP4_new.set(WLeptonic_new);
        AddFourMomenta addP4_JEC_up;
        addP4_JEC_up.set(WLeptonic_JEC_up);
        AddFourMomenta addP4_JEC_down;
        addP4_JEC_down.set(WLeptonic_JEC_down);
        AddFourMomenta addP4_JER_up;
        addP4_JER_up.set(WLeptonic_JER_up);
        AddFourMomenta addP4_JER_down;
        addP4_JER_down.set(WLeptonic_JER_down);
        gleptonicV_new.SetPtEtaPhiM(WLeptonic_new.pt(),WLeptonic_new.eta(),WLeptonic_new.phi(),WLeptonic_new.mass());
        gleptonicV_JEC_up.SetPtEtaPhiM(WLeptonic_JEC_up.pt(),WLeptonic_JEC_up.eta(),WLeptonic_JEC_up.phi(),WLeptonic_JEC_up.mass());
        gleptonicV_JEC_down.SetPtEtaPhiM(WLeptonic_JEC_down.pt(),WLeptonic_JEC_down.eta(),WLeptonic_JEC_down.phi(),WLeptonic_JEC_down.mass());
        gleptonicV_JER_down.SetPtEtaPhiM(WLeptonic_JER_down.pt(),WLeptonic_JER_down.eta(),WLeptonic_JER_down.phi(),WLeptonic_JER_down.mass());
        gleptonicV_JER_up.SetPtEtaPhiM(WLeptonic_JER_up.pt(),WLeptonic_JER_up.eta(),WLeptonic_JER_up.phi(),WLeptonic_JER_up.mass());
        
        ptVlepJEC_new    = WLeptonic_new.pt();
        yVlepJEC_new     = WLeptonic_new.eta();
        phiVlepJEC_new   = WLeptonic_new.phi();
        massVlepJEC_new  = WLeptonic_new.mass();
        mtVlepJEC_new    = WLeptonic_new.mt();
        //cout<<ptVlep<<" lep W "<<ptVlepJEC<<"   "<<yVlep<<" lep W "<<yVlepJEC<<"   "<<phiVlep<<" lep W "<<phiVlepJEC<<"   "<<massVlep<<" lep W "<<massVlepJEC<<"   "<<endl;
        //cout<<ptVlep<<" lep Wnew "<<ptVlepJEC_new<<"   "<<yVlep<<" lep W "<<yVlepJEC_new<<"   "<<phiVlep<<" lep W "<<phiVlepJEC_new<<"   "<<massVlep<<" lep W "<<massVlepJEC_new<<"   "<<endl;
        
        ptVlepJEC_JEC_up    = WLeptonic_JEC_up.pt();
        yVlepJEC_JEC_up     = WLeptonic_JEC_up.eta();
        phiVlepJEC_JEC_up   = WLeptonic_JEC_up.phi();
        massVlepJEC_JEC_up  = WLeptonic_JEC_up.mass();
        mtVlepJEC_JEC_up    = WLeptonic_JEC_up.mt();
        
        ptVlepJEC_JEC_down    = WLeptonic_JEC_down.pt();
        yVlepJEC_JEC_down     = WLeptonic_JEC_down.eta();
        phiVlepJEC_JEC_down   = WLeptonic_JEC_down.phi();
        massVlepJEC_JEC_down  = WLeptonic_JEC_down.mass();
        mtVlepJEC_JEC_down    = WLeptonic_JEC_down.mt();
        
        ptVlepJEC_JER_up    = WLeptonic_JER_up.pt();
        yVlepJEC_JER_up     = WLeptonic_JER_up.eta();
        phiVlepJEC_JER_up   = WLeptonic_JER_up.phi();
        massVlepJEC_JER_up  = WLeptonic_JER_up.mass();
        mtVlepJEC_JER_up    = WLeptonic_JER_up.mt();
        
        ptVlepJEC_JER_down    = WLeptonic_JER_down.pt();
        yVlepJEC_JER_down     = WLeptonic_JER_down.eta();
        phiVlepJEC_JER_down   = WLeptonic_JER_down.phi();
        massVlepJEC_JER_down  = WLeptonic_JER_down.mass();
        mtVlepJEC_JER_down    = WLeptonic_JER_down.mt();}
        //cout<<"mtVlepJEC"<<mtVlepJEC<<endl;
        ////////////////////////JEC for AK8/////////////////////////////////

        reco::Candidate::LorentzVector uncorrPrunedJet;

        bool doPuppi  = iEvent.getByToken(puppijetInputToken_, puppijets_ );

        if( doPuppi ){//1

            for(size_t ij=0; ij<puppijets_->size()&&ij<4;ij++){
                corr_AK8puppi[ij] = 1;
                corr_AK8puppiSD[ij] = 1;
                const pat::Jet& hadronicVa = puppijets_->at(ij);
                reco::Candidate::LorentzVector uncorrJet;
                if(not isJEC_) doCorrOnTheFly_ = false;
                if( doCorrOnTheFly_ ){
                    uncorrJet = hadronicVa.correctedP4(0);
                    jecAK8puppi_->setJetEta( uncorrJet.eta()          );
                    jecAK8puppi_->setJetPt ( uncorrJet.pt()           );
                    jecAK8puppi_->setJetE  ( uncorrJet.energy()       );
                    jecAK8puppi_->setRho   (fastJetRho);
                    jecAK8puppi_->setNPV   (nVtx);
                    jecAK8puppi_->setJetA  (hadronicVa.jetArea());
                    corr_AK8puppi[ij] = jecAK8puppi_->getCorrection();
                    jecAK8puppiGroomed_->setJetEta( uncorrJet.eta()          );
                    jecAK8puppiGroomed_->setJetPt ( uncorrJet.pt()           );
                    jecAK8puppiGroomed_->setJetE  ( uncorrJet.energy()       );
                    jecAK8puppiGroomed_->setRho   (fastJetRho);
                    jecAK8puppiGroomed_->setNPV   (nVtx);
                    jecAK8puppiGroomed_->setJetA  (hadronicVa.jetArea());
                    corr_AK8puppiSD[ij] = jecAK8puppiGroomed_->getCorrection();
                }
                else{uncorrJet = hadronicVa.p4();}

                if(ij<4){
                    jetAK8puppi_pt1[ij] = corr_AK8puppi[ij]*uncorrJet.pt();
                    jetAK8puppi_mass1[ij] = corr_AK8puppi[ij]*uncorrJet.mass();
                    jetAK8puppi_eta1[ij] = uncorrJet.eta();
                    jetAK8puppi_jec1[ij] = corr_AK8puppi[ij];
                    jetAK8puppiSD_jec1[ij] = corr_AK8puppiSD[ij];
                }
		jetAK8puppi_pt1_m[ij]=hadronicVa.p4().pt();
		//cout<<ij<<"   "<<jetAK8puppi_pt1[ij]<<" PF "<<(*puppijets_)[ij].isPFJet()<<endl;
        	if (RunOnMC_){ 
                jetAK8puppi_pt1_newnew[ij]=(*puppijets_)[ij].userFloat("SmearedPt_JEC_central");
                jetAK8puppi_pt1_new[ij]=(*puppijets_)[ij].userFloat("SmearedPt");
                jetAK8puppi_pt1_JEC_up[ij]=(*puppijets_)[ij].userFloat("SmearedPt_JEC_up");
                jetAK8puppi_pt1_JEC_down[ij]=(*puppijets_)[ij].userFloat("SmearedPt_JEC_down");
                jetAK8puppi_pt1_JER_up[ij]=(*puppijets_)[ij].userFloat("SmearedPt_JER_up");
                jetAK8puppi_pt1_JER_down[ij]=(*puppijets_)[ij].userFloat("SmearedPt_JER_down");
                jetAK8puppi_e1_new[ij]=(*puppijets_)[ij].userFloat("SmearedE");
                jetAK8puppi_e1_JEC_up[ij]=(*puppijets_)[ij].userFloat("SmearedE_JEC_up");
                jetAK8puppi_e1_JEC_down[ij]=(*puppijets_)[ij].userFloat("SmearedE_JEC_down");
                jetAK8puppi_e1_JER_up[ij]=(*puppijets_)[ij].userFloat("SmearedE_JER_up");
                jetAK8puppi_e1_JER_down[ij]=(*puppijets_)[ij].userFloat("SmearedE_JER_down");
                }

            }
            int usenumber3 = -1; double pt_larger=0;
            int numvhad = puppijets_->size();
            for( int inum = 0; inum< numvhad; inum++){
                const pat::Jet& Vpuppi = puppijets_->at(inum);
                if(tightJetIDpuppi(Vpuppi)<1) continue;
                if(jetAK8puppi_pt1[inum] > pt_larger && fabs(jetAK8puppi_eta1[inum])<2.4 && inum<4) {pt_larger = jetAK8puppi_pt1[inum]; usenumber3 = inum; continue;}
            }
            //cout<<"usenumber3"<<usenumber3<<endl;
            if (usenumber3>-1) {//2
                const pat::Jet& hadronicVpuppi = puppijets_->at(usenumber3);
                // DeepAK8
                jetAK8puppi_dnnTop       = hadronicVpuppi.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:TvsQCD"); 
                jetAK8puppi_dnnW         = hadronicVpuppi.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:WvsQCD"); 
                jetAK8puppi_dnnH4q       = hadronicVpuppi.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:H4qvsQCD"); 
                jetAK8puppi_dnnZ         = hadronicVpuppi.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:ZvsQCD"); 
                jetAK8puppi_dnnZbb       = hadronicVpuppi.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:ZbbvsQCD"); 
                jetAK8puppi_dnnHbb       = hadronicVpuppi.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:HbbvsQCD"); 
                jetAK8puppi_dnnqcd       = hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probQCDbb")+hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probQCDcc")+hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probQCDb")+hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probQCDc")+hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probQCDothers"); 
                jetAK8puppi_dnntop       = hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probTbcq")+hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probTbqq"); 
                jetAK8puppi_dnnw         = hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probWcq")+hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probWqq"); 
                jetAK8puppi_dnnz         = hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probZbb")+hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probZcc")+hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probZqq"); 
                jetAK8puppi_dnnzbb       = hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probZbb"); 
                jetAK8puppi_dnnhbb       = hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probHbb"); 
                jetAK8puppi_dnnh4q       = hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probHqqqq"); 
                // Decorrelated DeepAK8
                jetAK8puppi_dnnDecorrTop       = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:TvsQCD"); 
                jetAK8puppi_dnnDecorrW         = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:WvsQCD"); 
                jetAK8puppi_dnnDecorrH4q       = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:H4qvsQCD"); 
                jetAK8puppi_dnnDecorrZ         = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZvsQCD"); 
                jetAK8puppi_dnnDecorrZbb       = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZbbvsQCD"); 
                jetAK8puppi_dnnDecorrHbb       = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:HbbvsQCD"); 
                jetAK8puppi_dnnDecorrqcd       = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDbb")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDcc")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDb")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDc")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDothers"); 
                jetAK8puppi_dnnDecorrbb        = (hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDbb"))/(hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+jetAK8puppi_dnnDecorrqcd); 
                jetAK8puppi_dnnDecorrcc        = (hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDcc"))/(hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+jetAK8puppi_dnnDecorrqcd); 
                jetAK8puppi_dnnDecorrbbnog     = (hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb"))/(hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+jetAK8puppi_dnnDecorrqcd); 
                jetAK8puppi_dnnDecorrccnog     = (hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc"))/(hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+jetAK8puppi_dnnDecorrqcd); 
                jetAK8puppi_dnnDecorrtop       = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probTbcq")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probTbqq"); 
                jetAK8puppi_dnnDecorrw         = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probWcq")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probWqq"); 
                jetAK8puppi_dnnDecorrz         = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZqq"); 
                jetAK8puppi_dnnDecorrzbb       = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb"); 
                jetAK8puppi_dnnDecorrhbb       = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb"); 
                jetAK8puppi_dnnDecorrh4q       = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHqqqq"); 

                 // ------
                jetAK8puppi_ptJEC       = jetAK8puppi_pt1[usenumber3]; // unpruned corrected jet pt
                jetAK8puppi_ptJEC_m       = jetAK8puppi_pt1_m[usenumber3];
        	if (RunOnMC_){ 
                jetAK8puppi_ptJEC_new       = jetAK8puppi_pt1_new[usenumber3];
                jetAK8puppi_ptJEC_newnew       = jetAK8puppi_pt1_newnew[usenumber3];
                jetAK8puppi_ptJEC_JEC_up       = jetAK8puppi_pt1_JEC_up[usenumber3];
                jetAK8puppi_ptJEC_JEC_down       = jetAK8puppi_pt1_JEC_down[usenumber3];
                jetAK8puppi_ptJEC_JER_up       = jetAK8puppi_pt1_JER_up[usenumber3];
                jetAK8puppi_ptJEC_JER_down       = jetAK8puppi_pt1_JER_down[usenumber3];
                jetAK8puppi_e_new       = jetAK8puppi_e1_new[usenumber3];
                jetAK8puppi_e_JEC_up       = jetAK8puppi_e1_JEC_up[usenumber3];
                jetAK8puppi_e_JEC_down       = jetAK8puppi_e1_JEC_down[usenumber3];
                jetAK8puppi_e_JER_up       = jetAK8puppi_e1_JER_up[usenumber3];
                jetAK8puppi_e_JER_down       = jetAK8puppi_e1_JER_down[usenumber3];
		}
                jetAK8puppi_eta     = jetAK8puppi_eta1[usenumber3]; // unpruned (w/o jec) jet eta
                jetAK8puppi_phi      = hadronicVpuppi.phi(); // unpruned (w/o jec) jet phi
                jetAK8puppi_tau1         = hadronicVpuppi.userFloat("NjettinessAK8Puppi:tau1");
                jetAK8puppi_tau2         = hadronicVpuppi.userFloat("NjettinessAK8Puppi:tau2");
                jetAK8puppi_tau3         = hadronicVpuppi.userFloat("NjettinessAK8Puppi:tau3");
                jetAK8puppi_tau21        = jetAK8puppi_tau2/jetAK8puppi_tau1;
                jetAK8puppi_tau4         = hadronicVpuppi.userFloat("NjettinessAK8Puppi:tau3");
                jetAK8puppi_tau42        = jetAK8puppi_tau4/jetAK8puppi_tau2;
                jetAK8puppi_sd       =  hadronicVpuppi.userFloat("ak8PFJetsPuppiSoftDropMass"); // uncorrected pruned mass
                jetAK8puppi_sdJEC  =corr_AK8puppiSD[usenumber3]*jetAK8puppi_sd;
		//cout<<"jetAK8puppi_sd"<<jetAK8puppi_sd<<"  "<<jetAK8puppi_sdJEC<<endl;
                Double_t gencorrect=1.0;
                Double_t recocorrect_0eta1p3=1.0;
                Double_t recocorrect_1p3eta2p5=1.0;
                gencorrect=1.006-1.062*pow(jetAK8puppi_ptJEC*0.08,-1.2);
                recocorrect_0eta1p3=1.093-1.501e-04*jetAK8puppi_ptJEC+3.449e-07*pow(jetAK8puppi_ptJEC,2)-2.681e-10*pow(jetAK8puppi_ptJEC,3)+8.674e-14*pow(jetAK8puppi_ptJEC,4)-1.001e-17*pow(jetAK8puppi_ptJEC,5);
                recocorrect_1p3eta2p5=1.272-5.72e-04*jetAK8puppi_ptJEC+8.37e-07*pow(jetAK8puppi_ptJEC,2)-5.204e-10*pow(jetAK8puppi_ptJEC,3)+1.454e-13*pow(jetAK8puppi_ptJEC,4)-1.504e-17*pow(jetAK8puppi_ptJEC,5);
                if (fabs(jetAK8puppi_eta)<=1.3){jetAK8puppi_sdcorr=jetAK8puppi_sd*gencorrect*recocorrect_0eta1p3;}
                else if (fabs(jetAK8puppi_eta)<2.5 && fabs(jetAK8puppi_eta)>1.3){jetAK8puppi_sdcorr=jetAK8puppi_sd*gencorrect*recocorrect_1p3eta2p5;}
                IDLoose = tightJetIDpuppi(hadronicVpuppi);
                IDTight = tightJetIDpuppi(hadronicVpuppi);
            }
            
            int usenumber2 = -1; double pt_larger2=0;
            for( int inum = 0; inum< numvhad; inum++){
                const pat::Jet& Vpuppi = puppijets_->at(inum);
                if(tightJetIDpuppi(Vpuppi)<1) continue;
                if(jetAK8puppi_pt1[inum] > pt_larger2 && fabs(jetAK8puppi_eta1[inum])<2.4 && inum != usenumber3 && inum<4) {pt_larger2 = jetAK8puppi_pt1[inum]; usenumber2 = inum; continue;}
            }
            //cout<<"usenumber2"<<usenumber2<<endl;
            if(usenumber2>-1)  {
                const pat::Jet& hadronicVpuppi_2 = puppijets_->at(usenumber2);
                // DeepAK8
                jetAK8puppi_dnnTop_2       = hadronicVpuppi_2.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:TvsQCD"); 
                jetAK8puppi_dnnW_2         = hadronicVpuppi_2.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:WvsQCD"); 
                jetAK8puppi_dnnH4q_2       = hadronicVpuppi_2.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:H4qvsQCD"); 
                jetAK8puppi_dnnZ_2         = hadronicVpuppi_2.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:ZvsQCD"); 
                jetAK8puppi_dnnZbb_2       = hadronicVpuppi_2.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:ZbbvsQCD"); 
                jetAK8puppi_dnnHbb_2       = hadronicVpuppi_2.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:HbbvsQCD"); 
                jetAK8puppi_dnnqcd_2       = hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probQCDbb")+hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probQCDcc")+hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probQCDb")+hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probQCDc")+hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probQCDothers"); 
                jetAK8puppi_dnntop_2       = hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probTbcq")+hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probTbqq"); 
                jetAK8puppi_dnnw_2         = hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probWcq")+hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probWqq"); 
                jetAK8puppi_dnnz_2         = hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probZbb")+hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probZcc")+hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probZqq"); 
                jetAK8puppi_dnnzbb_2       = hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probZbb"); 
                jetAK8puppi_dnnhbb_2       = hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probHbb"); 
                jetAK8puppi_dnnh4q_2       = hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probHqqqq"); 
                // Decorrelated DeepAK8
                jetAK8puppi_dnnDecorrTop_2       = hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:TvsQCD"); 
                jetAK8puppi_dnnDecorrW_2         = hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:WvsQCD"); 
                jetAK8puppi_dnnDecorrH4q_2       = hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:H4qvsQCD"); 
                jetAK8puppi_dnnDecorrZ_2         = hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZvsQCD"); 
                jetAK8puppi_dnnDecorrZbb_2       = hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZbbvsQCD"); 
                jetAK8puppi_dnnDecorrHbb_2       = hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:HbbvsQCD"); 
                jetAK8puppi_dnnDecorrqcd_2       = hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDbb")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDcc")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDb")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDc")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDothers"); 
                jetAK8puppi_dnnDecorrbb_2        = (hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDbb"))/(hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+jetAK8puppi_dnnDecorrqcd_2); 
                jetAK8puppi_dnnDecorrcc_2        = (hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDcc"))/(hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+jetAK8puppi_dnnDecorrqcd_2); 
                jetAK8puppi_dnnDecorrbbnog_2     = (hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb"))/(hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+jetAK8puppi_dnnDecorrqcd_2); 
                jetAK8puppi_dnnDecorrccnog_2     = (hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc"))/(hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+jetAK8puppi_dnnDecorrqcd_2); 
                jetAK8puppi_dnnDecorrtop_2       = hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probTbcq")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probTbqq"); 
                jetAK8puppi_dnnDecorrw_2         = hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probWcq")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probWqq"); 
                jetAK8puppi_dnnDecorrz_2         = hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZqq"); 
                jetAK8puppi_dnnDecorrzbb_2       = hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb"); 
                jetAK8puppi_dnnDecorrhbb_2       = hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb"); 
                jetAK8puppi_dnnDecorrh4q_2       = hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHqqqq"); 

                 // ------
                jetAK8puppi_ptJEC_2       = jetAK8puppi_pt1[usenumber2]; // unpruned corrected jet pt
        	if (RunOnMC_){ 
                jetAK8puppi_ptJEC_2_new       = jetAK8puppi_pt1_new[usenumber2];
                jetAK8puppi_ptJEC_2_JEC_up       = jetAK8puppi_pt1_JEC_up[usenumber2];
                jetAK8puppi_ptJEC_2_JEC_down       = jetAK8puppi_pt1_JEC_down[usenumber2];
                jetAK8puppi_ptJEC_2_JER_up       = jetAK8puppi_pt1_JER_up[usenumber2];
                jetAK8puppi_ptJEC_2_JER_down       = jetAK8puppi_pt1_JER_down[usenumber2];
                jetAK8puppi_e_2_new       = jetAK8puppi_e1_new[usenumber2];
                jetAK8puppi_e_2_JEC_up       = jetAK8puppi_e1_JEC_up[usenumber2];
                jetAK8puppi_e_2_JEC_down       = jetAK8puppi_e1_JEC_down[usenumber2];
                jetAK8puppi_e_2_JER_up       = jetAK8puppi_e1_JER_up[usenumber2];
                jetAK8puppi_e_2_JER_down       = jetAK8puppi_e1_JER_down[usenumber2];
		}
                jetAK8puppi_eta_2     = jetAK8puppi_eta1[usenumber2]; // unpruned (w/o jec) jet eta
                jetAK8puppi_phi_2      = hadronicVpuppi_2.phi(); // unpruned (w/o jec) jet phi
                jetAK8puppi_tau1_2         = hadronicVpuppi_2.userFloat("NjettinessAK8Puppi:tau1");
                jetAK8puppi_tau2_2         = hadronicVpuppi_2.userFloat("NjettinessAK8Puppi:tau2");
                jetAK8puppi_tau3_2         = hadronicVpuppi_2.userFloat("NjettinessAK8Puppi:tau3");
                jetAK8puppi_tau21_2        = jetAK8puppi_tau2_2/jetAK8puppi_tau1_2;
                jetAK8puppi_tau4_2         = hadronicVpuppi_2.userFloat("NjettinessAK8Puppi:tau3");
                jetAK8puppi_tau42_2        = jetAK8puppi_tau4_2/jetAK8puppi_tau2_2;
                jetAK8puppi_sd_2       =  hadronicVpuppi_2.userFloat("ak8PFJetsPuppiSoftDropMass"); // uncorrected pruned mass
                jetAK8puppi_sdJEC_2  =corr_AK8puppiSD[usenumber2]*jetAK8puppi_sd_2;
                Double_t gencorrect=1.0;
                Double_t recocorrect_0eta1p3=1.0;
                Double_t recocorrect_1p3eta2p5=1.0;
                gencorrect=1.006-1.062*pow(jetAK8puppi_ptJEC_2*0.08,-1.2);
                recocorrect_0eta1p3=1.093-1.501e-04*jetAK8puppi_ptJEC_2+3.449e-07*pow(jetAK8puppi_ptJEC_2,2)-2.681e-10*pow(jetAK8puppi_ptJEC_2,3)+8.674e-14*pow(jetAK8puppi_ptJEC_2,4)-1.001e-17*pow(jetAK8puppi_ptJEC_2,5);
                    recocorrect_1p3eta2p5=1.272-5.72e-04*jetAK8puppi_ptJEC_2+8.37e-07*pow(jetAK8puppi_ptJEC_2,2)-5.204e-10*pow(jetAK8puppi_ptJEC_2,3)+1.454e-13*pow(jetAK8puppi_ptJEC_2,4)-1.504e-17*pow(jetAK8puppi_ptJEC_2,5);
                if (fabs(jetAK8puppi_eta_2)<=1.3){jetAK8puppi_sdcorr_2=jetAK8puppi_sd_2*gencorrect*recocorrect_0eta1p3;}
                else if (fabs(jetAK8puppi_eta_2)<2.5 && fabs(jetAK8puppi_eta_2)>1.3){jetAK8puppi_sdcorr_2=jetAK8puppi_sd_2*gencorrect*recocorrect_1p3eta2p5;}
                IDLoose_2 = tightJetIDpuppi(hadronicVpuppi_2);
                IDTight_2 = tightJetIDpuppi(hadronicVpuppi_2);
            }

            int usenumber1 = -1; double pt_larger1=0;
            for( int inum = 0; inum< numvhad; inum++){
                const pat::Jet& Vpuppi = puppijets_->at(inum);
                if(tightJetIDpuppi(Vpuppi)<1) continue;
                if(jetAK8puppi_pt1[inum] > pt_larger1 && fabs(jetAK8puppi_eta1[inum])<2.4 && inum != usenumber3 && inum != usenumber2 && inum<4) {pt_larger1 = jetAK8puppi_pt1[inum]; usenumber1 = inum; continue;}
            }
            //cout<<"usenumber1"<<usenumber1<<endl;
            if(usenumber1>-1)  {
                const pat::Jet& hadronicVpuppi_3 = puppijets_->at(usenumber1);
                // DeepAK8
                jetAK8puppi_dnnTop_3       = hadronicVpuppi_3.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:TvsQCD"); 
                jetAK8puppi_dnnW_3         = hadronicVpuppi_3.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:WvsQCD"); 
                jetAK8puppi_dnnH4q_3       = hadronicVpuppi_3.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:H4qvsQCD"); 
                jetAK8puppi_dnnZ_3         = hadronicVpuppi_3.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:ZvsQCD"); 
                jetAK8puppi_dnnZbb_3       = hadronicVpuppi_3.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:ZbbvsQCD"); 
                jetAK8puppi_dnnHbb_3       = hadronicVpuppi_3.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:HbbvsQCD"); 
                jetAK8puppi_dnnqcd_3       = hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probQCDbb")+hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probQCDcc")+hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probQCDb")+hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probQCDc")+hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probQCDothers"); 
                jetAK8puppi_dnntop_3       = hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probTbcq")+hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probTbqq"); 
                jetAK8puppi_dnnw_3         = hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probWcq")+hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probWqq"); 
                jetAK8puppi_dnnz_3         = hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probZbb")+hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probZcc")+hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probZqq"); 
                jetAK8puppi_dnnzbb_3       = hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probZbb"); 
                jetAK8puppi_dnnhbb_3       = hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probHbb"); 
                jetAK8puppi_dnnh4q_3       = hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probHqqqq"); 
                // Decorrelated DeepAK8
                jetAK8puppi_dnnDecorrTop_3       = hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:TvsQCD"); 
                jetAK8puppi_dnnDecorrW_3         = hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:WvsQCD"); 
                jetAK8puppi_dnnDecorrH4q_3       = hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:H4qvsQCD"); 
                jetAK8puppi_dnnDecorrZ_3         = hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZvsQCD"); 
                jetAK8puppi_dnnDecorrZbb_3       = hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZbbvsQCD"); 
                jetAK8puppi_dnnDecorrHbb_3       = hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:HbbvsQCD"); 
                jetAK8puppi_dnnDecorrqcd_3       = hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDbb")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDcc")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDb")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDc")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDothers"); 
                jetAK8puppi_dnnDecorrbb_3        = (hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDbb"))/(hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+jetAK8puppi_dnnDecorrqcd_3); 
                jetAK8puppi_dnnDecorrcc_3        = (hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDcc"))/(hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+jetAK8puppi_dnnDecorrqcd_3); 
                jetAK8puppi_dnnDecorrbbnog_3     = (hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb"))/(hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+jetAK8puppi_dnnDecorrqcd_3); 
                jetAK8puppi_dnnDecorrccnog_3     = (hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc"))/(hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+jetAK8puppi_dnnDecorrqcd_3); 
                jetAK8puppi_dnnDecorrtop_3       = hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probTbcq")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probTbqq"); 
                jetAK8puppi_dnnDecorrw_3         = hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probWcq")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probWqq"); 
                jetAK8puppi_dnnDecorrz_3         = hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZqq"); 
                jetAK8puppi_dnnDecorrzbb_3       = hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb"); 
                jetAK8puppi_dnnDecorrhbb_3       = hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb"); 
                jetAK8puppi_dnnDecorrh4q_3       = hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHqqqq"); 

                 // ------
                jetAK8puppi_ptJEC_3       = jetAK8puppi_pt1[usenumber1]; // unpruned corrected jet pt
        	if (RunOnMC_){ 
                jetAK8puppi_ptJEC_3_new       = jetAK8puppi_pt1_new[usenumber1];
                jetAK8puppi_ptJEC_3_JEC_up       = jetAK8puppi_pt1_JEC_up[usenumber1];
                jetAK8puppi_ptJEC_3_JEC_down       = jetAK8puppi_pt1_JEC_down[usenumber1];
                jetAK8puppi_ptJEC_3_JER_up       = jetAK8puppi_pt1_JER_up[usenumber1];
                jetAK8puppi_ptJEC_3_JER_down       = jetAK8puppi_pt1_JER_down[usenumber1];
                jetAK8puppi_e_3_new       = jetAK8puppi_e1_new[usenumber1];
                jetAK8puppi_e_3_JEC_up       = jetAK8puppi_e1_JEC_up[usenumber1];
                jetAK8puppi_e_3_JEC_down       = jetAK8puppi_e1_JEC_down[usenumber1];
                jetAK8puppi_e_3_JER_up       = jetAK8puppi_e1_JER_up[usenumber1];
                jetAK8puppi_e_3_JER_down       = jetAK8puppi_e1_JER_down[usenumber1];
		}
                jetAK8puppi_eta_3     = jetAK8puppi_eta1[usenumber1]; // unpruned (w/o jec) jet eta
                jetAK8puppi_phi_3      = hadronicVpuppi_3.phi(); // unpruned (w/o jec) jet phi
                jetAK8puppi_tau1_3         = hadronicVpuppi_3.userFloat("NjettinessAK8Puppi:tau1");
                jetAK8puppi_tau2_3         = hadronicVpuppi_3.userFloat("NjettinessAK8Puppi:tau2");
                jetAK8puppi_tau3_3         = hadronicVpuppi_3.userFloat("NjettinessAK8Puppi:tau3");
                jetAK8puppi_tau21_3        = jetAK8puppi_tau2_3/jetAK8puppi_tau1_3;
                jetAK8puppi_tau4_3         = hadronicVpuppi_3.userFloat("NjettinessAK8Puppi:tau3");
                jetAK8puppi_tau42_3        = jetAK8puppi_tau4_3/jetAK8puppi_tau2_3;
                jetAK8puppi_sd_3       =  hadronicVpuppi_3.userFloat("ak8PFJetsPuppiSoftDropMass"); // uncorrected pruned mass
                jetAK8puppi_sdJEC_3  =corr_AK8puppiSD[usenumber1]*jetAK8puppi_sd_3;
                Double_t gencorrect=1.0;
                Double_t recocorrect_0eta1p3=1.0;
                Double_t recocorrect_1p3eta2p5=1.0;
                gencorrect=1.006-1.062*pow(jetAK8puppi_ptJEC_3*0.08,-1.2);
                recocorrect_0eta1p3=1.093-1.501e-04*jetAK8puppi_ptJEC_3+3.449e-07*pow(jetAK8puppi_ptJEC_3,2)-2.681e-10*pow(jetAK8puppi_ptJEC_3,3)+8.674e-14*pow(jetAK8puppi_ptJEC_3,4)-1.001e-17*pow(jetAK8puppi_ptJEC_3,5);
                recocorrect_1p3eta2p5=1.272-5.72e-04*jetAK8puppi_ptJEC_3+8.37e-07*pow(jetAK8puppi_ptJEC_3,2)-5.204e-10*pow(jetAK8puppi_ptJEC_3,3)+1.454e-13*pow(jetAK8puppi_ptJEC_3,4)-1.504e-17*pow(jetAK8puppi_ptJEC_3,5);
                if (fabs(jetAK8puppi_eta_3)<=1.3){jetAK8puppi_sdcorr_3=jetAK8puppi_sd_3*gencorrect*recocorrect_0eta1p3;}
                else if (fabs(jetAK8puppi_eta_3)<2.5 && fabs(jetAK8puppi_eta_3)>1.3){jetAK8puppi_sdcorr_3=jetAK8puppi_sd_3*gencorrect*recocorrect_1p3eta2p5;}
                IDLoose_3 = tightJetIDpuppi(hadronicVpuppi_3);
                IDTight_3 = tightJetIDpuppi(hadronicVpuppi_3);
            }

            int nak4 = 0;
            double tj1=-10.0, tj2=-10.0;
 
            for (size_t ik=0; ik<ak4jets->size();ik++)
            {//3
                double corr = 1;
                reco::Candidate::LorentzVector uncorrJet;
                if( doCorrOnTheFly_ ){
                    uncorrJet = (*ak4jets)[ik].correctedP4(0);
                    jecAK4_->setJetEta( uncorrJet.eta() );
                    jecAK4_->setJetPt ( uncorrJet.pt() );
                    jecAK4_->setJetE ( uncorrJet.energy() );
                    jecAK4_->setRho ( fastJetRho );
                    jecAK4_->setNPV ( vertices->size() );
                    jecAK4_->setJetA ( (*ak4jets)[ik].jetArea() );
                    corr = jecAK4_->getCorrection();
                } else {uncorrJet = (*ak4jets)[ik].p4();}
    
                //if( (corr*uncorrJet.pt())>20 && (fabs((*ak4jets)[ik].eta()) < 5.0) && looseJetID((*ak4jets)[ik])>0 && dtemp>0.8 && nak4<8){
                if( (corr*uncorrJet.pt())>20 && (fabs((*ak4jets)[ik].eta()) < 5.0) && tightJetID((*ak4jets)[ik])>0 && nak4<8){
                    ak4jet_hf[nak4]=(*ak4jets)[ik].hadronFlavour();
                    ak4jet_pf[nak4]=(*ak4jets)[ik].partonFlavour();
                    ak4jet_pt[nak4] =  corr*uncorrJet.pt();
                    ak4jet_pt_uncorr[nak4] =  uncorrJet.pt();
                    ak4jet_eta[nak4] = (*ak4jets)[ik].eta();
                    ak4jet_phi[nak4] = (*ak4jets)[ik].phi();
                    ak4jet_e[nak4] =   corr*uncorrJet.energy();
                    ak4jet_csv[nak4] = (*ak4jets)[ik].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
                    ak4jet_icsv[nak4] = (*ak4jets)[ik].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
                        ak4jet_deepcsvudsg[nak4] = (*ak4jets)[ik].bDiscriminator("pfDeepCSVJetTags:probudsg");
                        ak4jet_deepcsvb[nak4] = (*ak4jets)[ik].bDiscriminator("pfDeepCSVJetTags:probb");
                        ak4jet_deepcsvc[nak4] = (*ak4jets)[ik].bDiscriminator("pfDeepCSVJetTags:probc");
                        ak4jet_deepcsvbb[nak4] = (*ak4jets)[ik].bDiscriminator("pfDeepCSVJetTags:probbb");
                        ak4jet_deepcsvcc[nak4] = (*ak4jets)[ik].bDiscriminator("pfDeepCSVJetTags:probcc");
                    ak4jet_IDLoose[nak4] = tightJetID((*ak4jets)[ik]);
                    ak4jet_IDTight[nak4] = tightJetID((*ak4jets)[ik]);
                    if(ak4jet_pt[nak4]>tj1 ) {
                        if(tj1>tj2) {tj2=tj1; nj2=nj1;}
                        tj1=ak4jet_pt[nak4]; nj1=nak4;
                    }
                    else if(ak4jet_pt[nak4]>tj2){
                        tj2=ak4jet_pt[nak4]; nj2=nak4;}
                    nak4 = nak4 + 1;
                }
            
            }//3
         

            if(nj1>-1 && nj2>-1 && ak4jet_pt[nj1]>30. && ak4jet_pt[nj2]>30.) {
                vbfeta=fabs(ak4jet_eta[nj1]-ak4jet_eta[nj2]);
                TLorentzVector vbfj1, vbfj2;
                vbfj1.SetPtEtaPhiE(ak4jet_pt[nj1], ak4jet_eta[nj1], ak4jet_phi[nj1], ak4jet_e[nj1]);
                vbfj2.SetPtEtaPhiE(ak4jet_pt[nj2], ak4jet_eta[nj2], ak4jet_phi[nj2], ak4jet_e[nj2]);
                vbfmjj=(vbfj1+vbfj2).Mag();
            }

            if(vbfeta>4.0 && vbfmjj>400) {vbftag=1;}
	    

            
            deltaRlepjet = deltaR(etalep1,philep1,jetAK8puppi_eta,jetAK8puppi_phi);
            deltaRlepjet_2 = deltaR(etalep1,philep1,jetAK8puppi_eta_2,jetAK8puppi_phi_2);
            TLorentzVector ghadronicVpuppi, gravitonpuppiJEC,ghadronicVpuppi_2, gravitonpuppiJEC_2;
            ghadronicVpuppi.SetPtEtaPhiM(jetAK8puppi_ptJEC, jetAK8puppi_eta, jetAK8puppi_phi, jetAK8puppi_sd);
            ghadronicVpuppi_2.SetPtEtaPhiM(jetAK8puppi_ptJEC_2, jetAK8puppi_eta_2, jetAK8puppi_phi_2, jetAK8puppi_sd_2);
            gravitonpuppiJEC = gleptonicV + ghadronicVpuppi+ ghadronicVpuppi_2;
            candMasspuppiJEC     = gravitonpuppiJEC.Mag();
            m_jlv     = (gleptonicV + ghadronicVpuppi).Mag();

            TLorentzVector lvw[3];
            if (RunOnMC_){ 
            TLorentzVector ghadronicVpuppi_new, gravitonpuppiJEC_new,ghadronicVpuppi_2_new;
            ghadronicVpuppi_new.SetPtEtaPhiM(jetAK8puppi_ptJEC_new, jetAK8puppi_eta, jetAK8puppi_phi, jetAK8puppi_sd);
            ghadronicVpuppi_2_new.SetPtEtaPhiM(jetAK8puppi_ptJEC_2_new, jetAK8puppi_eta_2, jetAK8puppi_phi_2, jetAK8puppi_sd_2);
            gravitonpuppiJEC_new = gleptonicV_new + ghadronicVpuppi_new+ ghadronicVpuppi_2_new;
            candMasspuppiJEC_new     = gravitonpuppiJEC_new.Mag();
            m_jlv_new     = (gleptonicV_new + ghadronicVpuppi_new).Mag();
	    //cout<<WLeptonic.pt()<<" old "<<WLeptonic.eta()<<"   "<<WLeptonic.phi()<<"   "<<WLeptonic.mass()<<endl;
	    //cout<<WLeptonic_new.pt()<<" new "<<WLeptonic_new.eta()<<"   "<<WLeptonic_new.phi()<<"   "<<WLeptonic_new.mass()<<endl;
            //cout<<jetAK8puppi_sd<<"  jetAK8puppi_sd  "<<ghadronicVpuppi.Mag()<<"   "<<ghadronicVpuppi_new.Mag()<<"   "<<ghadronicVpuppi.E()<<"   "<<ghadronicVpuppi_new.E()<<endl;    
            TLorentzVector ghadronicVpuppi_JEC_up, gravitonpuppiJEC_JEC_up,ghadronicVpuppi_2_JEC_up;
            ghadronicVpuppi_JEC_up.SetPtEtaPhiM(jetAK8puppi_ptJEC_JEC_up, jetAK8puppi_eta, jetAK8puppi_phi,jetAK8puppi_sd);
            ghadronicVpuppi_2_JEC_up.SetPtEtaPhiM(jetAK8puppi_ptJEC_2_JEC_up, jetAK8puppi_eta_2, jetAK8puppi_phi_2, jetAK8puppi_sd_2);
            gravitonpuppiJEC_JEC_up = gleptonicV_JEC_up + ghadronicVpuppi_JEC_up+ ghadronicVpuppi_2_JEC_up;
            candMasspuppiJEC_JEC_up     = gravitonpuppiJEC_JEC_up.Mag();
            m_jlv_JEC_up     = (gleptonicV_JEC_up + ghadronicVpuppi_JEC_up).Mag();
	    //cout<<ghadronicVpuppi_2_JEC_up.Pt()<<"  "<<candMasspuppiJEC_JEC_up<<endl;
            //cout<<jetAK8puppi_ptJEC_JEC_up<<"  "<<gleptonicV_JEC_up.Pt()<<"  "<<m_jlv_JEC_up<<endl;
            
            TLorentzVector ghadronicVpuppi_JEC_down, gravitonpuppiJEC_JEC_down,ghadronicVpuppi_2_JEC_down;
            ghadronicVpuppi_JEC_down.SetPtEtaPhiM(jetAK8puppi_ptJEC_JEC_down, jetAK8puppi_eta, jetAK8puppi_phi, jetAK8puppi_sd);
            ghadronicVpuppi_2_JEC_down.SetPtEtaPhiM(jetAK8puppi_ptJEC_2_JEC_down, jetAK8puppi_eta_2, jetAK8puppi_phi_2, jetAK8puppi_sd_2);
            gravitonpuppiJEC_JEC_down = gleptonicV_JEC_down + ghadronicVpuppi_JEC_down+ ghadronicVpuppi_2_JEC_down;
            candMasspuppiJEC_JEC_down     = gravitonpuppiJEC_JEC_down.Mag();
            m_jlv_JEC_down     = (gleptonicV_JEC_down + ghadronicVpuppi_JEC_down).Mag();
            
            TLorentzVector ghadronicVpuppi_JER_up, gravitonpuppiJEC_JER_up,ghadronicVpuppi_2_JER_up;
            ghadronicVpuppi_JER_up.SetPtEtaPhiM(jetAK8puppi_ptJEC_JER_up, jetAK8puppi_eta, jetAK8puppi_phi, jetAK8puppi_sd);
            ghadronicVpuppi_2_JER_up.SetPtEtaPhiM(jetAK8puppi_ptJEC_2_JER_up, jetAK8puppi_eta_2, jetAK8puppi_phi_2, jetAK8puppi_sd_2);
            gravitonpuppiJEC_JER_up = gleptonicV_JER_up + ghadronicVpuppi_JER_up+ ghadronicVpuppi_2_JER_up;
            candMasspuppiJEC_JER_up     = gravitonpuppiJEC_JER_up.Mag();
            m_jlv_JER_up     = (gleptonicV_JER_up + ghadronicVpuppi_JER_up).Mag();
            
            TLorentzVector ghadronicVpuppi_JER_down, gravitonpuppiJEC_JER_down,ghadronicVpuppi_2_JER_down;
            ghadronicVpuppi_JER_down.SetPtEtaPhiM(jetAK8puppi_ptJEC_JER_down, jetAK8puppi_eta, jetAK8puppi_phi, jetAK8puppi_sd);
            ghadronicVpuppi_2_JER_down.SetPtEtaPhiM(jetAK8puppi_ptJEC_2_JER_down, jetAK8puppi_eta_2, jetAK8puppi_phi_2,jetAK8puppi_sd_2);
            gravitonpuppiJEC_JER_down = gleptonicV_JER_down + ghadronicVpuppi_JER_down+ ghadronicVpuppi_2_JER_down;
            candMasspuppiJEC_JER_down     = gravitonpuppiJEC_JER_down.Mag();
            m_jlv_JER_down     = (gleptonicV_JER_down + ghadronicVpuppi_JER_down).Mag();

	    //swap var and var_new
            double jetAK8puppi_ptJEC_tmp = jetAK8puppi_ptJEC ; jetAK8puppi_ptJEC = jetAK8puppi_ptJEC_new ; jetAK8puppi_ptJEC_new = jetAK8puppi_ptJEC_tmp;
            double jetAK8puppi_ptJEC_2_tmp = jetAK8puppi_ptJEC_2 ; jetAK8puppi_ptJEC_2 = jetAK8puppi_ptJEC_2_new ; jetAK8puppi_ptJEC_2_new = jetAK8puppi_ptJEC_2_tmp;
            double jetAK8puppi_ptJEC_3_tmp = jetAK8puppi_ptJEC_3 ; jetAK8puppi_ptJEC_3 = jetAK8puppi_ptJEC_3_new ; jetAK8puppi_ptJEC_3_new = jetAK8puppi_ptJEC_3_tmp;
            double jetAK8puppi_e_tmp = jetAK8puppi_e ; jetAK8puppi_e = jetAK8puppi_e_new ; jetAK8puppi_e_new = jetAK8puppi_e_tmp;
            double jetAK8puppi_e_2_tmp = jetAK8puppi_e_2 ; jetAK8puppi_e_2 = jetAK8puppi_e_2_new ; jetAK8puppi_e_2_new = jetAK8puppi_e_2_tmp;
            double jetAK8puppi_e_3_tmp = jetAK8puppi_e_3 ; jetAK8puppi_e_3 = jetAK8puppi_e_3_new ; jetAK8puppi_e_3_new = jetAK8puppi_e_3_tmp;
            double ptVlepJEC_tmp = ptVlepJEC ; ptVlepJEC = ptVlepJEC_new ; ptVlepJEC_new = ptVlepJEC_tmp;
            double yVlepJEC_tmp = yVlepJEC ; yVlepJEC = yVlepJEC_new ; yVlepJEC_new = yVlepJEC_tmp;
            double phiVlepJEC_tmp = phiVlepJEC ; phiVlepJEC = phiVlepJEC_new ; phiVlepJEC_new = phiVlepJEC_tmp;
            double massVlepJEC_tmp = massVlepJEC ; massVlepJEC = massVlepJEC_new ; massVlepJEC_new = massVlepJEC_tmp;
            double mtVlepJEC_tmp = mtVlepJEC ; mtVlepJEC = mtVlepJEC_new ; mtVlepJEC_new = mtVlepJEC_tmp;
            double MET_et_tmp = MET_et ; MET_et = MET_et_new ; MET_et_new = MET_et_tmp;
            double MET_phi_tmp = MET_phi ; MET_phi = MET_phi_new ; MET_phi_new = MET_phi_tmp;
            double m_jlv_tmp = m_jlv ; m_jlv = m_jlv_new ; m_jlv_new = m_jlv_tmp;
            double candMasspuppiJEC_tmp = candMasspuppiJEC ; candMasspuppiJEC = candMasspuppiJEC_new ; candMasspuppiJEC_new = candMasspuppiJEC_tmp;

            lvw[0] = gleptonicV_new;
            lvw[1] = ghadronicVpuppi_new;
            lvw[2] = ghadronicVpuppi_2_new;
	    }
            delPhijetmet = deltaPhi(jetAK8puppi_phi, MET_phi);
            delPhijetlep = deltaPhi(jetAK8puppi_phi, phiVlepJEC);
            delPhijetmet_2 = deltaPhi(jetAK8puppi_phi_2, MET_phi);
            delPhijetlep_2 = deltaPhi(jetAK8puppi_phi_2, phiVlepJEC);
        
            delPhilepmet = deltaPhi(philep1, MET_phi);
            mtVlepJEC       =   sqrt(2*ptlep1*MET_et*(1.0-cos(philep1-MET_phi))); //WLeptonic.mt();

            if (!RunOnMC_){ 
            lvw[0] = gleptonicV;
            lvw[1] = ghadronicVpuppi;
            lvw[2] = ghadronicVpuppi_2;}
            Double_t Wpt[3];
            Wpt[0]=ptVlepJEC;
            Wpt[1]=jetAK8puppi_ptJEC;
            Wpt[2]=jetAK8puppi_ptJEC_2;
            Int_t *indexx=new Int_t[3];
            TMath::Sort(3,Wpt,indexx,1);
            //cout<<Wpt[indexx[0]]<<"   "<<Wpt[indexx[1]]<<"   "<<Wpt[indexx[2]]<<"   "<<endl;
            massww[0] = (lvw[indexx[0]]+lvw[indexx[1]]).Mag();
            massww[1] = (lvw[indexx[0]]+lvw[indexx[2]]).Mag();
            massww[2] = (lvw[indexx[1]]+lvw[indexx[2]]).Mag();

            masslvj1 = (lvw[0]+lvw[1]).Mag();
            masslvj2 = (lvw[0]+lvw[2]).Mag();
            massj1j2 = (lvw[1]+lvw[2]).Mag();
           
            edm::Handle<edm::View<pat::Jet> > jetsAK8;
            //jetsAK8Label_=puppijetInputToken_;
            iEvent.getByToken(jetsAK8Label_, jetsAK8);
        
            edm::View<pat::Jet>::const_iterator beginAK8 = jetsAK8->begin();
            edm::View<pat::Jet>::const_iterator endAK8 = jetsAK8->end();
            edm::View<pat::Jet>::const_iterator ijetAK8 = beginAK8;
    
            edm::View<pat::Jet>::const_iterator ijetAK8_j1,ijetAK8_j2;
            double drak8jetmatch1=10000.,drak8jetmatch2=10000.,tmpdrak8jet1=10000.,tmpdrak8jet2=10000.;
            // Loop over the "hard" jets
            for(ijetAK8 = beginAK8; ijetAK8 != endAK8; ++ijetAK8 ) {
                if(ijetAK8->pt()>0){
                    tmpdrak8jet1=deltaR(ijetAK8->eta(),ijetAK8->phi(),jetAK8puppi_eta,jetAK8puppi_phi);
                    //cout<<tmpdrak8jet1<<endl;
                    tmpdrak8jet2=deltaR(ijetAK8->eta(),ijetAK8->phi(),jetAK8puppi_eta_2,jetAK8puppi_phi_2);
                    if (tmpdrak8jet1<drak8jetmatch1) {drak8jetmatch1=tmpdrak8jet1; ijetAK8_j1=ijetAK8;}
                    if (tmpdrak8jet2<drak8jetmatch2) {drak8jetmatch2=tmpdrak8jet2; ijetAK8_j2=ijetAK8;}
                }
            }
            //if(jetAK8puppi_ptJEC>0) cout<<ijetAK8_j1->pt()<<"   "<<jetAK8puppi_ptJEC<<"  "<<drak8jetmatch1<<endl;
            //if(jetAK8puppi_ptJEC>0) cout<<ijetAK8_j1->eta()<<"   "<<ijetAK8_j1->phi()<<"   "<<jetAK8puppi_eta<<"   "<<jetAK8puppi_phi<<"   "<<endl;
            //if(jetAK8puppi_ptJEC_2>0) cout<<ijetAK8_j2->pt()<<"   "<<jetAK8puppi_ptJEC_2<<"  "<<drak8jetmatch2<<endl;
            //if(jetAK8puppi_ptJEC_2>0) cout<<ijetAK8_j2->eta()<<"   "<<ijetAK8_j2->phi()<<"   "<<jetAK8puppi_eta_2<<"   "<<jetAK8puppi_phi_2<<"   "<<endl;
            TLorentzVector aa,bb;
            puppi_softdropj1.SetPtEtaPhiM(0,0,0,0);
            puppi_softdropj2.SetPtEtaPhiM(0,0,0,0);

            auto const & sdSubjetsPuppi = ijetAK8_j1->subjets("SoftDropPuppi");
            //cout<<(ijetAK8_j1->subjets("SoftDropPuppi")).size()<<endl;
            //cout<<ijetAK8_j1->numberOfDaughters()<<endl;
            Int_t nsj=0;
            //cout<<sdSubjetsPuppi.at(1)->pt()<<endl;
            for ( auto const & puppiSDSJ : sdSubjetsPuppi ) {
                if(jetAK8puppi_ptJEC>0&&tmpdrak8jet1<1.0){
                    aa.SetPtEtaPhiM(puppiSDSJ->pt(),puppiSDSJ->eta(),puppiSDSJ->phi(),puppiSDSJ->mass());
                    //cout<<puppiSDSJ->jetArea()/3.14*3.14<<endl;
                    if (nsj==0)
                        ak8sj11.SetPtEtaPhiM(puppiSDSJ->pt(),puppiSDSJ->eta(),puppiSDSJ->phi(),puppiSDSJ->mass());
                    if (nsj==1)
                        ak8sj12.SetPtEtaPhiM(puppiSDSJ->pt(),puppiSDSJ->eta(),puppiSDSJ->phi(),puppiSDSJ->mass());
                    if (nsj==2)
                        ak8sj13.SetPtEtaPhiM(puppiSDSJ->pt(),puppiSDSJ->eta(),puppiSDSJ->phi(),puppiSDSJ->mass());
                    if (nsj==3)
                        ak8sj14.SetPtEtaPhiM(puppiSDSJ->pt(),puppiSDSJ->eta(),puppiSDSJ->phi(),puppiSDSJ->mass());
                    if (nsj==4)
                        ak8sj15.SetPtEtaPhiM(puppiSDSJ->pt(),puppiSDSJ->eta(),puppiSDSJ->phi(),puppiSDSJ->mass());
                    nsj++;
                    puppi_softdropj1+=aa;
                }
                //if(jetAK8puppi_ptJEC>0) cout<<nsj<<"   "<<puppiSDSJ->correctedP4(0).phi()<<"   "<<puppiSDSJ->phi()<<endl;
                //cout<<puppi_softdrop_subjet.Pt()<<"   "<<puppiSDSJ->pt()<<endl;
            }
            auto const & sdSubjetsPuppi_2 = ijetAK8_j2->subjets("SoftDropPuppi");
            Int_t nsj2=0;
            for ( auto const & puppiSDSJ_2 : sdSubjetsPuppi_2 ) {
                if(jetAK8puppi_ptJEC_2>0&&tmpdrak8jet2<1.0){
                    bb.SetPtEtaPhiM(puppiSDSJ_2->pt(),puppiSDSJ_2->eta(),puppiSDSJ_2->phi(),puppiSDSJ_2->mass());
                    if (nsj2==0)
                        ak8sj21.SetPtEtaPhiM(puppiSDSJ_2->pt(),puppiSDSJ_2->eta(),puppiSDSJ_2->phi(),puppiSDSJ_2->mass());
                    if (nsj2==1)
                        ak8sj22.SetPtEtaPhiM(puppiSDSJ_2->pt(),puppiSDSJ_2->eta(),puppiSDSJ_2->phi(),puppiSDSJ_2->mass());
                    if (nsj2==2)
                        ak8sj23.SetPtEtaPhiM(puppiSDSJ_2->pt(),puppiSDSJ_2->eta(),puppiSDSJ_2->phi(),puppiSDSJ_2->mass());
                    if (nsj2==3)
                        ak8sj24.SetPtEtaPhiM(puppiSDSJ_2->pt(),puppiSDSJ_2->eta(),puppiSDSJ_2->phi(),puppiSDSJ_2->mass());
                    if (nsj2==4)
                        ak8sj25.SetPtEtaPhiM(puppiSDSJ_2->pt(),puppiSDSJ_2->eta(),puppiSDSJ_2->phi(),puppiSDSJ_2->mass());
                    nsj2++;
                    puppi_softdropj2+=bb;
                }
            }
        }//1
        if((nLooseEle==1||nLooseMu==1)&&jetAK8puppi_ptJEC>200&&IDLoose==1&& fabs(jetAK8puppi_eta)<2.4 && IDLoose>0&&((jetAK8puppi_ptJEC_2>100)?(fabs(jetAK8puppi_eta_2)<2.4 && IDLoose_2>0):1)&&((jetAK8puppi_ptJEC_3>100)?(fabs(jetAK8puppi_eta_3)<2.4 && IDLoose_3>0):1)) outTree_->Fill();
    outTreew_->Fill();
    //outTree_->Fill();
	}
    
    else {
        outTreew_->Fill();
    //outTree_->Fill();
    }
//cout<< "test end3" <<endl;
}
//-------------------------------------------------------------------------------------------------------------------------------------//




// ------------ method called once each job just before starting event loop  ------------

// ------------ method called once each job just after ending the event loop  ------------


//define this as a plug-in
DEFINE_FWK_MODULE(VVVTreeMaker);
