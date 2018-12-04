// -*- C++ -*-
//
// Package:    DiphotonAnalyzer/TreeProducer
// Class:      TreeProducer
//
/**\class TreeProducer TreeProducer.cc DiphotonAnalyzer/TreeProducer/plugins/TreeProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Laurent Forthomme
//         Created:  Tue, 13 Sep 2016 03:57:43 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
//                               JW
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/DataFormats/interface/Jet.h"
//
#include "flashgg/DataFormats/interface/Met.h"

#include "DiphotonAnalyzer/TreeProducer/interface/HLTMatcher.h"
#include "DiphotonAnalyzer/TreeProducer/interface/SelectionUtils.h"
#include "DiphotonAnalyzer/TreeProducer/interface/FillNumberLUTHandler.h"
#include "DiphotonAnalyzer/TreeProducer/interface/TreeEvent.h"

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

class TreeProducer : public edm::one::EDAnalyzer<edm::one::WatchRuns,edm::one::SharedResources>
{
  public:
    explicit TreeProducer( const edm::ParameterSet& );
    ~TreeProducer();

    /// fill 'descriptions' with the allowed parameters for the module
    static void fillDescriptions( edm::ConfigurationDescriptions& );

  private:
    /// method called once each job just before starting event loop
    void beginJob() override;
    void beginRun( const edm::Run&, const edm::EventSetup& ) override;
    /// method called for each event
    void analyze( const edm::Event&, const edm::EventSetup& ) override;
    void endRun( const edm::Run&, const edm::EventSetup& ) override {}
    /// method called once each job just after ending the event loop
    void endJob() override {}

    // ----------member data ---------------------------

    /// get the trigger information from the event
    void analyzeTriggers( const edm::Event&, const edm::EventSetup& );

    edm::InputTag triggerResults_;

    edm::EDGetTokenT<edm::View<CTPPSLocalTrackLite> > totemRPTracksToken_;
    edm::EDGetTokenT<edm::View<flashgg::DiPhotonCandidate> > diphotonToken_;
    edm::EDGetTokenT<edm::View<flashgg::Met> > metToken_;
    edm::EDGetTokenT<edm::View<reco::Vertex> > vtxToken_;
    edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
//                                 JW
    edm::EDGetTokenT<edm::View<flashgg::Electron> > electronToken_;
    edm::EDGetTokenT<edm::View<flashgg::Muon> > muonToken_;
    edm::EDGetTokenT<edm::View<vector<flashgg::Jet> > > jetToken_;
//
    edm::EDGetTokenT<edm::View<reco::GenParticle> > genPartToken_;
    edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > genPhoToken_;
    edm::EDGetTokenT<edm::View<PileupSummaryInfo> > pileupToken_;
    edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;

    bool isData_;
    double sqrtS_;
    double singlePhotonMinPt_, singlePhotonMaxEta_, singlePhotonMinR9_;
    double photonPairMinMass_, photonPairMaxMass_;
    std::string filename_;
    double maxGenLevelDR_;
    edm::FileInPath puMCfile_, puDatafile_;
    std::string puMCpath_, puDatapath_;

    std::vector<std::string> triggersList_;
    HLTConfigProvider hlt_config_;
    HLTPrescaleProvider hlt_prescale_;
    HLTMatcher hlt_matcher_;

    std::unique_ptr<CTPPSAlCa::FillNumberLUTHandler> fillLUTHandler_;
    std::unique_ptr<edm::LumiReWeighting> lumiReweighter_;

    std::unique_ptr<TFile> file_;
    TTree* tree_;

    TreeEvent ev_;
};

TreeProducer::TreeProducer( const edm::ParameterSet& iConfig ) :
  triggerResults_( iConfig.getParameter<edm::InputTag>( "triggerResults" ) ),
  totemRPTracksToken_ ( consumes<edm::View<CTPPSLocalTrackLite> >       ( iConfig.getParameter<edm::InputTag>( "totemRPTracksLabel") ) ),
  diphotonToken_      ( consumes<edm::View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<edm::InputTag>( "diphotonLabel" ) ) ),
  metToken_           ( consumes<edm::View<flashgg::Met> >              ( iConfig.getParameter<edm::InputTag>( "metLabel") ) ),
  vtxToken_           ( consumes<edm::View<reco::Vertex> >              ( iConfig.getParameter<edm::InputTag>( "vertexLabel" ) ) ),
  beamSpotToken_      ( consumes<reco::BeamSpot>                        ( iConfig.getParameter<edm::InputTag>( "beamSpotLabel" ) ) ),
//                               JW
  electronToken_      ( consumes<edm::View<flashgg::Electron> >         ( iConfig.getParameter<edm::InputTag>( "electronLabel") ) ),
  muonToken_          ( consumes<edm::View<flashgg::Muon> >             ( iConfig.getParameter<edm::InputTag>( "muonLabel") ) ),
  jetToken_           ( consumes<edm::View<std::vector<flashgg::Jet> > >( iConfig.getParameter<edm::InputTag>( "jetLabel") ) ),
//
  genPartToken_       ( consumes<edm::View<reco::GenParticle> >         ( iConfig.getParameter<edm::InputTag>( "genPartLabel" ) ) ),
  genPhoToken_        ( consumes<edm::View<pat::PackedGenParticle> >    ( iConfig.getParameter<edm::InputTag>( "genPhoLabel" ) ) ),
  pileupToken_        ( consumes<edm::View<PileupSummaryInfo> >         ( iConfig.getUntrackedParameter<edm::InputTag>( "pileupLabel", edm::InputTag( "slimmedAddPileupInfo" ) ) ) ),
  triggerResultsToken_( consumes<edm::TriggerResults>                   ( triggerResults_ ) ),
  isData_             ( iConfig.getParameter<bool>       ( "isData" ) ),
  sqrtS_              ( iConfig.getParameter<double>     ( "sqrtS" ) ),
  singlePhotonMinPt_  ( iConfig.getParameter<double>     ( "minPtSinglePhoton" ) ),
  singlePhotonMaxEta_ ( iConfig.getParameter<double>     ( "maxEtaSinglePhoton" ) ),
  singlePhotonMinR9_  ( iConfig.getParameter<double>     ( "minR9SinglePhoton" ) ),
  photonPairMinMass_  ( iConfig.getParameter<double>     ( "minMassDiPhoton" ) ),
  photonPairMaxMass_  ( iConfig.getUntrackedParameter<double>( "maxMassDiPhoton", -1. ) ),
  filename_           ( iConfig.getParameter<std::string>( "outputFilename" ) ),
  maxGenLevelDR_      ( iConfig.getParameter<double>     ( "maxGenLevelDeltaR" ) ),
  puMCpath_           ( iConfig.getUntrackedParameter<std::string>( "pileupMCPath", "pileup" ) ),
  puDatapath_         ( iConfig.getUntrackedParameter<std::string>( "pileupDataPath", "pileup" ) ),
  triggersList_       ( iConfig.getParameter<std::vector<std::string> >( "triggersList" ) ),
  hlt_prescale_       ( iConfig, consumesCollector(), *this ),
  hlt_matcher_( triggersList_ ),
  file_( TFile::Open( filename_.c_str(), "recreate" ) ),
  tree_( nullptr )
{
  if ( !isData_ ) {
    puMCfile_ =iConfig.getParameter<edm::FileInPath>( "pileupMCFile" );
    puDatafile_ = iConfig.getParameter<edm::FileInPath>( "pileupDataFile" );
    std::cout << ">> Pileup reweighting will be used with files:\n\t" << puMCfile_.fullPath() << " for MC ;\n\t" << puDatafile_.fullPath() << " for data" << std::endl;
    lumiReweighter_ = std::make_unique<edm::LumiReWeighting>( puMCfile_.fullPath(), puDatafile_.fullPath(), puMCpath_, puDatapath_ );
  }
  else {
    fillLUTHandler_ = std::make_unique<CTPPSAlCa::FillNumberLUTHandler>( iConfig.getParameter<edm::FileInPath>( "fillNumLUTFile" ).fullPath().c_str() );
  }

  file_->cd();
  tree_ = new TTree( "ntp", "diphoton ntuple" );
}


TreeProducer::~TreeProducer()
{
  if ( file_ ) {
    file_->Write();
    file_->Close();
  }
  if ( tree_ ) delete tree_;
}

void
TreeProducer::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  #include <iostream> // for debugging purposes
  if ( isData_ && !iEvent.isRealData() )
    throw cms::Exception( "TreeProducer" ) << "Looks like you are trying to run on MC with the 'isData' flag set to True... Aborting!";

  analyzeTriggers( iEvent, iSetup );

  //std::cout << "---> passing the trigger!" << std::endl;

  ev_.clear();
  const reco::Candidate::Point orig( -999., -999, -999. );

  // Run and BX information
  ev_.bunch_crossing = iEvent.bunchCrossing();
  ev_.run_id = iEvent.id().run();
  ev_.lumisection = iEvent.luminosityBlock();
  ev_.event_number = iEvent.id().event();

  // get the fill number from the run id <-> fill number LUT
  if ( fillLUTHandler_ ) {
    ev_.fill_number = fillLUTHandler_->getFillNumber( iEvent.id().run() );
  }
  else ev_.fill_number = 1;

  //----- gen-level information -----

  edm::Handle<edm::View<pat::PackedGenParticle> > genPhotons;
  edm::Handle<edm::View<reco::GenParticle> > genParts;

  if ( !iEvent.isRealData() ) {
    iEvent.getByToken( genPhoToken_, genPhotons );
    for ( unsigned int i=0; i<genPhotons->size() && ev_.num_gen_photon<ev_.MAX_GEN_PHOTON; i++ ) {
      const edm::Ptr<pat::PackedGenParticle> genPho = genPhotons->ptrAt( i );

      //----- gen-level information

      if ( genPho->pdgId()!=22 ) continue; // only keep gen-level photons

      ev_.gen_photon_pt[ev_.num_gen_photon] = genPho->pt();
      ev_.gen_photon_eta[ev_.num_gen_photon] = genPho->eta();
      ev_.gen_photon_phi[ev_.num_gen_photon] = genPho->phi();
      ev_.gen_photon_energy[ev_.num_gen_photon] = genPho->energy();

      ev_.gen_photon_vertex_x[ev_.num_gen_photon] = genPho->vx();
      ev_.gen_photon_vertex_y[ev_.num_gen_photon] = genPho->vy();
      ev_.gen_photon_vertex_z[ev_.num_gen_photon] = genPho->vz();

      ev_.num_gen_photon++;
    }

    //----- smeared-level information

    iEvent.getByToken( genPartToken_, genParts );
    for ( unsigned int i=0; i<genParts->size() && ev_.num_gen_part<ev_.MAX_GEN_PART; i++ ) {
      const edm::Ptr<reco::GenParticle> genPart = genParts->ptrAt( i );

      ev_.gen_part_pdgid[i] = genPart->pdgId();
      ev_.gen_part_status[i] = genPart->status();

      ev_.gen_part_pt[i] = genPart->pt();
      ev_.gen_part_eta[i] = genPart->eta();
      ev_.gen_part_phi[i] = genPart->phi();
      ev_.gen_part_energy[i] = genPart->energy();

      ev_.gen_part_vertex_x[i] = genPart->vx();
      ev_.gen_part_vertex_y[i] = genPart->vy();
      ev_.gen_part_vertex_z[i] = genPart->vz();

      ev_.num_gen_part++;
    }

  }
  //std::cout << "---> passing the mc matching!" << std::endl;

  //----- diphoton candidates -----

  edm::Handle<edm::View<flashgg::DiPhotonCandidate> > diphotons;
  iEvent.getByToken( diphotonToken_, diphotons );

  edm::Handle<edm::View<std::vector<flashgg::Jet> > > jetsColls;
  iEvent.getByToken( jetToken_, jetsColls );

  ev_.num_diphoton = 0;
  ev_.num_jet = 0;
  edm::Ptr<reco::Vertex> diphoton_vtx[ev_.MAX_DIPHOTON];
  for ( unsigned int i=0; i<diphotons->size() && ev_.num_diphoton<ev_.MAX_DIPHOTON; i++ ) {
    const edm::Ptr<flashgg::DiPhotonCandidate> diphoton = diphotons->ptrAt( i );

    if ( diphoton->leadPhotonId()<-0.9 ) continue;
    if ( diphoton->subLeadPhotonId()<-0.9 ) continue;

    const flashgg::Photon* lead_pho = diphoton->leadingPhoton(),
                          *sublead_pho = diphoton->subLeadingPhoton();
    if ( !passSinglePhotonCuts( lead_pho ) ) continue;
    if ( !passSinglePhotonCuts( sublead_pho ) ) continue;

    if ( fabs( lead_pho->eta() )>=singlePhotonMaxEta_ || fabs( sublead_pho->eta() )>=singlePhotonMaxEta_ ) continue;
    if ( lead_pho->pt()<singlePhotonMinPt_ || sublead_pho->pt()<singlePhotonMinPt_ ) continue;
    if ( lead_pho->full5x5_r9()<singlePhotonMinR9_ || sublead_pho->full5x5_r9()<singlePhotonMinR9_ ) continue;

    if ( photonPairMaxMass_>0. && diphoton->mass()>photonPairMaxMass_ ) continue;
    if ( diphoton->mass()<photonPairMinMass_ ) continue;

    //ev_.diphoton_vertex_tracks[ev_.num_diphoton] = diphoton->vtx()->nTracks(); // 0 for miniAOD...
    ev_.diphoton_vertex_x[ev_.num_diphoton] = diphoton->vtx()->x();
    ev_.diphoton_vertex_y[ev_.num_diphoton] = diphoton->vtx()->y();
    ev_.diphoton_vertex_z[ev_.num_diphoton] = diphoton->vtx()->z();
    diphoton_vtx[ev_.num_diphoton] = diphoton->vtx();

    ev_.diphoton_pt1[ev_.num_diphoton] = lead_pho->pt();
    ev_.diphoton_eta1[ev_.num_diphoton] = lead_pho->eta();
    ev_.diphoton_phi1[ev_.num_diphoton] = lead_pho->phi();
    ev_.diphoton_energy1[ev_.num_diphoton] = lead_pho->energy();
    ev_.diphoton_r91[ev_.num_diphoton] = lead_pho->full5x5_r9();
    ev_.diphoton_id1[ev_.num_diphoton] = lead_pho->phoIdMvaDWrtVtx( diphoton->vtx() );
    ev_.diphoton_sigeove1[ev_.num_diphoton] = lead_pho->sigEOverE();

    ev_.diphoton_pt2[ev_.num_diphoton] = sublead_pho->pt();
    ev_.diphoton_eta2[ev_.num_diphoton] = sublead_pho->eta();
    ev_.diphoton_phi2[ev_.num_diphoton] = sublead_pho->phi();
    ev_.diphoton_energy2[ev_.num_diphoton] = sublead_pho->energy();
    ev_.diphoton_r92[ev_.num_diphoton] = sublead_pho->full5x5_r9();
    ev_.diphoton_id2[ev_.num_diphoton] = sublead_pho->phoIdMvaDWrtVtx( diphoton->vtx() );
    ev_.diphoton_sigeove2[ev_.num_diphoton] = sublead_pho->sigEOverE();

    ev_.diphoton_mass[ev_.num_diphoton] = diphoton->mass();
    ev_.diphoton_rapidity[ev_.num_diphoton] = diphoton->rapidity();
    ev_.diphoton_pt[ev_.num_diphoton] = diphoton->pt();

    float dphi = lead_pho->phi()-sublead_pho->phi();
    while ( dphi<-TMath::Pi() ) dphi += 2.*TMath::Pi();
    while ( dphi> TMath::Pi() ) dphi -= 2.*TMath::Pi();
    ev_.diphoton_dphi[ev_.num_diphoton] = dphi;

    //----- retrieve the SC information

    ev_.diphoton_supercluster_x1[ev_.num_diphoton] = lead_pho->superCluster()->x();
    ev_.diphoton_supercluster_y1[ev_.num_diphoton] = lead_pho->superCluster()->y();
    ev_.diphoton_supercluster_z1[ev_.num_diphoton] = lead_pho->superCluster()->z();
    ev_.diphoton_supercluster_x2[ev_.num_diphoton] = sublead_pho->superCluster()->x();
    ev_.diphoton_supercluster_y2[ev_.num_diphoton] = sublead_pho->superCluster()->y();
    ev_.diphoton_supercluster_z2[ev_.num_diphoton] = sublead_pho->superCluster()->z();

    //----- retrieve the gen-level information

    if ( !iEvent.isRealData() ) {
      reco::Candidate::Point vtx1( orig ), vtx2( orig ),
                             vtx1_smear( orig ), vtx2_smear( orig );
      if ( lead_pho->matchedGenPhoton() ) {
        ev_.diphoton_genpt1[ev_.num_diphoton] = lead_pho->matchedGenPhoton()->pt();
        ev_.diphoton_geneta1[ev_.num_diphoton] = lead_pho->matchedGenPhoton()->eta();
        ev_.diphoton_genphi1[ev_.num_diphoton] = lead_pho->matchedGenPhoton()->phi();
        ev_.diphoton_genenergy1[ev_.num_diphoton] = lead_pho->matchedGenPhoton()->energy();
        vtx1 = lead_pho->matchedGenPhoton()->vertex();
      }
      if ( sublead_pho->matchedGenPhoton() ) {
        ev_.diphoton_genpt2[ev_.num_diphoton] = sublead_pho->matchedGenPhoton()->pt();
        ev_.diphoton_geneta2[ev_.num_diphoton] = sublead_pho->matchedGenPhoton()->eta();
        ev_.diphoton_genphi2[ev_.num_diphoton] = sublead_pho->matchedGenPhoton()->phi();
        ev_.diphoton_genenergy2[ev_.num_diphoton] = sublead_pho->matchedGenPhoton()->energy();
        vtx2 = sublead_pho->matchedGenPhoton()->vertex();
      }
      if ( vtx2!=vtx1 ) {
        std::cerr << "-> 2 different gen-level vertices: " << std::endl
                  << "      leading photon vertex: (" << vtx1.x() << ", " << vtx1.y() << ", " << vtx1.z() << ")" << std::endl
                  << "   subleading photon vertex: (" << vtx2.x() << ", " << vtx2.y() << ", " << vtx2.z() << ")" << std::endl;
      }
      if ( vtx1!=orig ) {
        ev_.diphoton_genvertex_x = vtx1.x();
        ev_.diphoton_genvertex_y = vtx1.y();
        ev_.diphoton_genvertex_z = vtx1.z();
      }
      else {
        ev_.diphoton_genvertex_x = vtx2.x();
        ev_.diphoton_genvertex_y = vtx2.y();
        ev_.diphoton_genvertex_z = vtx2.z();
      }

      //----- retrieve the smeared-level generated vertex

      for ( unsigned int j=0; j<genParts->size(); j++ ) {
        const edm::Ptr<reco::GenParticle> genPart = genParts->ptrAt( j );
        if ( sqrt( pow( genPart->eta()-   lead_pho->eta(), 2 ) + pow( genPart->phi()-   lead_pho->phi(), 2 ) )<maxGenLevelDR_ ) { vtx1_smear = genPart->vertex(); continue; }
        if ( sqrt( pow( genPart->eta()-sublead_pho->eta(), 2 ) + pow( genPart->phi()-sublead_pho->phi(), 2 ) )<maxGenLevelDR_ ) { vtx2_smear = genPart->vertex(); continue; }
      }
      /*if ( vtx2_smear!=vtx1_smear ) {
        std::cerr << "-> 2 different smeared-level vertices: " << std::endl
                  << "      leading photon vertex: (" << vtx1_smear.x() << ", " << vtx1_smear.y() << ", " << vtx1_smear.z() << ")" << std::endl
                  << "   subleading photon vertex: (" << vtx2_smear.x() << ", " << vtx2_smear.y() << ", " << vtx2_smear.z() << ")" << std::endl;
      }*/
      if ( vtx1_smear!=orig ) {
        ev_.diphoton_genvertex_smeared_x = vtx1_smear.x();
        ev_.diphoton_genvertex_smeared_y = vtx1_smear.y();
        ev_.diphoton_genvertex_smeared_z = vtx1_smear.z();
      }
      else {
        ev_.diphoton_genvertex_smeared_x = vtx2_smear.x();
        ev_.diphoton_genvertex_smeared_y = vtx2_smear.y();
        ev_.diphoton_genvertex_smeared_z = vtx2_smear.z();
      }
    }

    //----- retrieve the associated jets

    if ( diphoton->jetCollectionIndex()<jetsColls->size() ) {
      for ( unsigned int j=0; j<jetsColls->at( diphoton->jetCollectionIndex() ).size(); j++ ) {
        if ( ev_.num_jet>=ev_.MAX_JET ) { std::cerr << ">> More jets than expected in this event (" << ev_.num_jet << ">MAX_JET=" << ev_.MAX_JET << "). Increase MAX_JET for safety" << std::endl; }

        const flashgg::Jet jet = jetsColls->at( diphoton->jetCollectionIndex() ).at( j );
        ev_.jet_pt[ev_.num_jet]   = jet.pt();
        ev_.jet_eta[ev_.num_jet]  = jet.eta();
        ev_.jet_phi[ev_.num_jet]  = jet.phi();
        ev_.jet_energy[ev_.num_jet]    = jet.energy();
        ev_.jet_mass[ev_.num_jet] = jet.mass();

        ev_.jet_vtx_x[ev_.num_jet] = jet.vertex().x();
        ev_.jet_vtx_y[ev_.num_jet] = jet.vertex().y();
        ev_.jet_vtx_z[ev_.num_jet] = jet.vertex().z();
        ev_.jet_dipho_match[ev_.num_jet] = i;
        ev_.num_jet++;
      }
    }

    ev_.num_diphoton++;
  }

  //----- only store events with > 0 diphoton candidate(s) -----

  if ( ev_.num_diphoton<1 ) return;

  //----- beam spot -----

  edm::Handle<reco::BeamSpot> beamSpot;
  iEvent.getByToken( beamSpotToken_, beamSpot );
  if ( beamSpot.isValid() ) {
    ev_.bs_x0 = beamSpot->x0();
    ev_.bs_y0 = beamSpot->y0();
    ev_.bs_z0 = beamSpot->z0();
    ev_.bs_sigma_z = beamSpot->sigmaZ();
    ev_.bs_dxdz = beamSpot->dxdz();
    ev_.bs_beam_width_x = beamSpot->BeamWidthX();
    ev_.bs_beam_width_y = beamSpot->BeamWidthY();
  }

  //----- forward RP tracks -----

  if ( iEvent.isRealData() ) {
    edm::Handle<edm::View<CTPPSLocalTrackLite> > rpLocalTracks;
    iEvent.getByToken( totemRPTracksToken_, rpLocalTracks );

    ev_.num_proton_track = 0;
    for ( const auto& trk : *rpLocalTracks ) {
      const CTPPSDetId detid( trk.getRPId() );

      ev_.proton_track_x[ev_.num_proton_track] = trk.getX() * 1.e-3; // store in m
      ev_.proton_track_y[ev_.num_proton_track] = trk.getY() * 1.e-3; // store in m
      ev_.proton_track_arm[ev_.num_proton_track] = detid.arm(); // 0 = left (45) ; 1 = right (56)
      ev_.proton_track_station[ev_.num_proton_track] = detid.station();
      ev_.proton_track_pot[ev_.num_proton_track] = detid.rp(); // 2 = 210n ; 3 = 210f

      //ev_.proton_track_chi2[ev_.num_proton_track] = trk.getChiSquared();
      //ev_.proton_track_normchi2[ev_.num_proton_track] = trk.getChiSquaredOverNDF();

      ev_.num_proton_track++;
    }
  }
  //                               JW

  //----- electrons collection -----

  edm::Handle<edm::View<flashgg::Electron> > electrons;
  iEvent.getByToken( electronToken_, electrons );

  ev_.num_electron = 0;
  for ( unsigned int i=0; i<electrons->size() && ev_.num_electron<ev_.MAX_ELECTRON; i++ ) {
    const edm::Ptr<flashgg::Electron> electron = electrons->ptrAt( i );
    if ( ev_.num_electron>=ev_.MAX_ELECTRON ) { std::cerr << ">> More electrons than expected in this event (" << ev_.num_electron << ">MAX_ELECTRON=" << ev_.MAX_ELECTRON << "). Increase MAX_ELECTRON for safety" << std::endl; }

    ev_.electron_pt[ev_.num_electron] = electron->pt();
    ev_.electron_eta[ev_.num_electron] = electron->eta();
    ev_.electron_phi[ev_.num_electron] = electron->phi();
    ev_.electron_energy[ev_.num_electron] = electron->energy();

    ev_.electron_vtx_x[ev_.num_electron] = electron->vertex().x();
    ev_.electron_vtx_y[ev_.num_electron] = electron->vertex().y();
    ev_.electron_vtx_z[ev_.num_electron] = electron->vertex().z();
    ev_.num_electron++;
  }

  //----- muons collection -----

  edm::Handle<edm::View<flashgg::Muon> > muons;
  iEvent.getByToken(muonToken_,muons);

  ev_.num_muon = 0;
  for ( unsigned int i=0; i<muons->size() && ev_.num_muon<ev_.MAX_MUON; i++ ) {
    const edm::Ptr<flashgg::Muon> muon = muons->ptrAt( i );
    if ( ev_.num_muon>=ev_.MAX_MUON ) { std::cerr << ">> More muons than expected in this event (" << ev_.num_muon << ">MAX_MUON=" << ev_.MAX_MUON << "). Increase MAX_MUON for safety" << std::endl; }

    ev_.muon_pt[ev_.num_muon] = muon->pt();
    ev_.muon_eta[ev_.num_muon] = muon->eta();
    ev_.muon_phi[ev_.num_muon] = muon->phi();
    ev_.muon_energy[ev_.num_muon] = muon->energy();

    ev_.muon_vtx_x[ev_.num_muon] = muon->vertex().x();
    ev_.muon_vtx_y[ev_.num_muon] = muon->vertex().y();
    ev_.muon_vtx_z[ev_.num_muon] = muon->vertex().z();
    ev_.num_muon++;
  }

  //----- event-wide information -----

  std::cout << "# found " << ev_.num_diphoton << " diphoton candidate(s) with " << ev_.num_proton_track << " proton track(s) (all pots)!" << std::endl;
  // retrieve the missing ET
  edm::Handle<edm::View<flashgg::Met> > mets;
  iEvent.getByToken( metToken_, mets );

  const edm::View<flashgg::Met>* metColl = mets.product();
  edm::View<flashgg::Met>::const_iterator met = metColl->begin();
  ev_.met = met->et();
  ev_.met_phi = met->phi();
  ev_.met_significance = met->significance();

  //std::cout << "passed MET" << std::endl;
  //--- pileup information ---
  if ( !iEvent.isRealData() ) {
    edm::Handle<edm::View<PileupSummaryInfo> > pileupInfo;
    iEvent.getByToken( pileupToken_, pileupInfo );
    int npv0true = 0;
    for ( unsigned int i=0; i<pileupInfo->size(); i++ ) {
      const edm::Ptr<PileupSummaryInfo> PVI = pileupInfo->ptrAt( i );
      const int beamXing = PVI->getBunchCrossing(), npvtrue = PVI->getTrueNumInteractions();
      if ( beamXing==0 ) { npv0true += npvtrue; }
    }
    ev_.pileup_weight = lumiReweighter_->weight( npv0true );
  }

  //std::cout << "passed PU" << std::endl;
  //----- vertexing information -----

  edm::Handle<edm::View<reco::Vertex> > vertices;
  iEvent.getByToken( vtxToken_, vertices );

  ev_.num_vertex = 0;
  for ( unsigned int i=0; i<vertices->size() && i<ev_.MAX_VERTEX; i++ ) {
    const edm::Ptr<reco::Vertex> vtx = vertices->ptrAt( i );
    if ( ev_.num_vertex>=ev_.MAX_VERTEX ) { std::cerr << ">> More vertices than expected in this event (" << ev_.num_vertex << ">MAX_VERTEX=" << ev_.MAX_VERTEX << "). Increase MAX_VERTEX for safety" << std::endl; }
    if ( !vtx->isValid() ) continue;

    // loop over all the diphoton candidates to find the closest vertices
    for ( unsigned int j=0; j<ev_.num_diphoton; j++ ) {
      if ( diphoton_vtx[j]->position()==vtx->position() ) continue; // found the diphoton vertex
      const float vtx_dist = sqrt( pow( diphoton_vtx[j]->x()-vtx->x(), 2 )+pow( diphoton_vtx[j]->y()-vtx->y(), 2 )+pow( diphoton_vtx[j]->z()-vtx->z(), 2 ) );
      if ( vtx_dist < ev_.diphoton_nearestvtxdist[j] ) ev_.diphoton_nearestvtxdist[j] = vtx_dist;
      if ( vtx_dist <= 0.1 ) ev_.diphoton_vertex_vtx1mmdist[j]++;
      if ( vtx_dist <= 0.2 ) ev_.diphoton_vertex_vtx2mmdist[j]++;
      if ( vtx_dist <= 0.5 ) ev_.diphoton_vertex_vtx5mmdist[j]++;
      if ( vtx_dist <= 1.0 ) ev_.diphoton_vertex_vtx1cmdist[j]++;
      ev_.diphoton_vertex_id[j] = ev_.num_vertex;
    }

    ev_.vertex_x[ev_.num_vertex] = vtx->x();
    ev_.vertex_y[ev_.num_vertex] = vtx->y();
    ev_.vertex_z[ev_.num_vertex] = vtx->z();

    // tracks content not stored in the miniAOD event format...
    /*ev_.vertex_tracks[ev_.num_vertex] = vtx->nTracks();
    ev_.vertex_tracksWght0p75[ev_.num_vertex] = vtx->nTracks( 0.75 );
    ev_.vertex_tracksWght0p9[ev_.num_vertex] = vtx->nTracks( 0.9 );
    ev_.vertex_tracksWght0p95[ev_.num_vertex] = vtx->nTracks( 0.95 );*/
    ev_.num_vertex++;
  }

  std::cout << ">> Filling the tree" << std::endl;
  tree_->Fill();
}

void
TreeProducer::analyzeTriggers( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  ev_.num_hlt = triggersList_.size();
  edm::Handle<edm::TriggerResults> hltResults;
  iEvent.getByToken( triggerResultsToken_, hltResults );
  const edm::TriggerNames& trigNames = iEvent.triggerNames( *hltResults );

  std::ostringstream os;
  os << "Prescale set: " << hlt_prescale_.prescaleSet(iEvent, iSetup) << "\n"
     << "Trigger names: " << std::endl;
  for ( unsigned int i=0; i<trigNames.size(); i++ ) {
    const int trigNum = hlt_matcher_.TriggerNum( trigNames.triggerNames().at( i ) );
    if ( trigNum<0 ) continue; // Trigger didn't match the interesting ones
    ev_.hlt_accept[trigNum] = hltResults->accept( i );

    // extract prescale value for this path
    ev_.hlt_prescale[trigNum] = 1;
    if ( ev_.hlt_accept[trigNum] == 1 && isData_ ) { //FIXME
      int prescale_set = hlt_prescale_.prescaleSet( iEvent, iSetup );
      ev_.hlt_prescale[trigNum] = ( prescale_set < 0 )
        ? 0.
        : hlt_config_.prescaleValue( prescale_set, trigNames.triggerNames().at( i ) );
    }
    os << "--> (fired? " << ev_.hlt_accept[trigNum] << ") " << trigNames.triggerNames().at( i ) << "\tprescale: " << ev_.hlt_prescale[trigNum] << std::endl;
  }
  LogDebug( "TreeProducer" ) << os.str();
  //#include <iostream>
  //std::cout << os.str() << std::endl;
}

void
TreeProducer::beginJob()
{
  ev_.create( tree_, isData_ );
}

void
TreeProducer::beginRun( const edm::Run& iRun, const edm::EventSetup& iSetup )
{
  bool changed = true;
  if ( !hlt_prescale_.init( iRun, iSetup, triggerResults_.process(), changed ) )
    throw cms::Exception("TreeProducer") << " prescales extraction failure with process name " << triggerResults_.process();
  hlt_config_ = hlt_prescale_.hltConfigProvider();
  if ( hlt_config_.size() == 0 )
    edm::LogError("TreeProducer") << "HLT config size error";
  if ( changed ) {
    // The HLT config has actually changed wrt the previous Run, hence rebook your
    // histograms or do anything else dependent on the revised HLT config
  }
}

void
TreeProducer::fillDescriptions( edm::ConfigurationDescriptions& descriptions )
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault( desc );
}

//define this as a plug-in
DEFINE_FWK_MODULE( TreeProducer );
