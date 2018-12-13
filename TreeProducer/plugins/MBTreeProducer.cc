// -*- C++ -*-
//
// Package:    DiphotonAnalyzer/TreeProducer
// Class:      MBTreeProducer
//
/**\class MBTreeProducer MBTreeProducer.cc DiphotonAnalyzer/TreeProducer/plugins/MBTreeProducer.cc

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

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"

#include "DiphotonAnalyzer/TreeProducer/interface/SelectionUtils.h"
#include "DiphotonAnalyzer/TreeProducer/interface/FillNumberLUTHandler.h"
#include "DiphotonAnalyzer/TreeProducer/interface/MBTreeEvent.h"

#include "TFile.h"
#include "TTree.h"

class MBTreeProducer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
  public:
    explicit MBTreeProducer( const edm::ParameterSet& );
    ~MBTreeProducer();

    static void fillDescriptions( edm::ConfigurationDescriptions& );


  private:
    virtual void analyze( const edm::Event&, const edm::EventSetup& ) override;

    // ----------member data ---------------------------

    edm::EDGetTokenT<edm::View<CTPPSLocalTrackLite> > protonTracksToken_;
    edm::EDGetTokenT<edm::View<reco::Vertex> > vtxToken_;
    edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;

    std::string filename_;

    std::unique_ptr<CTPPSAlCa::FillNumberLUTHandler> fillLUTHandler_;

    TFile* file_;
    TTree* tree_;
    MBTreeEvent ev_;
};

MBTreeProducer::MBTreeProducer( const edm::ParameterSet& iConfig ) :
  protonTracksToken_( consumes<edm::View<CTPPSLocalTrackLite> >( iConfig.getParameter<edm::InputTag>( "protonTracksLabel") ) ),
  vtxToken_         ( consumes<edm::View<reco::Vertex> >       ( iConfig.getParameter<edm::InputTag>( "vertexLabel" ) ) ),
  beamSpotToken_    ( consumes<reco::BeamSpot>                 ( iConfig.getParameter<edm::InputTag>( "beamSpotLabel" ) ) ),
  filename_         ( iConfig.getParameter<std::string>( "outputFilename" ) ),
  fillLUTHandler_ ( new CTPPSAlCa::FillNumberLUTHandler( iConfig.getParameter<edm::FileInPath>( "fillNumLUTFile" ).fullPath().c_str() ) ),
  file_( nullptr ), tree_( nullptr )
{
  file_ = new TFile( filename_.c_str(), "recreate" );
  file_->cd();

  tree_ = new TTree( "ntp", "mb ntuple" );
  ev_.create( tree_ );
}


MBTreeProducer::~MBTreeProducer()
{
  if ( file_ ) {
    file_->Write();
    file_->Close();
    delete file_;
  }
  if ( tree_ ) delete tree_;
}

void
MBTreeProducer::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  ev_.clear();

  // Run and BX information
  ev_.bunch_crossing = iEvent.bunchCrossing();
  ev_.run_id = iEvent.id().run();
  ev_.lumisection = iEvent.luminosityBlock();
  ev_.event_number = iEvent.id().event();

  // get the fill number from the run id <-> fill number LUT
  ev_.fill_number = ( fillLUTHandler_ ) ? fillLUTHandler_->getFillNumber( iEvent.id().run() ) : 1;

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

  //----- vertexing information -----

  edm::Handle<edm::View<reco::Vertex> > vertices;
  iEvent.getByToken( vtxToken_, vertices );

//std::cout << ">> " << vertices->size() << std::endl;
  ev_.num_vertex = 0;
  // at least one vertex is required
  if ( vertices->size() < 1 ) return;

  for ( unsigned int i = 0; i < vertices->size() && i < ev_.MAX_VERTEX; ++i ) {
    const edm::Ptr<reco::Vertex> vtx = vertices->ptrAt( i );
    if ( ev_.num_vertex >= ev_.MAX_VERTEX ) {
      edm::LogError("MBTreeProducer") << "More vertices than expected in this event (" << ev_.num_vertex << ">MAX_VERTEX=" << ev_.MAX_VERTEX << "). Increase MAX_VERTEX for safety";
    }
    if ( !vtx->isValid() ) continue;

    ev_.vertex_x[ev_.num_vertex] = vtx->x();
    ev_.vertex_y[ev_.num_vertex] = vtx->y();
    ev_.vertex_z[ev_.num_vertex] = vtx->z();

    ev_.vertex_tracks[ev_.num_vertex] = vtx->nTracks();
    ev_.vertex_tracks_wgt0p75[ev_.num_vertex] = vtx->nTracks( 0.75 );
    ev_.vertex_tracks_wgt0p90[ev_.num_vertex] = vtx->nTracks( 0.90 );
    ev_.vertex_tracks_wgt0p95[ev_.num_vertex] = vtx->nTracks( 0.95 );
    ev_.num_vertex++;
  }

  //----- forward RP tracks -----

  edm::Handle<edm::View<CTPPSLocalTrackLite> > protonTracks;
  iEvent.getByToken( protonTracksToken_, protonTracks );

  ev_.num_fwd_track = 0;
  for ( const auto& trk : *protonTracks) {
    const CTPPSDetId detid( trk.getRPId() );

    ev_.fwd_track_x[ev_.num_fwd_track] = trk.getX() * 1.e-3; // store in m
    ev_.fwd_track_y[ev_.num_fwd_track] = trk.getY() * 1.e-3; // store in m
    ev_.fwd_track_arm[ev_.num_fwd_track] = detid.arm(); // 0 = left (45) ; 1 = right (56)
    ev_.fwd_track_station[ev_.num_fwd_track] = detid.station();
    ev_.fwd_track_pot[ev_.num_fwd_track] = detid.rp(); // 2 = 210m ; 3 = 220m

    //ev_.fwd_track_chi2[ev_.num_fwd_track] = trk.getChiSquared();
    //ev_.fwd_track_normchi2[ev_.num_fwd_track] = trk.getChiSquaredOverNDF();

    ev_.num_fwd_track++;
  }
  if ( ev_.num_fwd_track == 0 ) return; // do not store in the tree if no valid tracks are found

  tree_->Fill();
}

void
MBTreeProducer::fillDescriptions( edm::ConfigurationDescriptions& descriptions )
{
  edm::ParameterSetDescription desc;

  //--- general parameters
  desc.add<std::string>( "outputFilename", "output.root" );
  //--- input collections
  desc.add<edm::InputTag>( "vertexLabel", edm::InputTag( "offlinePrimaryVertices" ) );
  desc.add<edm::InputTag>( "beamSpotLabel", edm::InputTag( "offlineBeamSpot" ) );
  //--- totem RP information extraction
  desc.add<edm::InputTag>( "protonTracksLabel", edm::InputTag( "ctppsLocalTrackLiteProducer" ) );
  desc.add<edm::FileInPath>( "fillNumLUTFile", edm::FileInPath( "DiphotonAnalyzer/TreeProducer/data/fill_run_lut_v2.dat" ) );

  descriptions.add( "mbTreeProducer", desc );
}

//define this as a plug-in
DEFINE_FWK_MODULE( MBTreeProducer );

