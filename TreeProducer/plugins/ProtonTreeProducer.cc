// -*- C++ -*-
//
// Package:    DiphotonAnalyzer/TreeProducer
// Class:      ProtonTreeProducer
//
/**\class ProtonTreeProducer ProtonTreeProducer.cc DiphotonAnalyzer/TreeProducer/plugins/ProtonTreeProducer.cc

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
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"

//#include "DiphotonAnalyzer/TreeProducer/interface/SelectionUtils.h"
#include "DiphotonAnalyzer/TreeProducer/interface/FillNumberLUTHandler.h"
#include "DiphotonAnalyzer/TreeProducer/interface/ProtonInfoEvent.h"

#include "TTree.h"

class ProtonTreeProducer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
  public:
    explicit ProtonTreeProducer( const edm::ParameterSet& );
    static void fillDescriptions( edm::ConfigurationDescriptions& );

  private:
    virtual void analyze( const edm::Event&, const edm::EventSetup& ) override;

    // ----------member data ---------------------------

    edm::EDGetTokenT<std::vector<CTPPSLocalTrackLite> > protonTracksToken_;
    std::string filename_;
    std::unique_ptr<CTPPSAlCa::FillNumberLUTHandler> fillLUTHandler_;

    TTree* tree_;
    ProtonInfoEvent ev_;
};

ProtonTreeProducer::ProtonTreeProducer( const edm::ParameterSet& iConfig ) :
  protonTracksToken_( consumes<std::vector<CTPPSLocalTrackLite> >( iConfig.getParameter<edm::InputTag>( "protonTracksLabel") ) ),
  filename_         ( iConfig.getParameter<std::string>( "outputFilename" ) ),
  fillLUTHandler_   ( new CTPPSAlCa::FillNumberLUTHandler( iConfig.getParameter<edm::FileInPath>( "fillNumLUTFile" ).fullPath().c_str() ) ),
  tree_( nullptr )
{
  // Book the output tree
  usesResource( "TFileService" );
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>( "ntp", "Protons ntuple" );
  ev_.create( tree_ );
}

void
ProtonTreeProducer::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  ev_.clear();

  // Run and BX information
  ev_.bunch_crossing = iEvent.bunchCrossing();
  ev_.run_id = iEvent.id().run();
  ev_.lumisection = iEvent.luminosityBlock();
  ev_.event_number = iEvent.id().event();
  // get the fill number from the run id <-> fill number LUT
  ev_.fill_number = ( fillLUTHandler_ ) ? fillLUTHandler_->getFillNumber( iEvent.id().run() ) : 0;

  //----- forward RP tracks -----

  edm::Handle<std::vector<CTPPSLocalTrackLite> > protonTracks;
  iEvent.getByToken( protonTracksToken_, protonTracks );

  ev_.num_fwd_track = 0;
  for ( const auto& trk : *protonTracks) {
    const CTPPSDetId detid( trk.getRPId() );

    ev_.fwd_track_x[ev_.num_fwd_track] = trk.getX() * 1.e-3; // store in m
    ev_.fwd_track_y[ev_.num_fwd_track] = trk.getY() * 1.e-3; // store in m
    ev_.fwd_track_x_unc[ev_.num_fwd_track] = trk.getXUnc() * 1.e-3; // store in m
    ev_.fwd_track_y_unc[ev_.num_fwd_track] = trk.getYUnc() * 1.e-3; // store in m
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
ProtonTreeProducer::fillDescriptions( edm::ConfigurationDescriptions& descriptions )
{
  edm::ParameterSetDescription desc;

  //--- general parameters
  desc.add<std::string>( "outputFilename", "output.root" );
  //--- RP input collections
  desc.add<edm::InputTag>( "protonTracksLabel", edm::InputTag( "ctppsLocalTrackLiteProducer" ) );
  desc.add<edm::FileInPath>( "fillNumLUTFile", edm::FileInPath( "DiphotonAnalyzer/TreeProducer/data/fill_run_lut_v2.dat" ) );

  descriptions.add( "protonTreeProducer", desc );
}

//define this as a plug-in
DEFINE_FWK_MODULE( ProtonTreeProducer );
