#include "pot_alignment.h"
#include "xi_reconstruction.h"

#include "DiphotonAnalyzer/TreeProducer/interface/TreeEvent.h"

#include "TFile.h"
#include "TTree.h"

void test_xireco()
{
  TFile f( "Samples/output_Run2016BCG_looseCuts_28jun.root" );
  TTree* tr = dynamic_cast<TTree*>( f.Get( "ntp" ) );

  pot_align::load_file( "TreeProducer/data/alignment_collection_v2.out" );
  xi_reco::load_optics_file( "TreeProducer/data/optics_17may22.root" );

  TreeEvent ev;
  ev.attach( tr, true );

  const unsigned long long num_events = tr->GetEntriesFast();
  for ( unsigned long long i = 0; i < num_events; ++i ) {
    tr->GetEntry( i );

    // first loop to identify the tracks and their respective pot

    auto align = pot_align::get_alignments( ev.fill_number );

    for ( unsigned short j = 0; j < ev.num_proton_track; ++j ) {
      const unsigned short pot_id = 100*ev.proton_track_side[j]+ev.proton_track_pot[j];
      const auto& al = align[pot_id];

      //----- reconstruct the kinematics
      double xi = 0., xi_err = 0.;
      xi_reco::reconstruct( ev.proton_track_x[j]+al.x, ev.proton_track_side[j], ev.proton_track_pot[j], xi, xi_err );
    }

  } // loop on events
}
