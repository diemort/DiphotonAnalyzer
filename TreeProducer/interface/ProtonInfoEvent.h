#ifndef DiphotonAnalyzer_TreeProducer_ProtonInfoEvent_h
#define DiphotonAnalyzer_TreeProducer_ProtonInfoEvent_h

#include "TTree.h"

struct ProtonInfoEvent
{
  static constexpr unsigned short MAX_PROTON_TRK = 50;

  void create( TTree* tree ) {
    if ( !tree ) return;

    tree->Branch( "run_id", &run_id, "run_id/i" );
    tree->Branch( "fill_number", &fill_number, "fill_number/i" );
    tree->Branch( "lumisection", &lumisection, "lumisection/i" );
    tree->Branch( "bunch_crossing", &bunch_crossing, "bunch_crossing/i" );
    tree->Branch( "event_number", &event_number, "event_number/l" );

    tree->Branch( "num_fwd_track", &num_fwd_track, "num_fwd_track/i" );
    tree->Branch( "fwd_track_x", fwd_track_x, "fwd_track_x[num_fwd_track]/F" );
    tree->Branch( "fwd_track_y", fwd_track_y, "fwd_track_y[num_fwd_track]/F" );
    tree->Branch( "fwd_track_x_unc", fwd_track_x_unc, "fwd_track_x_unc[num_fwd_track]/F" );
    tree->Branch( "fwd_track_y_unc", fwd_track_y_unc, "fwd_track_y_unc[num_fwd_track]/F" );
    tree->Branch( "fwd_track_arm", fwd_track_arm, "fwd_track_arm[num_fwd_track]/i" );
    tree->Branch( "fwd_track_station", fwd_track_station, "fwd_track_station[num_fwd_track]/i" );
    tree->Branch( "fwd_track_pot", fwd_track_pot, "fwd_track_pot[num_fwd_track]/i" );
    //tree->Branch( "fwd_track_chi2", fwd_track_chi2, "fwd_track_chi2[num_fwd_track]/F" );
    //tree->Branch( "fwd_track_normchi2", fwd_track_normchi2, "fwd_track_normchi2[num_fwd_track]/F" );
  }

  void attach( TTree* tree, const std::vector<std::string>& list = {} ) {
    if ( !tree ) return;

    tree->SetBranchAddress( "run_id", &run_id );
    tree->SetBranchAddress( "fill_number", &fill_number );
    tree->SetBranchAddress( "lumisection", &lumisection );
    tree->SetBranchAddress( "bunch_crossing", &bunch_crossing );
    tree->SetBranchAddress( "event_number", &event_number );

    tree->SetBranchAddress( "num_fwd_track", &num_fwd_track );
    tree->SetBranchAddress( "fwd_track_x", fwd_track_x );
    tree->SetBranchAddress( "fwd_track_y", fwd_track_y );
    tree->SetBranchAddress( "fwd_track_x_unc", fwd_track_x_unc );
    tree->SetBranchAddress( "fwd_track_y_unc", fwd_track_y_unc );
    tree->SetBranchAddress( "fwd_track_arm", fwd_track_arm );
    tree->SetBranchAddress( "fwd_track_station", fwd_track_station );
    tree->SetBranchAddress( "fwd_track_pot", fwd_track_pot );
    //tree->SetBranchAddress( "fwd_track_chi2", fwd_track_chi2 );
    //tree->SetBranchAddress( "fwd_track_normchi2", fwd_track_normchi2 );

    if ( list.empty() )
      return;
    tree->SetBranchStatus( "*", 0 );
    for ( const auto& br : list )
      tree->SetBranchStatus( br.c_str(), 1 );
  }

  void clear() {
    bunch_crossing = run_id = fill_number = lumisection = event_number = 0;
    num_fwd_track = 0;
    for ( unsigned int i = 0; i < MAX_PROTON_TRK; ++i ) {
      fwd_track_x[i] = fwd_track_y[i] = -1.;
      fwd_track_x_unc[i] = fwd_track_y_unc[i] = -1.;
      //fwd_track_chi2[i] = fwd_track_normchi2[i] = -1.;
      fwd_track_arm[i] = 2; //invalid
      fwd_track_station[i] = fwd_track_pot[i] = 0; //invalid
    }
  }

  // --- tree components ---

  unsigned int bunch_crossing, fill_number, run_id, lumisection;
  unsigned long long event_number;

  unsigned int num_fwd_track;
  float fwd_track_x[MAX_PROTON_TRK], fwd_track_y[MAX_PROTON_TRK];
  float fwd_track_x_unc[MAX_PROTON_TRK], fwd_track_y_unc[MAX_PROTON_TRK];
  //float fwd_track_chi2[MAX_PROTON_TRK], fwd_track_normchi2[MAX_PROTON_TRK];
  unsigned int fwd_track_arm[MAX_PROTON_TRK], fwd_track_station[MAX_PROTON_TRK], fwd_track_pot[MAX_PROTON_TRK];
};

#endif
