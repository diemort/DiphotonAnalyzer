#ifndef DiphotonAnalyzer_TreeProducer_MBTreeEvent_h
#define DiphotonAnalyzer_TreeProducer_MBTreeEvent_h

#include "TTree.h"

struct MBTreeEvent
{
  static constexpr unsigned short MAX_PROTON_TRK = 20;
  static constexpr unsigned short MAX_VERTEX = 100;

  void create( TTree* tree ) {
    if ( !tree ) return;

    tree->Branch( "run_id", &run_id, "run_id/i" );
    tree->Branch( "fill_number", &fill_number, "fill_number/i" );
    tree->Branch( "lumisection", &lumisection, "lumisection/i" );
    tree->Branch( "bunch_crossing", &bunch_crossing, "bunch_crossing/i" );
    tree->Branch( "event_number", &event_number, "event_number/l" );

    tree->Branch( "num_strips_track", &num_strips_track, "num_strips_track/i" );
    tree->Branch( "strips_track_x", strips_track_x, "strips_track_x[num_strips_track]/F" );
    tree->Branch( "strips_track_y", strips_track_y, "strips_track_y[num_strips_track]/F" );
    tree->Branch( "strips_track_tx", strips_track_tx, "strips_track_tx[num_strips_track]/F" );
    tree->Branch( "strips_track_ty", strips_track_ty, "strips_track_ty[num_strips_track]/F" );
    tree->Branch( "strips_track_arm", strips_track_arm, "strips_track_arm[num_strips_track]/i" );
    tree->Branch( "strips_track_station", strips_track_station, "strips_track_station[num_strips_track]/i" );
    tree->Branch( "strips_track_pot", strips_track_pot, "strips_track_pot[num_strips_track]/i" );
    tree->Branch( "strips_track_chi2", strips_track_chi2, "strips_track_chi2[num_strips_track]/F" );
    tree->Branch( "strips_track_normchi2", strips_track_normchi2, "strips_track_normchi2[num_strips_track]/F" );

    tree->Branch( "num_vertex", &num_vertex, "num_vertex/i" );
    tree->Branch( "vertex_x", vertex_x, "vertex_x[num_vertex]/F" );
    tree->Branch( "vertex_y", vertex_y, "vertex_y[num_vertex]/F" );
    tree->Branch( "vertex_z", vertex_z, "vertex_z[num_vertex]/F" );
    tree->Branch( "vertex_tracks", vertex_tracks, "vertex_tracks[num_vertex]/i" );
    tree->Branch( "vertex_tracks_wgt0p75", vertex_tracks_wgt0p75, "vertex_tracks_wgt0p75[num_vertex]/i" );
    tree->Branch( "vertex_tracks_wgt0p90", vertex_tracks_wgt0p90, "vertex_tracks_wgt0p90[num_vertex]/i" );
    tree->Branch( "vertex_tracks_wgt0p95", vertex_tracks_wgt0p95, "vertex_tracks_wgt0p95[num_vertex]/i" );

    tree->Branch( "bs_x0", &bs_x0, "bs_x0/F" );
    tree->Branch( "bs_y0", &bs_y0, "bs_y0/F" );
    tree->Branch( "bs_z0", &bs_z0, "bs_z0/F" );
    tree->Branch( "bs_sigma_z", &bs_sigma_z, "bs_sigma_z/F" );
    tree->Branch( "bs_dxdz", &bs_dxdz, "bs_dxdz/F" );
    tree->Branch( "bs_beam_width_x", &bs_beam_width_x, "bs_beam_width_x/F" );
    tree->Branch( "bs_beam_width_y", &bs_beam_width_y, "bs_beam_width_y/F" );
  }

  void attach( TTree* tree ) {
    if ( !tree ) return;

    tree->SetBranchAddress( "run_id", &run_id );
    tree->SetBranchAddress( "fill_number", &fill_number );
    tree->SetBranchAddress( "lumisection", &lumisection );
    tree->SetBranchAddress( "bunch_crossing", &bunch_crossing );
    tree->SetBranchAddress( "event_number", &event_number);

    tree->SetBranchAddress( "num_strips_track", &num_strips_track );
    tree->SetBranchAddress( "strips_track_x", strips_track_x );
    tree->SetBranchAddress( "strips_track_y", strips_track_y );
    tree->SetBranchAddress( "strips_track_tx", strips_track_tx );
    tree->SetBranchAddress( "strips_track_ty", strips_track_ty );
    tree->SetBranchAddress( "strips_track_arm", strips_track_arm );
    tree->SetBranchAddress( "strips_track_station", strips_track_station );
    tree->SetBranchAddress( "strips_track_pot", strips_track_pot );
    tree->SetBranchAddress( "strips_track_chi2", strips_track_chi2 );
    tree->SetBranchAddress( "strips_track_normchi2", strips_track_normchi2 );

    tree->SetBranchAddress( "num_vertex", &num_vertex );
    tree->SetBranchAddress( "vertex_x", vertex_x );
    tree->SetBranchAddress( "vertex_y", vertex_y );
    tree->SetBranchAddress( "vertex_z", vertex_z );
    tree->SetBranchAddress( "vertex_tracks", vertex_tracks );
    tree->SetBranchAddress( "vertex_tracks_wgt0p75", vertex_tracks_wgt0p75 );
    tree->SetBranchAddress( "vertex_tracks_wgt0p90", vertex_tracks_wgt0p90 );
    tree->SetBranchAddress( "vertex_tracks_wgt0p95", vertex_tracks_wgt0p95 );

    tree->SetBranchAddress( "bs_x0", &bs_x0 );
    tree->SetBranchAddress( "bs_y0", &bs_y0 );
    tree->SetBranchAddress( "bs_z0", &bs_z0 );
    tree->SetBranchAddress( "bs_sigma_z", &bs_sigma_z );
    tree->SetBranchAddress( "bs_dxdz", &bs_dxdz );
    tree->SetBranchAddress( "bs_beam_width_x", &bs_beam_width_x );
    tree->SetBranchAddress( "bs_beam_width_y", &bs_beam_width_y );
  }

  void clear() {
    bunch_crossing = run_id = fill_number = lumisection = event_number = 0;
    num_strips_track = 0;
    for ( unsigned int i = 0; i < MAX_PROTON_TRK; ++i ) {
      strips_track_x[i] = strips_track_y[i] = -1.;
      strips_track_tx[i] = strips_track_ty[i] = -1.;
      strips_track_chi2[i] = strips_track_normchi2[i] = -1.;
      strips_track_arm[i] = 2; //invalid
      strips_track_station[i] = strips_track_pot[i] = 0; //invalid
    }

    num_vertex = 0;
    for ( unsigned int i = 0; i < MAX_VERTEX; ++i ) {
      vertex_x[i] = vertex_y[i] = vertex_z[i] = -999.;
      vertex_tracks[i] = vertex_tracks_wgt0p75[i] = vertex_tracks_wgt0p90[i] = vertex_tracks_wgt0p95[i] = 0;
    }
    bs_x0 = bs_y0 = bs_z0 = bs_sigma_z = bs_dxdz = bs_beam_width_x = bs_beam_width_y = -999.;
  }

  // --- tree components ---

  unsigned int bunch_crossing, fill_number, run_id, lumisection;
  unsigned long long event_number;

  unsigned int num_strips_track;
  float strips_track_x[MAX_PROTON_TRK], strips_track_y[MAX_PROTON_TRK];
  float strips_track_tx[MAX_PROTON_TRK], strips_track_ty[MAX_PROTON_TRK];
  float strips_track_chi2[MAX_PROTON_TRK], strips_track_normchi2[MAX_PROTON_TRK];
  unsigned int strips_track_arm[MAX_PROTON_TRK], strips_track_station[MAX_PROTON_TRK], strips_track_pot[MAX_PROTON_TRK];

  unsigned int num_vertex;
  float vertex_x[MAX_VERTEX], vertex_y[MAX_VERTEX], vertex_z[MAX_VERTEX];
  unsigned int vertex_tracks[MAX_VERTEX];
  unsigned int vertex_tracks_wgt0p75[MAX_VERTEX], vertex_tracks_wgt0p90[MAX_VERTEX], vertex_tracks_wgt0p95[MAX_VERTEX];
  float bs_x0, bs_y0, bs_z0, bs_sigma_z, bs_dxdz, bs_beam_width_x, bs_beam_width_y;
};

#endif
