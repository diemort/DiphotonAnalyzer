#ifndef DiphotonAnalyzer_TreeProducer_TreeEvent_h
#define DiphotonAnalyzer_TreeProducer_TreeEvent_h

#include "TTree.h"
#include <string>

struct TreeEvent
{
  static constexpr unsigned short MAX_HLT = 10;
  static constexpr unsigned short MAX_PROTON_TRK = 20;
  static constexpr unsigned short MAX_DIPHOTON = 10;
  //                               JW
  static constexpr unsigned short MAX_ELECTRON = 20;
  static constexpr unsigned short MAX_MUON = 20;
  static constexpr unsigned short MAX_JET = 50;
  //
  static constexpr unsigned short MAX_VERTEX = 100;
  static constexpr unsigned short MAX_GEN_PHOTON = 10;
  static constexpr unsigned short MAX_GEN_PART = 20;

  typedef enum { all, ebeb, ebee, eeee, invalid, num_classes } classes_;
  static constexpr const char* classes[num_classes] = { "", "EBEB", "EBEE", "EEEE", "invalid" };

  void create( TTree* tree, bool data = false ) {
    if ( !tree ) return;

    tree->Branch( "run_id", &run_id, "run_id/i" );
    tree->Branch( "fill_number", &fill_number, "fill_number/i" );
    tree->Branch( "lumisection", &lumisection, "lumisection/i" );
    tree->Branch( "bunch_crossing", &bunch_crossing, "bunch_crossing/i" );
    tree->Branch( "event_number", &event_number, "event_number/l" );

    tree->Branch( "num_hlt", &num_hlt, "num_hlt/I" );
    tree->Branch( "hlt_accept", hlt_accept, "hlt_accept[num_hlt]/I" );
    tree->Branch( "hlt_prescale", hlt_prescale, "hlt_prescale[num_hlt]/I" );

    /*std::vector<std::string>* HLT_Name;
    tree->Branch( "hlt_name", &HLT_Name );
    *HLT_Name = triggersList_;*/

    if ( data ) {
      tree->Branch( "num_fwd_track", &num_fwd_track, "num_fwd_track/i" );
      tree->Branch( "fwd_track_x", fwd_track_x, "fwd_track_x[num_fwd_track]/F" );
      tree->Branch( "fwd_track_y", fwd_track_y, "fwd_track_y[num_fwd_track]/F" );
      tree->Branch( "fwd_track_arm", fwd_track_arm, "fwd_track_arm[num_fwd_track]/i" );
      tree->Branch( "fwd_track_station", fwd_track_station, "fwd_track_station[num_fwd_track]/i" );
      tree->Branch( "fwd_track_pot", fwd_track_pot, "fwd_track_pot[num_fwd_track]/i" );
      //tree->Branch( "fwd_track_chi2", fwd_track_chi2, "fwd_track_chi2[num_fwd_track]/F" );
      //tree->Branch( "fwd_track_normchi2", fwd_track_normchi2, "fwd_track_normchi2[num_fwd_track]/F" );
      tree->Branch( "num_proton_track", &num_proton_track, "num_proton_track/i" );
      tree->Branch( "proton_track_vx", proton_track_vx, "proton_track_vx[num_proton_track]/F" );
      tree->Branch( "proton_track_vy", proton_track_vy, "proton_track_vy[num_proton_track]/F" );
      tree->Branch( "proton_track_xi", proton_track_xi, "proton_track_xi[num_proton_track]/F" );
      tree->Branch( "proton_track_arm", proton_track_arm, "proton_track_arm[num_proton_track]/i" );
      tree->Branch( "proton_track_pot", proton_track_pot, "proton_track_pot[num_proton_track]/i" );
      tree->Branch( "proton_track_method", proton_track_method, "proton_track_method[num_proton_track]/i" );
      tree->Branch( "proton_track_chi2", proton_track_chi2, "proton_track_chi2[num_proton_track]/F" );
      tree->Branch( "proton_track_normchi2", proton_track_normchi2, "proton_track_normchi2[num_proton_track]/F" );
    }
    if ( !data ) {
      tree->Branch( "num_gen_photon", &num_gen_photon, "num_gen_photon/i" );
      tree->Branch( "gen_photon_pt", gen_photon_pt, "gen_photon_pt[num_gen_photon]/F" );
      tree->Branch( "gen_photon_eta", gen_photon_eta, "gen_photon_eta[num_gen_photon]/F" );
      tree->Branch( "gen_photon_phi", gen_photon_phi, "gen_photon_phi[num_gen_photon]/F" );
      tree->Branch( "gen_photon_energy", gen_photon_energy, "gen_photon_energy[num_gen_photon]/F" );
      tree->Branch( "gen_photon_vertex_x", gen_photon_vertex_x, "gen_photon_vertex_x[num_gen_photon]/F" );
      tree->Branch( "gen_photon_vertex_y", gen_photon_vertex_y, "gen_photon_vertex_y[num_gen_photon]/F" );
      tree->Branch( "gen_photon_vertex_z", gen_photon_vertex_z, "gen_photon_vertex_z[num_gen_photon]/F" );

      tree->Branch( "num_gen_part", &num_gen_part, "num_gen_part/i" );
      tree->Branch( "gen_part_pdgid", gen_part_pdgid, "gen_part_pdgid[num_gen_part]/I" );
      tree->Branch( "gen_part_status", gen_part_status, "gen_part_status[num_gen_part]/I" );
      tree->Branch( "gen_part_pt", gen_part_pt, "gen_part_pt[num_gen_part]/F" );
      tree->Branch( "gen_part_eta", gen_part_eta, "gen_part_eta[num_gen_part]/F" );
      tree->Branch( "gen_part_phi", gen_part_phi, "gen_part_phi[num_gen_part]/F" );
      tree->Branch( "gen_part_energy", gen_part_energy, "gen_part_energy[num_gen_part]/F" );
      tree->Branch( "gen_part_vertex_x", gen_part_vertex_x, "gen_part_vertex_x[num_gen_part]/F" );
      tree->Branch( "gen_part_vertex_y", gen_part_vertex_y, "gen_part_vertex_y[num_gen_part]/F" );
      tree->Branch( "gen_part_vertex_z", gen_part_vertex_z, "gen_part_vertex_z[num_gen_part]/F" );
    }

    tree->Branch( "num_diphoton", &num_diphoton, "num_diphoton/i" );
    tree->Branch( "diphoton_pt1", diphoton_pt1, "diphoton_pt1[num_diphoton]/F" );
    tree->Branch( "diphoton_pt2", diphoton_pt2, "diphoton_pt2[num_diphoton]/F" );
    tree->Branch( "diphoton_eta1", diphoton_eta1, "diphoton_eta1[num_diphoton]/F" );
    tree->Branch( "diphoton_eta2", diphoton_eta2, "diphoton_eta2[num_diphoton]/F" );
    tree->Branch( "diphoton_phi1", diphoton_phi1, "diphoton_phi1[num_diphoton]/F" );
    tree->Branch( "diphoton_phi2", diphoton_phi2, "diphoton_phi2[num_diphoton]/F" );
    tree->Branch( "diphoton_energy1", diphoton_energy1, "diphoton_energy1[num_diphoton]/F" );
    tree->Branch( "diphoton_energy2", diphoton_energy2, "diphoton_energy2[num_diphoton]/F" );
    tree->Branch( "diphoton_r91", diphoton_r91, "diphoton_r91[num_diphoton]/F" );
    tree->Branch( "diphoton_r92", diphoton_r92, "diphoton_r92[num_diphoton]/F" );
    tree->Branch( "diphoton_eleveto1", diphoton_ele_veto1, "diphoton_eleveto1[num_diphoton]/i" );
    tree->Branch( "diphoton_eleveto2", diphoton_ele_veto2, "diphoton_eleveto2[num_diphoton]/i" );
    tree->Branch( "diphoton_mass", diphoton_mass, "diphoton_mass[num_diphoton]/F" );
    tree->Branch( "diphoton_rapidity", diphoton_rapidity, "diphoton_rapidity[num_diphoton]/F" );
    tree->Branch( "diphoton_pt", diphoton_pt, "diphoton_pt[num_diphoton]/F" );
    tree->Branch( "diphoton_dphi", diphoton_dphi, "diphoton_dphi[num_diphoton]/F" );
    tree->Branch( "diphoton_id1", diphoton_id1, "diphoton_id1[num_diphoton]/F" );
    tree->Branch( "diphoton_id2", diphoton_id2, "diphoton_id2[num_diphoton]/F" );
    tree->Branch( "diphoton_sigeove1", diphoton_sigeove1, "diphoton_sigeove1[num_diphoton]/F" );
    tree->Branch( "diphoton_sigeove2", diphoton_sigeove2, "diphoton_sigeove2[num_diphoton]/F" );

    tree->Branch( "diphoton_supercluster_x1", diphoton_supercluster_x1, "diphoton_supercluster_x1[num_diphoton]/F" );
    tree->Branch( "diphoton_supercluster_y1", diphoton_supercluster_y1, "diphoton_supercluster_y1[num_diphoton]/F" );
    tree->Branch( "diphoton_supercluster_z1", diphoton_supercluster_z1, "diphoton_supercluster_z1[num_diphoton]/F" );
    tree->Branch( "diphoton_supercluster_x2", diphoton_supercluster_x2, "diphoton_supercluster_x2[num_diphoton]/F" );
    tree->Branch( "diphoton_supercluster_y2", diphoton_supercluster_y2, "diphoton_supercluster_y2[num_diphoton]/F" );
    tree->Branch( "diphoton_supercluster_z2", diphoton_supercluster_z2, "diphoton_supercluster_z2[num_diphoton]/F" );

    if ( !data ) {
      tree->Branch( "diphoton_genpt1", diphoton_genpt1, "diphoton_genpt1[num_diphoton]/F" );
      tree->Branch( "diphoton_genpt2", diphoton_genpt2, "diphoton_genpt2[num_diphoton]/F" );
      tree->Branch( "diphoton_geneta1", diphoton_geneta1, "diphoton_geneta1[num_diphoton]/F" );
      tree->Branch( "diphoton_geneta2", diphoton_geneta2, "diphoton_geneta2[num_diphoton]/F" );
      tree->Branch( "diphoton_genphi1", diphoton_genphi1, "diphoton_genphi1[num_diphoton]/F" );
      tree->Branch( "diphoton_genphi2", diphoton_genphi2, "diphoton_genphi2[num_diphoton]/F" );
      tree->Branch( "diphoton_genenergy1", diphoton_genenergy1, "diphoton_genenergy1[num_diphoton]/F" );
      tree->Branch( "diphoton_genenergy2", diphoton_genenergy2, "diphoton_genenergy2[num_diphoton]/F" );
    }

    //tree->Branch( "diphoton_vertex_tracks", diphoton_vertex_tracks, "diphoton_vertex_tracks[num_diphoton]/i" );
    tree->Branch( "diphoton_vertex_id", diphoton_vertex_id, "diphoton_vertex_id[num_diphoton]/I" );
    tree->Branch( "diphoton_vertex_x", diphoton_vertex_x, "diphoton_vertex_x[num_diphoton]/F" );
    tree->Branch( "diphoton_vertex_y", diphoton_vertex_y, "diphoton_vertex_y[num_diphoton]/F" );
    tree->Branch( "diphoton_vertex_z", diphoton_vertex_z, "diphoton_vertex_z[num_diphoton]/F" );
    tree->Branch( "diphoton_vertex_nearestvtxdist", diphoton_nearestvtxdist, "diphoton_vertex_nearestvtxdist[num_diphoton]/F" );
    tree->Branch( "diphoton_vertex_vtx1mmdist", diphoton_vertex_vtx1mmdist, "diphoton_vertex_vtx1mmdist[num_diphoton]/i" );
    tree->Branch( "diphoton_vertex_vtx2mmdist", diphoton_vertex_vtx2mmdist, "diphoton_vertex_vtx2mmdist[num_diphoton]/i" );
    tree->Branch( "diphoton_vertex_vtx5mmdist", diphoton_vertex_vtx5mmdist, "diphoton_vertex_vtx5mmdist[num_diphoton]/i" );
    tree->Branch( "diphoton_vertex_vtx1cmdist", diphoton_vertex_vtx1cmdist, "diphoton_vertex_vtx1cmdist[num_diphoton]/i" );

    if ( !data ) {
      tree->Branch( "diphoton_genvertex_x", &diphoton_genvertex_x, "diphoton_genvertex_x/F" );
      tree->Branch( "diphoton_genvertex_y", &diphoton_genvertex_y, "diphoton_genvertex_y/F" );
      tree->Branch( "diphoton_genvertex_z", &diphoton_genvertex_z, "diphoton_genvertex_z/F" );
      tree->Branch( "diphoton_genvertex_smeared_x", &diphoton_genvertex_smeared_x, "diphoton_genvertex_smeared_x/F" );
      tree->Branch( "diphoton_genvertex_smeared_y", &diphoton_genvertex_smeared_y, "diphoton_genvertex_smeared_y/F" );
      tree->Branch( "diphoton_genvertex_smeared_z", &diphoton_genvertex_smeared_z, "diphoton_genvertex_smeared_z/F" );
    }

    tree->Branch( "num_electron", &num_electron, "num_electron/i" );
    tree->Branch( "electron_pt", electron_pt, "electron_pt[num_electron]/F" );
    tree->Branch( "electron_eta", electron_eta, "electron_eta[num_electron]/F" );
    tree->Branch( "electron_phi", electron_phi, "electron_phi[num_electron]/F" );
    tree->Branch( "electron_energy", electron_energy, "electron_energy[num_electron]/F" );
    tree->Branch( "electron_vtx_x", electron_vtx_x, "electron_vtx_x[num_electron]/F" );
    tree->Branch( "electron_vtx_y", electron_vtx_y, "electron_vtx_y[num_electron]/F" );
    tree->Branch( "electron_vtx_z", electron_vtx_z, "electron_vtx_z[num_electron]/F" );

    tree->Branch( "num_muon", &num_muon, "num_muon/i" );
    tree->Branch( "muon_pt", muon_pt, "muon_pt[num_muon]/F" );
    tree->Branch( "muon_eta", muon_eta, "muon_eta[num_muon]/F" );
    tree->Branch( "muon_phi", muon_phi, "muon_phi[num_muon]/F" );
    tree->Branch( "muon_energy", muon_energy, "muon_energy[num_muon]/F" );
    tree->Branch( "muon_vtx_x", muon_vtx_x, "muon_vtx_x[num_muon]/F" );
    tree->Branch( "muon_vtx_y", muon_vtx_y, "muon_vtx_y[num_muon]/F" );
    tree->Branch( "muon_vtx_z", muon_vtx_z, "muon_vtx_z[num_muon]/F" );

    tree->Branch( "num_jet", &num_jet, "num_jet/i" );
    tree->Branch( "jet_pt", jet_pt, "jet_pt[num_jet]/F" );
    tree->Branch( "jet_eta", jet_eta, "jet_eta[num_jet]/F" );
    tree->Branch( "jet_phi", jet_phi, "jet_phi[num_jet]/F" );
    tree->Branch( "jet_energy", jet_energy, "jet_energy[num_jet]/F" );
    tree->Branch( "jet_mass", jet_mass, "jet_mass[num_jet]/F" );
    tree->Branch( "jet_vtx_x", jet_vtx_x, "jet_vtx_x[num_jet]/F" );
    tree->Branch( "jet_vtx_y", jet_vtx_y, "jet_vtx_y[num_jet]/F" );
    tree->Branch( "jet_vtx_z", jet_vtx_z, "jet_vtx_z[num_jet]/F" );
    tree->Branch( "jet_dipho_match", jet_dipho_match, "jet_dipho_match[num_jet]/I" );

    tree->Branch( "num_vertex", &num_vertex, "num_vertex/i" );
    tree->Branch( "vertex_x", vertex_x, "vertex_x[num_vertex]/F" );
    tree->Branch( "vertex_y", vertex_y, "vertex_y[num_vertex]/F" );
    tree->Branch( "vertex_z", vertex_z, "vertex_z[num_vertex]/F" );
    /*tree->Branch( "vertex_tracks", vertex_tracks, "vertex_tracks[num_vertex]/i" );
    tree->Branch( "vertex_tracks_weight0p75", vertex_tracksWght0p75, "vertex_tracks_weight0p75[num_vertex]/i" );
    tree->Branch( "vertex_tracks_weight0p9", vertex_tracksWght0p9, "vertex_tracks_weight0p9[num_vertex]/i" );
    tree->Branch( "vertex_tracks_weight0p95", vertex_tracksWght0p95, "vertex_tracks_weight0p95[num_vertex]/i" );*/

    tree->Branch( "met", &met, "met/F" );
    tree->Branch( "met_phi", &met_phi, "met_phi/F" );
    tree->Branch( "met_significance", &met_significance, "met_significance/F" );

    tree->Branch( "bs_x0", &bs_x0, "bs_x0/F" );
    tree->Branch( "bs_y0", &bs_y0, "bs_y0/F" );
    tree->Branch( "bs_z0", &bs_z0, "bs_z0/F" );
    tree->Branch( "bs_sigma_z", &bs_sigma_z, "bs_sigma_z/F" );
    tree->Branch( "bs_dxdz", &bs_dxdz, "bs_dxdz/F" );
    tree->Branch( "bs_beam_width_x", &bs_beam_width_x, "bs_beam_width_x/F" );
    tree->Branch( "bs_beam_width_y", &bs_beam_width_y, "bs_beam_width_y/F" );

    tree->Branch( "pileup_weight", &pileup_weight, "pileup_weight/F" );
  }

  void enableBranches( const std::vector<std::string>& list ) {
    if ( list.empty() )
      return;
    tree->SetBranchStatus( "*", 0 );
    for ( const auto& br : list )
      tree->SetBranchStatus( br, 1 );
  }

  void attach( TTree* tree, bool data = false ) {
    if ( !tree ) return;

    tree->SetBranchAddress( "run_id", &run_id );
    tree->SetBranchAddress( "fill_number", &fill_number );
    tree->SetBranchAddress( "lumisection", &lumisection );
    tree->SetBranchAddress( "bunch_crossing", &bunch_crossing );
    tree->SetBranchAddress( "event_number", &event_number );

    tree->SetBranchAddress( "num_hlt", &num_hlt );
    tree->SetBranchAddress( "hlt_accept", hlt_accept );
    tree->SetBranchAddress( "hlt_prescale", hlt_prescale );

    /*std::vector<std::string>* HLT_Name;
    tree->SetBranchAddress( "hlt_name", &HLT_Name );
    *HLT_Name = triggersList_;*/

    if ( data ) {
      tree->SetBranchAddress( "num_fwd_track", &num_fwd_track );
      tree->SetBranchAddress( "fwd_track_x", fwd_track_x );
      tree->SetBranchAddress( "fwd_track_y", fwd_track_y );
      tree->SetBranchAddress( "fwd_track_arm", fwd_track_arm );
      tree->SetBranchAddress( "fwd_track_station", fwd_track_station );
      tree->SetBranchAddress( "fwd_track_pot", fwd_track_pot );
      //tree->SetBranchAddress( "fwd_track_chi2", fwd_track_chi2 );
      //tree->SetBranchAddress( "fwd_track_normchi2", fwd_track_normchi2 );
      tree->SetBranchAddress( "num_proton_track", &num_proton_track );
      tree->SetBranchAddress( "proton_track_vx", proton_track_vx );
      tree->SetBranchAddress( "proton_track_vy", proton_track_vy );
      tree->SetBranchAddress( "proton_track_xi", proton_track_xi );
      tree->SetBranchAddress( "proton_track_arm", proton_track_arm );
      tree->SetBranchAddress( "proton_track_pot", proton_track_pot );
      tree->SetBranchAddress( "proton_track_method", proton_track_method );
      tree->SetBranchAddress( "proton_track_chi2", proton_track_chi2 );
      tree->SetBranchAddress( "proton_track_normchi2", proton_track_normchi2 );
    }
    if ( !data ) {
      tree->SetBranchAddress( "num_gen_photon", &num_gen_photon );
      tree->SetBranchAddress( "gen_photon_pt", gen_photon_pt );
      tree->SetBranchAddress( "gen_photon_eta", gen_photon_eta );
      tree->SetBranchAddress( "gen_photon_phi", gen_photon_phi );
      tree->SetBranchAddress( "gen_photon_energy", gen_photon_energy );
      tree->SetBranchAddress( "gen_photon_vertex_x", gen_photon_vertex_x );
      tree->SetBranchAddress( "gen_photon_vertex_y", gen_photon_vertex_y );
      tree->SetBranchAddress( "gen_photon_vertex_z", gen_photon_vertex_z );

      tree->SetBranchAddress( "num_gen_part", &num_gen_part );
      tree->SetBranchAddress( "gen_part_pdgid", gen_part_pdgid );
      tree->SetBranchAddress( "gen_part_status", gen_part_status );
      tree->SetBranchAddress( "gen_part_pt", gen_part_pt );
      tree->SetBranchAddress( "gen_part_eta", gen_part_eta );
      tree->SetBranchAddress( "gen_part_phi", gen_part_phi );
      tree->SetBranchAddress( "gen_part_energy", gen_part_energy );
      tree->SetBranchAddress( "gen_part_vertex_x", gen_part_vertex_x );
      tree->SetBranchAddress( "gen_part_vertex_y", gen_part_vertex_y );
      tree->SetBranchAddress( "gen_part_vertex_z", gen_part_vertex_z );
    }

    tree->SetBranchAddress( "num_diphoton", &num_diphoton );
    tree->SetBranchAddress( "diphoton_pt1", diphoton_pt1 );
    tree->SetBranchAddress( "diphoton_pt2", diphoton_pt2 );
    tree->SetBranchAddress( "diphoton_eta1", diphoton_eta1 );
    tree->SetBranchAddress( "diphoton_eta2", diphoton_eta2 );
    tree->SetBranchAddress( "diphoton_phi1", diphoton_phi1 );
    tree->SetBranchAddress( "diphoton_phi2", diphoton_phi2 );
    tree->SetBranchAddress( "diphoton_energy1", diphoton_energy1 );
    tree->SetBranchAddress( "diphoton_energy2", diphoton_energy2 );
    tree->SetBranchAddress( "diphoton_r91", diphoton_r91 );
    tree->SetBranchAddress( "diphoton_r92", diphoton_r92 );
    tree->SetBranchAddress( "diphoton_eleveto1", diphoton_ele_veto1 );
    tree->SetBranchAddress( "diphoton_eleveto2", diphoton_ele_veto2 );
    tree->SetBranchAddress( "diphoton_mass", diphoton_mass );
    tree->SetBranchAddress( "diphoton_rapidity", diphoton_rapidity );
    tree->SetBranchAddress( "diphoton_pt", diphoton_pt );
    tree->SetBranchAddress( "diphoton_dphi", diphoton_dphi );
    tree->SetBranchAddress( "diphoton_id1", diphoton_id1 );
    tree->SetBranchAddress( "diphoton_id2", diphoton_id2 );
    tree->SetBranchAddress( "diphoton_sigeove1", diphoton_sigeove1 );
    tree->SetBranchAddress( "diphoton_sigeove2", diphoton_sigeove2 );

    tree->SetBranchAddress( "diphoton_supercluster_x1", diphoton_supercluster_x1 );
    tree->SetBranchAddress( "diphoton_supercluster_y1", diphoton_supercluster_y1 );
    tree->SetBranchAddress( "diphoton_supercluster_z1", diphoton_supercluster_z1 );
    tree->SetBranchAddress( "diphoton_supercluster_x2", diphoton_supercluster_x2 );
    tree->SetBranchAddress( "diphoton_supercluster_y2", diphoton_supercluster_y2 );
    tree->SetBranchAddress( "diphoton_supercluster_z2", diphoton_supercluster_z2 );

    if ( !data ) {
      tree->SetBranchAddress( "diphoton_genpt1", diphoton_genpt1 );
      tree->SetBranchAddress( "diphoton_genpt2", diphoton_genpt2 );
      tree->SetBranchAddress( "diphoton_geneta1", diphoton_geneta1 );
      tree->SetBranchAddress( "diphoton_geneta2", diphoton_geneta2 );
      tree->SetBranchAddress( "diphoton_genphi1", diphoton_genphi1 );
      tree->SetBranchAddress( "diphoton_genphi2", diphoton_genphi2 );
      tree->SetBranchAddress( "diphoton_genenergy1", diphoton_genenergy1 );
      tree->SetBranchAddress( "diphoton_genenergy2", diphoton_genenergy2 );
    }

    tree->SetBranchAddress( "diphoton_vertex_id", diphoton_vertex_id );
    tree->SetBranchAddress( "diphoton_vertex_x", diphoton_vertex_x );
    tree->SetBranchAddress( "diphoton_vertex_y", diphoton_vertex_y );
    tree->SetBranchAddress( "diphoton_vertex_z", diphoton_vertex_z );
    tree->SetBranchAddress( "diphoton_vertex_nearestvtxdist", diphoton_nearestvtxdist );
    tree->SetBranchAddress( "diphoton_vertex_vtx1mmdist", diphoton_vertex_vtx1mmdist );
    tree->SetBranchAddress( "diphoton_vertex_vtx2mmdist", diphoton_vertex_vtx2mmdist );
    tree->SetBranchAddress( "diphoton_vertex_vtx5mmdist", diphoton_vertex_vtx5mmdist );
    tree->SetBranchAddress( "diphoton_vertex_vtx1cmdist", diphoton_vertex_vtx1cmdist );

    if ( !data ) {
      tree->SetBranchAddress( "diphoton_genvertex_x", &diphoton_genvertex_x );
      tree->SetBranchAddress( "diphoton_genvertex_y", &diphoton_genvertex_y );
      tree->SetBranchAddress( "diphoton_genvertex_z", &diphoton_genvertex_z );
      tree->SetBranchAddress( "diphoton_genvertex_smeared_x", &diphoton_genvertex_smeared_x );
      tree->SetBranchAddress( "diphoton_genvertex_smeared_y", &diphoton_genvertex_smeared_y );
      tree->SetBranchAddress( "diphoton_genvertex_smeared_z", &diphoton_genvertex_smeared_z );
    }

    tree->SetBranchAddress( "num_electron", &num_electron );
    tree->SetBranchAddress( "electron_pt", electron_pt );
    tree->SetBranchAddress( "electron_eta", electron_eta );
    tree->SetBranchAddress( "electron_phi", electron_phi );
    tree->SetBranchAddress( "electron_energy", electron_energy );
    tree->SetBranchAddress( "electron_vtx_x", electron_vtx_x );
    tree->SetBranchAddress( "electron_vtx_y", electron_vtx_y );
    tree->SetBranchAddress( "electron_vtx_z", electron_vtx_z );

    tree->SetBranchAddress( "num_muon", &num_muon );
    tree->SetBranchAddress( "muon_pt", muon_pt );
    tree->SetBranchAddress( "muon_eta", muon_eta );
    tree->SetBranchAddress( "muon_phi", muon_phi );
    tree->SetBranchAddress( "muon_energy", muon_energy );
    tree->SetBranchAddress( "muon_vtx_x", muon_vtx_x );
    tree->SetBranchAddress( "muon_vtx_y", muon_vtx_y );
    tree->SetBranchAddress( "muon_vtx_z", muon_vtx_z );

    tree->SetBranchAddress( "num_jet", &num_jet );
    tree->SetBranchAddress( "jet_pt", jet_pt );
    tree->SetBranchAddress( "jet_eta", jet_eta );
    tree->SetBranchAddress( "jet_phi", jet_phi );
    tree->SetBranchAddress( "jet_energy", jet_energy );
    tree->SetBranchAddress( "jet_mass", jet_mass );
    tree->SetBranchAddress( "jet_vtx_x", jet_vtx_x );
    tree->SetBranchAddress( "jet_vtx_y", jet_vtx_y );
    tree->SetBranchAddress( "jet_vtx_z", jet_vtx_z );
    tree->SetBranchAddress( "jet_dipho_match", jet_dipho_match );

    tree->SetBranchAddress( "num_vertex", &num_vertex );
    tree->SetBranchAddress( "vertex_x", vertex_x );
    tree->SetBranchAddress( "vertex_y", vertex_y );
    tree->SetBranchAddress( "vertex_z", vertex_z );

    tree->SetBranchAddress( "met", &met );
    tree->SetBranchAddress( "met_phi", &met_phi );
    tree->SetBranchAddress( "met_significance", &met_significance );

    tree->SetBranchAddress( "bs_x0", &bs_x0 );
    tree->SetBranchAddress( "bs_y0", &bs_y0 );
    tree->SetBranchAddress( "bs_z0", &bs_z0 );
    tree->SetBranchAddress( "bs_sigma_z", &bs_sigma_z );
    tree->SetBranchAddress( "bs_dxdz", &bs_dxdz );
    tree->SetBranchAddress( "bs_beam_width_x", &bs_beam_width_x );
    tree->SetBranchAddress( "bs_beam_width_y", &bs_beam_width_y );

    tree->SetBranchAddress( "pileup_weight", &pileup_weight );
  }

  void clear() {
    bunch_crossing = run_id = fill_number = lumisection = event_number = 0;
    num_hlt = 0;
    for ( unsigned int i=0; i<MAX_HLT; i++ ) {
      hlt_accept[i] = hlt_prescale[i] = -1;
    }
    num_fwd_track = num_proton_track = 0;
    for ( unsigned int i=0; i<MAX_PROTON_TRK; i++ ) {
      fwd_track_x[i] = fwd_track_y[i] = -1.;
      fwd_track_chi2[i] = fwd_track_normchi2[i] = -1.;
      fwd_track_arm[i] = 2; //invalid
      fwd_track_station[i] = fwd_track_pot[i] = 0;
      proton_track_vx[i] = proton_track_vy[i] = proton_track_xi[i] = -1.;
      proton_track_arm[i] = proton_track_method[i] = 2; //invalid
      proton_track_pot[i] = 0; //invalid
      proton_track_chi2[i] = proton_track_normchi2[i] = -1.;
    }

    num_diphoton = 0;
    for ( unsigned int i=0; i<MAX_DIPHOTON; i++ ) {
      diphoton_pt1[i] = diphoton_pt2[i] = -1.;
      diphoton_eta1[i] = diphoton_eta2[i] = -1.;
      diphoton_phi1[i] = diphoton_phi2[i] = -1.;
      diphoton_energy1[i] = diphoton_energy2[i] = -1.;
      diphoton_r91[i] = diphoton_r92[i] = -1.;
      diphoton_ele_veto1[i] = diphoton_ele_veto2[i] = 2; //invalid
      diphoton_id1[i] = diphoton_id2[i] = -1.;
      diphoton_sigeove1[i] = diphoton_sigeove2[i] = -1.;
      diphoton_mass[i] = diphoton_rapidity[i] = diphoton_pt[i] = diphoton_dphi[i] = -1.;
      //diphoton_vertex_tracks[i] = 0;
      diphoton_vertex_id[i] = -1;
      diphoton_vertex_vtx1mmdist[i] = diphoton_vertex_vtx2mmdist[i] = diphoton_vertex_vtx5mmdist[i] = diphoton_vertex_vtx1cmdist[i] = 0;
      diphoton_vertex_x[i] = diphoton_vertex_y[i] = diphoton_vertex_z[i] = -1.;
      diphoton_nearestvtxdist[i] = 999.;

      diphoton_genpt1[i] = diphoton_genpt2[i] = -1.;
      diphoton_geneta1[i] = diphoton_geneta2[i] = -1.;
      diphoton_genphi1[i] = diphoton_genphi2[i] = -1.;
      diphoton_genenergy1[i] = diphoton_genenergy2[i] = -1.;

      diphoton_supercluster_x1[i] = diphoton_supercluster_y1[i] = diphoton_supercluster_z1[i] = -999.;
      diphoton_supercluster_x2[i] = diphoton_supercluster_y2[i] = diphoton_supercluster_z2[i] = -999.;
    }
    diphoton_genvertex_x = diphoton_genvertex_y = diphoton_genvertex_z = -999.;
    diphoton_genvertex_smeared_x = diphoton_genvertex_smeared_y = diphoton_genvertex_smeared_z = -999.;

    num_electron = 0;
    for ( unsigned int i=0; i<MAX_ELECTRON; i++ ) {
      electron_pt[i] = electron_eta[i] = electron_phi[i] = electron_energy[i] = -1.;
      electron_vtx_x[i] = electron_vtx_y[i] = electron_vtx_z[i] = -1.;
    }

    num_muon = 0;
    for ( unsigned int i=0; i<MAX_MUON; i++ ) {
      muon_pt[i] = muon_eta[i] = muon_phi[i] = muon_energy[i] = -1.;
      muon_vtx_x[i] = muon_vtx_y[i] = muon_vtx_z[i] = -1.;
    }

    num_jet = 0;
    for ( unsigned int i=0; i<MAX_JET; i++ ) {
      jet_pt[i] = jet_eta[i] = jet_phi[i] = jet_energy[i] = jet_mass[i] = -1.;
      jet_vtx_x[i] = jet_vtx_y[i] = jet_vtx_z[i] = -1.;
    }

    num_gen_photon = 0;
    for ( unsigned int i=0; i<MAX_GEN_PHOTON; i++ ) {
      gen_photon_pt[i] = gen_photon_eta[i] = gen_photon_phi[i] = gen_photon_energy[i] = -1.;
      gen_photon_vertex_x[i] = gen_photon_vertex_y[i] = gen_photon_vertex_z[i] = -999.;
    }

    num_gen_part = 0;
    for ( unsigned int i=0; i<MAX_GEN_PART; i++ ) {
      gen_part_pdgid[i] = 0;
      gen_part_status[i] = -999;
      gen_part_pt[i] = gen_part_eta[i] = gen_part_phi[i] = gen_part_energy[i] = -1.;
      gen_part_vertex_x[i] = gen_part_vertex_y[i] = gen_part_vertex_z[i] = -999.;
    }

    met = met_phi = met_significance = -1.;

    num_vertex = 0;
    for ( unsigned int i=0; i<MAX_VERTEX; i++ ) {
      vertex_x[i] = vertex_y[i] = vertex_z[i] = -999.;
      //vertex_tracks[i] = vertex_tracksWght0p75[i] = vertex_tracksWght0p9[i] = vertex_tracksWght0p95[i] = 0;
    }
    bs_x0 = bs_y0 = bs_z0 = bs_sigma_z = bs_dxdz = bs_beam_width_x = bs_beam_width_y = -999.;

    pileup_weight = 1.;
  }
  // --- tree components ---

  unsigned int bunch_crossing, fill_number, run_id, lumisection;
  unsigned long long event_number;

  int num_hlt;
  int hlt_accept[MAX_HLT], hlt_prescale[MAX_HLT];

  unsigned int num_fwd_track;
  float fwd_track_x[MAX_PROTON_TRK], fwd_track_y[MAX_PROTON_TRK];
  float fwd_track_chi2[MAX_PROTON_TRK], fwd_track_normchi2[MAX_PROTON_TRK];
  unsigned int fwd_track_arm[MAX_PROTON_TRK], fwd_track_station[MAX_PROTON_TRK], fwd_track_pot[MAX_PROTON_TRK];
  unsigned int num_proton_track;
  float proton_track_vx[MAX_PROTON_TRK], proton_track_vy[MAX_PROTON_TRK];
  float proton_track_xi[MAX_PROTON_TRK];
  unsigned int proton_track_arm[MAX_PROTON_TRK], proton_track_pot[MAX_PROTON_TRK], proton_track_method[MAX_PROTON_TRK];
  float proton_track_chi2[MAX_PROTON_TRK], proton_track_normchi2[MAX_PROTON_TRK];

  unsigned int num_electron;
  float electron_pt[MAX_ELECTRON], electron_eta[MAX_ELECTRON], electron_phi[MAX_ELECTRON], electron_energy[MAX_ELECTRON];
  float electron_vtx_x[MAX_ELECTRON], electron_vtx_y[MAX_ELECTRON], electron_vtx_z[MAX_ELECTRON];

  unsigned int num_muon;
  float muon_pt[MAX_MUON], muon_eta[MAX_MUON], muon_phi[MAX_MUON], muon_energy[MAX_MUON];
  float muon_vtx_x[MAX_MUON], muon_vtx_y[MAX_MUON], muon_vtx_z[MAX_MUON];

  unsigned int num_jet;
  float jet_pt[MAX_JET], jet_eta[MAX_JET], jet_phi[MAX_JET], jet_energy[MAX_JET], jet_mass[MAX_JET];
  float jet_vtx_x[MAX_JET], jet_vtx_y[MAX_JET], jet_vtx_z[MAX_JET];

  unsigned int num_gen_photon;
  float gen_photon_pt[MAX_GEN_PHOTON], gen_photon_eta[MAX_GEN_PHOTON], gen_photon_phi[MAX_GEN_PHOTON], gen_photon_energy[MAX_GEN_PHOTON];
  float gen_photon_vertex_x[MAX_GEN_PHOTON], gen_photon_vertex_y[MAX_GEN_PHOTON], gen_photon_vertex_z[MAX_GEN_PHOTON];

  unsigned int num_gen_part;
  int gen_part_pdgid[MAX_GEN_PART], gen_part_status[MAX_GEN_PART];
  float gen_part_pt[MAX_GEN_PART], gen_part_eta[MAX_GEN_PART], gen_part_phi[MAX_GEN_PART], gen_part_energy[MAX_GEN_PART];
  float gen_part_vertex_x[MAX_GEN_PART], gen_part_vertex_y[MAX_GEN_PART], gen_part_vertex_z[MAX_GEN_PART];

  unsigned int num_diphoton;
  float diphoton_pt1[MAX_DIPHOTON], diphoton_pt2[MAX_DIPHOTON];
  float diphoton_eta1[MAX_DIPHOTON], diphoton_eta2[MAX_DIPHOTON];

  float diphoton_phi1[MAX_DIPHOTON], diphoton_phi2[MAX_DIPHOTON];
  float diphoton_energy1[MAX_DIPHOTON], diphoton_energy2[MAX_DIPHOTON];
  float diphoton_r91[MAX_DIPHOTON], diphoton_r92[MAX_DIPHOTON];
  unsigned int diphoton_ele_veto1[MAX_DIPHOTON], diphoton_ele_veto2[MAX_DIPHOTON];
  float diphoton_id1[MAX_DIPHOTON], diphoton_id2[MAX_DIPHOTON];
  float diphoton_sigeove1[MAX_DIPHOTON], diphoton_sigeove2[MAX_DIPHOTON];
  float diphoton_mass[MAX_DIPHOTON], diphoton_rapidity[MAX_DIPHOTON];
  float diphoton_pt[MAX_DIPHOTON], diphoton_dphi[MAX_DIPHOTON];

  float diphoton_supercluster_x1[MAX_DIPHOTON], diphoton_supercluster_y1[MAX_DIPHOTON], diphoton_supercluster_z1[MAX_DIPHOTON];
  float diphoton_supercluster_x2[MAX_DIPHOTON], diphoton_supercluster_y2[MAX_DIPHOTON], diphoton_supercluster_z2[MAX_DIPHOTON];

  float diphoton_genpt1[MAX_DIPHOTON], diphoton_genpt2[MAX_DIPHOTON];
  float diphoton_geneta1[MAX_DIPHOTON], diphoton_geneta2[MAX_DIPHOTON];
  float diphoton_genphi1[MAX_DIPHOTON], diphoton_genphi2[MAX_DIPHOTON];
  float diphoton_genenergy1[MAX_DIPHOTON], diphoton_genenergy2[MAX_DIPHOTON];

  int diphoton_vertex_id[MAX_DIPHOTON];
  unsigned int diphoton_vertexTracks[MAX_DIPHOTON];
  unsigned int diphoton_vertex_vtx1mmdist[MAX_DIPHOTON], diphoton_vertex_vtx2mmdist[MAX_DIPHOTON], diphoton_vertex_vtx5mmdist[MAX_DIPHOTON], diphoton_vertex_vtx1cmdist[MAX_DIPHOTON];
  float diphoton_vertex_x[MAX_DIPHOTON], diphoton_vertex_y[MAX_DIPHOTON], diphoton_vertex_z[MAX_DIPHOTON];
  float diphoton_nearestvtxdist[MAX_DIPHOTON];

  float diphoton_genvertex_x, diphoton_genvertex_y, diphoton_genvertex_z;
  float diphoton_genvertex_smeared_x, diphoton_genvertex_smeared_y, diphoton_genvertex_smeared_z;

  float met, met_phi, met_significance;

  unsigned int num_vertex;
  float vertex_x[MAX_VERTEX], vertex_y[MAX_VERTEX], vertex_z[MAX_VERTEX];
  int jet_dipho_match[MAX_VERTEX];
  //unsigned int vertex_Tracks[MAX_VERTEX], vertex_TracksWght0p75[MAX_VERTEX], vertex_TracksWght0p9[MAX_VERTEX], vertex_TracksWght0p95[MAX_VERTEX];
  float bs_x0, bs_y0, bs_z0, bs_sigma_z, bs_dxdz, bs_beam_width_x, bs_beam_width_y;

  float pileup_weight;
};

#endif

