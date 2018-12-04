#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TStyle.h"
#include "TSystem.h"

#include "Canvas.h"
#include "Plotter.h"
#include "EventsSelector.h"

#include "pot_alignment.h"
#include "xi_reconstruction.h"

#include "diproton_candidate.h"

#include <fstream>
#include <iostream>

#define out_path "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp"
//#define default_ntp_file "Samples/output_Run2016B_mgg-ov-500.root"
//#define default_ntp_file "Samples/output_Run2016BC_looseCuts.root"
//#define default_ntp_file "Samples/output_Run2016BCG_looseCuts_prelim.root"
//#define default_ntp_file "Samples/output_Run2016BC_looseCuts_7mar.root"
//#define default_ntp_file "Samples/output_Run2016BC_looseCuts_9mar_xifix.root"
//#define default_ntp_file "test.root"
//#define default_ntp_file "Samples/output_Run2016BCG_looseCuts_10mar_xifix.root"
//#define default_ntp_file "Samples/output_Run2016BCG_looseCuts_9may.root"
#define default_ntp_file "Samples/output_Run2016BCG_looseCuts_28jun.root"

float photon_rel_energy_scale( const float& pt, const float& eta, const float& r9 );
bool is_matched( int n_sigma, float xi_rp, float xi_cs, float err_xi_rp, float err_xi_cs );

enum Sector { sector45 = 0, sector56 = 1 };
enum Pot { nearPot = 2, farPot = 3 };

void tree_reader( TString file=default_ntp_file )
{
  gSystem->Load( "libEventFilterUtilities.so" );
  EventsSelector ev_selector( "/afs/cern.ch/user/l/lforthom/public/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_PPSruns_preTS2.txt" );

  xi_reco::load_file( "TreeProducer/data/optics_jun22.root" );
  pot_align::load_file( "TreeProducer/data/alignment_collection_v2.out" );

  TFile f( file );
  if ( !f.IsOpen() ) return;

  const float lumi_b = 5.055026851041,
              lumi_c = 1.470474265456,
              lumi_g = 7.487318307770;

  //const float lumi = ( file.Contains( "2016BCG" ) ) ? lumi_b+lumi_c+lumi_g : lumi_b+lumi_c;
  const float lumi = 9.412003739742; // pre-TS2 runs with horizontal RPs inserted
  const float rel_err_xigg = 0.028;

  //----- cuts -----

  const float max_dist_vtx_inclobjects = 2.0, // in cm
              min_pt_photon = 75.,
              min_r9_photon = 0.94,
              //min_r9_photon = 0.85, //FIXME
              max_eta_photon = 2.5,
              min_etaveto_photon = 1.4442,
              max_etaveto_photon = 1.566,
              min_mass_diphoton = 350.;
  const float max_pt_diphoton = 50.,
              max_acopl = 0.005;
  const float oth_obj_minpt = 150.;
  const float n_sigma_1d = 2.0, // xi matching sigma number
              n_sigma_2d = 2.0;

  /*210-N (beam 1): 0.0334862
    210-N (beam 2): 0.0496449
    210-F (beam 1): 0.023913
    210-F (beam 2): 0.0368263*/
  const float lim_45n = 0.033;
  const float lim_45f = 0.024;
  const float lim_56n = 0.050;
  const float lim_56f = 0.037;

  TH1D* h_pt_incl = new TH1D( "pt_incl", "", 100, 0., 2000. );

  //----- tree loading -----

  TTree* tr = dynamic_cast<TTree*>( f.Get( "ntp" ) );
  // general quantities
  unsigned int run_id, lumisection, fill_number;
  unsigned long long event_number;
  tr->SetBranchAddress( "run_id", &run_id );
  tr->SetBranchAddress( "lumisection", &lumisection );
  tr->SetBranchAddress( "event_number", &event_number );
  tr->SetBranchAddress( "fill_number", &fill_number );
  // diphoton quantities
  unsigned int num_diphoton;
  const unsigned short max_diph = 10;
  float diphoton_pt1[max_diph], diphoton_eta1[max_diph], diphoton_phi1[max_diph], diphoton_r9_1[max_diph],
        diphoton_pt2[max_diph], diphoton_eta2[max_diph], diphoton_phi2[max_diph], diphoton_r9_2[max_diph];
  float diphoton_pt[max_diph], diphoton_mass[max_diph], diphoton_rapidity[max_diph], diphoton_dphi[max_diph];
  float diphoton_vertex_x[max_diph], diphoton_vertex_y[max_diph], diphoton_vertex_z[max_diph];
  float diphoton_sigeove1[max_diph], diphoton_sigeove2[max_diph];
  tr->SetBranchAddress( "num_diphoton", &num_diphoton );
  tr->SetBranchAddress( "diphoton_pt1", diphoton_pt1 );
  tr->SetBranchAddress( "diphoton_eta1", diphoton_eta1 );
  tr->SetBranchAddress( "diphoton_phi1", diphoton_phi1 );
  tr->SetBranchAddress( "diphoton_r91", diphoton_r9_1 );
  tr->SetBranchAddress( "diphoton_pt2", diphoton_pt2 );
  tr->SetBranchAddress( "diphoton_eta2", diphoton_eta2 );
  tr->SetBranchAddress( "diphoton_phi2", diphoton_phi2 );
  tr->SetBranchAddress( "diphoton_r92", diphoton_r9_2 );
  tr->SetBranchAddress( "diphoton_pt", diphoton_pt );
  tr->SetBranchAddress( "diphoton_dphi", diphoton_dphi );
  tr->SetBranchAddress( "diphoton_mass", diphoton_mass );
  tr->SetBranchAddress( "diphoton_rapidity", diphoton_rapidity );
  tr->SetBranchAddress( "diphoton_sigeove1", diphoton_sigeove1 );
  tr->SetBranchAddress( "diphoton_sigeove2", diphoton_sigeove2 );
  unsigned int diphoton_vertex_vtx1mmdist[5], diphoton_vertex_vtx2mmdist[5], diphoton_vertex_vtx5mmdist[5], diphoton_vertex_vtx1cmdist[5];
  float diphoton_vertex_nearestvtxdist[5];
  unsigned int diphoton_vertex_tracks[5];
  tr->SetBranchAddress( "diphoton_vertex_x", diphoton_vertex_x );
  tr->SetBranchAddress( "diphoton_vertex_y", diphoton_vertex_y );
  tr->SetBranchAddress( "diphoton_vertex_z", diphoton_vertex_z );
  tr->SetBranchAddress( "diphoton_vertex_vtx1mmdist", diphoton_vertex_vtx1mmdist );
  tr->SetBranchAddress( "diphoton_vertex_vtx2mmdist", diphoton_vertex_vtx2mmdist );
  tr->SetBranchAddress( "diphoton_vertex_vtx5mmdist", diphoton_vertex_vtx5mmdist );
  tr->SetBranchAddress( "diphoton_vertex_vtx1cmdist", diphoton_vertex_vtx1cmdist );
  tr->SetBranchAddress( "diphoton_vertex_nearestvtxdist", diphoton_vertex_nearestvtxdist );
  tr->SetBranchAddress( "diphoton_vertex_tracks", diphoton_vertex_tracks );
  // proton quantities
  unsigned int num_proton;
  const unsigned short max_pr = 20;
  unsigned int proton_side[max_pr], proton_pot[max_pr];
  float proton_x[max_pr], proton_y[max_pr];
  tr->SetBranchAddress( "num_proton_track", &num_proton );
  tr->SetBranchAddress( "proton_track_side", proton_side );
  tr->SetBranchAddress( "proton_track_pot", proton_pot );
  tr->SetBranchAddress( "proton_track_x", proton_x );
  tr->SetBranchAddress( "proton_track_y", proton_y );
  // vertex quantities
  unsigned int num_vertex;
  const unsigned short max_vtx = 80;
  float vertex_x[max_vtx], vertex_y[max_vtx], vertex_z[max_vtx];
  tr->SetBranchAddress( "num_vertex", &num_vertex );
  tr->SetBranchAddress( "vertex_x", vertex_x );
  tr->SetBranchAddress( "vertex_y", vertex_y );
  tr->SetBranchAddress( "vertex_z", vertex_z );

  // other quantities
  float met, met_phi;
  tr->SetBranchAddress( "met", &met );
  tr->SetBranchAddress( "met_phi", &met_phi );

  //                JW
  // Electron Quantities
  unsigned int num_electron;
  const unsigned short max_ele = 20;
  float electron_pt[max_ele], electron_eta[max_ele], electron_phi[max_ele], electron_energy[max_ele];
  float electron_vtx_x[max_ele], electron_vtx_y[max_ele], electron_vtx_z[max_ele];
  tr->SetBranchAddress( "num_electron", &num_electron );
  tr->SetBranchAddress( "electron_pt", electron_pt );
  tr->SetBranchAddress( "electron_eta", electron_eta );
  tr->SetBranchAddress( "electron_phi", electron_phi );
  tr->SetBranchAddress( "electron_energy", electron_energy );
  tr->SetBranchAddress( "electron_vtx_x", electron_vtx_x );
  tr->SetBranchAddress( "electron_vtx_y", electron_vtx_y );
  tr->SetBranchAddress( "electron_vtx_z", electron_vtx_z );
  // Muon Quantities
  unsigned int num_muon;
  const unsigned short max_mu = 20;
  float muon_pt[max_mu], muon_eta[max_mu], muon_phi[max_mu], muon_energy[max_mu];
  float muon_vtx_x[max_mu], muon_vtx_y[max_mu], muon_vtx_z[max_mu];
  tr->SetBranchAddress( "num_muon",&num_muon );
  tr->SetBranchAddress( "muon_pt", muon_pt );
  tr->SetBranchAddress( "muon_eta",muon_eta );
  tr->SetBranchAddress( "muon_phi",muon_phi );
  tr->SetBranchAddress( "muon_energy", muon_energy );
  tr->SetBranchAddress( "muon_vtx_x", muon_vtx_x );
  tr->SetBranchAddress( "muon_vtx_y", muon_vtx_y );
  tr->SetBranchAddress( "muon_vtx_z", muon_vtx_z );
  //Jet Quantities
  unsigned int num_jet;
  const unsigned short max_jet = 50;
  float jet_pt[max_jet], jet_eta[max_jet], jet_phi[max_jet], jet_energy[max_jet];
  float jet_vtx_x[max_jet], jet_vtx_y[max_jet], jet_vtx_z[max_jet];
  tr->SetBranchAddress( "num_jet", &num_jet );
  tr->SetBranchAddress( "jet_pt",  jet_pt );
  tr->SetBranchAddress( "jet_eta", jet_eta );
  tr->SetBranchAddress( "jet_phi", jet_phi );
  tr->SetBranchAddress( "jet_energy", jet_energy );
  tr->SetBranchAddress( "jet_vtx_x", jet_vtx_x );
  tr->SetBranchAddress( "jet_vtx_y", jet_vtx_y );
  tr->SetBranchAddress( "jet_vtx_z", jet_vtx_z );
  //

  TH1D* h_2dmatch_mpp_over_mgg = new TH1D( "mpp_over_mgg", "m_{pp}^{missing} / m_{#gamma#gamma} for double-tag events@@Events@@?.2f", 30, -2., 4. ),
       *h_2dmatch_ypp_minus_ygg = new TH1D( "ypp_minus_ygg", "y_{pp}^{missing} - y_{#gamma#gamma} for double-tag events@@Events@@?.2f", 50, -2.5, 2.5 ),
       *h_2dmatch_mpair = new TH1D( "mpair_2dmatch", "Diphoton mass (2D matching)@@Events@@GeV", 17, 300., 2000. ),
       *h_2dmatch_ptpair = new TH1D( "ptpair_2dmatch", "Diphoton p_{T} (2D matching)@@Events@@GeV", 10, 0., 200. ),
       *h_2dmatch_r9single_leadpho = new TH1D( "r9single_2dmatch_leadpho", "Single photon r_{9}@@Events@@?.2f", 25, 0.75, 1. ),
       *h_2dmatch_r9single_subleadpho = (TH1D*)h_2dmatch_r9single_leadpho->Clone( "r9single_2dmatch_subleadpho" );
  TH1D* h_met = new TH1D( "met", "Missing E_{T}@@Events@@GeV?.0f", 50, 0., 200. ),
       *h_met_1tag = (TH1D*)h_met->Clone( "met_1tag" ),
       *h_met_2tag = (TH1D*)h_met->Clone( "met_2tag" );
  // diphoton only
  TH1D* h_num_diphoton = new TH1D( "num_diphoton", "Diphoton multiplicity@@Event", 5, 0., 5. ),
       *h_num_diphoton_1tag = (TH1D*)h_num_diphoton->Clone( "num_diphoton_1tag" ),
       *h_num_diphoton_2tag = (TH1D*)h_num_diphoton->Clone( "num_diphoton_2tag" );
  TH1D* h_diphoton_pt = new TH1D( "diphoton_pt", "Diphoton p_{T}@@Events@@GeV?.0f", 40, 0., 400. ),
       *h_diphoton_pt_1tag = (TH1D*)h_diphoton_pt->Clone( "diphoton_pt_1tag" ),
       *h_diphoton_pt_2tag = (TH1D*)h_diphoton_pt->Clone( "diphoton_pt_2tag" );
  TH1D* h_diphoton_pt_zoom = new TH1D( "diphoton_pt_zoom", "Diphoton p_{T}@@Events@@GeV?.0f", 10, 0., 10. ),
       *h_diphoton_pt_zoom_1tag = (TH1D*)h_diphoton_pt_zoom->Clone( "diphoton_pt_zoom_1tag" ),
       *h_diphoton_pt_zoom_2tag = (TH1D*)h_diphoton_pt_zoom->Clone( "diphoton_pt_zoom_2tag" );
  TH1D* h_diphoton_leadpt = new TH1D( "leadphoton_pt", "Leading photon p_{T}@@Events@@GeV?.0f", ( 750.-min_pt_photon )/20, min_pt_photon, 750. ),
       *h_diphoton_leadpt_1tag = (TH1D*)h_diphoton_leadpt->Clone( "leadphoton_pt_1tag" ),
       *h_diphoton_leadpt_2tag = (TH1D*)h_diphoton_leadpt->Clone( "leadphoton_pt_2tag" );
  TH1D* h_diphoton_subleadpt = new TH1D( "subleadphoton_pt", "Subleading photon p_{T}@@Events@@GeV?.0f", ( 750.-min_pt_photon )/20, min_pt_photon, 750. ),
       *h_diphoton_subleadpt_1tag = (TH1D*)h_diphoton_subleadpt->Clone( "subleadphoton_pt_1tag" ),
       *h_diphoton_subleadpt_2tag = (TH1D*)h_diphoton_subleadpt->Clone( "subleadphoton_pt_2tag" );
  TH1D* h_diphoton_leadeta = new TH1D( "leadphoton_eta", "Leading photon #eta@@Events@@?.3f", 40, -2.5, 2.5 ),
       *h_diphoton_leadeta_1tag = (TH1D*)h_diphoton_leadeta->Clone( "leadphoton_eta_1tag" ),
       *h_diphoton_leadeta_2tag = (TH1D*)h_diphoton_leadeta->Clone( "leadphoton_eta_2tag" );
  TH1D* h_diphoton_subleadeta = new TH1D( "subleadphoton_eta", "Subleading photon #eta@@Events@@?.3f", 40, -2.5, 2.5 ),
       *h_diphoton_subleadeta_1tag = (TH1D*)h_diphoton_subleadeta->Clone( "subleadphoton_eta_1tag" ),
       *h_diphoton_subleadeta_2tag = (TH1D*)h_diphoton_subleadeta->Clone( "subleadphoton_eta_2tag" );
  TH1D* h_diphoton_dphi = new TH1D( "diphoton_dphi", "Diphoton 1-|#Delta#phi/#pi|@@Events@@?.2f", 50, 0., 1. ),
       *h_diphoton_dphi_1tag = (TH1D*)h_diphoton_dphi->Clone( "diphoton_dphi_1tag" ),
       *h_diphoton_dphi_2tag = (TH1D*)h_diphoton_dphi->Clone( "diphoton_dphi_2tag" );
  TH1D* h_diphoton_dphi_zoom = new TH1D( "diphoton_dphi_zoom", "Diphoton 1-|#Delta#phi/#pi|@@Events@@?.2f", 20, 0., 0.1 ),
       *h_diphoton_dphi_zoom_1tag = (TH1D*)h_diphoton_dphi_zoom->Clone( "diphoton_dphi_zoom_1tag" ),
       *h_diphoton_dphi_zoom_2tag = (TH1D*)h_diphoton_dphi_zoom->Clone( "diphoton_dphi_zoom_2tag" );
  TH1D* h_diphoton_mass = new TH1D( "diphoton_mass", "Diphoton mass@@Events@@GeV?.0f", ( int )( ( 2000.-min_mass_diphoton )/25 ), min_mass_diphoton, 2000. ),
       *h_diphoton_mass_1tag = (TH1D*)h_diphoton_mass->Clone( "diphoton_mass_1tag" ),
       *h_diphoton_mass_2tag = (TH1D*)h_diphoton_mass->Clone( "diphoton_mass_2tag" );
  TH1D* h_diphoton_mass_withmet = new TH1D( "diphotonmet_mass", "Diphoton + missing E_{T} mass@@Events@@GeV?.0f", ( int )( ( 2500.-min_mass_diphoton )/40 ), min_mass_diphoton, 2500. ),
       *h_diphoton_mass_withmet_1tag = (TH1D*)h_diphoton_mass_withmet->Clone( "diphotonmet_mass_1tag" ),
       *h_diphoton_mass_withmet_2tag = (TH1D*)h_diphoton_mass_withmet->Clone( "diphotonmet_mass_2tag" );
  TH1D* h_diphoton_mass_incl = new TH1D( "diphotonincl_mass", "Diphoton + other objects mass@@Events@@GeV?.0f", ( int )( ( 3500.-min_mass_diphoton )/50 ), min_mass_diphoton, 3500. ),
       *h_diphoton_mass_incl_1tag = (TH1D*)h_diphoton_mass_incl->Clone( "diphotonincl_mass_1tag" ),
       *h_diphoton_mass_incl_2tag = (TH1D*)h_diphoton_mass_incl->Clone( "diphotonincl_mass_2tag" );
  TH1D* h_diphoton_rap = new TH1D( "diphoton_rap", "Diphoton rapidity@@Events@@?.1f", 40, -2.5, 2.5 ),
       *h_diphoton_rap_1tag = (TH1D*)h_diphoton_rap->Clone( "diphoton_rap_1tag" ),
       *h_diphoton_rap_2tag = (TH1D*)h_diphoton_rap->Clone( "diphoton_rap_2tag" );
  TH1D* h_diphoton_ntrk = new TH1D( "diphoton_ntrk", "Number of tracks on diphoton vertex@@Events@@?.1f", 20, 0., 20. ),
       *h_diphoton_ntrk_1tag = (TH1D*)h_diphoton_ntrk->Clone( "diphoton_ntrk_1tag" ),
       *h_diphoton_ntrk_2tag = (TH1D*)h_diphoton_ntrk->Clone( "diphoton_ntrk_2tag" );
  TH1D* h_diphoton_vtxz = new TH1D( "diphoton_vtxz", "Diphoton longitudinal vertex position@@Events@@cm?.1f", 40, -10., 10. ),
       *h_diphoton_vtxz_1tag = (TH1D*)h_diphoton_vtxz->Clone( "diphoton_vtxz_1tag" ),
       *h_diphoton_vtxz_2tag = (TH1D*)h_diphoton_vtxz->Clone( "diphoton_vtxz_2tag" );
  // vertexing study
  TH1D* h_num_vtx = new TH1D( "num_vtx", "Number of primary vertices in event@@Events", 40, 0., 40. ),
       *h_num_vtx_1tag = (TH1D*)h_num_vtx->Clone( "num_vtx_1tag" ),
       *h_num_vtx_2tag = (TH1D*)h_num_vtx->Clone( "num_vtx_2tag" );
  TH1D* h_vtxz = new TH1D( "vtxz", "Longitudinal primary vertices position@@Events@@cm", 30, -15., 15. ),
       *h_vtxz_1tag = (TH1D*)h_vtxz->Clone( "vtxz_1tag" ),
       *h_vtxz_2tag = (TH1D*)h_vtxz->Clone( "vtxz_2tag" );
  TH1D* h_num_vtx_1mm = new TH1D( "num_vtx_1mm", "Number of primary vertices near diphoton vertex@@Events", 12, 0., 12. ),
       *h_num_vtx_2mm = (TH1D*)h_num_vtx_1mm->Clone( "num_vtx_2mm" ),
       *h_num_vtx_5mm = (TH1D*)h_num_vtx_1mm->Clone( "num_vtx_5mm" ),
       *h_num_vtx_1cm = (TH1D*)h_num_vtx_1mm->Clone( "num_vtx_1cm" ),
       *h_num_vtx_1mm_excl = (TH1D*)h_num_vtx_1mm->Clone( "num_vtx_excl_1mm" ),
       *h_num_vtx_2mm_excl = (TH1D*)h_num_vtx_1mm->Clone( "num_vtx_excl_2mm" ),
       *h_num_vtx_5mm_excl = (TH1D*)h_num_vtx_1mm->Clone( "num_vtx_excl_5mm" ),
       *h_num_vtx_1cm_excl = (TH1D*)h_num_vtx_1mm->Clone( "num_vtx_excl_1cm" );
  TH1D* h_diphoton_closestvtx = new TH1D( "diphoton_closestvtx", "Distance diphoton/nearest vertex@@Events@@cm?.1f", 25, 0., 5.0 ),
       *h_diphoton_closestvtx_1tag = (TH1D*)h_diphoton_closestvtx->Clone( "diphoton_closestvtx_1tag" ),
       *h_diphoton_closestvtx_2tag = (TH1D*)h_diphoton_closestvtx->Clone( "diphoton_closestvtx_2tag" );
  TH1D* h_electron_pt_1mm = new TH1D( "electron_pt_1mm", "Close electron E_{T}@@Events@@GeV", 25, 0., 500. ),
       *h_electron_pt_2mm = (TH1D*)h_electron_pt_1mm->Clone( "electron_pt_2mm" ),
       *h_electron_pt_5mm = (TH1D*)h_electron_pt_1mm->Clone( "electron_pt_5mm" ),
       *h_electron_pt_1cm = (TH1D*)h_electron_pt_1mm->Clone( "electron_pt_1cm" );
  TH1D* h_muon_pt_1mm = new TH1D( "muon_pt_1mm", "Close muon p_{T}@@Events@@GeV", 25, 0., 250. ),
       *h_muon_pt_2mm = (TH1D*)h_muon_pt_1mm->Clone( "muon_pt_2mm" ),
       *h_muon_pt_5mm = (TH1D*)h_muon_pt_1mm->Clone( "muon_pt_5mm" ),
       *h_muon_pt_1cm = (TH1D*)h_muon_pt_1mm->Clone( "muon_pt_1cm" );
  TH1D* h_jet_pt_1mm = new TH1D( "jet_pt_1mm", "Close jet p_{T}@@Events@@GeV", 50, 0., 500. ),
       *h_jet_pt_2mm = (TH1D*)h_jet_pt_1mm->Clone( "jet_pt_2mm" ),
       *h_jet_pt_5mm = (TH1D*)h_jet_pt_1mm->Clone( "jet_pt_5mm" ),
       *h_jet_pt_1cm = (TH1D*)h_jet_pt_1mm->Clone( "jet_pt_1cm" );
  //  JW
  TH1D* h_electron_pt = new TH1D("electron_pt", "Close electrons' #Sigma(E_{T})@@Events@@GeV", 50, 0., 250.),
       *h_electron_pt_1tag = (TH1D*)h_electron_pt->Clone( "electron_pt_1tag" ),
       *h_electron_pt_2tag = (TH1D*)h_electron_pt->Clone( "electron_pt_2tag" );
  TH1D* h_muon_pt = new TH1D("muon_pt", "Close muons' #Sigma(p_{T})@@Events@@GeV", 30, 0., 150.),
       *h_muon_pt_1tag = (TH1D*)h_muon_pt->Clone( "muon_pt_1tag" ),
       *h_muon_pt_2tag = (TH1D*)h_muon_pt->Clone( "muon_pt_2tag" );
  TH1D* h_jet_pt = new TH1D("jet_pt", "Close jets' #Sigma(p_{T})@@Events@@GeV", 30, 0., 150.),
       *h_jet_pt_1tag = (TH1D*)h_jet_pt->Clone( "jet_pt_1tag" ),
       *h_jet_pt_2tag = (TH1D*)h_jet_pt->Clone( "jet_pt_2tag" );
  //
  TH2D* h_met_vs_pt = new TH2D( "met_vs_pt", "Missing E_{T} (GeV)@@Diphoton p_{T} (GeV)", 40, 0., 400., 40, 0., 400. ),
       *h_met_vs_pt_2tag = (TH2D*)h_met_vs_pt->Clone( "met_vs_pt_2tag" ),
       *h_metx_vs_mety = new TH2D( "metx_vs_mety", "#slash{E}_{T,x} (GeV)@@#slash{E}_{T,y} (GeV)", 50, -100., 100., 50, -100., 100. ),
       *h_metx_vs_mety_1tag = (TH2D*)h_metx_vs_mety->Clone( "metx_vs_mety_1tag" ),
       *h_metx_vs_mety_2tag = (TH2D*)h_metx_vs_mety->Clone( "metx_vs_mety_2tag" );
  TH2D* h_mggmet_vs_mgg = new TH2D( "mggmet_vs_mgg", "Diphoton + #slash{E}_{T} mass (GeV)@@Diphoton mass (GeV)", ( int )( ( 2000.-min_mass_diphoton )/50 ), min_mass_diphoton, 2000., ( int )( ( 2000.-min_mass_diphoton )/50 ), min_mass_diphoton, 2000. );
  // balances
  TGraphErrors h_ptgg_vs_mpp, h_ptgg_vs_mgg;
  TGraphErrors h_mgg_vs_mpp, h_mgg_vs_mpp_candm, h_mgg_vs_mpp_candy,
               h_mggmet_vs_mpp, h_mggmet_vs_mpp_candm, h_mggmet_vs_mpp_candy,
               h_mggincl_vs_mpp, h_mggincl_vs_mpp_candm, h_mggincl_vs_mpp_candy;
  TGraphErrors h_ygg_vs_ypp, h_ygg_vs_ypp_candm, h_ygg_vs_ypp_candy,
               h_yggmet_vs_ypp, h_yggmet_vs_ypp_candm, h_yggmet_vs_ypp_candy,
               h_yggincl_vs_ypp, h_yggincl_vs_ypp_candm, h_yggincl_vs_ypp_candy;
  TGraphErrors h_2dmatch_withmet_metvspt;
  // proton reco study
  TH1D* h_num_proton = new TH1D( "num_proton", "Forward track multiplicity@@Events", 6, 0., 6. ),
       *h_num_proton_45 = (TH1D*)h_num_proton->Clone( "num_proton_45" ),
       *h_num_proton_56 = (TH1D*)h_num_proton->Clone( "num_proton_56" );
  TH1D* h_diproton_mass = new TH1D( "diproton_mass", "Diproton mass@@Events@@GeV", 40, 300., 2300. );

  const unsigned short num_bins = 50;
  const float max_bin_dist = 5.;
  float dist_bin[num_bins];
  unsigned short num_ele_bin[num_bins], num_ele_ptcut_bin[num_bins],
                 num_mu_bin[num_bins], num_mu_ptcut_bin[num_bins],
                 num_jet_bin[num_bins], num_jet_ptcut_bin[num_bins];
  for ( unsigned short i=0; i<num_bins; i++ ) {
    dist_bin[i] = max_bin_dist*i/( num_bins+1 );
  }

  TH1D* gr_veto_ele = new TH1D( "gr_veto_ele", "Veto on distance to diphoton vertex (cm)@@#LT# nearby objects#GT / event / diphoton cand.", num_bins-1, dist_bin ),
       *gr_veto_mu = ( TH1D* )gr_veto_ele->Clone( "gr_veto_mu" ),
       *gr_veto_jet = ( TH1D* )gr_veto_ele->Clone( "gr_veto_jet" ),
       *gr_veto_ele_ptcut = ( TH1D* )gr_veto_ele->Clone( "gr_veto_ele_ptcut" ),
       *gr_veto_mu_ptcut = ( TH1D* )gr_veto_ele->Clone( "gr_veto_mu_ptcut" ),
       *gr_veto_jet_ptcut = ( TH1D* )gr_veto_ele->Clone( "gr_veto_jet_ptcut" );

  ofstream events_list( "events_list.txt" );

  unsigned int num_evts_notag = 0, num_evts_with_tag = 0;
  TLorentzVector pho1, pho2, electron, muon, jet;

  unsigned int num_match_diph = 0, num_match_withmet = 0, num_match_incl = 0;

  // tree readout stage
  for ( unsigned int i=0; i<tr->GetEntries(); i++ ) {
    tr->GetEntry( i );

    //if ( run_id>280385 ) { continue; } // skip post-TS2 runs
    if ( !ev_selector.isSelected( run_id, lumisection, event_number ) ) continue;

    auto align = pot_align::get_alignments( fill_number );

    //----- diphotons retrieval part -----

    unsigned short num_diphoton_cand = 0, num_diphoton_cand_1tag = 0, num_diphoton_cand_2tag = 0;
    unsigned short num_proton_45 = 0, num_proton_56 = 0;
    typedef pair<float,float> ximeas_t;
    vector<ximeas_t> fwdtrk_45n, fwdtrk_45f, fwdtrk_56n, fwdtrk_56f;

    // first fill the proton hitmaps outside the diphoton loop
    for ( unsigned int j=0; j<num_proton; j++ ) {
      const unsigned short pot_id = 100*proton_side[j]+proton_pot[j];
      const auto& al = align[pot_id];
      double xi, xi_err;
      xi_reco::reconstruct( proton_x[j]+al.x, proton_side[j], proton_pot[j], xi, xi_err );

      if ( proton_side[j]==0 && proton_pot[j]==2 ) {
        fwdtrk_45n.emplace_back( xi, xi_err );
        num_proton_45++;
      }
      if ( proton_side[j]==0 && proton_pot[j]==3 ) {
        fwdtrk_45f.emplace_back( xi, xi_err );
        num_proton_45++;
      }
      if ( proton_side[j]==1 && proton_pot[j]==2 ) {
        fwdtrk_56n.emplace_back( xi, xi_err );
        num_proton_56++;
      }
      if ( proton_side[j]==1 && proton_pot[j]==3 ) {
        fwdtrk_56f.emplace_back( xi, xi_err );
        num_proton_56++;
      }
    }

    for ( unsigned short i=0; i<num_bins; i++ ) {
      num_ele_bin[i] = num_mu_bin[i] = num_jet_bin[i] = 0;
      num_ele_ptcut_bin[i] = num_mu_ptcut_bin[i] = num_jet_ptcut_bin[i] = 0;
    }

    //----- list all possible 45/56 diproton matches

    vector<diproton_candidate_t> candidates;
    for ( const auto& trk45 : fwdtrk_45n ) {
      for ( const auto& trk56 : fwdtrk_56n ) {
        candidates.emplace_back( trk45.first, trk45.second, trk56.first, trk56.second );
      }
      for ( const auto& trk56 : fwdtrk_56f ) {
        candidates.emplace_back( trk45.first, trk45.second, trk56.first, trk56.second );
      }
    }
    for ( const auto& trk45 : fwdtrk_45f ) {
      for ( const auto& trk56 : fwdtrk_56n ) {
        candidates.emplace_back( trk45.first, trk45.second, trk56.first, trk56.second );
      }
      for ( const auto& trk56 : fwdtrk_56f ) {
        candidates.emplace_back( trk45.first, trk45.second, trk56.first, trk56.second );
      }
    }

    //----- BEFORE THE LOOP ON DIPHOTON -----

    for ( unsigned int j=0; j<num_diphoton; j++ ) {

      const float acopl = 1-fabs( diphoton_dphi[j]/TMath::Pi() );

      /*const float energy_corr_pho1 = photon_rel_energy_scale( diphoton_pt1[j], diphoton_eta1[j], diphoton_r9_1[j] ) * diphoton_pt1[j],
	energy_corr_pho2 = photon_rel_energy_scale( diphoton_pt2[j], diphoton_eta2[j], diphoton_r9_2[j] ) * diphoton_pt2[j];*/
//      const float energy_corr_pho1 = 3.5, energy_corr_pho2 = 3.5; //FIXME
      const float energy_corr_pho1 = diphoton_sigeove1[j] * diphoton_pt1[j],
                  energy_corr_pho2 = diphoton_sigeove2[j] * diphoton_pt2[j];
      //cout << ">>> " << diphoton_sigeove1[j]*diphoton_pt1[j] << " :: " << photon_rel_energy_scale( diphoton_pt1[j], diphoton_eta1[j], diphoton_r9_1[j] )*diphoton_pt1[j] << "\t" << diphoton_sigeove2[j]*diphoton_pt2[j] << " :: " << photon_rel_energy_scale( diphoton_pt2[j], diphoton_eta2[j], diphoton_r9_2[j] )*diphoton_pt2[j] << endl;

      const float xim_reco = ( diphoton_pt1[j] * exp( -diphoton_eta1[j] ) + diphoton_pt2[j] * exp( -diphoton_eta2[j] ) )/sqrt_s,
                  xip_reco = ( diphoton_pt1[j] * exp(  diphoton_eta1[j] ) + diphoton_pt2[j] * exp(  diphoton_eta2[j] ) )/sqrt_s,
      /*err_xim_reco = sqrt( pow( energy_corr_pho1, 2 ) * exp( -2.*diphoton_eta1[j] ) + pow( energy_corr_pho2, 2 ) * exp( -2.*diphoton_eta2[j] ) )/sqrt_s,
      err_xip_reco = sqrt( pow( energy_corr_pho1, 2 ) * exp(  2.*diphoton_eta1[j] ) + pow( energy_corr_pho2, 2 ) * exp(  2.*diphoton_eta2[j] ) )/sqrt_s; //FIXME*/
	/*err_xim_reco = sqrt( pow( energy_corr_pho1, 2 ) * exp( -2.*diphoton_eta1[j] ) + pow( energy_corr_pho2, 2 ) * exp( -2.*diphoton_eta2[j] ) )/sqrt_s,
	  err_xip_reco = sqrt( pow( energy_corr_pho1, 2 ) * exp(  2.*diphoton_eta1[j] ) + pow( energy_corr_pho2, 2 ) * exp(  2.*diphoton_eta2[j] ) )/sqrt_s; //FIXME*/
      err_xim_reco = xim_reco*rel_err_xigg,
      err_xip_reco = xip_reco*rel_err_xigg;

      //cout << xim_reco << " +- " << err_xim_reco << " // " << xip_reco << " +- " << err_xip_reco << endl;
      //cout << ">>>> " << err_xim_reco << "\t" << err_xip_reco << endl;
      //cout << ">>>> " << err_xim_reco/xim_reco << "\t" << err_xip_reco/xip_reco << endl;

      const float xim_reco_withmet = xim_reco + met/sqrt_s, err_xim_reco_withmet = xim_reco_withmet*rel_err_xigg,
                  xip_reco_withmet = xip_reco + met/sqrt_s, err_xip_reco_withmet = xip_reco_withmet*rel_err_xigg;

      //----- PRESELECTION ------

      //----- quality cuts for the photons -----

      if ( diphoton_pt1[j]<min_pt_photon || diphoton_pt2[j]<min_pt_photon ) continue;
      if ( fabs( diphoton_eta1[j] )>max_eta_photon || fabs( diphoton_eta2[j] )>max_eta_photon ) continue;
      if ( ( fabs( diphoton_eta1[j] )>=min_etaveto_photon && fabs( diphoton_eta1[j] )<=max_etaveto_photon ) ||
           ( fabs( diphoton_eta2[j] )>=min_etaveto_photon && fabs( diphoton_eta2[j] )<=max_etaveto_photon ) ) continue;
      if ( diphoton_r9_1[j]<min_r9_photon || diphoton_r9_2[j]<min_r9_photon ) continue;
      if ( diphoton_mass[j]<min_mass_diphoton ) continue;

      //----- build the event's kinematics -----

      pho1.SetPtEtaPhiM( diphoton_pt1[j], diphoton_eta1[j], diphoton_phi1[j], 0. );
      pho2.SetPtEtaPhiM( diphoton_pt2[j], diphoton_eta2[j], diphoton_phi2[j], 0. );

      const TVector3 diph_vtx( diphoton_vertex_x[j], diphoton_vertex_y[j], diphoton_vertex_z[j] );
      float min_dist_vtx = 999.999;
      for ( unsigned int k=0; k<num_vertex; k++ ) {
        const TVector3 oth_vtx( vertex_x[k], vertex_y[k], vertex_z[k] );
        float dist = ( oth_vtx-diph_vtx ).Mag();
        if ( dist>0. && dist<min_dist_vtx ) {
          min_dist_vtx = dist;
        }
      }
      h_diphoton_closestvtx->Fill( min_dist_vtx );

      //----- leptons+jets retrieval part -----

      bool has_ele = false, has_muon = false, has_jet = false;
      TLorentzVector electrons, muons, jets;
      vector< pair<float,TLorentzVector> > electrons_list, muons_list, jets_list;

      unsigned short num_ele_included = 0,
                     num_mu_included = 0,
                     num_jet_included = 0;

      for ( unsigned int k=0; k<num_electron; k++ ) {
        TLorentzVector ele; ele.SetPtEtaPhiE( electron_pt[k], electron_eta[k], electron_phi[k], electron_energy[k] );
        const TVector3 ele_vtx( electron_vtx_x[k], electron_vtx_y[k], electron_vtx_z[k] );
        const float dist = ( ele_vtx-diph_vtx ).Mag();
        if ( dist<0.1 ) { h_electron_pt_1mm->Fill( ele.Pt() ); }
        if ( dist<0.2 ) { h_electron_pt_2mm->Fill( ele.Pt() ); }
        if ( dist<0.5 ) { h_electron_pt_5mm->Fill( ele.Pt() ); }
        if ( dist<1.0 ) { h_electron_pt_1cm->Fill( ele.Pt() ); }
        for ( unsigned short l=0; l<num_bins; l++ ) {
          if ( dist<dist_bin[l] ) num_ele_bin[l]++;
          if ( ele.Pt()>oth_obj_minpt && dist<dist_bin[l] ) num_ele_ptcut_bin[l]++;
        }
        if ( dist<max_dist_vtx_inclobjects ) {
          electrons += ele;
          electrons_list.push_back( make_pair( dist, ele ) );
          has_ele = true;
          num_ele_included++;
        }
      }

      for ( unsigned int k=0; k<num_muon; k++ ) {
        TLorentzVector mu; mu.SetPtEtaPhiE( muon_pt[k], muon_eta[k], muon_phi[k], muon_energy[k] );
        const TVector3 mu_vtx( muon_vtx_x[k], muon_vtx_y[k], muon_vtx_z[k] );
        const float dist = ( mu_vtx-diph_vtx ).Mag();
        if ( dist<0.1 ) { h_muon_pt_1mm->Fill( mu.Pt() ); }
        if ( dist<0.2 ) { h_muon_pt_2mm->Fill( mu.Pt() ); }
        if ( dist<0.5 ) { h_muon_pt_5mm->Fill( mu.Pt() ); }
        if ( dist<1.0 ) { h_muon_pt_1cm->Fill( mu.Pt() ); }
        for ( unsigned short l=0; l<num_bins; l++ ) {
          if ( dist<dist_bin[l] ) num_mu_bin[l]++;
          if ( mu.Pt()>oth_obj_minpt && dist<dist_bin[l] ) num_mu_ptcut_bin[l]++;
        }
        if ( dist<max_dist_vtx_inclobjects ) {
          muons += mu;
          muons_list.push_back( make_pair( dist, mu ) );
          has_muon = true;
          num_mu_included++;
        }
      }

      for ( unsigned int k=0; k<num_jet; k++ ) {
        TLorentzVector jet; jet.SetPtEtaPhiE( jet_pt[k], jet_eta[k], jet_phi[k], jet_energy[k] );
        const TVector3 jet_vtx( jet_vtx_x[k], jet_vtx_y[k], jet_vtx_z[k] );
        const float dist = ( jet_vtx-diph_vtx ).Mag();
        if ( dist<0.1 ) { h_jet_pt_1mm->Fill( jet.Pt() ); }
        if ( dist<0.2 ) { h_jet_pt_2mm->Fill( jet.Pt() ); }
        if ( dist<0.5 ) { h_jet_pt_5mm->Fill( jet.Pt() ); }
        if ( dist<1.0 ) { h_jet_pt_1cm->Fill( jet.Pt() ); }
        for ( unsigned short l=0; l<num_bins; l++ ) {
          if ( dist<dist_bin[l] ) num_jet_bin[l]++;
          if ( jet.Pt()>oth_obj_minpt && dist<dist_bin[l] ) num_jet_ptcut_bin[l]++;
        }
        /*if ( dist<max_dist_vtx_inclobjects ) {
          jets += jet;
          jets_list.push_back( make_pair( dist, jet ) );
          has_jet = true;
          num_jet_included++;
	  }*/
      }

      for ( unsigned short k=0; k<num_bins; k++ ) {
        gr_veto_ele->Fill( dist_bin[k], num_ele_bin[k]/num_diphoton );
        gr_veto_mu->Fill( dist_bin[k], num_mu_bin[k]/num_diphoton );
        gr_veto_jet->Fill( dist_bin[k], num_jet_bin[k]/num_diphoton );
        gr_veto_ele_ptcut->Fill( dist_bin[k], num_ele_ptcut_bin[k]/num_diphoton );
        gr_veto_mu_ptcut->Fill( dist_bin[k], num_mu_ptcut_bin[k]/num_diphoton );
        gr_veto_jet_ptcut->Fill( dist_bin[k], num_jet_ptcut_bin[k]/num_diphoton );
      }

      const float met_x = met*cos( met_phi ),
                  met_y = met*sin( met_phi );

      const TLorentzVector lv_met( met_x, met_y, 0., met ),
                           dipho_met = pho1+pho2+lv_met,
                           //dipho_incl = pho1+pho2+electrons+muons+jets+lv_met;
                           dipho_incl = pho1+pho2+electrons+muons+jets;
      //cout << diphoton_mass[j] << "\t" << dipho_met.M() << "\t" << dipho_incl.M() << endl;
      const float diphoton_plus_met_mass = dipho_met.M(),
                  diphoton_plus_met_rap = dipho_met.Rapidity(),
                  diphoton_incl_mass = dipho_incl.M(),
                  diphoton_incl_rap = dipho_incl.Rapidity();

      //----- N-1 plots -----

      h_diphoton_pt->Fill( diphoton_pt[j] );
      h_diphoton_pt_zoom->Fill( diphoton_pt[j] );
      h_diphoton_mass->Fill( diphoton_mass[j] );
      h_diphoton_mass_withmet->Fill( diphoton_plus_met_mass );
      h_diphoton_mass_incl->Fill( diphoton_incl_mass );
      h_diphoton_rap->Fill( diphoton_rapidity[j] );
      h_diphoton_ntrk->Fill( diphoton_vertex_tracks[j] );
      h_diphoton_dphi->Fill( acopl );
      h_diphoton_dphi_zoom->Fill( acopl );
      h_diphoton_leadpt->Fill( diphoton_pt1[j] );
      h_diphoton_subleadpt->Fill( diphoton_pt2[j] );
      h_diphoton_leadeta->Fill( diphoton_eta1[j] );
      h_diphoton_subleadeta->Fill( diphoton_eta2[j] );
      h_diphoton_vtxz->Fill( diphoton_vertex_z[j] );

      h_electron_pt->Fill( electrons.Pt() );
      h_muon_pt->Fill( muons.Pt() );
      h_jet_pt->Fill( jets.Pt() );

      h_met_vs_pt->Fill( met, diphoton_pt[j] );
      h_metx_vs_mety->Fill( met_x, met_y );
      h_mggmet_vs_mgg->Fill( diphoton_plus_met_mass, diphoton_mass[j] );

      // rough check of the proton tagging information for N-1 plots
      // (to extract the tagged arms only ; not yet at a Roman pot granularity)

      bool has_singletag_45 = false, has_singletag_56 = false;

      for ( unsigned short k=0; k<num_proton; k++ ) {
        if ( proton_side[k]==sector45 ) {
          if ( proton_pot[k]==nearPot ) { has_singletag_45 = true; }
          if ( proton_pot[k]==farPot ) { has_singletag_45 = true; }
        }
        if ( proton_side[k]==sector56 ) {
          if ( proton_pot[k]==nearPot ) { has_singletag_56 = true; }
          if ( proton_pot[k]==farPot ) { has_singletag_56 = true; }
        }
      }

      bool has_singletag = ( has_singletag_45 || has_singletag_56 ),
           has_doubletag = ( has_singletag_45 && has_singletag_56 );

      if ( has_singletag ) {
        h_diphoton_pt_1tag->Fill( diphoton_pt[j] );
        h_diphoton_pt_zoom_1tag->Fill( diphoton_pt[j] );
        h_diphoton_mass_1tag->Fill( diphoton_mass[j] );
        h_diphoton_mass_withmet_1tag->Fill( diphoton_plus_met_mass );
        h_diphoton_mass_incl_1tag->Fill( diphoton_incl_mass );
        h_diphoton_rap_1tag->Fill( diphoton_rapidity[j] );
        h_diphoton_closestvtx_1tag->Fill( min_dist_vtx );
        h_diphoton_ntrk_1tag->Fill( diphoton_vertex_tracks[j] );
        h_diphoton_dphi_1tag->Fill( acopl );
        h_diphoton_dphi_zoom_1tag->Fill( acopl );
        h_diphoton_leadpt_1tag->Fill( diphoton_pt1[j] );
        h_diphoton_subleadpt_1tag->Fill( diphoton_pt2[j] );
        h_diphoton_leadeta_1tag->Fill( diphoton_eta1[j] );
        h_diphoton_subleadeta_1tag->Fill( diphoton_eta2[j] );
        h_diphoton_vtxz_1tag->Fill( diphoton_vertex_z[j] );

        h_electron_pt_1tag->Fill( electrons.Pt() );
        h_muon_pt_1tag->Fill( muons.Pt() );
        h_jet_pt_1tag->Fill( jets.Pt() );

        h_metx_vs_mety_1tag->Fill( met_x, met_y );

        num_diphoton_cand_1tag++;
      }
      if ( has_doubletag ) {
        //cout << "maximal diproton mass: " << max_diproton_mass << " +- " << max_diproton_mass_error << " (rapidity=" << max_diproton_mass_rap << ")" << endl;
        h_diphoton_pt_2tag->Fill( diphoton_pt[j] );
        h_diphoton_pt_zoom_2tag->Fill( diphoton_pt[j] );
        h_diphoton_mass_2tag->Fill( diphoton_mass[j] );
        h_diphoton_mass_withmet_2tag->Fill( diphoton_plus_met_mass );
        h_diphoton_mass_incl_2tag->Fill( diphoton_incl_mass );
        h_diphoton_rap_2tag->Fill( diphoton_rapidity[j] );
        h_diphoton_closestvtx_2tag->Fill( min_dist_vtx );
        h_diphoton_ntrk_2tag->Fill( diphoton_vertex_tracks[j] );
        h_diphoton_dphi_2tag->Fill( acopl );
        h_diphoton_dphi_zoom_2tag->Fill( acopl );
        h_diphoton_leadpt_2tag->Fill( diphoton_pt1[j] );
        h_diphoton_subleadpt_2tag->Fill( diphoton_pt2[j] );
        h_diphoton_leadeta_2tag->Fill( diphoton_eta1[j] );
        h_diphoton_subleadeta_2tag->Fill( diphoton_eta2[j] );
        h_diphoton_vtxz_2tag->Fill( diphoton_vertex_z[j] );

        h_ptgg_vs_mgg.SetPoint( h_ptgg_vs_mgg.GetN(), diphoton_pt[j], diphoton_mass[j] );

        h_electron_pt_2tag->Fill( electrons.Pt() );
        h_muon_pt_2tag->Fill( muons.Pt() );
        h_jet_pt_2tag->Fill( jets.Pt() );

        h_met_vs_pt_2tag->Fill( met, diphoton_pt[j] );
        h_metx_vs_mety_2tag->Fill( met_x, met_y );

        num_evts_with_tag++;
        num_diphoton_cand_2tag++;
      }

      h_num_vtx_1mm->Fill( diphoton_vertex_vtx1mmdist[j] );
      h_num_vtx_2mm->Fill( diphoton_vertex_vtx2mmdist[j] );
      h_num_vtx_5mm->Fill( diphoton_vertex_vtx5mmdist[j] );
      h_num_vtx_1cm->Fill( diphoton_vertex_vtx1cmdist[j] );

      num_evts_notag++;
      num_diphoton_cand++;

      //----- exclusivity cuts -----

      //if ( has_ele || has_muon || has_jet ) continue; //FIXME FIXME FIXME FIXME FIXME
      //if ( diphoton_pt[j]>max_pt_diphoton ) continue;
      if ( acopl>max_acopl ) continue;
      //if ( diphoton_pt1[j]/diphoton_pt2[j]<0.95 ) continue; // matches Fichet, Royon et al
      //if ( has_ele || has_muon || has_jet ) continue; //FIXME FIXME FIXME FIXME FIXME
      //if ( diphoton_pt[j]>max_pt_diphoton ) continue;

      //----- FROM THAT POINT ON, EXCLUSIVE DIPHOTON CANDIDATE -----

      h_num_vtx_1mm_excl->Fill( diphoton_vertex_vtx1mmdist[j] );
      h_num_vtx_2mm_excl->Fill( diphoton_vertex_vtx2mmdist[j] );
      h_num_vtx_5mm_excl->Fill( diphoton_vertex_vtx5mmdist[j] );
      h_num_vtx_1cm_excl->Fill( diphoton_vertex_vtx1cmdist[j] );

      //----- loop over diproton candidates

      unsigned short cand_2tag_id = 0;
      for ( const auto& cand : candidates ) { // double-tag candidates
        h_diproton_mass->Fill( cand.mass() );

        double diphoton_mass_err = diphoton_mass[j]*0.03;
        double diphoton_rapidity_err = diphoton_rapidity[j]*0.0; //FIXME

        bool mass_match = is_matched( n_sigma_2d, diphoton_mass[j], cand.mass(), diphoton_mass_err, cand.mass_error() );
        bool rap_match = is_matched( n_sigma_2d, diphoton_rapidity[j], cand.rapidity(), diphoton_rapidity_err, cand.rapidity_error() );
        bool mass_match_withmet = is_matched( n_sigma_2d, diphoton_plus_met_mass, cand.mass(), 0., cand.mass_error() );
        bool rap_match_withmet = is_matched( n_sigma_2d, diphoton_plus_met_rap, cand.rapidity(), 0., cand.rapidity_error() );
        bool mass_match_incl = is_matched( n_sigma_2d, diphoton_incl_mass, cand.mass(), 0., cand.mass_error() );
        bool rap_match_incl = is_matched( n_sigma_2d, diphoton_incl_rap, cand.rapidity(), 0., cand.rapidity_error() );

        //----- matching (diphoton only) -----
        h_mgg_vs_mpp.SetPoint( cand_2tag_id, diphoton_mass[j], cand.mass() );
        h_mgg_vs_mpp.SetPointError( cand_2tag_id, diphoton_mass_err, cand.mass_error() );
        h_ygg_vs_ypp.SetPoint( cand_2tag_id, diphoton_rapidity[j], cand.rapidity() );
        h_ygg_vs_ypp.SetPointError( cand_2tag_id, diphoton_rapidity_err, cand.rapidity_error() );

        if ( mass_match ) {
          const unsigned short id_match = h_mgg_vs_mpp_candm.GetN();
          h_mgg_vs_mpp_candm.SetPoint( id_match, diphoton_mass[j], cand.mass() );
          h_ygg_vs_ypp_candm.SetPoint( id_match, diphoton_rapidity[j], cand.rapidity() );
        }
        if ( rap_match ) {
          const unsigned short id_match = h_mgg_vs_mpp_candy.GetN();
          h_mgg_vs_mpp_candy.SetPoint( id_match, diphoton_mass[j], cand.mass() );
          h_ygg_vs_ypp_candy.SetPoint( id_match, diphoton_rapidity[j], cand.rapidity() );
        }

        if ( mass_match && rap_match ) {
          num_match_diph++;
          cout << "--------> matching!!!!!" << endl;
          cout << " >>> " << diphoton_mass[j] << "/" << cand.mass() << " && " << diphoton_rapidity[j] << "/" << cand.rapidity() << endl;
          cout << "     " << run_id << ":" << lumisection << ":" << event_number << "\t" << "pt=" << diphoton_pt[j] << endl;
          cout << "  diphoton vertex: " << diphoton_vertex_x[j] << ", " << diphoton_vertex_y[j] << ", " << diphoton_vertex_z[j] << endl;
          h_2dmatch_mpair->Fill( diphoton_mass[j] );
          h_2dmatch_ptpair->Fill( diphoton_pt[j] );
          h_2dmatch_mpp_over_mgg->Fill( cand.mass()/diphoton_mass[j] );
          h_2dmatch_ypp_minus_ygg->Fill( cand.rapidity() - diphoton_rapidity[j] );
          h_2dmatch_r9single_leadpho->Fill( diphoton_r9_1[j] );
          h_2dmatch_r9single_subleadpho->Fill( diphoton_r9_2[j] );
          //events_list << "2\t" << diphoton_mass[j] << "\t" << diphoton_rapidity[j] << endl;
        }

        //----- matching (diphoton+MET) -----
        h_mggmet_vs_mpp.SetPoint( cand_2tag_id, diphoton_plus_met_mass, cand.mass() );
        h_mggmet_vs_mpp.SetPointError( cand_2tag_id, 0., cand.mass_error() );
        h_yggmet_vs_ypp.SetPoint( cand_2tag_id, diphoton_plus_met_rap, cand.rapidity() );
        h_yggmet_vs_ypp.SetPointError( cand_2tag_id, 0., cand.rapidity_error() );
        if ( mass_match_withmet ) {
          const unsigned short id_match = h_mggmet_vs_mpp_candm.GetN();
          h_mggmet_vs_mpp_candm.SetPoint( id_match, diphoton_plus_met_mass, cand.mass() );
          h_yggmet_vs_ypp_candm.SetPoint( id_match, diphoton_plus_met_rap, cand.rapidity() );
        }
        if ( rap_match_withmet ) {
          const unsigned short id_match = h_mggmet_vs_mpp_candy.GetN();
          h_mggmet_vs_mpp_candy.SetPoint( id_match, diphoton_plus_met_mass, cand.mass() );
          h_yggmet_vs_ypp_candy.SetPoint( id_match, diphoton_plus_met_rap, cand.rapidity() );
        }

        if ( mass_match_withmet && rap_match_withmet ) {
          cout << "--------> matching with met!!!!!" << endl;
          cout << "     " << run_id << ":" << lumisection << ":" << event_number << "\t" << "pt=" << diphoton_pt[j] << ", " << "MET=" << met << endl;
          cout << "      r9: " << diphoton_r9_1[j] << " / " << diphoton_r9_2[j] << endl;
          h_2dmatch_withmet_metvspt.SetPoint( h_2dmatch_withmet_metvspt.GetN(), diphoton_pt[j], met );
          num_match_withmet++;
        }

        //----- matching (diphoton+other objects) -----
        h_mggincl_vs_mpp.SetPoint( cand_2tag_id, diphoton_incl_mass, cand.mass() );
        h_mggincl_vs_mpp.SetPointError( cand_2tag_id, 0., cand.mass_error() );
        h_yggincl_vs_ypp.SetPoint( cand_2tag_id, diphoton_incl_rap, cand.rapidity() );
        h_yggincl_vs_ypp.SetPointError( cand_2tag_id, 0., cand.rapidity_error() );
        if ( mass_match_incl ) {
          const unsigned short id_match = h_mggincl_vs_mpp_candm.GetN();
          h_mggincl_vs_mpp_candm.SetPoint( id_match, diphoton_incl_mass, cand.mass() );
          h_yggincl_vs_ypp_candm.SetPoint( id_match, diphoton_incl_rap, cand.rapidity() );
          //if ( fabs( diphoton_incl_rap-cand.rapidity() )<1.0*cand.rapidity_error() ) {
          //cout << run_id << ":" << lumisection << ":" << event_number << "\t" << "pt=" << dipho_incl.Pt() << "\t" << diphoton_incl_mass << "\t" << diphoton_mass[j] << "\t" << num_ele_included << "\t" << num_mu_included << "\t" << num_jet_included << endl;
          //cout << "--> " << diphoton_eta1[j] << "/" << diphoton_eta2[j] << " ::: " << diphoton_pt1[j] << "\t" << diphoton_pt2[j] << endl;
          //cout << "--> " << diphoton_incl_rap << "\t" << cand.rapidity() << endl;
          //h_pt_incl->Fill( diphoton_pt[j] );
          //}
        }
        if ( rap_match_incl ) {
          const unsigned short id_match = h_mggincl_vs_mpp_candy.GetN();
          h_mggincl_vs_mpp_candy.SetPoint( id_match, diphoton_incl_mass, cand.mass() );
          h_yggincl_vs_ypp_candy.SetPoint( id_match, diphoton_incl_rap, cand.rapidity() );
        }

        if ( mass_match_incl && rap_match_incl ) {
          num_match_incl++;
          cout << "--------> matching with other objects!!!!!" << endl;
          /*cout << " >>> " << diphoton_incl_mass << "/" << cand.mass() << " && " << diphoton_incl_rap << "/" << cand.rapidity() << endl;
          cout << "     " << run_id << ":" << lumisection << ":" << event_number << "\t" << "pt=" << dipho_incl.Pt() << endl;
          cout << " computed with:" << endl;
          cout << num_ele_included << " electron(s)" << endl;
          for ( size_t k=0; k<electrons_list.size(); k++ ) {
            cout << "  ---> " << "electron " << k << ": pt=" << electrons_list[k].second.Pt() << ", eta=" << electrons_list[k].second.Eta() << ", vtx distance(diph vtx)=" << electrons_list[k].first << " cm" << endl;
          }
          cout << num_mu_included << " muon(s)" << endl;
          for ( size_t k=0; k<muons_list.size(); k++ ) {
            cout << "  ---> " << "muon " << k << ": pt=" << muons_list[k].second.Pt() << ", eta=" << muons_list[k].second.Eta() << ", vtx distance(diph vtx)=" << muons_list[k].first << " cm" << endl;
          }
          cout << num_jet_included << " jet(s)" << endl;
          for ( size_t k=0; k<jets_list.size(); k++ ) {
            cout << "  ---> " << "jet " << k << ": pt=" << jets_list[k].second.Pt() << ", eta=" << jets_list[k].second.Eta() << ", vtx distance(diph vtx)=" << jets_list[k].first << " cm" << endl;
          }*/
	}
/*          }
        }*/

	cand_2tag_id++;
      } // loop over diproton candidates

    } // loop over diphotons

    //----- AFTER THE LOOP ON DIPHOTON -----

    h_num_vtx->Fill( num_vertex );
    h_met->Fill( met );
    if ( num_proton_45>0 || num_proton_56>0 ) {
      h_num_vtx_1tag->Fill( num_vertex );
      h_met_1tag->Fill( met );
    }
    if ( num_proton_45>0 && num_proton_56>0 ) {
      h_num_vtx_2tag->Fill( num_vertex );
      h_met_2tag->Fill( met );
    }

    for ( unsigned short j=0; j<num_vertex; j++ ) {
      h_vtxz->Fill( vertex_z[j] );
      if ( num_proton_45>0 || num_proton_56>0 ) { // double tag
        h_vtxz_1tag->Fill( vertex_z[j] );
      }
      if ( num_proton_45>0 && num_proton_56>0 ) { // single tag
        h_vtxz_2tag->Fill( vertex_z[j] );
      }
    }

    h_num_diphoton->Fill( num_diphoton_cand );
    //h_num_diphoton->Fill( num_diphoton );
    h_num_diphoton_1tag->Fill( num_diphoton_cand_1tag );
    h_num_diphoton_2tag->Fill( num_diphoton_cand_2tag );

    h_num_proton->Fill( num_proton );
    h_num_proton_45->Fill( num_proton_45 );
    h_num_proton_56->Fill( num_proton_56 );

  }

  cout << "===== SUMMARY =====" << endl;
  cout << " --> match with diphoton only: " << num_match_diph << endl;
  cout << " --> match with diphoton+MET:  " << num_match_withmet << endl;
  cout << " --> match with diphoton+obj:  " << num_match_incl << endl;


  // plotting stage
  cout << "events: " << num_evts_notag << ", with tag: " << num_evts_with_tag << endl;

  gStyle->SetOptStat( 0 );

  /*const string top_label_str = ( lumi<10. )
    ? Form( "CMS+CTPPS Preliminary 2016, #sqrt{s} = 13 TeV, L = %.2f fb^{-1}", lumi )
    : Form( "CMS+CTPPS Preliminary 2016, #sqrt{s} = 13 TeV, L = %.1f fb^{-1}", lumi );*/
  const string top_label_str = Form( "CMS+CTPPS Preliminary 2016, #sqrt{s} = 13 TeV, L = %.1f fb^{-1}", lumi );
  const char* top_label = top_label_str.c_str();

  const Plotter plt( out_path, top_label );

  {
    //----- preselection level -----

    plt.plot_3hists( "presel_diphoton_mult", h_num_diphoton, h_num_diphoton_1tag, h_num_diphoton_2tag );

    plt.plot_3hists( "presel_diphoton_mass", h_diphoton_mass, h_diphoton_mass_1tag, h_diphoton_mass_2tag );
    plt.plot_3hists( "presel_diphoton_mass_withmet", h_diphoton_mass_withmet, h_diphoton_mass_withmet_1tag, h_diphoton_mass_withmet_2tag );
    plt.plot_3hists( "presel_diphoton_mass_incl", h_diphoton_mass_incl, h_diphoton_mass_incl_1tag, h_diphoton_mass_incl_2tag );
    plt.plot_3hists( "presel_diphoton_pt", h_diphoton_pt, h_diphoton_pt_1tag, h_diphoton_pt_2tag );
    plt.plot_3hists( "presel_diphoton_pt_zoom", h_diphoton_pt_zoom, h_diphoton_pt_zoom_1tag, h_diphoton_pt_zoom_2tag, false );
    plt.plot_3hists( "presel_diphoton_lead_pt", h_diphoton_leadpt, h_diphoton_leadpt_1tag, h_diphoton_leadpt_2tag );
    plt.plot_3hists( "presel_diphoton_sublead_pt", h_diphoton_subleadpt, h_diphoton_subleadpt_1tag, h_diphoton_subleadpt_2tag );
    plt.plot_3hists( "presel_diphoton_vtxz", h_diphoton_vtxz, h_diphoton_vtxz_1tag, h_diphoton_vtxz_2tag );

    /*h_diphoton_leadeta->GetYaxis()->SetRangeUser(0., 55.);
    h_diphoton_subleadeta->GetYaxis()->SetRangeUser(0., 55.);*/
    h_diphoton_leadeta->SetMinimum( 0. );
    h_diphoton_subleadeta->SetMinimum( 0. );

    plt.plot_3hists( "presel_diphoton_lead_eta", h_diphoton_leadeta, h_diphoton_leadeta_1tag, h_diphoton_leadeta_2tag, false );
    plt.plot_3hists( "presel_diphoton_sublead_eta", h_diphoton_subleadeta, h_diphoton_subleadeta_1tag, h_diphoton_subleadeta_2tag, false );
    plt.plot_3hists( "presel_diphoton_dphi", h_diphoton_dphi, h_diphoton_dphi_1tag, h_diphoton_dphi_2tag );
    plt.plot_3hists( "presel_diphoton_dphi_zoom", h_diphoton_dphi_zoom, h_diphoton_dphi_zoom_1tag, h_diphoton_dphi_zoom_2tag, false );
    plt.plot_3hists( "presel_diphoton_rapidity", h_diphoton_rap, h_diphoton_rap_1tag, h_diphoton_rap_2tag );
    plt.plot_3hists( "presel_diphoton_closest_vtx", h_diphoton_closestvtx, h_diphoton_closestvtx_1tag, h_diphoton_closestvtx_2tag );
    plt.plot_3hists( "presel_diphoton_vtx_numtracks", h_diphoton_ntrk, h_diphoton_ntrk_1tag, h_diphoton_ntrk_2tag );
    plt.plot_3hists( "presel_num_vertex", h_num_vtx, h_num_vtx_1tag, h_num_vtx_2tag );
    plt.plot_3hists( "presel_vertexz_allvtx", h_vtxz, h_vtxz_1tag, h_vtxz_2tag );
    plt.plot_3hists( "presel_event_met", h_met, h_met_1tag, h_met_2tag );
    //              JW
    plt.plot_3hists( "electron_pt", h_electron_pt, h_electron_pt_1tag, h_electron_pt_2tag );
    plt.plot_3hists( "muon_pt", h_muon_pt, h_muon_pt_1tag, h_muon_pt_2tag );
    plt.plot_3hists( "jet_pt", h_jet_pt, h_jet_pt_1tag, h_jet_pt_2tag );
    //

    //plt.plot_balances( "presel_diphoton_pt_vs_diproton_mass", "Diphoton p_{T} (GeV)@@Diproton mass (GeV)", &h_ptgg_vs_mpp );

    {
      Plotter::HistsMap hm;
      hm.push_back( std::make_pair( "Sector 45", h_num_proton_45 ) );
      hm.push_back( std::make_pair( "Sector 56", h_num_proton_56 ) );
      hm.push_back( std::make_pair( "Both sides", h_num_proton ) );
      plt.draw_multiplot( "presel_num_proton", hm );
    }
    {
      Plotter::HistsMap hm;
      hm.push_back( std::make_pair( "at 1 mm distance", h_electron_pt_1mm ) );
      hm.push_back( std::make_pair( "at 2 mm distance", h_electron_pt_2mm ) );
      hm.push_back( std::make_pair( "at 5 mm distance", h_electron_pt_5mm ) );
      hm.push_back( std::make_pair( "at 1 cm distance", h_electron_pt_1cm ) );
      plt.draw_multiplot( "presel_electron_pt_close_vertex", hm, true, true );
    }
    {
      Plotter::HistsMap hm;
      hm.push_back( std::make_pair( "at 1 mm distance", h_muon_pt_1mm ) );
      hm.push_back( std::make_pair( "at 2 mm distance", h_muon_pt_2mm ) );
      hm.push_back( std::make_pair( "at 5 mm distance", h_muon_pt_5mm ) );
      hm.push_back( std::make_pair( "at 1 cm distance", h_muon_pt_1cm ) );
      plt.draw_multiplot( "presel_muon_pt_close_vertex", hm, true, true );
    }
    {
      Plotter::HistsMap hm;
      hm.push_back( std::make_pair( "at 1 mm distance", h_jet_pt_1mm ) );
      hm.push_back( std::make_pair( "at 2 mm distance", h_jet_pt_2mm ) );
      hm.push_back( std::make_pair( "at 5 mm distance", h_jet_pt_5mm ) );
      hm.push_back( std::make_pair( "at 1 cm distance", h_jet_pt_1cm ) );
      plt.draw_multiplot( "presel_jet_pt_close_vertex", hm, true, true );
    }
    {
      Plotter::HistsMap hm;
      hm.push_back( std::make_pair( "at 1 mm distance", h_num_vtx_1mm ) );
      hm.push_back( std::make_pair( "at 2 mm distance", h_num_vtx_2mm ) );
      hm.push_back( std::make_pair( "at 5 mm distance", h_num_vtx_5mm ) );
      hm.push_back( std::make_pair( "at 1 cm distance", h_num_vtx_1cm ) );
      plt.draw_multiplot( "presel_num_close_vertex", hm, true, true );
    }

    {
      gr_veto_ele->Scale( 1./tr->GetEntries() );
      gr_veto_mu->Scale( 1./tr->GetEntries() );
      gr_veto_jet->Scale( 1./tr->GetEntries() );
      Plotter::HistsMap gm;
      gm.push_back( std::make_pair( "all electrons", gr_veto_ele ) );
      gm.push_back( std::make_pair( "all muons", gr_veto_mu ) );
      gm.push_back( std::make_pair( "all jets", gr_veto_jet ) );
      plt.draw_multiplot( "presel_veto_eff_dist", gm );
    }
    {
      gr_veto_ele_ptcut->Scale( 1./tr->GetEntries() );
      gr_veto_mu_ptcut->Scale( 1./tr->GetEntries() );
      gr_veto_jet_ptcut->Scale( 1./tr->GetEntries() );
      Plotter::HistsMap gm;
      gm.push_back( std::make_pair( Form( "electrons, p_{T} > %.f GeV", oth_obj_minpt ), gr_veto_ele_ptcut ) );
      gm.push_back( std::make_pair( Form( "muons, p_{T} > %.f GeV", oth_obj_minpt ), gr_veto_mu_ptcut ) );
      gm.push_back( std::make_pair( Form( "jets, p_{T} > %.f GeV", oth_obj_minpt ), gr_veto_jet_ptcut ) );
      plt.draw_multiplot( "presel_veto-ptcut_eff_dist", gm );
    }

    cout << "total candidates: " << h_mgg_vs_mpp.Integral() << endl
         << " -> with mass matching: " << h_mgg_vs_mpp_candm.Integral() << endl
         << " -> with rapidity matching: " << h_mgg_vs_mpp_candy.Integral() << endl;

    {
      Plotter::GraphsMap gm;
      gm.push_back( make_pair( "Exclusive diphotons", &h_mgg_vs_mpp ) );
      gm.push_back( make_pair( Form( "Mass match (%.0f#sigma)", n_sigma_2d ), &h_mgg_vs_mpp_candm ) );
      gm.push_back( make_pair( Form( "Rapidity match (%.0f#sigma)", n_sigma_2d ), &h_mgg_vs_mpp_candy ) );
      plt.plot_balances( "presel_mass_balance", "Diphoton mass (GeV)@@Diproton missing mass (GeV)", gm, min_mass_diphoton, 2500. );
    }
    {
      Plotter::GraphsMap gm;
      gm.push_back( make_pair( "Exclusive diphotons", &h_mggmet_vs_mpp ) );
      gm.push_back( make_pair( Form( "Mass match (%.0f#sigma)", n_sigma_2d ), &h_mggmet_vs_mpp_candm ) );
      gm.push_back( make_pair( Form( "Rapidity match (%.0f#sigma)", n_sigma_2d ), &h_mggmet_vs_mpp_candy ) );
      plt.plot_balances( "presel_mass_balance_withmet", "Diphoton + #slash{E}_{T} mass (GeV)@@Diproton missing mass (GeV)", gm, min_mass_diphoton, 3000. );
    }
    {
      Plotter::GraphsMap gm;
      gm.push_back( make_pair( "Exclusive diphotons", &h_mggincl_vs_mpp ) );
      gm.push_back( make_pair( Form( "Mass match (%.0f#sigma)", n_sigma_2d ), &h_mggincl_vs_mpp_candm ) );
      gm.push_back( make_pair( Form( "Rapidity match (%.0f#sigma)", n_sigma_2d ), &h_mggincl_vs_mpp_candy ) );
      plt.plot_balances( "presel_mass_balance_incl", "Diphoton + surrounding objects mass (GeV)@@Diproton missing mass (GeV)", gm, min_mass_diphoton, 3000. );
    }
    {
      Plotter::GraphsMap gm;
      gm.push_back( make_pair( "Exclusive diphotons", &h_ygg_vs_ypp ) );
      gm.push_back( make_pair( Form( "Mass match (%.0f#sigma)", n_sigma_2d ), &h_ygg_vs_ypp_candm ) );
      gm.push_back( make_pair( Form( "Rapidity match (%.0f#sigma)", n_sigma_2d ), &h_ygg_vs_ypp_candy ) );
      plt.plot_balances( "presel_rapidity_balance", "Diphoton rapidity@@Diproton rapidity", gm, -3., 3. );
    }
    {
      Plotter::GraphsMap gm;
      gm.push_back( make_pair( "Exclusive diphotons", &h_yggmet_vs_ypp ) );
      gm.push_back( make_pair( Form( "Mass match (%.0f#sigma)", n_sigma_2d ), &h_yggmet_vs_ypp_candm ) );
      gm.push_back( make_pair( Form( "Rapidity match (%.0f#sigma)", n_sigma_2d ), &h_yggmet_vs_ypp_candy ) );
      plt.plot_balances( "presel_rapidity_balance_withmet", "Diphoton + #slash{E}_{T} rapidity@@Diproton rapidity", gm, -1., 1. );
    }
    {
      Plotter::GraphsMap gm;
      gm.push_back( make_pair( "Exclusive diphotons", &h_yggincl_vs_ypp ) );
      gm.push_back( make_pair( Form( "Mass match (%.0f#sigma)", n_sigma_2d ), &h_yggincl_vs_ypp_candm ) );
      gm.push_back( make_pair( Form( "Rapidity match (%.0f#sigma)", n_sigma_2d ), &h_yggincl_vs_ypp_candy ) );
      plt.plot_balances( "presel_rapidity_balance_incl", "Diphoton + surrounding objects rapidity@@Diproton rapidity", gm, -1.5, 1.5 );
      plt.plot_balances( "presel_rapidity_balance_incl_zoom", "Diphoton + surrounding objects rapidity@@Diproton rapidity", gm, -0.5, 0.5 );
    }

    //----- EXCLUSIVE DIPHOTONS SELECTION -----

    //----- matching plots with diphoton only -----
    {
      Plotter::HistsMap hm;
      hm.push_back( std::make_pair( "at 1 mm distance", h_num_vtx_1mm_excl ) );
      hm.push_back( std::make_pair( "at 2 mm distance", h_num_vtx_2mm_excl ) );
      hm.push_back( std::make_pair( "at 5 mm distance", h_num_vtx_5mm_excl ) );
      hm.push_back( std::make_pair( "at 1 cm distance", h_num_vtx_1cm_excl ) );
      plt.draw_multiplot( "excl_num_close_vertex", hm, true );
    }
  }

  {
    Canvas c( "met_x_vs_y", top_label );
    h_metx_vs_mety->Draw( "colz" );
    h_metx_vs_mety_1tag->Draw( "p same" );
    h_metx_vs_mety_1tag->SetMarkerStyle( 26 );
    h_metx_vs_mety_1tag->SetMarkerColor( kBlack );
    h_metx_vs_mety_2tag->Draw( "p same" );
    h_metx_vs_mety_2tag->SetMarkerStyle( 24 );
    h_metx_vs_mety_2tag->SetMarkerColor( kRed );
    c.SetLegendX1( 0.15 );
    c.SetLegendY1( 0.15 );
    c.AddLegendEntry( h_metx_vs_mety_1tag, "#geq 1 proton tags", "p" );
    c.AddLegendEntry( h_metx_vs_mety_2tag, "#geq 2 proton tags", "p" );
    c.Prettify( h_metx_vs_mety );
    c.Save( "pdf,png", out_path );
  }
  {
    Canvas c( "diphoton_mass_vs diphotonmet_mass", top_label );
    h_mggmet_vs_mgg->Draw( "colz" );
    //c.DrawDiagonal( h_mggmet_vs_mgg->GetXaxis()->GetXmin(), h_mggmet_vs_mgg->GetXaxis()->GetXmax() );
    h_mggmet_vs_mgg->Draw( "colz same" );
    c.Prettify( h_mggmet_vs_mgg );
    c.Save( "pdf,png", out_path );
  }
  {
    Canvas c( "diphoton_pt_vs_met", top_label );
    h_met_vs_pt->Draw( "colz" );
    h_met_vs_pt_2tag->Draw( "p same" );
    h_met_vs_pt_2tag->SetMarkerStyle( 24 );
    h_met_vs_pt_2tag->SetMarkerColor( kRed );
    c.AddLegendEntry( h_met_vs_pt_2tag, "#geq 2 proton tags", "p" );
    //c.DrawDiagonal( h_met_vs_pt->GetXaxis()->GetXmin(), h_met_vs_pt->GetXaxis()->GetXmax() );
    c.Prettify( h_met_vs_pt );
    c.Save( "pdf,png", out_path );
  }
  /*{
    Canvas c( "diphoton_pt_vs_diproton_mass", top_label );
    h_ptgg_vs_mpp->Draw( "colz" );
    c.DrawDiagonal( h_ptgg_vs_mpp->GetXaxis()->GetXmin(), h_ptgg_vs_mpp->GetXaxis()->GetXmax() );
    c.Prettify( h_ptgg_vs_mpp );
    c.Save( "pdf", out_path );
    c.Save( "png", out_path );
  }*/
  {
    Canvas c( "match2d_singlephoton_r9", top_label );
    h_2dmatch_r9single_leadpho->Sumw2();
    h_2dmatch_r9single_leadpho->Draw();
    h_2dmatch_r9single_leadpho->SetMarkerStyle( 20 );
    h_2dmatch_r9single_leadpho->SetLineColor( kBlack );
    h_2dmatch_r9single_subleadpho->Sumw2();
    h_2dmatch_r9single_subleadpho->Draw( "same" );
    h_2dmatch_r9single_subleadpho->SetMarkerStyle( 24 );
    c.AddLegendEntry( h_2dmatch_r9single_leadpho, "Leading photon" );
    c.AddLegendEntry( h_2dmatch_r9single_subleadpho, "Subleading photon" );
    h_2dmatch_r9single_subleadpho->SetLineColor( kBlack );
    c.Prettify( h_2dmatch_r9single_leadpho );
    c.Save( "pdf,png", out_path );
  }
  {
    Canvas c( "match2d_diphoton_pt", top_label );
    h_2dmatch_ptpair->Sumw2();
    h_2dmatch_ptpair->Draw();
    h_2dmatch_ptpair->SetMarkerStyle( 20 );
    h_2dmatch_ptpair->SetLineColor( kBlack );
    c.Prettify( h_2dmatch_ptpair );
    c.Save( "pdf,png", out_path );
  }
  {
    Canvas c( "match2d_diphoton_mass", top_label );
    h_2dmatch_mpair->Sumw2();
    h_2dmatch_mpair->Draw();
    h_2dmatch_mpair->SetMarkerStyle( 20 );
    h_2dmatch_mpair->SetLineColor( kBlack );
    c.Prettify( h_2dmatch_mpair );
    c.Save( "pdf,png", out_path );
  }
  {
    Canvas c( "match2d_mass_ratio", top_label );
    h_2dmatch_mpp_over_mgg->Sumw2();
    h_2dmatch_mpp_over_mgg->Draw();
    h_2dmatch_mpp_over_mgg->SetMarkerStyle( 20 );
    h_2dmatch_mpp_over_mgg->SetLineColor( kBlack );
    c.Prettify( h_2dmatch_mpp_over_mgg );
    c.Save( "pdf,png", out_path );
  }
  {
    Canvas c( "match2d_rapidity_difference", top_label );
    h_2dmatch_ypp_minus_ygg->Sumw2();
    h_2dmatch_ypp_minus_ygg->Draw();
    h_2dmatch_ypp_minus_ygg->SetMarkerStyle( 20 );
    h_2dmatch_ypp_minus_ygg->SetLineColor( kBlack );
    c.Prettify( h_2dmatch_ypp_minus_ygg );
    c.Save( "pdf,png", out_path );
  }
  {
    Plotter::GraphsMap gm; gm.push_back( make_pair( "Matching with #slash{E}_{T}", &h_2dmatch_withmet_metvspt ) );
    //plt.plot_balances( "match2d_withmet_ptvsmetbalance", "Diphoton p_{T} (GeV)@@#slash{E}_{T} (GeV)", &h_2dmatch_withmet_metvspt, 0., 150. );
    plt.plot_balances( "match2d_withmet_ptvsmetbalance", "Diphoton p_{T} (GeV)@@#slash{E}_{T} (GeV)", gm, 0., 150. );
  }
  {
    Canvas c( "match2d_inclusive_ptpair", top_label );
    h_pt_incl->Draw();
    c.Prettify( h_pt_incl );
    c.Save( "pdf,png", out_path );
  }
}

float photon_rel_energy_scale( const float& pt, const float& eta, const float& r9 )
{
  if ( r9<0.94 ) return 0.004;
  if ( fabs( eta )<1.4442 ) return 0.002;
  if ( fabs( eta )>1.566 ) return 0.005;
  return 1.;
}

