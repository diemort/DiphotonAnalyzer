#include "Canvas.h"
#include "pot_alignment.h"
#include "xi_reconstruction.h"
#include "diproton_candidate.h"

#include "DiphotonAnalyzer/TreeProducer/interface/TreeEvent.h"

#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH2.h"

#include <fstream>

#define MAX_PROTONS 30
#define MAX_DIPH 20
#define MAX_ELE 50
#define MAX_MUON 50
#define MAX_JET 200

TGraphAsymmErrors* asym_error_bars( const TH1D* hist ) {
  const double alpha = 1 - 0.6827;
  auto g = new TGraphAsymmErrors( hist );
  for ( int i = 0; i < g->GetN(); ++i ) {
    int N = g->GetY()[i];
    if ( N == 0 ) continue; //FIXME skip the empty bins??
    double L = ( N == 0 ) ? 0. : ( ROOT::Math::gamma_quantile( 0.5*alpha, N, 1. ) );
    double U = ( N == 0 ) ? ( ROOT::Math::gamma_quantile_c( alpha, N+1, 1 ) ) : ( ROOT::Math::gamma_quantile_c( 0.5*alpha, N+1, 1 ) );
    g->SetPointEXlow( i, 0. ); //FIXME
    g->SetPointEXhigh( i, 0. ); //FIXME
    g->SetPointEYlow( i, N-L );
    g->SetPointEYhigh( i, U-N );
  }
  return g;
}

const string up_label = "9.4 fb^{-1} (13 TeV)";
//const string up_label = "2016B - 4.3 fb^{-1} (13 TeV)";

const float max_xi = 0.25;
struct track_t {
  track_t() : xi( 0. ), err_xi( 0. ), x( 0. ), y( 0. ) {}
  track_t( float xi, float err_xi, float x, float y ) : xi( xi ), err_xi( err_xi ), x( x ), y( y ) {}
  float xi, err_xi, x, y;
};

void plot_matching( double num_sigma, const char* name, TGraphErrors&, TGraphErrors&, TGraphErrors&, TGraphErrors&, double min, double max, bool right_leg = true );
vector<pair<float,float> > merge_nearfar( const vector<track_t>& near_tracks, const vector<track_t>& far_tracks, float xdiff_cut=0.01 );

map<unsigned short,float> pots_accept = { { 2, 0.033 }, { 3, 0.024 }, { 102, 0.050 }, { 103, 0.037 } };
//map<unsigned short,float> pots_accept = { { 2, 0.067 }, { 3, 0.066 }, { 102, 0.070 }, { 103, 0.067 } }; // less than 10% radiation damage
map<unsigned short,float> pots_accept_match = { { 2, 0.034 }, { 3, 0.023 }, { 102, 0.042 }, { 103, 0.032 } };
map<unsigned short,float> pots_accept_90pc = { { 2, 0.067 }, { 3, 0.066 }, { 102, 0.070 }, { 103, 0.067 } }; // less than 10% radiation damage
map<unsigned short,const char*> pots_names = { { 2, "45N" }, { 3, "45F" }, { 102, "56N" }, { 103, "56F" } };

//pots_accept = pots_accept_90pc;

void massrap_matcher( double num_sigma = 2., const char* sample = "/eos/cms/store/group/phys_pps/diphoton/DoubleEG/ntuple-Run2016BCG_94Xrereco_v1.root" )
{
  TFile f( sample );
  TTree* tr = dynamic_cast<TTree*>( f.Get( "ntp" ) );

  xi_reco::load_optics_file( "TreeProducer/data/optics_17may22.root" );
  pot_align::load_file( "TreeProducer/data/alignment_collection_v2.out" );

  //const float rel_err_xi_gg = 0.039;
  const double rel_err_mass = 0.02, rel_err_rap = 0.074;

  gggg::TreeEvent ev;
  ev.attach( tr, true );

  TGraphErrors gr_mass_massrapmatch, gr_mass_massmatch, gr_mass_rapmatch, gr_mass_nomatch;
  TGraphErrors gr_rap_massrapmatch, gr_rap_massmatch, gr_rap_rapmatch, gr_rap_nomatch;

  unsigned int num_massmatch = 0, num_rapmatch = 0, num_massrapmatch = 0, num_nomatch = 0;

  TH1D* h_mass_all = new TH1D( "mass_all", "Diproton missing mass@@Events@@GeV", 12, 200., 2000. );
  TH1D* h_rap_all = new TH1D( "rap_all", "Diproton rapidity@@Events", 20, -1., 1. );
  map<unsigned short,TH1D*> m_h_xi;
  for ( const auto& pid : pots_names )
    m_h_xi[pid.first] = new TH1D( Form( "h_xi_%d", pid.first ), Form( "#xi (%s)@@Events", pid.second ), 40, 0., 0.24 );
  TH1D* h_num_45 = new TH1D( "num_45", "Forward tracks multiplicity@@Events", 3, -0.5, 2.5 );
  TH1D* h_num_56 = (TH1D*)h_num_45->Clone( "num_56" );
  TH1D* h_num_sect = new TH1D( "num_sect", ";;Events", 4, -0.5, 3.5 );
  /*auto h_massratio = new TH1D( "mass_ratio", "m_{pp}/m_{#gamma#gamma}@@Events@@?.g", 20, -2., 4. );
  auto h_rapdiff = new TH1D( "rap_diff", "y_{pp}-y_{#gamma#gamma}@@Events@@?.g", 20, -2.5, 2.5 );*/
  auto h_massratio = new TH1D( "mass_ratio", "(m_{#gamma#gamma}-m_{pp})/#sigma(m_{#gamma#gamma}-m_{pp})@@Events@@?.g", 20, -20., 20. );
  auto h_rapdiff = new TH1D( "rap_diff", "(y_{#gamma#gamma}-y_{pp})/#sigma(y_{#gamma#gamma}-y_{pp})@@Events@@?.g", 20, -20., 20. );

  const unsigned long long num_events = tr->GetEntriesFast();
  for ( unsigned long long i = 0; i < num_events; ++i ) {
    tr->GetEntry( i );

    if ( ev.hlt_accept[0] == 0 ) continue;

    //cout << "event " << i << ": " << ev.num_fwd_track << " proton tracks, " << ev.num_diphoton << " diphoton candidates" << endl;

    // first loop to identify the tracks and their respective pot

    auto align = pot_align::get_alignments( ev.fill_number );

    vector<track_t> xi_45n, xi_45f, xi_56n, xi_56f;
    vector<pair<float,float> > xi_45, xi_56;

    for ( unsigned short j = 0; j < ev.num_fwd_track; ++j ) {
      //if ( ev.proton_track_method[j] != 1 ) continue; // only use multi-pot reconstruction
      const unsigned short pot_id = 100*ev.fwd_track_arm[j]+ev.fwd_track_pot[j];
      const auto& al = align[pot_id];

      double xi, xi_err;
      xi_reco::reconstruct( ev.fwd_track_x[j]+al.x, ev.fwd_track_arm[j], ev.fwd_track_pot[j], xi, xi_err );

      //----- reconstruct the kinematics
      //if ( ev.proton_track_xi[j] < pots_accept[pot_id] )
      //if ( xi < pots_accept[pot_id] )
      if ( xi < pots_accept_90pc[pot_id] || xi > 0.15 ) continue; //FIXME FIXME FIXME

      //----- associate each track to a RP
      if      ( ev.fwd_track_arm[j] == 0 && ev.fwd_track_pot[j] == 2 ) xi_45n.emplace_back( xi, xi_err, ev.fwd_track_x[j]+al.x, ev.fwd_track_y[j]-al.y );
      else if ( ev.fwd_track_arm[j] == 0 && ev.fwd_track_pot[j] == 3 ) xi_45f.emplace_back( xi, xi_err, ev.fwd_track_x[j]+al.x, ev.fwd_track_y[j]-al.y );
      else if ( ev.fwd_track_arm[j] == 1 && ev.fwd_track_pot[j] == 2 ) xi_56n.emplace_back( xi, xi_err, ev.fwd_track_x[j]+al.x, ev.fwd_track_y[j]-al.y );
      else if ( ev.fwd_track_arm[j] == 1 && ev.fwd_track_pot[j] == 3 ) xi_56f.emplace_back( xi, xi_err, ev.fwd_track_x[j]+al.x, ev.fwd_track_y[j]-al.y );
//NEW      if ( ev.proton_track_arm[j] == 0 ) xi_45.emplace_back( ev.proton_track_xi[j], ev.proton_track_xi[j]*0.1 ); //FIXME FIXME FIXME
//NEW      else if ( ev.proton_track_arm[j] == 1 ) xi_56.emplace_back( ev.proton_track_xi[j], ev.proton_track_xi[j]*0.1 ); //FIXME FIXME FIXME
    }

    //----- merge 2 tracks in one if N-F pot content is similar

    const float xdiff_cut = 0.01;

    //--- sector 45

    xi_45 = merge_nearfar( xi_45n, xi_45f, xdiff_cut );
    xi_56 = merge_nearfar( xi_56n, xi_56f, xdiff_cut );

    /*//FIXME FIXME
    for ( const auto& trk : xi_45n ) xi_45.emplace_back( trk.xi, trk.err_xi );
    for ( const auto& trk : xi_45f ) xi_45.emplace_back( trk.xi, trk.err_xi );
    for ( const auto& trk : xi_56n ) xi_56.emplace_back( trk.xi, trk.err_xi );
    for ( const auto& trk : xi_56f ) xi_56.emplace_back( trk.xi, trk.err_xi );
    //FIXME FIXME*/

    //---- identify the diproton candidates

    vector<diproton_candidate_t> candidates;
    for ( const auto trk45 : xi_45 )
      for ( const auto trk56 : xi_56 )
        candidates.emplace_back( trk45.first, trk45.second, trk56.first, trk56.second );
    //cout << candidates.size() << " diproton candidate(s) in total!" << endl;

    const float eta_cut = 2.5, min_etaveto = 1.4442, max_etaveto = 1.566;

    //----- identify the diphoton candidates

    bool has_diph_cand = false;

    for ( unsigned short j = 0; j < ev.num_diphoton; ++j ) {

      //----- photon quality cuts

      // EB: 0 < |eta| < 1.4442
      // EE: |eta| > 1.566
      unsigned short ev_class = gggg::TreeEvent::invalid;
      /*if ( ev.diphoton_eta1[j] < min_etaveto && ev.diphoton_eta2[j] > max_etaveto ) ev_class = gggg::TreeEvent::ebee;
      if ( ev.diphoton_eta1[j] > max_etaveto && ev.diphoton_eta2[j] < min_etaveto ) ev_class = gggg::TreeEvent::ebee;
      else if ( ev.diphoton_eta1[j] < min_etaveto && ev.diphoton_eta2[j] < min_etaveto ) ev_class = gggg::TreeEvent::ebeb;
      else if ( ev.diphoton_eta1[j] > max_etaveto && ev.diphoton_eta2[j] > max_etaveto ) ev_class = gggg::TreeEvent::eeee;*/
      if ( ev.diphoton_ele_veto1[j] != 1 ) continue;
      if ( ev.diphoton_ele_veto2[j] != 1 ) continue;
      if ( fabs( ev.diphoton_eta1[j] ) <= eta_cut && fabs( ev.diphoton_eta2[j] ) <= eta_cut ) {
        if ( fabs( ev.diphoton_eta1[j] ) < min_etaveto && fabs( ev.diphoton_eta2[j] ) > max_etaveto ) ev_class = gggg::TreeEvent::ebee;
        else if ( fabs( ev.diphoton_eta2[j] ) < min_etaveto && fabs( ev.diphoton_eta1[j] ) > max_etaveto ) ev_class = gggg::TreeEvent::ebee;
        else if ( fabs( ev.diphoton_eta1[j] ) < min_etaveto && fabs( ev.diphoton_eta2[j] ) < min_etaveto ) ev_class = gggg::TreeEvent::ebeb;
        else if ( fabs( ev.diphoton_eta1[j] ) > max_etaveto && fabs( ev.diphoton_eta2[j] ) > max_etaveto ) ev_class = gggg::TreeEvent::eeee;
      }
      //----- only keep EBEE and EBEB diphoton events
      if ( ev_class == gggg::TreeEvent::invalid ) continue;
      if ( ev_class == gggg::TreeEvent::eeee ) continue; //FIXME FIXME

      if ( ev.diphoton_pt1[j] < 75. ) continue;
      if ( ev.diphoton_pt2[j] < 75. ) continue;
      if ( fabs( ev.diphoton_eta1[j] ) > eta_cut || ( fabs( ev.diphoton_eta1[j] ) > min_etaveto && fabs( ev.diphoton_eta1[j] ) < max_etaveto ) ) continue;
      if ( fabs( ev.diphoton_eta2[j] ) > eta_cut || ( fabs( ev.diphoton_eta2[j] ) > min_etaveto && fabs( ev.diphoton_eta2[j] ) < max_etaveto ) ) continue;
      if ( ev.diphoton_r91[j] < 0.94 ) continue;
      if ( ev.diphoton_r92[j] < 0.94 ) continue;
      if ( ev.diphoton_mass[j] < 350. ) continue;

      //----- back-to-back photons

      if ( 1.-fabs( ev.diphoton_dphi[j] )/M_PI > 0.005 ) continue;
      //if ( 1.-fabs( ev.diphoton_dphi[j] )/M_PI > 0.1 ) continue;

      has_diph_cand = true;

      const float xip = ( ev.diphoton_pt1[j]*exp( +ev.diphoton_eta1[j] ) + ev.diphoton_pt2[j]*exp( +ev.diphoton_eta2[j] ) ) / sqrt_s,
                  xim = ( ev.diphoton_pt1[j]*exp( -ev.diphoton_eta1[j] ) + ev.diphoton_pt2[j]*exp( -ev.diphoton_eta2[j] ) ) / sqrt_s;

      //if ( xi < pots_accept_90pc[pot_id] || xi > 0.15 ) continue; //FIXME FIXME FIXME
      if ( ( xim < pots_accept_90pc[102] && xim < pots_accept_90pc[103] ) || xim > 0.15 ) continue; //FIXME FIXME
      if ( ( xip < pots_accept_90pc[2] && xip < pots_accept_90pc[3] ) || xip > 0.15 ) continue;

      //----- search for associated leptons

      /*const double lepton_mindist = 5.0; // in centimeters
      TVector3 diph_vtx( ev.diphoton_vertex_x[j], ev.diphoton_vertex_y[j], ev.diphoton_vertex_z[j] );
      unsigned short num_close_ele = 0, num_close_muon = 0;
      for ( unsigned short k = 0; k < ev.num_electron; ++k ) {
        TVector3 ele_vtx( electron_vtx_x[k], electron_vtx_y[k], electron_vtx_z[k] );
        const float ele_dist = ( ele_vtx-diph_vtx ).Mag();
        if ( ele_dist < lepton_mindist ) num_close_ele++;
      }
      for ( unsigned short k = 0; k < ev.num_muon; ++k ) {
        TVector3 mu_vtx( muon_vtx_x[k], muon_vtx_y[k], muon_vtx_z[k] );
        const float mu_dist = ( mu_vtx-diph_vtx ).Mag();
        if ( mu_dist < lepton_mindist ) num_close_muon++;
      }*/

      //----- search for associated jets

      TLorentzVector pho1, pho2;
      pho1.SetPtEtaPhiM( ev.diphoton_pt1[j], ev.diphoton_eta1[j], ev.diphoton_phi1[j], 0. );
      pho2.SetPtEtaPhiM( ev.diphoton_pt2[j], ev.diphoton_eta2[j], ev.diphoton_phi2[j], 0. );
      TLorentzVector cms = pho1+pho2;

      //const double min_extrajet_pt = 500.;
      /*unsigned short num_associated_jet = 0;
      for ( unsigned short k = 0; k < ev.num_jet; ++k ) {
        if ( jet_dipho_match[k] != j ) continue;
        //if ( jet_pt[k] > min_extrajet_pt ) num_associated_jet++;
        TLorentzVector jet;
        jet.SetPtEtaPhiE( jet_pt[k], jet_eta[k], jet_phi[k], jet_energy[k] );
        cms += jet;
      }*/

      //----- another batch of "exclusivity" cuts

      /*if ( num_close_ele > 0 ) continue;
      if ( num_close_muon > 0 ) continue;*/
      //if ( num_associated_jet > 0 ) continue;

      //----- reconstruct the energy loss from central system

      float diphoton_mass_error = ev.diphoton_mass[j]*rel_err_mass;
      float diphoton_rapidity_error = fabs( ev.diphoton_rapidity[j] )*rel_err_rap;

      for ( const auto cand : candidates ) {
        h_mass_all->Fill( cand.mass() );
        h_rap_all->Fill( cand.rapidity() );
        bool mass_match = is_matched( num_sigma, cms.M(), cand.mass(), diphoton_mass_error, cand.mass_error() );
        bool rap_match = is_matched( num_sigma, cms.Rapidity(), cand.rapidity(), diphoton_rapidity_error, cand.rapidity_error() );
        /*h_massratio->Fill( cand.mass()/cms.M() );
        h_rapdiff->Fill( cand.rapidity()-cms.Rapidity() );*/
        h_massratio->Fill( ( cms.M()-cand.mass() )/std::hypot( diphoton_mass_error, cand.mass_error() ) );
        h_rapdiff->Fill( ( cms.Rapidity()-cand.rapidity() )/std::hypot( diphoton_rapidity_error, cand.rapidity_error() ) );
        if ( mass_match && rap_match ) {
          cout << "@@@ DOUBLE TAGGING" << endl;
          cout << "event:" << ev.run_id << ":" << ev.lumisection << ":" << ev.event_number << endl;
          cout << "masses: central system: " << cms.M() << ", diphoton: " << ev.diphoton_mass[j] << " +/- " << diphoton_mass_error << ", diproton: " << cand.mass() << " +/- " << cand.mass_error() << endl;
          cout << "rapidities: central system: " << cms.Rapidity() << ", diphoton: " << ev.diphoton_rapidity[j] << " +/- " << diphoton_rapidity_error << ", diproton: " << cand.rapidity() << " +/- " << cand.rapidity_error() << endl;
/*          cout << "xi-per-pot: " << "45-near: " << ( xi_45n.size() > 0 ? xi_45n[0].xi : -1. )
               << ", 45-far: " << ( xi_45f.size() > 0 ? xi_45f[0].xi : -1. )
               << ", 56-near: " << ( xi_56n.size() > 0 ? xi_56n[0].xi : -1. )
               << ", 56-far: " << ( xi_56f.size() > 0 ? xi_56f[0].xi : -1. ) << endl;*/
          gr_mass_massrapmatch.SetPoint( num_massrapmatch, cand.mass(), cms.M() );
          gr_mass_massrapmatch.SetPointError( num_massrapmatch, cand.mass_error(), diphoton_mass_error );
          gr_rap_massrapmatch.SetPoint( num_massrapmatch, cand.rapidity(), cms.Rapidity() );
          gr_rap_massrapmatch.SetPointError( num_massrapmatch, cand.rapidity_error(), diphoton_rapidity_error );
          num_massrapmatch++;
        }
        else if ( mass_match ) { // only match in mass
          gr_mass_massmatch.SetPoint( num_massmatch, cand.mass(), cms.M() );
          gr_mass_massmatch.SetPointError( num_massmatch, cand.mass_error(), diphoton_mass_error );
          gr_rap_massmatch.SetPoint( num_massmatch, cand.rapidity(), cms.Rapidity() );
          gr_rap_massmatch.SetPointError( num_massmatch, cand.rapidity_error(), diphoton_rapidity_error );
          num_massmatch++;
        }
        else if ( rap_match ) { // only match in rapidity
          gr_mass_rapmatch.SetPoint( num_rapmatch, cand.mass(), cms.M() );
          gr_mass_rapmatch.SetPointError( num_rapmatch, cand.mass_error(), diphoton_mass_error );
          gr_rap_rapmatch.SetPoint( num_rapmatch, cand.rapidity(), cms.Rapidity() );
          gr_rap_rapmatch.SetPointError( num_rapmatch, cand.rapidity_error(), diphoton_rapidity_error );
          num_rapmatch++;
        }
        else { // no matching at all
          gr_mass_nomatch.SetPoint( num_nomatch, cand.mass(), cms.M() );
          gr_mass_nomatch.SetPointError( num_nomatch, cand.mass_error(), diphoton_mass_error );
          gr_rap_nomatch.SetPoint( num_nomatch, cand.rapidity(), cms.Rapidity() );
          gr_rap_nomatch.SetPointError( num_nomatch, cand.rapidity_error(), diphoton_rapidity_error );
          num_nomatch++;
        }
        cout << "matching: " << mass_match << "\t" << rap_match << endl;
      }
    } // loop on diphotons
    if ( has_diph_cand ) {
      unsigned short num_45 = 0, num_56 = 0;
      if ( !xi_45n.empty() )
        for ( const auto& m : xi_45n ) {
          m_h_xi[2]->Fill( m.xi );
          num_45++;
        }
      if ( !xi_45f.empty() )
        for ( const auto& m : xi_45f ) {
          m_h_xi[3]->Fill( m.xi );
          num_45++;
        }
      if ( !xi_56n.empty() )
        for ( const auto& m : xi_56n ) {
          m_h_xi[102]->Fill( m.xi );
          num_56++;
        }
      if ( !xi_56f.empty() )
        for ( const auto& m : xi_56f ) {
          m_h_xi[103]->Fill( m.xi );
          num_56++;
        }
      h_num_45->Fill( num_45 );
      h_num_56->Fill( num_56 );
      if ( num_45 > 0 && num_56 > 0 )
        h_num_sect->Fill( 3 );
      else if ( num_45 > 0 )
        h_num_sect->Fill( 1 );
      else if ( num_56 > 0 )
        h_num_sect->Fill( 2 );
      else
        h_num_sect->Fill( 0 );
    }
  } // loop on events
cout << "in plot:\n\t" << "not matching: " << num_nomatch << "\n\tmass match: " << num_massmatch << "\n\trap match: " << num_rapmatch << "\n\tboth match: " << num_massrapmatch << endl;

  //----- plotting part

  gr_mass_massrapmatch.SetTitle( "m_{pp} (GeV)@@m_{#gamma#gamma} (GeV)" );
  gr_rap_massrapmatch.SetTitle( "y_{pp}@@y_{#gamma#gamma}" );

  plot_matching( num_sigma, "2d_massmatch", gr_mass_nomatch, gr_mass_rapmatch, gr_mass_massmatch, gr_mass_massrapmatch, 300., 1800. );
  //plot_matching( num_sigma, "2d_rapmatch", gr_rap_nomatch, gr_rap_rapmatch, gr_rap_massmatch, gr_rap_massrapmatch, -3., 3., 0.15 );
  plot_matching( num_sigma, "2d_rapmatch", gr_rap_nomatch, gr_rap_rapmatch, gr_rap_massmatch, gr_rap_massrapmatch, -2., 2., false );

  for ( auto& nh : map<const char*,TH1D*>{ { "1d_massmatch", h_mass_all }, { "1d_rapmatch", h_rap_all } } ) {
    Canvas c( nh.first, up_label.c_str(), "Preliminary" );
    nh.second->Sumw2();
    nh.second->Draw();
    c.Prettify( nh.second );
    c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
  }
  {
    gStyle->SetOptStat( 0 );
    Canvas c( "elastic_numsect", up_label.c_str(), "Preliminary" );
    h_num_sect->Sumw2();
    h_num_sect->Draw( "p" );
    h_num_sect->SetLineColor( kBlack );
    h_num_sect->SetLineWidth( 2 );
    h_num_sect->SetMarkerStyle( 20 );
    auto axis = h_num_sect->GetXaxis();
    axis->SetBinLabel( 1, "No forward track" );
    axis->SetBinLabel( 2, "Sector 45 only" );
    axis->SetBinLabel( 3, "Sector 56 only" );
    axis->SetBinLabel( 4, "Both sectors" );
    h_num_sect->SetMinimum( 0. );
    c.Prettify( h_num_sect );
    c.SetGrid( 0, 1 );
    PaveText::topLabel( "Elastic selection" );
    c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
  }
  {
    Canvas c( "elastic_numtrks", up_label.c_str(), "Preliminary" );
    THStack hs;
    h_num_45->SetLineColor( kBlue+1 );
    h_num_45->SetLineWidth( 2 );
    h_num_45->SetMarkerStyle( 24 );
    h_num_45->SetMarkerColor( kBlue+1 );
    h_num_56->SetLineColor( kRed+1 );
    h_num_56->SetLineWidth( 2 );
    h_num_56->SetMarkerStyle( 25 );
    h_num_56->SetMarkerColor( kRed+1 );
    c.AddLegendEntry( h_num_45, "Sector 45" );
    c.AddLegendEntry( h_num_56, "Sector 56" );
    hs.Add( h_num_45 );
    hs.Add( h_num_56 );
    hs.Draw( "e,nostack" );
    hs.GetHistogram()->SetTitle( h_num_45->GetTitle() );
    c.Prettify( hs.GetHistogram() );
    hs.GetHistogram()->GetXaxis()->SetNdivisions( 5 );
    c.SetGrid( 0, 1 );
    PaveText::topLabel( "Elastic selection" );
    c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
  }
  for ( auto& nh : map<const char*,TH1D*>{ { "mass_ratio", h_massratio }, { "rapidity_difference", h_rapdiff } } ) {
    gStyle->SetOptStat( 0 );
    Canvas c( nh.first, up_label.c_str(), "Preliminary" );
    nh.second->Sumw2();
    nh.second->Draw();
    nh.second->SetMarkerStyle( 24 );
    nh.second->SetLineColor( kBlack );
    TLine lim;
    lim.SetLineColor( kRed+1 );
    lim.SetLineStyle( 2 );
    lim.SetLineWidth( 3 );
    //const double min = nh.second->GetMinimum(), max = nh.second->GetMaximum()*1.55;
    //const double min = nh.second->GetYaxis()->GetXmin(), max = nh.second->GetYaxis()->GetXmax();
    //const double min = nh.second->GetMinimum(), max = nh.second->GetMaximum()*1.445;
    const double min = 0., max = 9.5;
    nh.second->GetYaxis()->SetRangeUser( min, max );
    auto l_2sigma = lim.DrawLine( -2., min, -2., max );
    lim.DrawLine( 2., min, 2., max );
    lim.SetLineStyle( 1 );
    auto l_3sigma = lim.DrawLine( -3., min, -3., max );
    lim.DrawLine( 3., min, 3., max );
    nh.second->Draw( "same" );
    c.SetLegendX1( 0.65 );
    c.AddLegendEntry( l_2sigma, "#pm 2#sigma", "l" );
    c.AddLegendEntry( l_3sigma, "#pm 3#sigma", "l" );
    c.Prettify( nh.second );
    PaveText::topLabel( "Elastic selection" );
    c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
  }
  for ( const auto& p : pots_names ) {
    Canvas c( Form( "elastic_xi_%u", p.first ), up_label.c_str(), "Preliminary" );
    TGraphAsymmErrors* p_xi = asym_error_bars( m_h_xi[p.first] );
    p_xi->Draw( "ap" );
    p_xi->SetMarkerStyle( 24 );
    p_xi->SetLineWidth( 2 );
    p_xi->SetLineColor( kBlack );
    p_xi->GetHistogram()->SetMaximum( p_xi->GetHistogram()->GetMaximum()*1.1 );
    //p_xi->SetTitle( "" );
    //c.Prettify( m_h_xi[p.first] );
    c.Prettify( p_xi->GetHistogram() );
    auto acc_exp = new TLine( pots_accept[p.first], p_xi->GetHistogram()->GetMinimum(), pots_accept[p.first], p_xi->GetHistogram()->GetMaximum() );
    acc_exp->SetLineStyle( 2 );
    acc_exp->SetLineWidth( 3 );
    acc_exp->SetLineColor( kGreen+2 );
    acc_exp->Draw();
    auto acc_exp_match = new TLine( pots_accept_match[p.first], p_xi->GetHistogram()->GetMinimum(), pots_accept_match[p.first], p_xi->GetHistogram()->GetMaximum() );
    acc_exp_match->SetLineStyle( 1 );
    acc_exp_match->SetLineWidth( 3 );
    acc_exp_match->SetLineColor( kGreen+2 );
    acc_exp_match->Draw();
    auto acc_exp_90pc = new TLine( pots_accept_90pc[p.first], p_xi->GetHistogram()->GetMinimum(), pots_accept_90pc[p.first], p_xi->GetHistogram()->GetMaximum() );
    acc_exp_90pc->SetLineStyle( 1 );
    acc_exp_90pc->SetLineColor( kRed+1 );
    acc_exp_90pc->SetLineWidth( 3 );
    acc_exp_90pc->Draw();
    auto arr_90pc = new TArrow( pots_accept_90pc[p.first], p_xi->GetHistogram()->GetMaximum()*0.9,
                                pots_accept_90pc[p.first]+0.025, p_xi->GetHistogram()->GetMaximum()*0.9, 0.015, ">" );
    arr_90pc->SetLineColor( kRed+1 );
    arr_90pc->SetLineWidth( 3 );
    arr_90pc->Draw();
    c.SetLegendX1( 0.45 );
    c.AddLegendEntry( acc_exp, "Expected acceptance", "l" );
    c.AddLegendEntry( acc_exp_match, "Observed acceptance", "l" );
    c.AddLegendEntry( acc_exp_90pc, "< 10% rad.damage ineff.", "l" );
    PaveText::topLabel( "Elastic selection" );
    c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
  }
}

void plot_matching( double num_sigma, const char* name, TGraphErrors& gr_nomatch, TGraphErrors& gr_rapmatch, TGraphErrors& gr_massmatch, TGraphErrors& gr_massrapmatch, double min, double max, bool right_leg )
{
  Canvas c( name, up_label.c_str(), "Preliminary" );
  /*auto tmp = new TH2D( Form( "tmp_%s", name ), gr_massrapmatch.GetTitle(), 2, min, max, 2, min, max );
  tmp->Draw();*/
  auto diag = new TF1( "diag", "x", min, max );
  diag->SetLineColor( kGray );

  c.SetLegendX1( right_leg ? 0.5 : 0.15 );
  c.SetLegendY1( right_leg ? 0.7 : 0.65 );
  TMultiGraph mg;
  gr_nomatch.SetMarkerStyle( 24 );
  gr_rapmatch.SetMarkerStyle( 22 );
  gr_rapmatch.SetMarkerSize( 1.2 );
  gr_massmatch.SetMarkerStyle( 20 );
  gr_massmatch.SetMarkerSize( 1.0 );
  gr_massrapmatch.SetMarkerStyle( 21 );
  gr_massrapmatch.SetMarkerColor( kBlue+1 );
  gr_massrapmatch.SetMarkerSize( 1.3 );
  gr_massmatch.SetMarkerColor( kRed+1 );
  gr_rapmatch.SetMarkerColor( kGreen+2 );
  gr_massrapmatch.SetLineWidth( 2 );
  gr_massmatch.SetLineWidth( 2 );
  gr_rapmatch.SetLineWidth( 2 );
  gr_nomatch.SetLineWidth( 2 );
  gr_massrapmatch.SetFillStyle( 0 );
  gr_massmatch.SetFillStyle( 0 );
  gr_rapmatch.SetFillStyle( 0 );
  gr_nomatch.SetFillStyle( 0 );
  if ( gr_massrapmatch.GetN() > 0 )
    c.AddLegendEntry( &gr_massrapmatch, "2D matching", "pf" );
  c.AddLegendEntry( &gr_massmatch, "Mass matching", "pf" );
  c.AddLegendEntry( &gr_rapmatch, "Rapidity matching", "pf" );
  c.AddLegendEntry( &gr_nomatch, "No matching", "pf" );
  mg.Add( &gr_nomatch );
  mg.Add( &gr_massmatch );
  mg.Add( &gr_rapmatch );
  mg.Add( &gr_massrapmatch );

  diag->Draw( "l" );
  mg.Draw( "p" );

  //PaveText pt( 0.15, 0.8, 0.4, 0.95 );
  //PaveText pt( right_leg ? 0.5 : 0.15, 0.8, ( right_leg ? 0.5 : 0.15 )+0.25, right_leg ? 0.95 :  );
  PaveText pt( right_leg ? 0.5 : 0.15, right_leg ? 0.8 : 0.735 );
  pt.SetTextAlign( kHAlignLeft+kVAlignBottom );
  pt.AddText( "1-|#Delta#phi_{#gamma#gamma}/#pi| < 0.005" );
  pt.Draw();

  PaveText::topLabel( Form( "%g#sigma matching", num_sigma ) );

  diag->GetHistogram()->SetTitle( gr_massrapmatch.GetTitle() );
  c.Prettify( diag->GetHistogram() );
  //c.Prettify( tmp );
  //mg.GetXaxis()->SetLimits( min, max );
  diag->GetHistogram()->GetYaxis()->SetRangeUser( min, max );
  c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
}


/*void
plot_matching( const char* name, TGraphErrors& gr_unmatch, TGraphErrors& gr_match, TGraphErrors& gr_ooa, double limits )
{
  Canvas c( name, up_label.c_str(), "Preliminary" );

  TF1 lim( "lim", "x", 0., max_xi );
  lim.SetTitle( "#xi(RP)@@#xi(#gamma#gamma)" );
  lim.Draw("lr+");
  lim.SetLineColor(4);
  lim.SetLineWidth(2);

  int colour = kAzure+10;
  TBox box( 0., 0., limits, max_xi );
  box.SetFillColor( colour );
  box.SetLineColor( kBlack );
  box.SetLineWidth( 1 );
  box.Draw( "l" );

  lim.Draw( "lr,same" );

  gr_unmatch.SetLineWidth(2);
  gr_unmatch.SetMarkerStyle(22);
  gr_unmatch.SetMarkerSize(1.25);
  gr_unmatch.Draw("p same");

  gr_match.SetLineColor(2);
  gr_match.SetLineWidth(2);
  gr_match.SetMarkerColor(2);
  gr_match.SetMarkerStyle(20);
  gr_match.SetMarkerSize(1.25);
  gr_match.Draw( "p same" );

  gr_ooa.SetLineWidth(2);
  gr_ooa.SetMarkerStyle(25);
  gr_ooa.SetMarkerSize(1.25);
  gr_ooa.Draw( "p same" );

  c.SetLegendX1( 0.4 );
  c.AddLegendEntry( &gr_match, "Matching events", "lp" );
  c.AddLegendEntry( &gr_unmatch, "Non-matching events", "lp" );
  c.AddLegendEntry( &gr_ooa, "Out of acceptance events", "lp" );
  c.AddLegendEntry( &box, "No acceptance for RP", "f" );

  PaveText lab( 0.8, 0.2, 0.85, 0.25 );
  lab.SetTextSize( 0.075 );
  lab.SetFillStyle( 0 );
  lab.SetLineWidth( 0 );
  lab.AddText( gr_match.GetTitle() );
  lab.Draw( "same" );

  lim.GetYaxis()->SetRangeUser( 0., max_xi );

  TH2D *axiiis = new TH2D( Form( "axiiis_%s", name ), "", 10, 0, max_xi, 10, 0, max_xi );
  axiiis->Draw( "sameaxis" );
  c.GetLegend()->SetFillStyle( 0 );
  c.GetLegend()->SetLineWidth( 0 );
  c.Prettify( lim.GetHistogram() );

  c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
}

void
plot_xispectrum( const char* name, TH1D* spec, double limit )
{
  Canvas c( name, "9.4 fb^{-1} (13 TeV)", "Preliminary" );
  gStyle->SetStatX( 0.88 );
  gStyle->SetStatY( 0.92 );
  gStyle->SetOptStat( "erm" );
  const double max_y = spec->GetMaximum()*1.1;
  spec->Sumw2();
  spec->Draw( "hist" );
  spec->GetYaxis()->SetRangeUser( 0., max_y );
  spec->SetLineColor( kBlack );
  //spec->SetLineWidth( 2 );
  //spec->SetFillColor( kBlack );
  TBox lim( 0., 0., limit, max_y );
  spec->Draw( "p same" );
  lim.SetFillStyle( 3003 );
  lim.SetFillColor( kBlack );
  lim.SetLineColor( kBlack );
  lim.SetLineWidth( 1 );
  lim.Draw( "l" );
  c.Prettify( spec );
  c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
}

void
plot_generic( const char* name, TH1* plot, const char* plot_style, bool logy )
{
  Canvas c( name, "9.4 fb^{-1} (13 TeV)", "Preliminary" );
  plot->Sumw2();
  plot->Draw( plot_style );
  plot->SetMarkerStyle( 20 );
  plot->SetLineColor( kBlack );
  if ( logy ) c.SetLogy();
  c.Prettify( plot );
  c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
}
*/

vector<pair<float,float> > merge_nearfar( const vector<track_t>& near_tracks, const vector<track_t>& far_tracks, float xdiff_cut ) {
  vector<pair<float,float> > out;
  set<unsigned short> matched_far_ids;

  //--- first loop to extract near tracks with matching
  for ( const auto& near : near_tracks ) {
    float min_xidiff = 999.;
    int matched_far_id = -1;
    pair<track_t,track_t> assoc;
    unsigned short id_far = 0;
    for ( const auto& far : far_tracks ) {
      if ( fabs( near.x-far.x ) < xdiff_cut ) {
        //--- near-far track association candidate
        float xidiff = fabs( near.xi-far.xi );
        if ( xidiff < min_xidiff ) {
          assoc = make_pair( near, far );
          min_xidiff = xidiff;
          matched_far_id = id_far;
        }
      }
      id_far++;
    }
    if ( min_xidiff < 999. ) {
      out.emplace_back( assoc.second.xi, assoc.second.err_xi ); // store the far tracks info
      matched_far_ids.insert( matched_far_id );
    }
    //--- store the near track if no matching found with far track
    else out.emplace_back( near.xi, near.err_xi );
  }
  //--- second loop to add the remaining, unmatched far tracks
  unsigned short id_far = 0;
  for ( const auto& far : far_tracks ) {
    if ( matched_far_ids.count( id_far ) == 0 ) {
      //--- discard far track if a mapping is already found
      out.emplace_back( far.xi, far.err_xi );
    }
  }
  return out;
}
