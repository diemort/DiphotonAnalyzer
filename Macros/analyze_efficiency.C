#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TTree.h"
#include "TArrow.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFitResult.h"

#include <map>
#include <iostream>

#include "Canvas.h"
#include "xi_reconstruction.h"
#include "pot_alignment.h"
#include "DiphotonAnalyzer/TreeProducer/interface/MBTreeEvent.h"

#define loc_www "/afs/cern.ch/user/l/lforthom/www/private/ctpps/efficiency"

void analyze_efficiency()
{
  //const unsigned int ref_fill = 4985;
  const unsigned int ref_fill = 4828;
  //const unsigned int ref_fill = 4964;
  const double minimum_threshold = 0.95;

  pot_align::load_file( "TreeProducer/data/alignment_collection_v2.out" );
  xi_reco::load_file( "TreeProducer/data/optics_jun22.root" );

  map<unsigned short,const char*> pot_names = { { 2, "45N" }, { 3, "45F" }, { 102, "56N" }, { 103, "56F" } };
  map<unsigned short,pair<double,double> > pot_fit_limits = { { 2, { 9., 14. } }, { 3, { 12., 14. } }, { 102, { 7., 12. }  }, { 103, { 6., 12. } } };
  map<unsigned short,double> pot_x_mineff = { { 2, 7.6 }, { 3, 7.0 }, { 102, 5.5 }, { 103, 4.8 } }; // x such as eff(x) > 95%
  map<unsigned short,double> erf_fit_min_xi = { { 2, 0.045 }, { 3, 0.042 }, { 102, 0.045 }, { 103, 0.038 } };
  map<unsigned short,TH1D*> h_num_x, h_denom_x, h_num_y, h_denom_y, h_num_y_win, h_denom_y_win, h_num_xi, h_denom_xi, h_num_xi_optbins, h_denom_xi_optbins;
  map<unsigned short,TH2D*> h2_num_xy, h2_denom_xy;
  double y_bins[] = {
   -10., -9., -8., -7., -6.,
   -5., -4.5, -4., -3.5, -3., -2.5, -2.25,
   -2., -1.75, -1.5, -1.375, -1.25, -1.125,
   -1., -0.875, -0.75, -0.625, -0.5, -0.375, -0.25, -0.125,
    0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875,
    1., 1.125, 1.25, 1.375, 1.5, 1.75,
    2., 2.25, 2.5, 3., 3.5, 4., 4.5,
    5., 6., 7., 8., 9., 10. };
  double xi_bins[] = {
    0.020, 0.030,
    0.040, 0.045, 0.0475, 0.050, 0.0525, 0.055, 0.0575,
    0.060, 0.0625, 0.065, 0.070, 0.075,
    0.080, 0.090,
    0.100, 0.105, 0.110,
    0.120, 0.135, 0.150/*, 0.200*/ };
  double x_bins_2d[] = {
    0., 1., 1.5, 2., 2.25, 2.5, 2.75, 3., 3.25, 3.5, 3.75, 4., 4.25, 4.5, 4.75, 5., 5.25, 5.5, 5.75, 6., 6.25, 6.5, 6.75, 7., 7.25, 7.5, 7.75, 8., 8.25, 8.5, 9., 9.5, 10., 10.5, 11., 11.5, 12. };
  double y_bins_2d[] = {
   -8., -6., -5., -4., -3.5, -3., -2.5, -2., -1.5, -1.25, -1., -0.75, -0.5, -0.25, -0.125,
    0., 0.125, 0.25, 0.5, 0.75, 1., 1.25, 1.5, 2., 2.5, 3., 3.5, 4., 5., 6., 8. };
  map<unsigned short,map<unsigned short,TH1D*> > m_h_num_x, m_h_num_xi;
  for ( const auto& p : pot_names ) {
    h_num_x[p.first] = new TH1D( Form( "h_num_x_%d", p.first ), Form( "Track x (%s)@@Entries@@mm?.1f", p.second ), 40, 0., 20. );
//    h_num_x[p.first] = new TH1D( Form( "h_num_x_%d", p.first ), Form( "Track x (%s)@@Entries@@mm?.2f", p.second ), 100, 0., 20. );
    h_denom_x[p.first] = dynamic_cast<TH1D*>( h_num_x[p.first]->Clone( Form( "h_denom_x_%d", p.first ) ) );
    //h_num_y[p.first] = new TH1D( Form( "h_num_y_%d", p.first ), Form( "Track y (%s)@@Entries@@mm?.2f", p.second ), 100, -12.5, 12.5 );
    h_num_y[p.first] = new TH1D( Form( "h_num_y_%d", p.first ), Form( "Track y (%s)@@Entries@@mm?.2f", p.second ), sizeof( y_bins )/sizeof( double )-1, y_bins );
    h_denom_y[p.first] = dynamic_cast<TH1D*>( h_num_y[p.first]->Clone( Form( "h_denom_y_%d", p.first ) ) );
    h_num_y_win[p.first] = dynamic_cast<TH1D*>( h_num_y[p.first]->Clone( Form( "h_num_y_win_%d", p.first ) ) );
    h_num_y_win[p.first]->SetTitle( Form( "Track y (with cut in x) (%s)@@Entries@@mm?.2f", p.second ) );
    h_denom_y_win[p.first] = dynamic_cast<TH1D*>( h_num_y_win[p.first]->Clone( Form( "h_denom_y_win_%d", p.first ) ) );
    //h_num_xi[p.first] = new TH1D( Form( "h_num_xi_%d", p.first ), Form( "Track #xi (%s)@@Entries@@?.3f", p.second ), 26, 0.02, 0.15 );
    h_num_xi[p.first] = new TH1D( Form( "h_num_xi_%d", p.first ), Form( "Track #xi (%s)@@Entries@@?.3f", p.second ), 52, 0.02, 0.15 );
    h_denom_xi[p.first] = dynamic_cast<TH1D*>( h_num_xi[p.first]->Clone( Form( "h_denom_xi_%d", p.first ) ) );
    h_num_xi_optbins[p.first] = new TH1D( Form( "h_num_xi_optbins_%d", p.first ), Form( "Track #xi (%s)@@Entries@@?.3f", p.second ), sizeof( xi_bins )/sizeof( double )-1, xi_bins );
    h_denom_xi_optbins[p.first] = dynamic_cast<TH1D*>( h_num_xi_optbins[p.first]->Clone( Form( "h_denom_xi_optbins_%d", p.first ) ) );
    h2_num_xy[p.first] = new TH2D( Form( "h2_num_xy_%d", p.first ), Form( "Track x (%s) (mm)@@Track y (%s) (mm)", p.second, p.second ), 48, 0., 12., 64, -8., 8. );
//    h2_num_xy[p.first] = new TH2D( Form( "h2_num_xy_%d", p.first ), Form( "Track x (%s) (mm)@@Track y (%s) (mm)", p.second, p.second ), sizeof( x_bins_2d )/sizeof( double )-1, x_bins_2d, sizeof( y_bins_2d )/sizeof( double )-1, y_bins_2d );
    h2_denom_xy[p.first] = dynamic_cast<TH2D*>( h2_num_xy[p.first]->Clone( Form( "h2_denom_xy_%d", p.first ) ) );
  }

  //map<unsigned int,pot_align::align_t> align_el = { { 2, { 1.7773, 0., 0.5484, 0. } }, { 3, { -1.25905, 0., -0.533, 0. } }, { 102, { -0.0385, 0., 0.77165, 0. } }, { 103, { -0.55638, 0., 0.5982, 0. } } };

  MBTreeEvent mb_ev;
  auto mb_tree = dynamic_cast<TTree*>( TFile::Open( "Samples/output_alignmentrun_mbtree.root" )->Get( "ntp" ) );
  mb_ev.attach( mb_tree );
  for ( long long i = 0; i < mb_tree->GetEntriesFast(); ++i ) {
    mb_tree->GetEntry( i );
    for ( unsigned int j = 0; j < mb_ev.num_strips_track; ++j ) {
      const unsigned short pid = 100 * mb_ev.strips_track_arm[j] + mb_ev.strips_track_pot[j];
      /*const double trk_x = mb_ev.strips_track_x[j] * 1.e3 - align_el[pid].x;
      const double trk_y = mb_ev.strips_track_y[j] * 1.e3 + align_el[pid].y;*/
      const double trk_x = mb_ev.strips_track_x[j] * 1.e3;
      const double trk_y = mb_ev.strips_track_y[j] * 1.e3;
      double xi = 0., xi_err = 0.;
      xi_reco::reconstruct( trk_x*1.e-3, mb_ev.strips_track_arm[j], mb_ev.strips_track_pot[j], xi, xi_err );
      h_denom_x[pid]->Fill( trk_x );
      h_denom_y[pid]->Fill( trk_y );
      h_denom_xi[pid]->Fill( xi );
      h_denom_xi_optbins[pid]->Fill( xi );
      if ( trk_x > pot_x_mineff[pid] ) h_denom_y_win[pid]->Fill( trk_y );
      h2_denom_xy[pid]->Fill( trk_x, trk_y );
    }
  }

  auto f = TFile::Open( "/eos/cms/store/user/lforthom/ctpps/efficiency_study/merged_eff_outputBCG.root" );
  auto tree = dynamic_cast<TTree*>( f->Get( "ntp" ) );
  tree->SetBranchStatus( "*", 0 );
  unsigned int fill_number;
  tree->SetBranchStatus( "fill_number", 1 ); tree->SetBranchAddress( "fill_number", &fill_number );
  const unsigned short max_tracks = 20;
  unsigned int num_proton_track, proton_track_side[max_tracks], proton_track_pot[max_tracks];
  float proton_track_x[max_tracks], proton_track_y[max_tracks];
  tree->SetBranchStatus( "num_proton_track", 1 ); tree->SetBranchAddress( "num_proton_track", &num_proton_track );
  tree->SetBranchStatus( "proton_track_x", 1 ); tree->SetBranchAddress( "proton_track_x", proton_track_x );
  tree->SetBranchStatus( "proton_track_y", 1 ); tree->SetBranchAddress( "proton_track_y", proton_track_y );
  tree->SetBranchStatus( "proton_track_side", 1 ); tree->SetBranchAddress( "proton_track_side", proton_track_side );
  tree->SetBranchStatus( "proton_track_pot", 1 ); tree->SetBranchAddress( "proton_track_pot", proton_track_pot );

  const unsigned long long num_entries = tree->GetEntriesFast();
  for ( unsigned long long i = 0; i < num_entries; ++i ) {
    tree->GetEntry( i );
    if ( i % 1000000 == 0 ) cout << "-- event " << i << "/" << num_entries << endl;
    //const bool is_ref_fill = ( fill_number == ref_fill );
    //if ( fill_number < ref_fill ) continue; // FIXME
    if ( fill_number < 4964 ) continue; //FIXME skipping fills with margin

    auto align = pot_align::get_alignments( fill_number );

    for ( unsigned int j = 0; j < num_proton_track; ++j ) {
      const unsigned short pid = 100 * proton_track_side[j] + proton_track_pot[j];
      if ( m_h_num_xi[pid].count( fill_number ) == 0 )
        m_h_num_xi[pid][fill_number] = dynamic_cast<TH1D*>( h_num_xi[pid]->Clone( Form( "h_num_xi_%d_%d", pid, fill_number ) ) );
      if ( m_h_num_x[pid].count( fill_number ) == 0 )
        m_h_num_x[pid][fill_number] = dynamic_cast<TH1D*>( h_num_x[pid]->Clone( Form( "h_num_x_%d_%d", pid, fill_number ) ) );

      double xi = 0., xi_err = 0.;
      xi_reco::reconstruct( proton_track_x[j]+align[pid].x, proton_track_side[j], proton_track_pot[j], xi, xi_err );
      const double trk_x = ( proton_track_x[j]+align[pid].x )*1.e3;
      const double trk_y = ( proton_track_y[j]-align[pid].y )*1.e3;
      h_num_x[pid]->Fill( trk_x );
      h_num_y[pid]->Fill( trk_y );
      h_num_xi[pid]->Fill( xi );
      h_num_xi_optbins[pid]->Fill( xi );
      m_h_num_x[pid][fill_number]->Fill( trk_x );
      m_h_num_xi[pid][fill_number]->Fill( xi );
      if ( trk_x > pot_x_mineff[pid] ) h_num_y_win[pid]->Fill( trk_y );
      h2_num_xy[pid]->Fill( trk_x, trk_y, 1./num_entries );
    }
  }

  // rescaling for a given x range

  for ( const auto& p : pot_names ) {
    // first start by scaling the two histograms to 1
    h_num_x[p.first]->Sumw2();
    h_denom_x[p.first]->Sumw2();
    h_num_y[p.first]->Sumw2();
    h_denom_y[p.first]->Sumw2();
    h_num_xi[p.first]->Sumw2();
    h_denom_xi[p.first]->Sumw2();
    h_num_xi_optbins[p.first]->Sumw2();
    h_denom_xi_optbins[p.first]->Sumw2();
    h_num_y_win[p.first]->Sumw2();
    h_denom_y_win[p.first]->Sumw2();
    h_num_x[p.first]->Scale( 1./h_num_x[p.first]->Integral() );
    h_denom_x[p.first]->Scale( 1./h_denom_x[p.first]->Integral() );
    const auto limits = pot_fit_limits[p.first];
    const double weight_num = h_num_x[p.first]->Integral( h_num_x[p.first]->GetXaxis()->FindBin( limits.first ), h_num_x[p.first]->GetXaxis()->FindBin( limits.second ) );
    const double weight_denom = h_denom_x[p.first]->Integral( h_denom_x[p.first]->GetXaxis()->FindBin( limits.first ), h_denom_x[p.first]->GetXaxis()->FindBin( limits.second ) );
    //const double weight_denom2 = h_denom_x[p.first]->Integral( h_denom_x[p.first]->GetXaxis()->FindBin( ( p.first / 100 == 0 ) ? 9. : 7.5 ), h_denom_x[p.first]->GetXaxis()->FindBin( 13. ) );
    const double norm = weight_denom/weight_num;
    cout << p.second << "::" << weight_num << "|" << weight_denom << "|ratio=" << norm << endl;
    h_num_x[p.first]->Scale( norm );
    h_num_y[p.first]->Scale( norm/h_num_y[p.first]->Integral() );
    h_denom_y[p.first]->Scale( 1./h_denom_y[p.first]->Integral() );
    h_num_y_win[p.first]->Scale( norm/h_num_y_win[p.first]->Integral() );
    h_denom_y_win[p.first]->Scale( 1./h_denom_y_win[p.first]->Integral() );
    h_num_xi[p.first]->Scale( norm/h_num_xi[p.first]->Integral() );
    h_denom_xi_optbins[p.first]->Scale( 1./h_denom_xi[p.first]->Integral() );
    h_num_xi_optbins[p.first]->Scale( norm/h_num_xi_optbins[p.first]->Integral() );
    h_denom_xi[p.first]->Scale( 1./h_denom_xi[p.first]->Integral() );
    h2_num_xy[p.first]->Scale( norm/h2_num_xy[p.first]->Integral( "width" ) );
    h2_denom_xy[p.first]->Scale( 1./h2_denom_xy[p.first]->Integral( "width" ) );
    for ( auto& h_fill : m_h_num_x[p.first] ) {
      h_fill.second->Sumw2();
      h_fill.second->Scale( 1./h_fill.second->Integral() );
      const double weight_num_fill = h_fill.second->Integral( h_fill.second->GetXaxis()->FindBin( limits.first ), h_fill.second->GetXaxis()->FindBin( limits.second ) );
      //const double weight_num_fill = h_fill.second->Integral( h_fill.second->GetXaxis()->FindBin( ( p.first / 100 == 0 ) ? 9. : 7.5 ), h_fill.second->GetXaxis()->FindBin( 13. ) );
      const double norm_fill = weight_denom/weight_num_fill;
      h_fill.second->Scale( norm_fill );
      m_h_num_xi[p.first][h_fill.first]->Scale( norm_fill/m_h_num_xi[p.first][h_fill.first]->Integral() );
    }
  }

  // plotting part

  const string top_title = "CMS-TOTEM Preliminary 2016, #sqrt{s} = 13 TeV, L = 9.4 fb^{-1}";

  auto text = new TText();
  //text->SetTextColor( kGray+3 );
  text->SetTextColor( kRed+1 );
  text->SetTextFont( 42 );
  text->SetTextAlign( 21 );
  text->SetTextSize( 0.035 );

  auto effErf = [](double* x, double* p) { return ( TMath::Erf( ( x[0] - p[0] )/p[1] )+1. )*0.5*p[2]; };

  map<string,pair<map<unsigned short,TH1D*>,map<unsigned short,TH1D*> > > hists = { { "x", { h_num_x, h_denom_x } }, { "y", { h_num_y, h_denom_y } }, { "xi", { h_num_xi, h_denom_xi } }, { "xi_optbins", { h_num_xi_optbins, h_denom_xi_optbins } }, { "y_win", { h_num_y_win, h_denom_y_win } } };
  map<string,pair<map<unsigned short,map<unsigned short,TH1D*> >,pair<map<unsigned short,TH1D*>,map<unsigned short,TH1D*> > > > hists_perfill = { { "x", { m_h_num_x, { h_num_x, h_denom_x } } }, { "xi", { m_h_num_xi, { h_num_xi, h_denom_xi } } } };
  unsigned short i = 0;
  for ( auto& h_map : hists ) {
    string distrib = h_map.first;
    auto hist = h_map.second;
    for ( const auto& p : pot_names ) {
      auto limits = pot_fit_limits[p.first];
      auto range = new TArrow( limits.first, 0.015, limits.second, 0.015, 0.015, "<>" );
      //range->SetLineColor( kGray+3 );
      range->SetLineColor( kRed+1 );
      range->SetLineWidth( 3 );
      { // plot both the numerator and denominator
        Canvas c( Form( "dist_%s_%s", distrib.c_str(), p.second ), top_title.c_str() );
        if ( distrib == "y" || distrib == "y_win" ) c.SetLegendX1( 0.15 );
        else c.SetLegendX1( 0.475 );
        THStack hs;
        hs.Add( hist.first[p.first], "hist" );
        c.AddLegendEntry( hist.first[p.first], "All runs in 2016B/C/G", "f" );
        hist.first[p.first]->SetLineColor( kBlack );
        hist.first[p.first]->SetFillStyle( 3004 );
        hist.first[p.first]->SetFillColor( kBlack );
        hs.Add( hist.second[p.first], "e" );
        c.AddLegendEntry( hist.second[p.first], Form( "Reference fill (%d)", ref_fill ), "lp" );
        hist.second[p.first]->SetMarkerStyle( 24 );
        //hist.second[p.first]->SetMarkerSize( 0.8 );
        hist.second[p.first]->SetLineColor( kBlack );
        //hist.second[p.first]->SetLineWidth( 2 );
        hs.Draw( "nostack" );
        hs.GetHistogram()->SetTitle( hist.first[p.first]->GetTitle() );
        c.Prettify( hs.GetHistogram() );
        if ( distrib == "xi" ) hs.SetMaximum( 0.052 );
        else hs.SetMaximum( 0.1 );
        if ( distrib == "x" ) {
          range->Draw();
          text->DrawText( 0.5*( limits.first+limits.second ), 0.0175, "Norm. range" );
        }
        c.Save( "pdf,png", loc_www );
      }
      { // plot the efficiencies
        Canvas c( Form( "ratio_%s_%s", distrib.c_str(), p.second ), top_title.c_str() );
        gStyle->SetOptStat( 0 );
        auto ratio = dynamic_cast<TH1D*>( hist.first[p.first]->Clone() ), den = dynamic_cast<TH1D*>( hist.second[p.first]->Clone() );
        ratio->SetTitle( TString( ratio->GetTitle() ).ReplaceAll( "Entries", "Efficiency" ) );
        ratio->Divide( den );
        ratio->Draw( "e0" );
        /*auto ratio = new TGraphAsymmErrors( (TH1D*)hist.first[p.first]->Clone(), (TH1D*)hist.second[p.first]->Clone(), "pois" );
        ratio->Draw( "ap" );
        ratio->SetTitle( TString( hist.first[p.first]->GetTitle() ).ReplaceAll( "Entries", "Efficiency" ) );
        c.Prettify( ratio->GetHistogram() );*/
        c.Prettify( ratio );
        ratio->SetMarkerStyle( 24 );
        ratio->SetLineColor( kBlack );
        ratio->GetYaxis()->SetRangeUser( 0.01, 1.49 );
        if ( distrib == "x" || distrib == "xi" || distrib == "xi_optbins" ) {
          const double min_x = ( distrib == "x" )
            ? ( ( p.first / 100 == 0 ) ? 3. : 2. )
            : 0.03;
          cout << p.first << "::" << min_x << endl;
          c.SetLegendX1( 0.45 );
          c.SetLegendY1( 0.18 );
          c.SetLegendSizeY( 0.2 );
          double range_min = 0., range_max = 0.;
          if ( distrib == "x" ) {
            range_min = ( p.first/100 == 0 ) ? 5.2 : 3.2;
            range_max = 12.;
          }
          else {
            range_min = erf_fit_min_xi[p.first];
            range_max = 0.15;
          }
          auto erf = new TF1( "myErf", effErf, range_min, range_max, 3 );
          erf->SetParameter( 0, 0.05 );
          erf->SetParameter( 1, 5. );
          erf->SetParameter( 2, 1. );
          auto fit_res = ratio->Fit( erf, "rsn" );
          erf->Draw( "same" );
          erf->SetLineColor( kRed+1 );
          erf->SetLineWidth( 2 );
          if ( (int)fit_res == 0 ) {
            c.AddLegendEntry( erf, Form( "Fit (%.3f, %.3f, %.2f),", fit_res->Parameter( 0 ), fit_res->Parameter( 1 ), fit_res->Parameter( 2 ) ), "l" );
            c.AddLegendEntry( 0, Form( "#chi^{2}/ndf = %.2f/%d", fit_res->Chi2(), fit_res->Ndf() ), "" );
            cout << distrib << "::" << p.first << "::" << erf->GetX( 0.85 ) << "|" << erf->GetX( 0.9 ) << "|" << erf->GetX( 0.95 ) << endl;
            auto line90 = new TLine( erf->GetX( 0.9 ), 0.01, erf->GetX( 0.9 ), 0.9 );
            line90->Draw( "same" );
            line90->SetLineStyle( 2 );
            line90->SetLineWidth( 3 );
            line90->SetLineColor( kGray+2 );
            c.AddLegendEntry( line90, Form( "90%% threshold (%.3f)", erf->GetX( 0.9 ) ), "l" );
            auto line95 = new TLine( erf->GetX( 0.95 ), 0.01, erf->GetX( 0.95 ), 0.95 );
            line95->Draw( "same" );
            line95->SetLineStyle( 2 );
            line95->SetLineWidth( 3 );
            line95->SetLineColor( kGray+3 );
            c.AddLegendEntry( line95, Form( "95%% threshold (%.3f)", erf->GetX( 0.95 ) ), "l" );
          }
        }
        if ( distrib == "x" ) range->Draw( "same" );
        if ( distrib == "y" ) {
          auto ratio2 = dynamic_cast<TH1D*>( h_num_y_win[p.first]->Clone() );
          ratio2->Divide( h_denom_y_win[p.first] );
          ratio2->Draw( "e0,same" );
          ratio2->SetMarkerStyle( 25 );
          ratio2->SetMarkerColor( kRed+1 );
          ratio2->SetLineColor( kRed+1 );
          c.SetLegendX1( 0.15 );
          c.SetLegendY1( 0.15 );
          c.AddLegendEntry( ratio, "Full distribution", "p" );
          c.AddLegendEntry( ratio2, "With x cut", "p" );
        }
        c.SetGrid( 1, 1 );
        c.Save( "pdf,png", loc_www );
      }
    }
    i++;
  }
  for ( const auto& p : pot_names ) {
    {
      Canvas c( Form( "ratio2d_xy_%s", p.second ), top_title.c_str() );
      auto ratio = dynamic_cast<TH2D*>( h2_num_xy[p.first]->Clone() );
      ratio->Divide( h2_denom_xy[p.first] );
      ratio->Draw( "colz" );
      //gStyle->SetPalette( kRainBow );
      c.Pad()->SetRightMargin( 0.175 );
      c.Prettify( ratio );
      //ratio->Scale( 1./ratio->Integral() );
      ratio->GetZaxis()->SetRangeUser( 0., 1. );
      ratio->GetZaxis()->SetTitle( "All fills / reference fill" );
      ratio->GetZaxis()->SetTitleOffset( 1.4 );
      c.Save( "pdf,png", loc_www );
    }
    for ( auto& dist : hists_perfill ) {
      auto fills_map = dist.second.first[p.first];
      auto num = dist.second.second.first[p.first], denom = dist.second.second.second[p.first];
      {
        Canvas c( Form( "dist_%s_perfill_%s", dist.first.c_str(), p.second ), top_title.c_str() );
        THStack hs;
        c.SetLegendX1( 0.57 );
        c.AddLegendEntry( denom, "Reference fill", "l" );
        unsigned short j = 0;
        for ( auto& h : fills_map ) {
          hs.Add( h.second );
          h.second->SetLineColor( kGray+1 );
          if ( j == 0 ) c.AddLegendEntry( h.second, "Indiv. fills", "l" );
          j++;
        }
        c.AddLegendEntry( num, "#Sigma fills", "l" );
        hs.Add( denom );
        denom->SetLineColor( kRed+1 );
        denom->SetLineWidth( 3 );
        denom->SetLineStyle( 2 );
        hs.Add( num );
        num->SetLineColor( kBlack );
        num->SetLineWidth( 3 );
        num->SetFillStyle( 0 );
        hs.Draw( "hist,nostack" );
        //hs.SetMaximum( hs.GetMaximum()*1.1 );
        hs.GetHistogram()->SetTitle( num->GetTitle() );
        c.Prettify( hs.GetHistogram() );
        hs.SetMaximum( denom->GetMaximum()*1.25 );
        c.Save( "pdf,png", loc_www );
      }
      {
        Canvas c( Form( "ratio_%s_perfill_%s", dist.first.c_str(), p.second ), top_title.c_str() );
        THStack hs;
        for ( auto& h : fills_map ) {
          auto ratio = dynamic_cast<TH1D*>( h.second->Clone() );
          ratio->Divide( denom );
          //ratio->Sumw2( false );
          hs.Add( ratio, "hist" );
          ratio->SetLineColor( kGray+1 );
        }
        auto ratio_tot = dynamic_cast<TH1D*>( num->Clone() );
        ratio_tot->Divide( denom );
        hs.Add( ratio_tot, "e" );
        ratio_tot->SetLineColor( kBlack );
        ratio_tot->SetLineWidth( 3 );
        ratio_tot->SetFillStyle( 0 );
        hs.Draw( "nostack" );
        hs.GetHistogram()->SetTitle( TString( num->GetTitle() ).ReplaceAll( "Entries", "Efficiency" ) );
        hs.SetMinimum( 0.45 );
        hs.SetMaximum( 1.35 );
        c.SetGrid( 1, 1 );
        c.Prettify( hs.GetHistogram() );
        auto one = new TF1( "one", "1", 0.02, 0.15 );
        one->SetLineColor( kRed+1 );
        one->SetLineWidth( 2 );
        one->Draw( "same" );
        c.Save( "pdf,png", loc_www );
      }
      {
        Canvas c( Form( "ratio_ratio_%s_perfill_%s", dist.first.c_str(), p.second ), top_title.c_str() );
        THStack hs;
        auto ratio_tot = dynamic_cast<TH1D*>( num->Clone() );
        ratio_tot->Divide( denom );
        for ( auto& h : fills_map ) {
          auto ratio = dynamic_cast<TH1D*>( h.second->Clone() );
          ratio->Divide( denom );
          ratio->Divide( ratio_tot );
          //ratio->Sumw2( false );
          hs.Add( ratio, "hist" );
          ratio->SetLineColor( kGray+1 );
        }
        hs.Draw( "nostack" );
        hs.GetHistogram()->SetTitle( TString( num->GetTitle() ).ReplaceAll( "Entries", "Efficiency (fill) / Efficiency (total)" ) );
        hs.SetMinimum( 0.35 );
        hs.SetMaximum( 1.45 );
        c.SetGrid( 1, 1 );
        c.Prettify( hs.GetHistogram() );
        auto one = new TF1( "one", "1", 0.02, 0.15 );
        one->SetLineColor( kRed+1 );
        one->SetLineWidth( 2 );
        one->Draw( "same" );
        c.Save( "pdf,png", loc_www );
      }
    }
  }
}

