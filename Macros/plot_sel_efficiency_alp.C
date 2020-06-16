#include "Canvas.h"
#include "mass_templates.h"
#include "DiphotonAnalyzer/TreeProducer/interface/TreeEvent.h"

#include <fstream>
#include <iostream>
#include <unordered_map>

#include "TFile.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH1.h"
#include "THStack.h"
#include "TString.h"

#define OUT_PATH "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp/limits"

void plot_sel_efficiency_alp()
{
  struct sample_t
  {
    double m, k;
    string path;
    double xsec, xsec_bare;
    unsigned int num_events;
  };
  //map<string,float> pots_accept_tight = { { "45N", 0.067 }, { "45F", 0.066 }, { "56N", 0.070 }, { "56F", 0.061 } };
  map<string,float> pots_accept_tight = { { "45N", 0.068 }, { "45F", 0.064 }, { "56N", 0.069 }, { "56F", 0.060 } }; //may20
  TH1D h_binned_mass_data( "data_obs", ";m (GeV);Entries", alp::mass_bins.size()-1, alp::mass_bins.data() );
  const string combine_path = "/afs/cern.ch/work/l/lforthom/private/twophoton/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit";
  const string aux_combine_file = "sum_bckg_mass.root";
  TH1D* h_sig = nullptr, *h_incl = nullptr;
  {
    auto f_in = TFile::Open( aux_combine_file.c_str() );
    h_sig = dynamic_cast<TH1D*>( f_in->Get( "aaaa" )->Clone() );
    h_incl = dynamic_cast<TH1D*>( f_in->Get( "incl" )->Clone() );
//    delete f_in;
  }
  const float the_lumi = ( 4.291753570355+1.438424638452+3.625791461729 ) * 1.e3; //FIXME FIXME FIXME
  //const double probed_coupl = 0.5;
  const double probed_coupl = 0.1;
  //const double probed_coupl = 0.0316;

  //retrieved from /afs/cern.ch/user/j/juwillia/public/forLaurent/AC_microAOD
  vector<sample_t> samples = {
    { 500.,  5.,     "samples/ntuple-alp-m500_f5e0.root",    -1., -1., 5000 },
    { 500.,  1.,     "samples/ntuple-alp-m500_f1e0.root",    -1., -1., 5000 },
    { 500.,  0.6,    "samples/ntuple-alp-m500_f6e-1.root",    8.343e-2, 2.004e-1, 5000 },
    { 500.,  0.5,    "samples/ntuple-alp-m500_f5e-1.root",    5.884e-2, 1.504e-1, 5000 },
    { 500.,  0.25,   "samples/ntuple-alp-m500_f2p5e-1.root",  1.448e-2, 3.521e-2, 5000 },
    { 500.,  0.15,   "samples/ntuple-alp-m500_f1p5e-1.root",  5.120e-3, 1.340e-2, 5000 },
    { 500.,  0.1,    "samples/ntuple-alp-m500_f1e-1.root",    2.310e-3, 5.912e-3, 5000 },
    { 500.,  0.06,   "samples/ntuple-alp-m500_f6e-2.root",    8.147e-4, 1.867e-3, 5000 },
    { 500.,  0.0316, "samples/ntuple-alp-m500_f1e-1p5.root",  2.297e-4, 8.623e-4, 5000 },
    { 750.,  5.,     "samples/ntuple-alp-m750_f5e0.root",    -1., -1., 5000 },
    { 750.,  1.,     "samples/ntuple-alp-m750_f1e0.root",    -1., -1., 5000 },
    { 750.,  0.5,    "samples/ntuple-alp-m750_f5e-1.root",    5.015e-2, 7.633e-2, 5000 },
    { 750.,  0.1,    "samples/ntuple-alp-m750_f1e-1.root",    1.978e-3, 3.116e-3, 5000 },
    { 750.,  0.0316, "samples/ntuple-alp-m750_f1e-1p5.root",  2.001e-4, 3.809e-4, 5000 },
    { 1000., 5.,     "samples/ntuple-alp-m1000_f5e0.root",    -1., -1., 5000 },
    { 1000., 1.,     "samples/ntuple-alp-m1000_f1e0.root",    -1., -1., 5000 },
    { 1000., 0.6,    "samples/ntuple-alp-m1000_f6e-1.root",   3.393e-2, 6.693e-2, 5000 },
    { 1000., 0.5,    "samples/ntuple-alp-m1000_f5e-1.root",   2.382e-2, 4.377e-2, 5000 },
    { 1000., 0.25,   "samples/ntuple-alp-m1000_f2p5e-1.root", 6.159e-3, 1.136e-2, 5000 },
    { 1000., 0.15,   "samples/ntuple-alp-m1000_f1p5e-1.root", 2.200e-3, 3.933e-3, 5000 },
    { 1000., 0.1,    "samples/ntuple-alp-m1000_f1e-1.root",   9.534e-4, 1.745e-3, 5000 },
    { 1000., 0.06,   "samples/ntuple-alp-m1000_f6e-2.root",   3.355e-4, 6.796e-4, 5000 },
    { 1000., 0.0316, "samples/ntuple-alp-m1000_f1e-1p5.root", 9.595e-5, 1.466e-4, 5000 },
    { 1250., 5.,     "samples/ntuple-alp-m1250_f5e0.root",    -1., -1., 5000 },
    { 1250., 1.,     "samples/ntuple-alp-m1250_f1e0.root",    -1., -1., 5000 },
    { 1250., 0.5,    "samples/ntuple-alp-m1250_f5e-1.root",   1.089e-2, 3.062e-2, 5000 },
    { 1250., 0.1,    "samples/ntuple-alp-m1250_f1e-1.root",   4.624e-4, 1.282e-3, 5000 },
    { 1250., 0.0316, "samples/ntuple-alp-m1250_f1e-1p5.root", 4.654e-5, 1.163e-4, 5000 },
    { 1500., 5.,     "samples/ntuple-alp-m1500_f5e0.root",    -1., -1., 5000 },
    { 1500., 1.,     "samples/ntuple-alp-m1500_f1e0.root",    -1., -1., 5000 },
    { 1500., 0.6,    "samples/ntuple-alp-m1500_f6e-1.root",   6.680e-3, 2.722e-2, 5000 },
    { 1500., 0.5,    "samples/ntuple-alp-m1500_f5e-1.root",   4.823e-3, 2.015e-2, 5000 },
    { 1500., 0.25,   "samples/ntuple-alp-m1500_f2p5e-1.root", 1.305e-3, 4.919e-3, 5000 },
    { 1500., 0.15,   "samples/ntuple-alp-m1500_f1p5e-1.root", 4.444e-4, 1.700e-3, 5000 },
    { 1500., 0.1,    "samples/ntuple-alp-m1500_f1e-1.root",   1.968e-4, 7.331e-4, 5000 },
    { 1500., 0.06,   "samples/ntuple-alp-m1500_f6e-2.root",   6.940e-5, 2.729e-4, 5000 },
    { 1500., 0.0316, "samples/ntuple-alp-m1500_f1e-1p5.root", 1.968e-5, 7.123e-5, 5000 },
    { 1750., 5.,     "samples/ntuple-alp-m1750_f5e0.root",    -1., -1., 5000 },
    { 1750., 1.,     "samples/ntuple-alp-m1750_f1e0.root",    -1., -1., 5000 },
    { 1750., 0.5,    "samples/ntuple-alp-m1750_f5e-1.root",   1.504e-3, 1.214e-2, 5000 },
    { 1750., 0.1,    "samples/ntuple-alp-m1750_f1e-1.root",   5.819e-5, 5.556e-4, 5000 },
    { 1750., 0.0316, "samples/ntuple-alp-m1750_f1e-1p5.root", 5.799e-6, 5.155e-5, 5000 },
    { 2000., 5.,     "samples/ntuple-alp-m2000_f5e0.root",    -1., -1., 5000 },
    { 2000., 1.,     "samples/ntuple-alp-m2000_f1e0.root",    -1., -1., 5000 },
    { 2000., 0.5,    "samples/ntuple-alp-m2000_f5e-1.root",   2.162e-4, 8.872e-3, 5000 },
    { 2000., 0.1,    "samples/ntuple-alp-m2000_f1e-1.root",   3.777e-7, 3.344e-4, 5000 },
    { 2000., 0.0316, "samples/ntuple-alp-m2000_f1e-1p5.root", 3.552e-9, 3.346e-5, 5000 },
    { 2250., 0.1,    "samples/ntuple-alp-m2250_f1e-1.root",   1.082e-7, 2.440e-4, 5000 },
  };
  const size_t num_samples = samples.size();
  vector<TH1D> v_h_binned_mass;
  vector<pair<double,double> > v_params;
  vector<TH1D> v_h_mass( num_samples ), v_h_ptpair( num_samples ), v_h_ptlead( num_samples ), v_h_acop( num_samples );
  TGraph g_el_vals;
  map<double,TGraphErrors> g_el_eff_m, g_el_eff_k;
  TGraph2D g_el_eff, g_el_acc, g_xs;
  gggg::TreeEvent ev;
  const float sqrt_s = 13.e3;
  const TString ylabel = "Events";
  //const TString ylabel = "Events fraction";
  unsigned short i = 0;
  for ( const auto& s : samples ) {
    unique_ptr<TFile> file( TFile::Open( s.path.c_str() ) );
    auto tree = dynamic_cast<TTree*>( file->Get( "ntp" ) );
    v_h_mass[i] = TH1D( Form( "mass_%d", i ), Form( "m_{#gamma#gamma}@@%s@@GeV", ylabel.Data() ), 37, 250., 2100. );
    if ( s.k == probed_coupl ) {
      v_h_binned_mass.emplace_back( Form( "alp_%zu", v_h_binned_mass.size() ), ";m (GeV);Entries", alp::mass_bins.size()-1, alp::mass_bins.data() );
      v_params.emplace_back( s.m, s.k );
    }
    v_h_ptpair[i] = TH1D( Form( "ptpair_%d", i ), Form( "p_{T}^{#gamma#gamma}@@%s@@GeV", ylabel.Data() ), 25, 0., 100. );
    v_h_ptlead[i] = TH1D( Form( "leadpt_%d", i ), Form( "p_{T}^{#gamma} (leading gamma)@@%s@@GeV", ylabel.Data() ), 50, 50., 1050. );
    v_h_acop[i] = TH1D( Form( "acop_%d", i ), Form( "1-|#Delta#phi_{#gamma#gamma}/#pi| (#times 10^{-3})@@%s@@?.g", ylabel.Data() ), 50, 0., 5. );
    ev.attach( tree );
    //cout << s.xsec << endl;
    unsigned short num_passing = 0;
    const double weight = s.xsec/s.num_events*the_lumi;
    //const double weight = 1./s.num_events;
    for ( long long j = 0; j < tree->GetEntriesFast(); ++j ) {
      tree->GetEntry( j );
      if ( ev.hlt_accept[0] == 0 ) continue;
      bool has_diphoton = false;
      for ( unsigned short k = 0; k < ev.num_diphoton; ++k ) {
        if ( ev.diphoton_mass[k] < 350. ) continue;
        if ( ev.diphoton_pt1[k] < 75. ) continue;
        if ( ev.diphoton_pt2[k] < 75. ) continue;
        if ( fabs( ev.diphoton_eta1[k] ) > 1.4442 && fabs( ev.diphoton_eta1[k] < 1.566 ) ) continue;
        if ( fabs( ev.diphoton_eta2[k] ) > 1.4442 && fabs( ev.diphoton_eta2[k] < 1.566 ) ) continue;
        if ( fabs( ev.diphoton_eta1[k] ) > 2.5 || fabs( ev.diphoton_eta2[k] ) > 2.5 ) continue;
        if ( fabs( ev.diphoton_eta1[k] ) > 1.566 && fabs( ev.diphoton_eta2[k] ) > 1.566 ) continue;
        if ( ev.diphoton_r91[k] < 0.94 ) continue;
        if ( ev.diphoton_r92[k] < 0.94 ) continue;
        if ( ev.diphoton_ele_veto1[k] == 0 ) continue;
        if ( ev.diphoton_ele_veto2[k] == 0 ) continue;
        const float acop = 1.-fabs( ev.diphoton_dphi[k]/M_PI );
        const float xip = ( ev.diphoton_pt1[k]*exp( +ev.diphoton_eta1[k] ) + ev.diphoton_pt2[k]*exp( +ev.diphoton_eta2[k] ) ) / sqrt_s,
                    xim = ( ev.diphoton_pt1[k]*exp( -ev.diphoton_eta1[k] ) + ev.diphoton_pt2[k]*exp( -ev.diphoton_eta2[k] ) ) / sqrt_s;
        v_h_acop[i].Fill( acop*1.e3, weight );
        if ( acop > 0.005 ) continue;
        v_h_mass[i].Fill( ev.diphoton_mass[k], weight );
        if ( s.k == probed_coupl )
          v_h_binned_mass.rbegin()->Fill( ev.diphoton_mass[k], weight );
        v_h_ptpair[i].Fill( ev.diphoton_pt[k], weight );
        v_h_ptlead[i].Fill( ev.diphoton_pt1[k], weight );
        if ( ( ( xim > pots_accept_tight["56N"] || xim > pots_accept_tight["56F"] ) && xim < 0.15 )
          && ( ( xip > pots_accept_tight["45N"] || xip > pots_accept_tight["45F"] ) && xip < 0.15 ) )
          has_diphoton = true;
      }
      if ( has_diphoton ) num_passing++;
    }
    const double el_eff = num_passing*1./s.num_events, el_eff_unc = el_eff*sqrt( 1./num_passing+1./s.num_events );
    const double el_acc = s.xsec / s.xsec_bare;
    {
      unsigned short k = g_el_eff_k[s.m].GetN();
      g_el_eff_k[s.m].SetPoint( k, s.k, el_eff );
      g_el_eff_k[s.m].SetPointError( k, 0., el_eff_unc );
    }
    {
      unsigned short k = g_el_eff_m[s.k].GetN();
      g_el_eff_m[s.k].SetPoint( k, s.m, el_eff );
      g_el_eff_m[s.k].SetPointError( k, 0., el_eff_unc );
    }
    const double log_k = log10( s.k );
    g_el_eff.SetPoint( i, s.m, log_k, el_eff );
    g_el_acc.SetPoint( i, s.m, log_k, el_acc );
    //g_xs.SetPoint( i, s.m, log_k, s.xsec_bare );
    g_xs.SetPoint( i, s.m, log_k, s.xsec*el_eff );
    //g_xs.SetPoint( i, s.m, log_k, s.xsec );
    g_el_vals.SetPoint( i, s.m, log_k );
    cout << s.path << ": " << el_eff << "(" << el_eff_unc/el_eff << ")|" << el_acc << "|" << s.xsec << endl;
    ++i;
  }
  g_el_vals.SetPoint( g_el_vals.GetN(), 0., 0. );
  //gStyle->SetPalette( kSunset );
  gStyle->SetPalette( kBeach );
  //--- prepare the data cards for limits setting
  {
    Canvas c( "alp_mass_templates", "FPMC #gamma#gamma #rightarrow #gamma#gamma + ALP", "Simulation preliminary", false, Canvas::Align::right );
    THStack hs;
    unsigned short idx = 0;
    for ( auto& h_bm : v_h_binned_mass ) {
      ifstream in_tmpl( combine_path+"/counting_experiment_protontag_tmpl.txt" );
      ofstream out_tmpl( combine_path+string( Form( "/counting_experiment_alp%d.txt", idx ) ) );
      string line;
      while ( getline( in_tmpl, line ) ) {
        TString ln( line.c_str() );
        out_tmpl
          << ln.ReplaceAll( "XIX", Form( "%d", idx ) )
               //.ReplaceAll( "XRX", Form( "%g", el_eff*s.xsec*the_lumi ) )
               //.ReplaceAll( "XRX", Form( "%g", s.xsec*the_lumi ) )
               .ReplaceAll( "XRX", Form( "%g", h_bm.Integral( 1, alp::mass_bins.size() ) ) )
               .ReplaceAll( "XLBLX", Form( "%g", h_sig->Integral( 1, alp::mass_bins.size() ) ) )
               .ReplaceAll( "XINCX", Form( "%g", h_incl->Integral( 1, alp::mass_bins.size() ) ) )
               .ReplaceAll( "XMX", Form( "%g", v_params[idx].first ) )
               .ReplaceAll( "XFM1X", Form( "%g", v_params[idx].second ) )
          << endl;
      }
      auto hist = (TH1D*)h_bm.Clone();
      hist->SetLineColor( Canvas::colour_pool[idx] );
      hist->SetLineWidth( 3 ) ;
      hs.Add( hist );
      ++idx;
    }
    hs.Draw( "hist,nostack" );
    hs.GetHistogram()->SetTitle( ";m_{#gamma#gamma};Events/bin" );
    c.Prettify( hs.GetHistogram() );
    PaveText::topLabel( Form( "f^{-1} = %g GeV^{-1}", probed_coupl ) );
    c.Save( "pdf,png", OUT_PATH );
  }
  //const char* top_label = "Elastic selection";
  const char* top_label = "#xi^{#pm}_{#gamma#gamma} in acc., #epsilon(#xi^{#pm}_{#gamma#gamma}) > 90%";
  {
    Canvas c( "scan_xsec_alp", "FPMC #gamma#gamma #rightarrow #gamma#gamma + ALP", "Simulation preliminary" );
    c.SetRightMargin( 0.15 );
    g_xs.SetTitle( ";m_{a} (GeV);log(f^{-1}/GeV^{-1});#sigma_{obs} (pb)" );
    //g_xs.SetTitle( ";m_{a} (GeV);log(f^{-1}/GeV^{-1});#sigma_{gen} (pb)" );
    g_el_vals.SetMarkerStyle( 20 );
    g_xs.Draw( "colz0" );
    //g_xs.Draw( "cont4z" );
    g_el_vals.Draw( "p same" );
    g_el_vals.SetMarkerStyle( 3 );
    g_xs.GetZaxis()->SetTitleOffset( 1.4 );
    c.SetLogz();
    c.Prettify( (TH1*)g_el_vals.GetHistogram() );
    c.Save( "pdf,png", OUT_PATH );
  }
  {
    Canvas c( "elastic_acc_scan_alp", "FPMC #gamma#gamma #rightarrow #gamma#gamma + ALP", "Simulation preliminary" );
    c.SetRightMargin( 0.15 );
    g_el_acc.SetTitle( ";m_{a} (GeV);log(f^{-1}/GeV^{-1});Signal acceptance" );
    g_el_vals.SetMarkerStyle( 20 );
    g_el_acc.Draw( "colz0" );
    g_el_vals.Draw( "p same" );
    g_el_vals.SetMarkerStyle( 3 );
    g_el_acc.GetZaxis()->SetTitleOffset( 1.4 );
    c.Prettify( (TH1*)g_el_vals.GetHistogram() );
    c.Save( "pdf,png", OUT_PATH );
  }
  {
    Canvas c( "elastic_eff_scan_alp", "FPMC #gamma#gamma #rightarrow #gamma#gamma + ALP", "Simulation preliminary" );
    c.SetRightMargin( 0.15 );
    g_el_eff.SetTitle( ";m_{a} (GeV);log(f^{-1}/GeV^{-1});Sel. efficiency" );
    g_el_vals.SetMarkerStyle( 20 );
    //g_el_vals.Draw( "a*" );
    //g_el_eff.Draw( "arr colz" );
    g_el_eff.Draw( "colz0" );
    //g_el_eff.SetMinimum( 0.5 ); g_el_eff.SetMaximum( 0.6 );
    //g_el_eff.SetMinimum( 0.35 ); g_el_eff.SetMaximum( 0.45 );
    //g_el_eff.SetMinimum( 0. ); g_el_eff.SetMaximum( 0.65 );
    g_el_eff.SetMinimum( 0. ); g_el_eff.SetMaximum( 0.75 );
    //c.SetLogx(); c.SetLogy();
    g_el_vals.Draw( "p same" );
    g_el_vals.SetMarkerStyle( 3 );
    g_el_eff.GetZaxis()->SetTitleOffset( 1.4 );
    PaveText::topLabel( top_label );
    c.Prettify( (TH1*)g_el_vals.GetHistogram() );
    c.Save( "pdf,png", OUT_PATH );
  }
  {
    Canvas c( "elastic_eff_scan_alp_1d_coupl", "FPMC #gamma#gamma #rightarrow #gamma#gamma + ALP", "Simulation preliminary" );
    c.SetLegendY1( 0.7 );
    c.SetLegendSizeY( 0.225 );
    TMultiGraph mg;
    unsigned short i = 0;
    TF1 f_fit( "f_fit", "pol2" );
    vector<TF1*> fits_res;
    for ( auto& g : g_el_eff_k ) {
      mg.Add( &g.second );
      g.second.SetLineColor( kBlack );
      g.second.SetMarkerStyle( Canvas::marker_pool[i] );
      g.second.SetMarkerColor( Canvas::colour_pool[i] );
      c.AddLegendEntry( &g.second, Form( "m = %g GeV", g.first ), "ep" );
      ++i;
    }
    mg.Draw( "ap" );
    for ( const auto& f : fits_res )
      f->DrawF1( 500., 2000., "same" );
    mg.SetMinimum( 0. );
    //mg.SetMaximum( 1.05 );
    mg.SetMaximum( 0.85 );
    mg.GetHistogram()->SetTitle( ";f^{-1} (GeV^{-1});Selection efficiency" );
    PaveText::topLabel( top_label );
    c.Prettify( mg.GetHistogram() );
    c.Save( "pdf,png", OUT_PATH );
  }
  {
    Canvas c( "elastic_eff_scan_alp_1d_mass", "FPMC #gamma#gamma #rightarrow #gamma#gamma + ALP", "Simulation preliminary" );
    c.SetLegendY1( 0.55 );
    c.SetLegendX1( 0.15 );
    c.SetLegendSizeY( 0.275 );
    TMultiGraph mg;
    unsigned short i = 0;
    TF1 f_fit( "f_fit", "pol2" );
    vector<TF1*> fits_res;
    for ( auto& g : g_el_eff_m ) {
      mg.Add( &g.second );
      g.second.SetLineColor( kBlack );
      g.second.SetMarkerStyle( Canvas::marker_pool[i] );
      g.second.SetMarkerColor( Canvas::colour_pool[i] );
/*      TFitResultPtr fit_res = g.second.Fit( &f_fit, "s0+" );
      if ( (int)fit_res ) {
        f_fit.SetLineStyle( i+1 );
        fits_res.emplace_back( (TF1*)f_fit.Clone() );
      }*/
      c.AddLegendEntry( &g.second, Form( "f^{-1} = %g GeV^{-1}", g.first ), "ep" );
      ++i;
    }
    mg.Draw( "ap" );
    for ( const auto& f : fits_res )
      f->DrawF1( 500., 2000., "same" );
    mg.SetMinimum( 0. );
    //mg.SetMaximum( 1.05 );
    mg.SetMaximum( 0.85 );
    mg.GetHistogram()->SetTitle( ";m_{a} (GeV);Selection efficiency" );
    PaveText::topLabel( top_label );
    c.Prettify( mg.GetHistogram() );
    c.Save( "pdf,png", OUT_PATH );
  }

  const unordered_map<string,vector<TH1D>&> plots = {
    { "anom_alp_elastic_ptpair", v_h_ptpair },
    { "anom_alp_elastic_acop", v_h_acop },
    { "anom_alp_elastic_ptlead", v_h_ptlead },
    { "anom_alp_elastic_mpair", v_h_mass }
  };
  i = 0;
  for ( auto& p : plots ) {
    Canvas c( p.first.c_str(), Form( "FPMC #gamma#gamma #rightarrow #gamma#gamma (ALP), yields for %.1f fb^{-1} (%g TeV)", the_lumi/1.e3, sqrt_s/1.e3 ), "Simulation preliminary" );
    //Canvas c( p.first.c_str(), Form( "FPMC #gamma#gamma #rightarrow #gamma#gamma (%g TeV)", sqrt_s/1.e3 ), "Simulation preliminary" );
    THStack hs;
    //if ( i > 1 )
      //c.SetLegendY1( 0.695 );
    //else {
      //c.SetLegendY1( 0.61 );
      c.SetLegendY1( 0.62 );
      c.SetLegendX1( 0.17 );
    //  c.SetLegendX1( 0.18 );
    //}
    c.SetLegendSizeX( 0.7 );
    c.SetLegendSizeY( 0.225 );
    unsigned short k = 0;
    for ( unsigned short j = 0; j < num_samples; ++j ) {
      if ( samples[j].k != 0.1 && samples[j].k != 0.5 && samples[j].k != 0.0316 )
        continue;
      //TString label = Form( "#zeta_{1} = %g, #zeta_{2} = %g", samples[j].zeta1, samples[j].zeta2 );
      TString label = Form( "(%g, %g)", samples[j].m, samples[j].k );
      /*p.second[j].SetLineColor( Canvas::colour_pool[j] );
      p.second[j].SetFillColorAlpha( Canvas::colour_pool[j], 0.25 );*/
      p.second[j].SetLineColor( Canvas::colour_pool[k/3] );
      p.second[j].SetFillColorAlpha( Canvas::colour_pool[k/3], 0.15+(k%3)*0.2 );
      p.second[j].SetLineStyle( 1+(k%3) );
      p.second[j].SetLineWidth( 2 );
      hs.Add( &p.second[j] );
      c.AddLegendEntry( &p.second[j], label, "f" );
      ++k;
    }
    c.GetLegend()->SetNColumns( 3 );
    hs.Draw( "hist,nostack" );
    hs.SetMaximum( hs.GetHistogram()->GetMaximum()*25. );
    hs.GetHistogram()->SetTitle( p.second[0].GetTitle() );
    //if ( i == 0 ) { // mass plot
      hs.SetMinimum( 1.e-2 );
      c.SetLogy();
    //}
    c.Prettify( hs.GetHistogram() );
    c.Save( "pdf,png", OUT_PATH );
    ++i;
  }
  {
    {
      TFile f( "input_malp.root", "recreate" );
      h_binned_mass_data.Write();
      for ( const auto& h : v_h_binned_mass )
        h.Write();
      if ( h_sig )  h_sig->Write( "aaaa" );
      if ( h_incl ) h_incl->Write( "incl" );
    }
  }
}
