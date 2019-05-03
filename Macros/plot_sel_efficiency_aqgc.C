#include "Canvas.h"
#include "Plotter.h"
#include "DiphotonAnalyzer/TreeProducer/interface/TreeEvent.h"

#include <iostream>
#include <unordered_map>

#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TH1.h"
#include "THStack.h"

void plot_sel_efficiency_aqgc()
{
  struct sample_t
  {
    double zeta1, zeta2;
    string path;
    double xsec;
    unsigned int num_events;
  };
  map<string,float> pots_accept_tight = { { "45N", 0.067 }, { "45F", 0.066 }, { "56N", 0.070 }, { "56F", 0.061 } };

  //retrieved from /afs/cern.ch/user/j/juwillia/public/forLaurent/AC_microAOD
  vector<sample_t> samples = {
    //{ 0., 0., "/eos/cms/store/user/lforthom/twophoton/samples_80x/output_GammaGammaToGammaGamma_fpmc_justin_sm.root", 1.3988e-05, 100000 },
    //{ 0., 0., "samples/ntuple-sm_gggg.root", 1.3988e-05, 5000 },
    { 0., 1.e-13, "samples/ntuple-aqgc-zeta1_0_zeta2_1e-13.root", 1.3068e-4, 5000 },
    { 1.e-13, 0., "samples/ntuple-aqgc-zeta1_1e-13_zeta2_0.root", 5.7066e-4, 5000 },
    { 1.e-13, 1.e-13, "samples/ntuple-aqgc-zeta1_1e-13_zeta2_1e-13.root", 1.1770e-3, 5000 },
    { 1.e-13, 1.e-12, "samples/ntuple-aqgc-zeta1_1e-13_zeta2_1e-12.root", 1.8399e-2, 5000 },
    { 5.e-13, 0., "samples/ntuple-aqgc-zeta1_5e-13_zeta2_0.root", 1.4263e-2, 4900 },
    { 5.e-13, 5.e-13, "samples/ntuple-aqgc-zeta1_5e-13_zeta2_5e-13.root", 2.9418e-2, 5000 },
    { 5.e-13, 1.e-13, "samples/ntuple-aqgc-zeta1_5e-13_zeta2_1e-13.root", 1.6771e-2, 5000 },
    { 1.e-12, 0., "samples/ntuple-aqgc-zeta1_1e-12_zeta2_0.root", 5.7055e-2, 5000 },
    { 1.e-12, 1.e-12, "samples/ntuple-aqgc-zeta1_1e-12_zeta2_1e-12.root", 0.1177, 5000 },
    { 5.e-12, 0., "samples/ntuples-aqgc-zeta1_5e-12_zeta2_0.root", 5.7055e-2, 5000 },
    // negative
    { 1.e-12, -1.e-12, "samples/ntuple-aqgc-zeta1_1e-12_zeta2_-1e-12.root", 2.2584e-2, 4200 },
    { -1.e-12, -1.e-12, "samples/ntuple-aqgc-zeta1_-1e-12_zeta2_-1e-12.root", 0.1177, 5000 },
    { -1.e-12, 1.e-12, "samples/ntuple-aqgc-zeta1_-1e-12_zeta2_1e-12.root", 2.2585e-2, 5000 },
    { -1.e-12, 0., "samples/ntuple-aqgc-zeta1_-1e-12_zeta2_0.root", 5.7644e-2, 5000 },
    { -5.e-12, -1.e-12, "samples/ntuple-aqgc-zeta1_-5e-12_zeta2_-1e-12.root", 0., 5000 },
    { -5.e-12, 0., "samples/ntuple-aqgc-zeta1_-5e-12_zeta2_0.root", 1.4410, 5000 },
    { -5.e-12, 1.e-12, "samples/ntuple-aqgc-zeta1_-5e-12_zeta2_1e-12.root", 1.2141, 5000 },
    //{ 0., -5.e-12, "samples/ntuple-aqgc-zeta1_0_zeta2_-5e-12.root", 0.3302, 5000 },
    { 1.e-12, 5.e-13, "samples/ntuple-aqgc-zeta1_1e-12_zeta2_5e-13.root", 8.4941e-2, 5000 },
    { 5.e-12, -1.e-12, "samples/ntuple-aqgc-zeta1_5e-12_zeta2_-1e-12.root", 1.2140, 5000 },
    { 5.e-12, 1.e-12, "samples/ntuple-aqgc-zeta1_5e-12_zeta2_1e-12.root", 1.6944, 5000 },
  };
  const size_t num_samples = samples.size();
  vector<TH1D> v_h_mass( num_samples ), v_h_ptpair( num_samples ), v_h_ptlead( num_samples ), v_h_acop( num_samples );
  TGraphErrors g_el_eff_z1, g_el_eff_z2;
  TGraph g_el_eff_vals;
  TGraph2D g_el_eff;
  gggg::TreeEvent ev;
  const float sqrt_s = 13.e3, lumi = 9.36e3;
  const TString ylabel = "Events";
  //const TString ylabel = "Events fraction";
  unsigned short i = 0;
  for ( const auto& s : samples ) {
    unique_ptr<TFile> file( TFile::Open( s.path.c_str() ) );
    auto tree = dynamic_cast<TTree*>( file->Get( "ntp" ) );
    v_h_mass[i] = TH1D( Form( "mass_%d", i ), Form( "m_{#gamma#gamma}@@%s@@GeV", ylabel.Data() ), 18, 200., 2000. );
    v_h_ptpair[i] = TH1D( Form( "ptpair_%d", i ), Form( "p_{T}^{#gamma#gamma}@@%s@@GeV", ylabel.Data() ), 25, 0., 100. );
    v_h_ptlead[i] = TH1D( Form( "leadpt_%d", i ), Form( "p_{T}^{#gamma} (leading gamma)@@%s@@GeV", ylabel.Data() ), 50, 50., 1050. );
    v_h_acop[i] = TH1D( Form( "acop_%d", i ), Form( "1-|#Delta#phi_{#gamma#gamma}/#pi|@@%s@@?.g", ylabel.Data() ), 50, 0., 0.005 );
    ev.attach( tree );
    //cout << s.xsec << endl;
    unsigned short num_passing = 0;
    const double weight = s.xsec/s.num_events*lumi;
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
        v_h_acop[i].Fill( acop, weight );
        if ( acop > 0.005 ) continue;
        v_h_mass[i].Fill( ev.diphoton_mass[k], weight );
        v_h_ptpair[i].Fill( ev.diphoton_pt[k], weight );
        v_h_ptlead[i].Fill( ev.diphoton_pt1[k], weight );
        if ( ( ( xim > pots_accept_tight["56N"] || xim > pots_accept_tight["56F"] ) && xim < 0.15 )
          && ( ( xip > pots_accept_tight["45N"] || xip > pots_accept_tight["45F"] ) && xip < 0.15 ) )
          has_diphoton = true;
      }
      if ( has_diphoton ) num_passing++;
    }
    const double el_eff = num_passing*1./s.num_events, el_eff_unc = el_eff*sqrt( 1./num_passing+1./s.num_events );
    if ( s.zeta1 == 0. ) {
      unsigned short k = g_el_eff_z2.GetN();
      g_el_eff_z2.SetPoint( k, s.zeta2, el_eff );
      g_el_eff_z2.SetPointError( k, 0., el_eff_unc );
    }
    if ( s.zeta2 == 0. ) {
      unsigned short k = g_el_eff_z1.GetN();
      g_el_eff_z1.SetPoint( k, s.zeta1, el_eff );
      g_el_eff_z1.SetPointError( k, 0., el_eff_unc );
    }
    g_el_eff.SetPoint( i, s.zeta1*1.e12, s.zeta2*1.e12, el_eff );
    g_el_eff_vals.SetPoint( i, s.zeta1*1.e12, s.zeta2*1.e12 );
    cout << s.path << ": " << el_eff << endl;
    ++i;
  }
  g_el_eff_vals.SetPoint( g_el_eff_vals.GetN(), 0., 0. );
  //gStyle->SetPalette( kSunset );
  gStyle->SetPalette( kBeach );
  //const char* top_label = "Elastic selection";
  const char* top_label = "#xi^{#pm}_{#gamma#gamma} in acc., #epsilon(#xi^{#pm}_{#gamma#gamma}) > 90%";
  {
    Canvas c( "elastic_eff_scan_aqgc", "FPMC #gamma#gamma #rightarrow #gamma#gamma, AQGC", "Simulation preliminary" );
    c.SetRightMargin( 0.15 );
    g_el_eff.SetTitle( ";#zeta_{1} (#times 10^{-12});#zeta_{2} (#times 10^{-12});#varepsilon" );
    g_el_eff_vals.SetMarkerStyle( 20 );
    //g_el_eff_vals.Draw( "a*" );
    //g_el_eff.Draw( "arr colz" );
    g_el_eff.Draw( "colz" );
    //g_el_eff.SetMinimum( 0.5 ); g_el_eff.SetMaximum( 0.6 );
    //g_el_eff.SetMinimum( 0.35 ); g_el_eff.SetMaximum( 0.45 );
    g_el_eff.SetMinimum( 0.4 ); g_el_eff.SetMaximum( 0.455 );
    //c.SetLogx(); c.SetLogy();
    g_el_eff_vals.Draw( "p same" );
    g_el_eff_vals.SetMarkerStyle( 3 );
    PaveText::topLabel( top_label );
    c.Prettify( (TH1*)g_el_eff_vals.GetHistogram() );
    c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
  }
  {
    Canvas c( "elastic_eff_scan_aqgc_1d", "FPMC #gamma#gamma #rightarrow #gamma#gamma, AQGC", "Simulation preliminary" );
    g_el_eff_z1.SetTitle( ";#zeta_{1,2};Selection efficiency" );
    g_el_eff_z1.Draw( "ap" );
    g_el_eff_z1.SetMarkerStyle( 24 );
    g_el_eff_z1.SetMinimum( 0. );
    g_el_eff_z1.SetMaximum( 1.05 );
    g_el_eff_z2.Draw( "p,same" );
    g_el_eff_z2.SetMarkerStyle( 20 );
    c.AddLegendEntry( &g_el_eff_z1, "#zeta_{1} (#zeta_{2} = 0)", "ep" );
    c.AddLegendEntry( &g_el_eff_z2, "#zeta_{2} (#zeta_{1} = 0)", "ep" );
    PaveText::topLabel( top_label );
    c.Prettify( g_el_eff_z1.GetHistogram() );
    c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
  }

  const unordered_map<string,vector<TH1D>&> plots = {
    { "anom_aqgc_elastic_ptpair", v_h_ptpair },
    { "anom_aqgc_elastic_acop", v_h_acop },
    { "anom_aqgc_elastic_ptlead", v_h_ptlead },
    { "anom_aqgc_elastic_mpair", v_h_mass }
  };
  i = 0;
  for ( auto& p : plots ) {
    Canvas c( p.first.c_str(), Form( "FPMC #gamma#gamma #rightarrow #gamma#gamma, yields @ %.1f fb^{-1} (%g TeV)", lumi/1.e3, sqrt_s/1.e3 ), "Simulation preliminary" );
    //Canvas c( p.first.c_str(), Form( "FPMC #gamma#gamma #rightarrow #gamma#gamma (%g TeV)", sqrt_s/1.e3 ), "Simulation preliminary" );
    THStack hs;
    /*if ( i++ > 2 )
      c.SetLegendY1( 0.695 );
    else {
      c.SetLegendY1( 0.61 );
      c.SetLegendX1( 0.18 );
    }
    c.SetLegendSizeY( 0.2 );*/
    c.SetLegendY1( 0.68 );
    c.SetLegendX1( 0.18 );
    c.SetLegendSizeX( 0.7 );
    c.SetLegendSizeY( 0.15 );
    for ( unsigned short j = 0; j < num_samples; ++j ) {
      //TString label = Form( "#zeta_{1} = %g, #zeta_{2} = %g", samples[j].zeta1, samples[j].zeta2 );
      TString label = Form( "(%g, %g)", samples[j].zeta1, samples[j].zeta2 );
      if ( samples[j].zeta1 == 0. && samples[j].zeta2 == 0. )
        label += " (SM)";
      p.second[j].SetLineColor( Plotter::colour_pool[j] );
      p.second[j].SetFillColorAlpha( Plotter::colour_pool[j], 0.25 );
      p.second[j].SetLineWidth( 2 );
      hs.Add( &p.second[j] );
      c.AddLegendEntry( &p.second[j], label, "f" );
    }
    c.GetLegend()->SetNColumns( 3 );
    hs.Draw( "hist,nostack" );
    hs.SetMaximum( hs.GetHistogram()->GetMaximum()*1.5 );
    hs.GetHistogram()->SetTitle( p.second[0].GetTitle() );
    c.Prettify( hs.GetHistogram() );
    c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
  }
}
