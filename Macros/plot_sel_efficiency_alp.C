#include "Canvas.h"
#include "Plotter.h"
#include "DiphotonAnalyzer/TreeProducer/interface/TreeEvent.h"

#include <iostream>
#include <unordered_map>

#include "TFile.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH1.h"
#include "THStack.h"

void plot_sel_efficiency_alp()
{
  struct sample_t
  {
    double m, k;
    string path;
    double xsec;
    unsigned int num_events;
  };
  map<string,float> pots_accept_tight = { { "45N", 0.067 }, { "45F", 0.066 }, { "56N", 0.070 }, { "56F", 0.061 } };

  //retrieved from /afs/cern.ch/user/j/juwillia/public/forLaurent/AC_microAOD
  double k2 = 0.031622777; // 1.e-1.5
  vector<sample_t> samples = {
    { 500.,  -1., "samples/ntuple-alp-m500_f1e-1.root", 2.310e-3, 5000 },
    { 1000., -1., "samples/ntuple-alp-m1000_f1e-1.root", 9.534e-4, 5000 },
    { 1500., -1., "samples/ntuple-alp-m1500_f1e-1.root", 1.968e-4, 5000 },
    { 2000., -1., "samples/ntuple-alp-m2000_f1e-1.root", 3.7767e-7, 5000 },
    { 500.,  -1.5, "samples/ntuple-alp-m500_f1e-1p5.root", 2.297e-4, 5000 },
    { 1000., -1.5, "samples/ntuple-alp-m1000_f1e-1p5.root", 9.595e-5, 5000 },
    { 1500., -1.5, "samples/ntuple-alp-m1500_f1e-1p5.root", 1.968e-5, 5000 },
    //{ 2000., -1.5, "samples/ntuple-alp-m2000_f1e-1p5.root", 3.5524-9, 5000 },
    { 500.,  -0.301, "samples/ntuple-alp-m500_f5e-1.root", 5.884e-2, 5000 },
    { 1000., -0.301, "samples/ntuple-alp-m1000_f5e-1.root", 2.382e-2, 5000 },
    { 1500., -0.301, "samples/ntuple-alp-m1500_f5e-1.root", 4.823e-3, 5000 },
    { 2000., -0.301, "samples/ntuple-alp-m2000_f5e-1.root", 2.1620e-4, 5000 },
  };
  const size_t num_samples = samples.size();
  vector<TH1D> v_h_mass( num_samples ), v_h_ptpair( num_samples ), v_h_ptlead( num_samples ), v_h_acop( num_samples );
  TGraph g_el_eff_vals;
  map<double,TGraphErrors> g_el_eff_m;
  TGraph2D g_el_eff;
  gggg::TreeEvent ev;
  const float sqrt_s = 13.e3, lumi = 9.36e3;
  const TString ylabel = "Events";
  //const TString ylabel = "Events fraction";
  unsigned short i = 0;
  for ( const auto& s : samples ) {
    unique_ptr<TFile> file( TFile::Open( s.path.c_str() ) );
    auto tree = dynamic_cast<TTree*>( file->Get( "ntp" ) );
    v_h_mass[i] = TH1D( Form( "mass_%d", i ), Form( "m_{#gamma#gamma}@@%s@@GeV", ylabel.Data() ), 50, 200., 2200. );
    v_h_ptpair[i] = TH1D( Form( "ptpair_%d", i ), Form( "p_{T}^{#gamma#gamma}@@%s@@GeV", ylabel.Data() ), 25, 0., 100. );
    v_h_ptlead[i] = TH1D( Form( "leadpt_%d", i ), Form( "p_{T}^{#gamma} (leading gamma)@@%s@@GeV", ylabel.Data() ), 50, 50., 1050. );
    v_h_acop[i] = TH1D( Form( "acop_%d", i ), Form( "1-|#Delta#phi_{#gamma#gamma}/#pi| (#times 10^{-3})@@%s@@?.g", ylabel.Data() ), 50, 0., 5. );
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
        v_h_acop[i].Fill( acop*1.e3, weight );
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
    {
      unsigned short k = g_el_eff_m[s.k].GetN();
      g_el_eff_m[s.k].SetPoint( k, s.m, el_eff );
      g_el_eff_m[s.k].SetPointError( k, 0., el_eff_unc );
    }
    g_el_eff.SetPoint( i, s.m, s.k, el_eff );
    g_el_eff_vals.SetPoint( i, s.m, s.k );
    cout << s.path << ": " << el_eff << endl;
    ++i;
  }
  g_el_eff_vals.SetPoint( g_el_eff_vals.GetN(), 0., 0. );
  //gStyle->SetPalette( kSunset );
  gStyle->SetPalette( kBeach );
  //const char* top_label = "Elastic selection";
  const char* top_label = "#xi^{#pm}_{#gamma#gamma} in acc., #epsilon(#xi^{#pm}_{#gamma#gamma}) > 90%";
  {
    Canvas c( "elastic_eff_scan_alp", "FPMC #gamma#gamma #rightarrow #gamma#gamma + ALP", "Simulation preliminary" );
    c.SetRightMargin( 0.15 );
    g_el_eff.SetTitle( ";m_{a} (GeV);f^{-1} (XXX^{-1});#varepsilon" );
    g_el_eff_vals.SetMarkerStyle( 20 );
    //g_el_eff_vals.Draw( "a*" );
    //g_el_eff.Draw( "arr colz" );
    g_el_eff.Draw( "colz" );
    //g_el_eff.SetMinimum( 0.5 ); g_el_eff.SetMaximum( 0.6 );
    //g_el_eff.SetMinimum( 0.35 ); g_el_eff.SetMaximum( 0.45 );
    g_el_eff.SetMinimum( 0. ); g_el_eff.SetMaximum( 0.65 );
    //c.SetLogx(); c.SetLogy();
    g_el_eff_vals.Draw( "p same" );
    g_el_eff_vals.SetMarkerStyle( 3 );
    PaveText::topLabel( top_label );
    c.Prettify( (TH1*)g_el_eff_vals.GetHistogram() );
    c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
  }
  {
    Canvas c( "elastic_eff_scan_alp_1d", "FPMC #gamma#gamma #rightarrow #gamma#gamma + ALP", "Simulation preliminary" );
    TMultiGraph mg;
    unsigned short i = 0;
    TF1 f_fit( "f_fit", "pol2" );
    vector<TF1*> fits_res;
    for ( auto& g : g_el_eff_m ) {
      mg.Add( &g.second );
      g.second.SetMarkerStyle( 24+i );
      TFitResultPtr fit_res = g.second.Fit( &f_fit, "ls0+" );
      if ( (int)fit_res ) {
        f_fit.SetLineStyle( i+1 );
        fits_res.emplace_back( (TF1*)f_fit.Clone() );
      }
      c.AddLegendEntry( &g.second, Form( "k = %g", g.first ), "ep" );
      ++i;
    }
    mg.Draw( "ap" );
    for ( const auto& f : fits_res )
      f->DrawF1( 500., 2000., "same" );
    mg.SetMinimum( 0. );
    mg.SetMaximum( 1.05 );
    mg.GetHistogram()->SetTitle( ";m_{a} (GeV);Selection efficiency" );
    PaveText::topLabel( top_label );
    c.Prettify( mg.GetHistogram() );
    c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
  }

  const unordered_map<string,vector<TH1D>&> plots = {
    { "anom_alp_elastic_ptpair", v_h_ptpair },
    { "anom_alp_elastic_acop", v_h_acop },
    { "anom_alp_elastic_ptlead", v_h_ptlead },
    { "anom_alp_elastic_mpair", v_h_mass }
  };
  i = 0;
  for ( auto& p : plots ) {
    Canvas c( p.first.c_str(), Form( "FPMC #gamma#gamma #rightarrow #gamma#gamma (ALP), yields for %.1f fb^{-1} (%g TeV)", lumi/1.e3, sqrt_s/1.e3 ), "Simulation preliminary" );
    //Canvas c( p.first.c_str(), Form( "FPMC #gamma#gamma #rightarrow #gamma#gamma (%g TeV)", sqrt_s/1.e3 ), "Simulation preliminary" );
    THStack hs;
    //if ( i > 1 )
      //c.SetLegendY1( 0.695 );
    //else {
      //c.SetLegendY1( 0.61 );
      c.SetLegendY1( 0.68 );
      c.SetLegendX1( 0.18 );
    //  c.SetLegendX1( 0.18 );
    //}
    c.SetLegendSizeX( 0.7 );
    c.SetLegendSizeY( 0.15 );
    for ( unsigned short j = 0; j < num_samples; ++j ) {
      //TString label = Form( "#zeta_{1} = %g, #zeta_{2} = %g", samples[j].zeta1, samples[j].zeta2 );
      TString label = Form( "(%g, %g)", samples[j].m, samples[j].k );
      p.second[j].SetLineColor( Plotter::colour_pool[j] );
      p.second[j].SetFillColorAlpha( Plotter::colour_pool[j], 0.25 );
      p.second[j].SetLineWidth( 2 );
      hs.Add( &p.second[j] );
      c.AddLegendEntry( &p.second[j], label, "f" );
    }
    c.GetLegend()->SetNColumns( 3 );
    hs.Draw( "hist,nostack" );
    hs.SetMaximum( hs.GetHistogram()->GetMaximum()*20. );
    hs.GetHistogram()->SetTitle( p.second[0].GetTitle() );
    //if ( i == 0 ) { // mass plot
      hs.SetMinimum( 1.e-2 );
      c.SetLogy();
    //}
    c.Prettify( hs.GetHistogram() );
    c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
    ++i;
  }
}
