#include "Canvas.h"
#include "DiphotonAnalyzer/TreeProducer/interface/TreeEvent.h"
#include "DiphotonAnalyzer/Macros/xi_acceptance.h"

#include <iostream>
#include <unordered_map>

#include "TFile.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH1.h"
#include "THStack.h"

#define OUT_PATH "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp/limits"

void plot_sel_efficiency_aqgc()
{
  struct sample_t
  {
    double zeta1, zeta2;
    string path;
    double xsec, xsec_bare;
    unsigned int num_events, num_passing;
  };
  //map<string,float> pots_accept_tight = { { "45N", 0.067 }, { "45F", 0.066 }, { "56N", 0.070 }, { "56F", 0.061 } };
  map<string,float> pots_accept_tight = { { "45N", 0.068 }, { "45F", 0.064 }, { "56N", 0.069 }, { "56F", 0.060 } }; //may20
  xi_accept::load_file( "output_efficiencies.root" );
  xi_accept::vary_uncertainty( -1. ); //nsigma

  //retrieved from /afs/cern.ch/user/j/juwillia/public/forLaurent/AC_microAOD
  vector<sample_t> samples = {
    //{ 0., 0., "/eos/cms/store/user/lforthom/twophoton/samples_80x/output_GammaGammaToGammaGamma_fpmc_justin_sm.root", 1.3988e-05, 100000 },
    //{ 0.,     0.,      "samples/ntuple-sm_gggg.root", 1.3988e-05, 0.1398, 5000 },
    { 0.,     1.e-13,  "samples/ntuple-aqgc-zeta1_0_zeta2_1e-13.root",        1.3068e-4, 9.3727e-3, 5000 },
    { 1.e-13, 0.,      "samples/ntuple-aqgc-zeta1_1e-13_zeta2_0.root",        5.7066e-4, 4.0904e-2, 5000 },
    { 1.e-13, 1.e-13,  "samples/ntuple-aqgc-zeta1_1e-13_zeta2_1e-13.root",    1.1770e-3, 8.4366e-2, 5000 },
    { 1.e-12, 0.,      "samples/ntuple-aqgc-zeta1_1e-12_zeta2_0.root",        5.7055e-2, 4.0906,    5000 },
    { 1.e-12, 1.e-12,  "samples/ntuple-aqgc-zeta1_1e-12_zeta2_1e-12.root",    0.1177,    8.4370,    5000 },
    //{ 5.e-12, 0.,      "samples/ntuples-aqgc-zeta1_5e-12_zeta2_0.root",       5.7055e-2, 102.2665,  5000 },
    { 5.e-12, 0.,      "samples/ntuples-aqgc-zeta1_5e-12_zeta2_0.root",       1.4410,    102.2665,  5000 },
    // negative
    { 1.e-12, -1.e-12,  "samples/ntuple-aqgc-zeta1_1e-12_zeta2_-1e-12.root",  2.2584e-2, 1.6192,    4200 },
    { -1.e-12, -1.e-12, "samples/ntuple-aqgc-zeta1_-1e-12_zeta2_-1e-12.root", 0.1177,    8.4370,    5000 },
    { -1.e-12, 1.e-12,  "samples/ntuple-aqgc-zeta1_-1e-12_zeta2_1e-12.root",  2.2585e-2, 1.6192,    5000 },
    { -1.e-12, 0.,      "samples/ntuple-aqgc-zeta1_-1e-12_zeta2_0.root",      5.7644e-2, 4.0907,    5000 },
    { -5.e-12, -1.e-12, "samples/ntuple-aqgc-zeta1_-5e-12_zeta2_-1e-12.root", 0.,        120.2486,  5000 },
    { -5.e-12, 0.,      "samples/ntuple-aqgc-zeta1_-5e-12_zeta2_0.root",      1.4410,    102.2667,  5000 },
    { -5.e-12, 1.e-12,  "samples/ntuple-aqgc-zeta1_-5e-12_zeta2_1e-12.root",  1.2141,    86.1597,   5000 },
    //{ 0., -5.e-12,      "samples/ntuple-aqgc-zeta1_0_zeta2_-5e-12.root",      0.3302,    23.4362,   5000 },
    { 1.e-12, 5.e-13,   "samples/ntuple-aqgc-zeta1_1e-12_zeta2_5e-13.root",   8.4941e-2, 6.0294,    5000 },
    { 5.e-12, -1.e-12,  "samples/ntuple-aqgc-zeta1_5e-12_zeta2_-1e-12.root",  1.2140,    86.1595,   5000 },
    { 5.e-12, 1.e-12,   "samples/ntuple-aqgc-zeta1_5e-12_zeta2_1e-12.root",   1.6944,    120.2484,  5000 },
    // to check?
    { 1.e-13, 1.e-12,  "samples/ntuple-aqgc-zeta1_1e-13_zeta2_1e-12.root",    1.8399e-2, 1.3192,    5000 },
    { 5.e-13, 0.,      "samples/ntuple-aqgc-zeta1_5e-13_zeta2_0.root",        1.4263e-2, 1.0227,    4900 },
    { 5.e-13, 5.e-13,  "samples/ntuple-aqgc-zeta1_5e-13_zeta2_5e-13.root",    2.9418e-2, 2.1092,    5000 },
    { 5.e-13, 1.e-13,  "samples/ntuple-aqgc-zeta1_5e-13_zeta2_1e-13.root",    1.6771e-2, 1.2025,    5000 },
  };
  /*vector<sample_t> samples = {
    ///// NEW SAMPLES /////
    { 0., 0., "samples/ntuple-sm_gggg_tightxi.root", 1.92274e-7, 0., 10000 },
    //{ 0.,     0.,      "samples/ntuple-sm_gggg.root", 1.3988e-05, 0.1398, 5000 }, //FIXME
    { 1.e-12, 0., "samples/ntuple-aqgc-zeta1_1e-12_zeta2_0_tightxi.root", 3.90508e-2, 0., 10000 },
    //{ 1.e-12, 0.,      "samples/ntuple-aqgc-zeta1_1e-12_zeta2_0.root",        5.7055e-2, 4.0906,    5000 }, //FIXME
    { 0., 1.e-12, "samples/ntuple-aqgc-zeta1_0_zeta2_1e-12_tightxi.root", 8.94720e-3, 0., 10000 },
    { 1.e-13, 1.e-13, "samples/ntuple-aqgc-zeta1_1e-13_zeta2_1e-13_tightxi.root", 8.04886e-4, 0., 10000 },
    //{ 1.e-13, 1.e-13,  "samples/ntuple-aqgc-zeta1_1e-13_zeta2_1e-13.root",    1.1770e-3, 8.4366e-2, 5000 }, //FIXME
    { 1.e-13, 0., "samples/ntuple-aqgc-zeta1_1e-13_zeta2_0_tightxi.root", 3.90226e-4, 0., 10000 },
    //{ 1.e-13, 0.,      "samples/ntuple-aqgc-zeta1_1e-13_zeta2_0.root",        5.7066e-4, 4.0904e-2, 5000 }, //FIXME
  };*/
  const sample_t sm_sample = { 0., 0., "samples/ntuple-sm_gggg_tightxi.root", 1.92274e-7, 0., 10000 };
  //samples.emplace_back( sm_sample );
  const size_t num_samples = samples.size();
  vector<TH1D> v_h_mass( num_samples ), v_h_ptpair( num_samples ), v_h_ptlead( num_samples ), v_h_acop( num_samples ), v_h_xi1( num_samples ), v_h_xi2( num_samples );
  TGraph g_el_vals;
  TGraphErrors g_el_eff_z1, g_el_eff_z2;
  TGraph g_xs_z1, g_xs_z2;
  TGraph2D g_xs, g_el_eff, g_el_acc;
  gggg::TreeEvent ev;
  const float sqrt_s = 13.e3, lumi = 9.36e3;
  const TString ylabel = "Events";
  //const TString ylabel = "Events fraction";
  TLorentzVector ip1, ip2, op1, op2, pho1, pho2;
  ip1.SetXYZM( 0., 0., 6500., 0.938 );
  ip2.SetXYZM( 0., 0., -6500., 0.938 );
  unsigned short i = 0;
  for ( auto& s : samples ) {
    unique_ptr<TFile> file( TFile::Open( s.path.c_str() ) );
    auto tree = dynamic_cast<TTree*>( file->Get( "ntp" ) );
    v_h_mass[i] = TH1D( Form( "mass_%d", i ), Form( "m_{#gamma#gamma}@@%s@@GeV", ylabel.Data() ), 18, 200., 2000. );
    v_h_ptpair[i] = TH1D( Form( "ptpair_%d", i ), Form( "p_{T}^{#gamma#gamma}@@%s@@GeV", ylabel.Data() ), 25, 0., 100. );
    v_h_ptlead[i] = TH1D( Form( "leadpt_%d", i ), Form( "p_{T}^{#gamma} (leading gamma)@@%s@@GeV", ylabel.Data() ), 50, 50., 1050. );
    v_h_acop[i] = TH1D( Form( "acop_%d", i ), Form( "1-|#Delta#phi_{#gamma#gamma}/#pi| (#times 10^{-3})@@%s@@?.g", ylabel.Data() ), 50, 0., 5. );
    v_h_xi1[i] = TH1D( Form( "xi1_%d", i ), Form( "#xi_{1}@@%s@@?.g", ylabel.Data() ), 20, 0.06, 0.16 );
    v_h_xi2[i] = TH1D( Form( "xi2_%d", i ), Form( "#xi_{2}@@%s@@?.g", ylabel.Data() ), 20, 0.06, 0.16 );
    ev.attach( tree );
    //cout << s.xsec << endl;
    double wgt_passing = 0;
    //const double weight = s.xsec/s.num_events*lumi;
    //const double weight = 1./s.num_events;
    s.num_passing = 0;
    for ( long long j = 0; j < tree->GetEntriesFast(); ++j ) {
      tree->GetEntry( j );
      bool has_diphoton = false;
      pho1.SetPtEtaPhiE( ev.gen_part_pt[0], ev.gen_part_eta[0], ev.gen_part_phi[0], ev.gen_part_energy[0] );
      pho2.SetPtEtaPhiE( ev.gen_part_pt[1], ev.gen_part_eta[1], ev.gen_part_phi[1], ev.gen_part_energy[1] );
      op1 = ip1-pho1, op2 = ip2-pho2;
      const double xi1 = 1.-op1.P()/6500., xi2 = 1.-op2.P()/6500.;
      v_h_xi1[i].Fill( xi1 );
      v_h_xi2[i].Fill( xi2 );
      if ( xi1 >= 0.07 && xi1 <= 0.111
        && xi2 >= 0.07 && xi2 <= 0.138 )
        ++s.num_passing;
      else continue;
      //else cout << "j="<<j<<":"<<xi1<<":"<<xi2<<endl;
      if ( ev.hlt_accept[0] == 0 ) continue;
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
        v_h_acop[i].Fill( acop*1.e3 );
        if ( acop > 0.005 ) continue;
        v_h_mass[i].Fill( ev.diphoton_mass[k] );
        v_h_ptpair[i].Fill( ev.diphoton_pt[k] );
        v_h_ptlead[i].Fill( ev.diphoton_pt1[k] );
        // 45 = 0, 56 = 1
        if ( ( ( xim > pots_accept_tight["56N"] || xim > pots_accept_tight["56F"] ) && xim < 0.138 )
          && ( ( xip > pots_accept_tight["45N"] || xip > pots_accept_tight["45F"] ) && xip < 0.111 ) )
          has_diphoton = true;
      }
      if ( has_diphoton ) {
        //const double wgt = xi_accept::weight( xi1, xi2 );
        const double wgt = 1.;
        //cout << xi1 << ":" << xi2 << "|" << wgt << endl;
        wgt_passing += wgt;
      }
    }
    const double el_eff = wgt_passing/s.num_passing, el_eff_unc = el_eff*sqrt( 1./wgt_passing+1./s.num_passing );
    const double el_acc = s.xsec/s.xsec_bare;
    if ( s.zeta1 == 0. ) {
      unsigned short k = g_el_eff_z2.GetN();
      g_el_eff_z2.SetPoint( k, s.zeta2, el_eff );
      g_el_eff_z2.SetPointError( k, 0., el_eff_unc );
      g_xs_z2.SetPoint( k, s.zeta2, s.xsec );
    }
    if ( s.zeta2 == 0. ) {
      unsigned short k = g_el_eff_z1.GetN();
      g_el_eff_z1.SetPoint( k, s.zeta1, el_eff );
      g_el_eff_z1.SetPointError( k, 0., el_eff_unc );
      g_xs_z1.SetPoint( k, s.zeta1, s.xsec );
    }
    //if ( fabs( s.zeta1 ) <= 1.e-12 )
      g_xs.SetPoint( i, s.zeta1*1.e12, s.zeta2*1.e12, s.xsec_bare );
    g_el_eff.SetPoint( i, s.zeta1*1.e12, s.zeta2*1.e12, el_eff );
    if ( el_acc > 0.01 ) //FIXME
      g_el_acc.SetPoint( g_el_acc.GetN(), s.zeta1*1.e12, s.zeta2*1.e12, el_acc );
    g_el_vals.SetPoint( i, s.zeta1*1.e12, s.zeta2*1.e12 );
    cout << s.path << "(" << s.num_events << ":" << s.num_passing << "): " << el_eff << "|" << el_acc << endl;
    ++i;
  }
  g_el_vals.SetPoint( g_el_vals.GetN(), 0., 0. );
  //gStyle->SetPalette( kSunset );
  gStyle->SetPalette( kBeach );
  //const char* top_label = "Elastic selection";
  const char* top_label = "#xi^{#pm}_{#gamma#gamma} in acc., #epsilon(#xi^{#pm}_{#gamma#gamma}) > 90%";
  {
    Canvas c( "elastic_xsec_scan_aqgc", "FPMC #gamma#gamma #rightarrow #gamma#gamma, AQGC", "CMS", "Simulation preliminary" );
    c.SetRightMargin( 0.15 );
    g_xs.SetTitle( ";#zeta_{1} (#times 10^{-12});#zeta_{2} (#times 10^{-12});Signal cross section (pb)" );
    g_el_vals.SetMarkerStyle( 20 );
    g_xs.Draw( "colz" );
//    g_xs.SetMaximum( 20. );
    g_el_vals.Draw( "p same" );
    g_el_vals.SetMarkerStyle( 3 );
    g_xs.GetZaxis()->SetTitleOffset( 1.4 );
//    PaveText::topLabel( top_label );
    c.Prettify( (TH1*)g_el_vals.GetHistogram() );
//    g_xs.GetYaxis()->SetLimits( -5.e-12, 5.e-12 );
//    g_el_vals.GetYaxis()->SetRangeUser( -5.e-12, 5.e-12 );
//    g_xs.GetYaxis()->SetRangeUser( -5.e-12, 5.e-12 );
    c.SetLogz();
    c.Save( "pdf,png", OUT_PATH );
  }
  {
    Canvas c( "elastic_acc_scan_aqgc", "FPMC #gamma#gamma #rightarrow #gamma#gamma, AQGC", "CMS", "Simulation preliminary" );
    c.SetRightMargin( 0.15 );
    g_el_acc.SetTitle( ";#zeta_{1} (#times 10^{-12});#zeta_{2} (#times 10^{-12});Signal acceptance" );
    g_el_vals.SetMarkerStyle( 20 );
    g_el_acc.Draw( "colz" );
    //g_el_acc.SetMinimum( 0.4 ); g_el_acc.SetMaximum( 0.455 );
    //g_el_acc.SetMinimum( 0.01 );
    //g_el_acc.SetMaximum( 0.02 );
    g_el_vals.Draw( "p same" );
    g_el_vals.SetMarkerStyle( 3 );
    g_el_acc.GetZaxis()->SetTitleOffset( 1.4 );
    c.SetLogz();
    c.Prettify( (TH1*)g_el_vals.GetHistogram() );
    c.Save( "pdf,png", OUT_PATH );
  }
  {
    Canvas c( "elastic_eff_scan_aqgc", "FPMC #gamma#gamma #rightarrow #gamma#gamma, AQGC", "CMS", "Simulation preliminary" );
    c.SetRightMargin( 0.15 );
    g_el_eff.SetTitle( ";#zeta_{1} (#times 10^{-12});#zeta_{2} (#times 10^{-12});Sel.efficiency" );
    g_el_vals.SetMarkerStyle( 20 );
    //g_el_vals.Draw( "a*" );
    //g_el_eff.Draw( "arr colz" );
    g_el_eff.Draw( "colz" );
    //g_el_eff.SetMinimum( 0.5 ); g_el_eff.SetMaximum( 0.6 );
    //g_el_eff.SetMinimum( 0.35 ); g_el_eff.SetMaximum( 0.45 );
    //g_el_eff.SetMinimum( 0.4 ); g_el_eff.SetMaximum( 0.455 );
    //g_el_eff.SetMinimum( 0.587 ); g_el_eff.SetMaximum( 0.612 );
    g_el_eff.SetMinimum( 0.62 ); g_el_eff.SetMaximum( 0.68 );
    //c.SetLogx(); c.SetLogy();
    g_el_vals.Draw( "p same" );
    g_el_vals.SetMarkerStyle( 3 );
    g_el_eff.GetZaxis()->SetTitleOffset( 1.45 );
    PaveText::topLabel( top_label );
    c.Prettify( (TH1*)g_el_vals.GetHistogram() );
    c.Save( "pdf,png", OUT_PATH );
  }
  {
    TF1 f_pol2( "f_pol2", "[0]+[1]*x+[2]*x**2" );
    Canvas c( "elastic_xsec_scan_aqgc_1d", "FPMC #gamma#gamma #rightarrow #gamma#gamma, AQGC", "CMS", "Simulation preliminary" );
    c.SetLegendX1( 0.68 );
    g_xs_z1.SetTitle( ";#zeta_{1,2};Visible cross section (pb)" );
    g_xs_z1.Draw( "ap" );
    g_xs_z1.SetMarkerStyle( 24 );
    //g_xs_z1.Fit( "pol2" );
    //g_xs_z1.Fit( &f_pol2, "s" );
    g_xs_z2.Draw( "p,same" );
    g_xs_z2.SetMarkerStyle( 20 );
    //g_xs_z2.Fit( "pol2" );
    //g_xs_z2.Fit( &f_pol2, "s+" );
    c.AddLegendEntry( &g_xs_z1, "#zeta_{2} = 0", "ep" );
    c.AddLegendEntry( &g_xs_z2, "#zeta_{1} = 0", "ep" );
    PaveText::topLabel( top_label );
    c.Prettify( g_xs_z1.GetHistogram() );
    c.Save( "pdf,png", OUT_PATH );
  }
  {
    Canvas c( "elastic_eff_scan_aqgc_1d", "FPMC #gamma#gamma #rightarrow #gamma#gamma, AQGC", "CMS", "Simulation preliminary" );
    gStyle->SetOptFit( 111 );
    c.SetLegendX1( 0.68 );
    auto pol01 = new TF1( "pol01", "pol0" ), pol02 = (TF1*)pol01->Clone( "pol02" );
    g_el_eff_z1.SetTitle( ";#zeta_{1,2};Selection efficiency" );
    g_el_eff_z1.Draw( "ap" );
    g_el_eff_z1.SetMarkerStyle( 20 );
    g_el_eff_z1.SetMarkerColor( kRed+1 );
    g_el_eff_z1.SetMinimum( 0. );
    g_el_eff_z1.SetMaximum( 1.05 );
    g_el_eff_z1.Fit( pol01, "s" );
    pol01->SetLineColor( kRed+1 );
    pol01->SetLineWidth( 2 );
    pol01->DrawF1( g_el_eff_z1.GetXaxis()->GetXmin(), g_el_eff_z1.GetXaxis()->GetXmax(), "same" );
    g_el_eff_z2.Draw( "p,same" );
    g_el_eff_z2.SetMarkerStyle( 21 );
    g_el_eff_z2.SetMarkerColor( kBlue-2 );
    g_el_eff_z2.Fit( pol02, "s" );
    pol02->SetLineColor( kBlue-2 );
    pol02->SetLineStyle( 2 );
    pol02->SetLineWidth( 4 );
    pol02->DrawF1( g_el_eff_z1.GetXaxis()->GetXmin(), g_el_eff_z2.GetXaxis()->GetXmax(), "same" );
    c.AddLegendEntry( &g_el_eff_z1, "#zeta_{2} = 0", "ep" );
    c.AddLegendEntry( &g_el_eff_z2, "#zeta_{1} = 0", "ep" );
    gPad->Update();
    auto stat1 = (TPaveStats*)g_el_eff_z1.GetListOfFunctions()->FindObject( "stats" );
    stat1->SetTextSize( 0.03 );
    stat1->SetY1NDC( 0.33 );
    stat1->SetY2NDC( 0.48 );
    stat1->SetTextColor( kRed+1 );
    stat1->Draw();
    auto stat2 = (TPaveStats*)g_el_eff_z2.GetListOfFunctions()->FindObject( "stats" );
    stat2->SetY1NDC( 0.17 );
    stat2->SetY2NDC( 0.32 );
    stat2->SetTextColor( kBlue-2 );
    stat2->SetTextSize( 0.03 );
    stat2->Draw();
    PaveText::topLabel( top_label );
    c.Prettify( g_el_eff_z1.GetHistogram() );
    c.Save( "pdf,png", OUT_PATH );
  }

  const unordered_map<string,vector<TH1D>&> plots = {
    { "anom_aqgc_elastic_ptpair", v_h_ptpair },
    { "anom_aqgc_elastic_acop", v_h_acop },
    { "anom_aqgc_elastic_ptlead", v_h_ptlead },
    { "anom_aqgc_elastic_mpair", v_h_mass },
    { "anom_aqgc_elastic_xi1", v_h_xi1 }, { "anom_aqgc_elastic_xi2", v_h_xi2 },
  };
  i = 0;
  for ( auto& p : plots ) {
    Canvas c( p.first.c_str(), Form( "FPMC #gamma#gamma #rightarrow #gamma#gamma (AQGC), yields for %.1f fb^{-1} (%g TeV)", lumi/1.e3, sqrt_s/1.e3 ), "CMS", "Simulation preliminary" );
    //Canvas c( p.first.c_str(), Form( "FPMC #gamma#gamma #rightarrow #gamma#gamma (%g TeV)", sqrt_s/1.e3 ), "Simulation preliminary" );
    THStack hs;
    /*if ( i++ > 2 )
      c.SetLegendY1( 0.695 );
    else {
      c.SetLegendY1( 0.61 );
      c.SetLegendX1( 0.18 );
    }
    c.SetLegendSizeY( 0.2 );*/
    /*c.SetLegendY1( 0.68 );
    c.SetLegendX1( 0.18 );
    c.SetLegendSizeX( 0.7 );
    c.SetLegendSizeY( 0.15 );*/
    c.SetLegendY1( 0.62 );
    c.SetLegendX1( 0.17 );
    c.SetLegendSizeX( 0.7 );
    c.SetLegendSizeY( 0.225 );
    for ( unsigned short j = 0; j < num_samples; ++j ) {
      //TString label = Form( "#zeta_{1} = %g, #zeta_{2} = %g", samples[j].zeta1, samples[j].zeta2 );
      TString label = Form( "(%g, %g)", samples[j].zeta1, samples[j].zeta2 );
      if ( samples[j].zeta1 == 0. && samples[j].zeta2 == 0. ) {
        label = "SM";
        p.second[j].SetLineColor( 1 );
        p.second[j].SetLineStyle( 2 );
        p.second[j].SetFillColor( 0 );
      }
      else {
        p.second[j].SetLineColor( Canvas::colour_pool[(j-1)/3] );
        p.second[j].SetLineStyle( 1+(j-1)%3 );
        p.second[j].SetFillColorAlpha( Canvas::colour_pool[(j-1)/3], 0.25 );
      }
      p.second[j].Scale( samples[j].xsec/samples[j].num_passing*lumi );
      p.second[j].SetLineWidth( 2 );
      hs.Add( &p.second[j] );
      c.AddLegendEntry( &p.second[j], label, "f" );
    }
    c.GetLegend()->SetNColumns( 3 );
    hs.Draw( "hist,nostack" );
    hs.SetMaximum( hs.GetHistogram()->GetMaximum()*1.5 );
    hs.GetHistogram()->SetTitle( p.second[0].GetTitle() );
    c.Prettify( hs.GetHistogram() );
    if ( p.first == "anom_aqgc_elastic_ptpair"
      || p.first == "anom_aqgc_elastic_mpair" ) {
      c.SetLogy();
      hs.SetMaximum( hs.GetHistogram()->GetMaximum()*25. );
      hs.SetMinimum( 5.e-2 );
    }
    c.Save( "pdf,png", OUT_PATH );
  }
}
