#include "DatasetHandler.h"
#include "EventsSelector.h"
#include "Plotter.h"
#include "Canvas.h"
#include "mass_templates.h"
#include "PhotonScalesParser.h"

#include "DiphotonAnalyzer/TreeProducer/interface/TreeEvent.h"

#include "TGaxis.h"
#include "THStack.h"
#include "TH1.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "TStyle.h"

#include "TMinuit.h"
#include "TVirtualFitter.h"

#include <string.h>
#include <iomanip>

using namespace gggg;

void logarithmicBins( TAxis* axis );
void plot_2ddiscrim( const char* name, TH2D* h2[], bool logx );

typedef enum {
  the_data = 0,
  mc_inclusive = 1,
  mc_signal = 2,
  num_types
} sample_types;

double loosemvaeff( double r9, double eta_sc ) {
  if ( fabs( eta_sc ) < 1.5 ) {
    if ( r9 < 0.85 ) return 0.9999;
    return 1.0003;
  }
  if ( r9 < 0.9 ) return 1.0003;
  return 1.0004;
}
const float sqrt_s = 13.e3;
//const float lumi = ( 5.055026851041 + 1.470474265456 + 7.487318307770 ) * 1.e3;
//const float lumi = ( 5.060574552184 + 1.485056762163 + 7.500305757473 ) * 1.e3; // computed from processedLumis.json
//const float the_lumi = ( 9.412003739742+5.08 ) * 1.e3; // pre+post-TS2 runs with pots inserted
//const float the_lumi = 5.081503324822e3; //post-TS2 runs with pots inserted
//const float the_lumi = 4.291753570355 * 1.e3; // B, pre-TS2 runs with pots inserted
//const float the_lumi = 9.412003739742 * 1.e3; // pre-TS2 runs with pots inserted
const float the_lumi = ( 4.291753570355+1.438424638452+3.625791461729 ) * 1.e3; //FIXME FIXME FIXME
//const float the_lumi = ( 9.732907106 ) * 1.e3;

map<string,float> pots_accept = { { "45N", 0.033 }, { "45F", 0.024 }, { "56N", 0.050 }, { "56F", 0.037 } };
map<string,float> pots_accept_tight = { { "45N", 0.067 }, { "45F", 0.066 }, { "56N", 0.070 }, { "56F", 0.061 } };

void
plotter()
{
  //-----------------------------------------------------------------------------------------------
  const float pt_cut = 75.;
  const float eta_cut = 2.5, min_etaveto = 1.4442, max_etaveto = 1.566;
  const float r9_cut = 0.94;
  const float mass_cut = 350.;
  //-----------------------------------------------------------------------------------------------
  const float scaling_signal = 5.e3;
  const float rescaling_background = 0.91;
  //const float rescaling_background = 1.;
  const bool rescale_background_parametric = false;
  const double alpha = 0.4; // plotting transparency
  //-----------------------------------------------------------------------------------------------

  TF1* incl_fit_func = nullptr;
  if ( rescale_background_parametric )
    incl_fit_func = (TF1*)( TFile::Open( "incl_fit_result.root" )->Get( "fit_func" ) );

  // turn off "Info in TCanvas::Print: file XXX has been created" messages
  //gErrorIgnoreLevel = kWarning;
  gErrorIgnoreLevel = kBreak;
  gStyle->SetEndErrorSize( 0. ); // remove ticks at end of error bars

  gSystem->Load( "libEventFilterUtilities.so" );
  DatasetHandler dsh( "datasets_list.json" );
  EventsSelector ev_selector( "/afs/cern.ch/user/l/lforthom/public/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_PPSruns_preTS2.txt" );

  const unsigned short num_samples = dsh.size();

  typedef enum {
    nosel = 0, presel,
    elastic, xicomp, xitight,
    qcd, incl,
    presel_noneleveto, singlevtx,
    num_regions
  } regions;
  const vector<string> kin_regions = {
    "nosel", "presel",
    "elastic", "xicomp", "xitight",
    "qcd", "incl",
    "presel_noneleveto", "singlevtx"
  };
  const vector<string> kin_regions_label = {
    "After HLT", "Preselection",
    "Elastic selection", "#xi^{#pm}_{#gamma#gamma} in acceptance", "#xi^{#pm}_{#gamma#gamma} in acc., #epsilon(#xi^{#pm}_{#gamma#gamma}) > 90%",
    "QCD selection", "Inclusive selection",
    "Preselection (no electron veto)", "Single vertex in event"
  };
  const vector<string> kin_regions_label_small = {
    "Preselect.", "Presel.",
    "Elastic", "Accept. #xi^{#pm}_{#gamma#gamma}", "Tight #xi^{#pm}_{#gamma#gamma}",
    "QCD sel.", "Inclusive",
    "Presel. (no e-veto)", "Single vertex"
  };
  const vector<string> types_names = { "data", "backgr", "signal" };
  const vector<regions> regions_in_cutflow = { nosel, incl, elastic, xicomp, xitight };

  //const PhotonScalesParser pho_scales( "/afs/cern.ch/work/l/lforthom/private/twophoton/CMSSW_8_0_26_patch1/src/DiphotonAnalyzer/egammaEffi.txt_EGM2D.root" );
  //const PhotonScalesParser pho_scales( "Fall17V2_2016_MVAwp80_photons.root" );
  const PhotonScalesParser pho_scales( "effLooseMvaBins.root", "h_eff" );

  TH1D* h_mass[num_regions][num_samples],
    *h_mass_rebin[num_regions][num_samples],
    *h_ptpair[num_regions][num_samples],
    *h_dphi[num_regions][num_samples],
    *h_dphi_varbins[num_regions][num_samples],
    *h_dphi_logbins[num_regions][num_samples],
    *h_dphi_tightbins[num_regions][num_samples],
    *h_met[num_regions][num_samples],
    *h_ptlead[num_regions][num_samples],
    *h_ptsublead[num_regions][num_samples],
    *h_pt[num_regions][num_samples],
    *h_pt_varbins[num_regions][num_samples],
    *h_etalead[num_regions][num_samples],
    *h_etasublead[num_regions][num_samples],
    *h_eta[num_regions][num_samples],
    *h_r9lead[num_regions][num_samples],
    *h_r9sublead[num_regions][num_samples],
    *h_r9[num_regions][num_samples],
    *h_xip[num_regions][num_samples],
    *h_xim[num_regions][num_samples],
    *h_diph_vtxz[num_regions][num_samples],
    *h_ndiph[num_regions][num_samples],
    *h_nvtx[num_regions][num_samples],
    *h_nvtx_unc[num_regions][num_samples],
    *h_diph_numjets[num_regions][num_samples],
    *h_diph_numleptons[num_regions][num_samples];
  TH2D* h2_excl_sel[num_types],
    *h2_excl_acop_dpt[num_types],
    *h2_excl_jet[num_types],
    *h2_eleveto[num_samples];
  TH1D* h_diph_nvtxaround[4];
  //TH1D* h_fwdtrk_x[num_regions][num_samples], h_fwdtrk_y[num_regions][num_samples];
  TH1D* h_cutflow[num_samples];

  for ( unsigned short i = 0; i < num_types; ++i ) {
    h2_excl_sel[i] = new TH2D( Form( "exclusivity_cuts_%d", i ), "Distance to the closest lepton vertex (cm)@@Acoplanarity 1-#||{#Delta#phi/#pi}", 10, -2.5, 0.8, 10, -5., 0. );
    h2_excl_acop_dpt[i] = new TH2D( Form( "excl_cuts_acop_dpt_%d", i ), "#Deltap_{T}(#gamma,#gamma) (GeV)@@Acoplanarity 1-#||{#Delta#phi/#pi}", 10, -2., 2.3, 10, -5., 0. );
    h2_excl_jet[i] = new TH2D( Form( "excl_cuts_jet_%d", i ), "High-p_{T} jets associated to the #gamma#gamma vertex@@Acoplanarity 1-#||{#Delta#phi/#pi}", 5, 0., 5., 10, -5., 0. );
    logarithmicBins( h2_excl_sel[i]->GetXaxis() );
    logarithmicBins( h2_excl_sel[i]->GetYaxis() );
    logarithmicBins( h2_excl_acop_dpt[i]->GetXaxis() );
    logarithmicBins( h2_excl_acop_dpt[i]->GetYaxis() );
    logarithmicBins( h2_excl_jet[i]->GetYaxis() );
  }
  for ( unsigned short i = 0; i < num_samples; ++i )
    h2_eleveto[i] = new TH2D( Form( "ele_veto_%d", i ), "Electron veto, leading #gamma@@Electron veto, subleading #gamma", 2, 0., 2., 2, 0., 2. );
  for ( unsigned short i = 0; i < 4; ++i )
    h_diph_nvtxaround[i] = new TH1D( Form( "h_diph_nvtxaround_%d", i ), "Number of vertices surrounding diphoton vertex@@Events", 10, 0., 10. );

  TF1 f_expo( "my_expo", "[0]*exp(-[1]*x)", 0., 0.25 );
  //TF1 f_expo( "my_expo", "[0]+[1]/x", 0., 0.25 );
  //TF1 f_expo( "my_expo", "[0]+[1]/x**2", 0., 0.25 );
  //TF1 f_expo( "my_expo", "[0]*exp(-[1]*x**2)", 0., 0.25 );
  f_expo.SetParName( 0, "Intercept" );
  f_expo.SetParName( 1, "Slope" );

  Plotter::HistsMap
    hm_ndiph[num_regions][3],
    hm_nvtx[num_regions][3], hm_nvtx_unc[num_regions][3],
    hm_mass[num_regions][3], hm_mass_rebin[num_regions][3],
    hm_ptpair[num_regions][3],
    hm_dphi[num_regions][3], hm_dphi_varbins[num_regions][3], hm_dphi_logbins[num_regions][3], hm_dphi_tightbins[num_regions][3],
    hm_met[num_regions][3],
    hm_ptlead[num_regions][3], hm_ptsublead[num_regions][3], hm_pt[num_regions][3], hm_pt_varbins[num_regions][3],
    hm_etalead[num_regions][3], hm_etasublead[num_regions][3], hm_eta[num_regions][3],
    hm_r9lead[num_regions][3], hm_r9sublead[num_regions][3], hm_r9[num_regions][3],
    hm_xip[num_regions][3], hm_xim[num_regions][3],
    hm_diph_vtxz[num_regions][3], hm_diph_nvtxaround[num_regions][1],
    hm_diph_numjets[num_regions][3], hm_diph_numleptons[num_regions][3],
    hm_fwdtrk_x[num_regions][3], hm_fwdtrk_y[num_regions][3],
    hm_cutflow[3];

  map<string,float> in_pot_accept;
  for ( const auto& pot : pots_accept )
    in_pot_accept[pot.first] = 0.; // first initialise the counter

  double evts_sel[num_samples][num_regions], num_after[num_regions][num_samples];
  bool rescaled[num_samples];

  for ( unsigned short i = 0; i < num_samples; ++i ) {
    Sample s = dsh.sample( i );
    cerr << ">>> Processing " << s.name() << endl;
    float sample_weight = s.type() == Sample::kData
      ? 1. : s.cross_section()/s.num_events()*the_lumi;

    for ( unsigned short j = 0; j < num_regions; ++j )
      evts_sel[i][j] = num_after[j][i] = 0.;

    sample_types sample_type = the_data;
    switch ( s.type() ) {
      case Sample::kData: default:
        sample_type = the_data;
        break;
      case Sample::kBackground:
        sample_type = mc_inclusive;
        break;
      case Sample::kSignal:
        sample_type = mc_signal;
        break;
    }

    const vector<double> pt_bins_arr = { 75., 100., 150., 250., 500. };
    //const vector<double> acop_bins_arr = { 0., 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1. };
    const vector<double> acop_bins_arr = { 0., 0.001, 0.002, 0.005, 0.01, 0.05, 0.25, 1. };
    const vector<double> lacop_bins_arr = { -6., -4., -2.5, -2., -1.5, -1.25, -1., -0.75, -0.5, -0.25, 0. };
    { // same plots features for all regions
      unsigned short j = 0;
      for ( const auto& reg : kin_regions ) {
        const char* region = reg.c_str();
        h_ndiph[j][i] = new TH1D( Form( "h_%s_ndiph_%d", region, i ), "Number of diphoton candidates@@Events", 7, 1., 8. );
        h_nvtx[j][i] = new TH1D( Form( "h_%s_nvtx_%d", region, i ), "Number of primary vertices@@Events", 50, 0., 50. );
        h_pt_varbins[j][i] = new TH1D( Form( "h_%s_pt_varbins_%d", region, i ), "p_{T}^{#gamma}@@Events@@GeV", pt_bins_arr.size()-1, pt_bins_arr.data() );
        h_dphi_varbins[j][i] = new TH1D( Form( "h_%s_dphi_varbins_%d", region, i ), "1-|#Delta#phi_{#gamma#gamma}/#pi|@@Events@@?.g", acop_bins_arr.size()-1, acop_bins_arr.data() );
        //h_dphi_logbins[j][i] = new TH1D( Form( "h_%s_dphi_logbins_%d", region, i ), "log_{10}(1-|#Delta#phi_{#gamma#gamma}/#pi|)@@Events@@?.g", 20, -4., 0. );
        h_dphi_logbins[j][i] = new TH1D( Form( "h_%s_dphi_logbins_%d", region, i ), "log_{10}(1-|#Delta#phi_{#gamma#gamma}/#pi|)@@Events@@?.g", lacop_bins_arr.size()-1, lacop_bins_arr.data() );
        h_dphi_tightbins[j][i] = new TH1D( Form( "h_%s_dphi_tightbins_%d", region, i ), "1-|#Delta#phi_{#gamma#gamma}/#pi| (#times 10^{-3})@@Events@@?.g", 25, 0., 5. );
        h_diph_numjets[j][i] = new TH1D( Form( "h_%s_diph_numjets_%d", region, i ), "Number of jets surrounding diphoton vertex@@Events", 10, 0., 10. );
        h_diph_numleptons[j][i] = new TH1D( Form( "h_%s_diph_numleptons_%d", region, i ), "Number of leptons surrounding diphoton vertex@@Events", 10, 0., 10. );
        //h_fwdtrk_x[j][i] = new TH1D( Form( "h_%s_fwdtrk_x_%d", region, i ), "Forward track x@@Events@@mm", 100, 0., 50. );
        h_mass_rebin[j][i] = new TH1D( Form( "h_%s_mass_rebin_%d", region, i ), ";m (GeV);Entries", alp::mass_bins.size()-1, alp::mass_bins.data() );
        ++j;
      }
      h_cutflow[i] = new TH1D( Form( "h_cutflow_%d", i ), ".@@Events", regions_in_cutflow.size(), 0., regions_in_cutflow.size() );
      //h_cutflow[i]->SetLabelSize( 0.02, "X" );
      //h_cutflow[i] = new TH1D( Form( "h_cutflow_%d", i ), ".@@Events", kin_regions.size(), 0., kin_regions.size() );
      //-------------------------------------------------------------------------------------------
      //h_mass[nosel][i] = new TH1D( Form( "h_%s_mass_%d", kin_regions[nosel].c_str(), i ), "m_{#gamma#gamma}@@Events@@GeV", 17, 350., 2050. );
      //h_mass[nosel][i] = new TH1D( Form( "h_%s_mass_%d", kin_regions[nosel].c_str(), i ), "m_{#gamma#gamma}@@Events@@GeV", 28, 350., 1050. );
      h_mass[nosel][i] = new TH1D( Form( "h_%s_mass_%d", kin_regions[nosel].c_str(), i ), "m_{#gamma#gamma}@@Events@@GeV", 34, 300., 2000. );
      h_mass[qcd][i] = new TH1D( Form( "h_%s_mass_%d", kin_regions[qcd].c_str(), i ), "m_{#gamma#gamma}@@Events@@GeV", 8, 300., 900. );
      h_mass[presel][i] = new TH1D( Form( "h_%s_mass_%d", kin_regions[presel].c_str(), i ), "m_{#gamma#gamma}@@Events@@GeV", 24, 350., 2750. );
      h_mass[presel_noneleveto][i] = (TH1D*)h_mass[presel][i]->Clone();
      h_mass[singlevtx][i] = (TH1D*)h_mass[presel][i]->Clone();
      h_mass[incl][i] = new TH1D( Form( "h_%s_mass_%d", kin_regions[incl].c_str(), i ), "m_{#gamma#gamma}@@Events@@GeV", 20, 250., 2750. );
      h_mass[elastic][i] = new TH1D( Form( "h_%s_mass_%d", kin_regions[elastic].c_str(), i ), "m_{#gamma#gamma}@@Events@@GeV", 16, 350., 1950. );
      h_mass[xicomp][i] = new TH1D( Form( "h_%s_mass_%d", kin_regions[xicomp].c_str(), i ), "m_{#gamma#gamma}@@Events@@GeV", 16, 350., 1950. );
      h_mass[xitight][i] = new TH1D( Form( "h_%s_mass_%d", kin_regions[xitight].c_str(), i ), "m_{#gamma#gamma}@@Events@@GeV", 8, 500., 2100. );
      //-------------------------------------------------------------------------------------------
      h_ptpair[nosel][i] = new TH1D( Form( "h_%s_ptpair_%d", kin_regions[nosel].c_str(), i ), "p_{T}^{#gamma#gamma}@@Events@@GeV", 25, 0., 500. );
      h_ptpair[qcd][i] = new TH1D( Form( "h_%s_ptpair_%d", kin_regions[qcd].c_str(), i ), "p_{T}^{#gamma#gamma}@@Events@@GeV", 11, 0., 440. );
      h_ptpair[presel][i] = new TH1D( Form( "h_%s_ptpair_%d", kin_regions[presel].c_str(), i ), "p_{T}^{#gamma#gamma}@@Events@@GeV", 25, 0., 500. );
      h_ptpair[presel_noneleveto][i] = (TH1D*)h_ptpair[presel][i]->Clone();
      h_ptpair[singlevtx][i] = (TH1D*)h_ptpair[presel][i]->Clone();
      h_ptpair[incl][i] = new TH1D( Form( "h_%s_ptpair_%d", kin_regions[incl].c_str(), i ), "p_{T}^{#gamma#gamma}@@Events@@GeV", 31, 0., 775. );
      h_ptpair[elastic][i] = new TH1D( Form( "h_%s_ptpair_%d", kin_regions[elastic].c_str(), i ), "p_{T}^{#gamma#gamma}@@Events@@GeV", 15, 0., 150. );
      h_ptpair[xicomp][i] = new TH1D( Form( "h_%s_ptpair_%d", kin_regions[xicomp].c_str(), i ), "p_{T}^{#gamma#gamma}@@Events@@GeV", 10, 0., 100. );
      h_ptpair[xitight][i] = new TH1D( Form( "h_%s_ptpair_%d", kin_regions[xitight].c_str(), i ), "p_{T}^{#gamma#gamma}@@Events@@GeV", 5, 0., 50. );
      //-------------------------------------------------------------------------------------------
      h_pt[nosel][i] = new TH1D( Form( "h_%s_pt_%d", kin_regions[nosel].c_str(), i ), "p_{T}^{#gamma}@@Events@@GeV", 40, 0., 600. );
      h_pt[qcd][i] = new TH1D( Form( "h_%s_pt_%d", kin_regions[qcd].c_str(), i ), "p_{T}^{#gamma}@@Events@@GeV", 24, 0., 200. );
      h_pt[presel][i] = new TH1D( Form( "h_%s_pt_%d", kin_regions[presel].c_str(), i ), "p_{T}^{#gamma}@@Events@@GeV", 40, 0., 600. );
      h_pt[presel_noneleveto][i] = (TH1D*)h_pt[presel][i]->Clone();
      h_pt[singlevtx][i] = (TH1D*)h_pt[presel][i]->Clone();
      h_pt[incl][i] = new TH1D( Form( "h_%s_pt_%d", kin_regions[incl].c_str(), i ), "p_{T}^{#gamma}@@Events@@GeV", 40, 150., 750. );
      h_pt[elastic][i] = new TH1D( Form( "h_%s_pt_%d", kin_regions[elastic].c_str(), i ), "p_{T}^{#gamma}@@Events@@GeV", 20, 0., 600. );
      h_pt[xicomp][i] = new TH1D( Form( "h_%s_pt_%d", kin_regions[xicomp].c_str(), i ), "p_{T}^{#gamma}@@Events@@GeV", 20, 0., 600. );
      h_pt[xitight][i] = new TH1D( Form( "h_%s_pt_%d", kin_regions[xitight].c_str(), i ), "p_{T}^{#gamma}@@Events@@GeV", 6, 50., 800. );
      //-------------------------------------------------------------------------------------------
      h_eta[nosel][i] = new TH1D( Form( "h_%s_eta_%d", kin_regions[nosel].c_str(), i ), "#eta^{#gamma}@@Events@@?g", 30, -3., 3. );
      h_eta[qcd][i] = new TH1D( Form( "h_%s_eta_%d", kin_regions[qcd].c_str(), i ), "#eta^{#gamma}@@Events@@?g", 20, -3., 3. );
      h_eta[presel][i] = new TH1D( Form( "h_%s_eta_%d", kin_regions[presel].c_str(), i ), "#eta^{#gamma}@@Events@@?g", 30, -3., 3. );
      h_eta[presel_noneleveto][i] = (TH1D*)h_eta[presel][i]->Clone();
      h_eta[singlevtx][i] = (TH1D*)h_eta[presel][i]->Clone();
      h_eta[incl][i] = new TH1D( Form( "h_%s_eta_%d", kin_regions[incl].c_str(), i ), "#eta^{#gamma}@@Events@@?g", 30, -3., 3. );
      h_eta[elastic][i] = new TH1D( Form( "h_%s_eta_%d", kin_regions[elastic].c_str(), i ), "#eta^{#gamma}@@Events@@?g", 30, -3., 3. );
      h_eta[xicomp][i] = new TH1D( Form( "h_%s_eta_%d", kin_regions[xicomp].c_str(), i ), "#eta^{#gamma}@@Events@@?g", 30, -3., 3. );
      h_eta[xitight][i] = new TH1D( Form( "h_%s_eta_%d", kin_regions[xitight].c_str(), i ), "#eta^{#gamma}@@Events@@?g", 30, -3., 3. );
      //-------------------------------------------------------------------------------------------
      h_r9[nosel][i] = new TH1D( Form( "h_%s_r9_%d", kin_regions[nosel].c_str(), i ), "R_{9}^{#gamma}@@Events@@?g", 50, 0.5, 1. );
      h_r9[qcd][i] = new TH1D( Form( "h_%s_r9_%d", kin_regions[qcd].c_str(), i ), "R_{9}^{#gamma}@@Events@@?g", 10, 0.5, 1. );
      h_r9[presel][i] = new TH1D( Form( "h_%s_r9_%d", kin_regions[presel].c_str(), i ), "R_{9}^{#gamma}@@Events@@?g", 30, 0.94, 1. );
      h_r9[presel_noneleveto][i] = (TH1D*)h_r9[presel][i]->Clone();
      h_r9[incl][i] = (TH1D*)h_r9[nosel][i]->Clone();
      h_r9[singlevtx][i] = (TH1D*)h_r9[nosel][i]->Clone();
      h_r9[elastic][i] = new TH1D( Form( "h_%s_r9_%d", kin_regions[elastic].c_str(), i ), "R_{9}^{#gamma}@@Events@@?g", 15, 0.94, 1. );
      h_r9[xicomp][i] = new TH1D( Form( "h_%s_r9_%d", kin_regions[xicomp].c_str(), i ), "R_{9}^{#gamma}@@Events@@?g", 6, 0.94, 1. );
      h_r9[xitight][i] = new TH1D( Form( "h_%s_r9_%d", kin_regions[xitight].c_str(), i ), "R_{9}^{#gamma}@@Events@@?g", 6, 0.94, 1. );
      //-------------------------------------------------------------------------------------------
      h_xip[nosel][i] = new TH1D( Form( "h_%s_xip_%d", kin_regions[nosel].c_str(), i ), "#xi_{#gamma#gamma}^{+}@@Events@@?g", 50, 0., 0.5 );
      h_xip[qcd][i] = new TH1D( Form( "h_%s_xip_%d", kin_regions[qcd].c_str(), i ), "#xi_{#gamma#gamma}^{+}@@Events@@?g", 25, 0., 0.5 );
      h_xip[presel][i] = new TH1D( Form( "h_%s_xip_%d", kin_regions[presel].c_str(), i ), "#xi_{#gamma#gamma}^{+}@@Events@@?g", 50, 0., 0.5 );
      h_xip[presel_noneleveto][i] = (TH1D*)h_xip[presel][i]->Clone();
      h_xip[singlevtx][i] = (TH1D*)h_xip[presel][i]->Clone();
      //h_xip[incl][i] = new TH1D( Form( "h_%s_xip_%d", kin_regions[incl].c_str(), i ), "#xi_{#gamma#gamma}^{+}@@Events@@?g", 25, 0., 0.5 );
      //////h_xip[incl][i] = new TH1D( Form( "h_%s_xip_%d", kin_regions[incl].c_str(), i ), "#xi_{#gamma#gamma}^{+}@@Events@@?g", 15, 0., 0.3 );
      h_xip[incl][i] = new TH1D( Form( "h_%s_xip_%d", kin_regions[incl].c_str(), i ), "#xi_{#gamma#gamma}^{+}@@Events@@?g", 25, 0., 0.2 );
      //h_xip[incl][i] = new TH1D( Form( "h_%s_xip_%d", kin_regions[incl].c_str(), i ), "#xi_{#gamma#gamma}^{+}@@Events@@?g", 30, 0., 0.3 );
      h_xip[elastic][i] = new TH1D( Form( "h_%s_xip_%d", kin_regions[elastic].c_str(), i ), "#xi_{#gamma#gamma}^{+}@@Events@@?g", 20, 0., 0.4 );
      h_xip[xicomp][i] = new TH1D( Form( "h_%s_xip_%d", kin_regions[xicomp].c_str(), i ), "#xi_{#gamma#gamma}^{+}@@Events@@?g", 13, 0.02, 0.15 );
      h_xip[xitight][i] = new TH1D( Form( "h_%s_xip_%d", kin_regions[xitight].c_str(), i ), "#xi_{#gamma#gamma}^{+}@@Events@@?g", 5, 0.05, 0.15 );
      //-------------------------------------------------------------------------------------------
      for ( const auto& reg : { nosel, qcd, presel, presel_noneleveto, incl, elastic, xicomp, xitight, singlevtx } ) {
        h_ptlead[reg][i] = (TH1D*)h_pt[reg][i]->Clone( Form( "h_%s_ptlead_%d", kin_regions[reg].c_str(), i ) );
        h_ptlead[reg][i]->SetTitle( "p_{T}^{#gamma} (leading #gamma)@@Events@@GeV" );
        h_ptsublead[reg][i] = (TH1D*)h_pt[reg][i]->Clone( Form( "h_%s_ptsublead_%d", kin_regions[reg].c_str(), i ) );
        h_ptsublead[reg][i]->SetTitle( "p_{T}^{#gamma} (subleading #gamma)@@Events@@GeV" );
        //-----------------------------------------------------------------------------------------
        h_etalead[reg][i] = (TH1D*)h_eta[reg][i]->Clone( Form( "h_%s_etalead_%d", kin_regions[reg].c_str(), i ) );
        h_etalead[reg][i]->SetTitle( "#eta^{#gamma} (leading #gamma)@@Events@@?g" );
        h_etasublead[reg][i] = (TH1D*)h_eta[reg][i]->Clone( Form( "h_%s_etasublead_%d", kin_regions[reg].c_str(), i ) );
        h_etasublead[reg][i]->SetTitle( "#eta^{#gamma} (subleading #gamma)@@Events@@?g" );
        //-----------------------------------------------------------------------------------------
        h_r9lead[reg][i] = (TH1D*)h_r9[reg][i]->Clone( Form( "h_%s_ptlead_%d", kin_regions[reg].c_str(), i ) );
        h_r9lead[reg][i]->SetTitle( "R_{9} (leading #gamma)@@Events@@?g" );
        h_r9sublead[reg][i] = (TH1D*)h_r9[reg][i]->Clone( Form( "h_%s_ptsublead_%d", kin_regions[reg].c_str(), i ) );
        h_r9sublead[reg][i]->SetTitle( "R_{9} (subleading #gamma)@@Events@@?g" );
        //-----------------------------------------------------------------------------------------
        h_xim[reg][i] = (TH1D*)h_xip[reg][i]->Clone( Form( "h_%s_xim_%d", kin_regions[reg].c_str(), i ) );
        h_xim[reg][i]->SetTitle( "#xi_{#gamma#gamma}^{-}@@Events@@?g" );
        //-----------------------------------------------------------------------------------------
        h_nvtx_unc[reg][i] = (TH1D*)h_nvtx[reg][i]->Clone( Form( "h_%s_nvtx_unc_%d", kin_regions[reg].c_str(), i ) );
        h_nvtx_unc[reg][i]->SetTitle( "Unreweighted number of primary vertices @@Events" );
      }
      //-------------------------------------------------------------------------------------------
      h_dphi[nosel][i] = new TH1D( Form( "h_%s_dphi_%d", kin_regions[nosel].c_str(), i ), "1-|#Delta#phi_{#gamma#gamma}/#pi|@@Events@@?.3f", 20, 0., 1. );
      //h_dphi[nosel][i] = new TH1D( Form( "h_%s_dphi_%d", kin_regions[nosel].c_str(), i ), "1-|#Delta#phi_{#gamma#gamma}/#pi|@@Events@@?.3f", 20, 0., 0.5 );
      h_dphi[qcd][i] = new TH1D( Form( "h_%s_dphi_%d", kin_regions[qcd].c_str(), i ), "1-|#Delta#phi_{#gamma#gamma}/#pi|@@Events@@?.3f", 20, 0., 1. );
      //h_dphi[presel][i] = new TH1D( Form( "h_%s_dphi_%d", kin_regions[presel].c_str(), i ), "1-|#Delta#phi_{#gamma#gamma}/#pi|@@Events@@?.3f", 40, 0., 1. );
      h_dphi[presel][i] = new TH1D( Form( "h_%s_dphi_%d", kin_regions[presel].c_str(), i ), "1-|#Delta#phi_{#gamma#gamma}/#pi|@@Events@@?.3f", 40, 0., 0.25 );
      h_dphi[presel_noneleveto][i] = (TH1D*)h_dphi[presel][i]->Clone();
      h_dphi[incl][i] = (TH1D*)h_dphi[presel][i]->Clone();
      h_dphi[singlevtx][i] = (TH1D*)h_dphi[presel][i]->Clone();
      h_dphi[elastic][i] = new TH1D( Form( "h_%s_dphi_%d", kin_regions[elastic].c_str(), i ), "1-|#Delta#phi_{#gamma#gamma}/#pi| (#times 10^{-3})@@Events@@?.4f", 10, 0., 5. );
      h_dphi[xicomp][i] = new TH1D( Form( "h_%s_dphi_%d", kin_regions[xicomp].c_str(), i ), "1-|#Delta#phi_{#gamma#gamma}/#pi| (#times 10^{-3})@@Events@@?.3f", 5, 0., 5. );
      h_dphi[xitight][i] = new TH1D( Form( "h_%s_dphi_%d", kin_regions[xitight].c_str(), i ), "1-|#Delta#phi_{#gamma#gamma}/#pi| (#times 10^{-3})@@Events@@?.3f", 5, 0., 5. );
      //-------------------------------------------------------------------------------------------
      h_met[nosel][i] = new TH1D( Form( "h_%s_met_%d", kin_regions[nosel].c_str(), i ), "p_{T}^{miss}@@Events@@GeV", 50, 0., 150. );
      h_met[qcd][i] = new TH1D( Form( "h_%s_met_%d", kin_regions[qcd].c_str(), i ), "p_{T}^{miss}@@Events@@GeV", 15, 0., 150. );
      h_met[presel][i] = new TH1D( Form( "h_%s_met_%d", kin_regions[presel].c_str(), i ), "p_{T}^{miss}@@Events@@GeV", 50, 0., 150. );
      h_met[presel_noneleveto][i] = (TH1D*)h_met[presel][i]->Clone();
      h_met[incl][i] = (TH1D*)h_met[presel][i]->Clone();
      h_met[singlevtx][i] = (TH1D*)h_met[presel][i]->Clone();
      h_met[elastic][i] = new TH1D( Form( "h_%s_met_%d", kin_regions[elastic].c_str(), i ), "p_{T}^{miss}@@Events@@GeV", 25, 0., 150. );
      h_met[xicomp][i] = new TH1D( Form( "h_%s_met_%d", kin_regions[xicomp].c_str(), i ), "p_{T}^{miss}@@Events@@GeV", 25, 0., 150. );
      h_met[xitight][i] = new TH1D( Form( "h_%s_met_%d", kin_regions[xitight].c_str(), i ), "p_{T}^{miss}@@Events@@GeV", 10, 0., 100. );
      //-------------------------------------------------------------------------------------------
      h_diph_vtxz[nosel][i] = new TH1D( Form( "h_%s_diph_vtxz_%d", kin_regions[nosel].c_str(), i ), "Diphoton vertex z@@Events@@cm", 40, -20., 20. );
      h_diph_vtxz[qcd][i] = new TH1D( Form( "h_%s_diph_vtxz_%d", kin_regions[qcd].c_str(), i ), "Diphoton vertex z@@Events@@cm", 40, -20., 20. );
      h_diph_vtxz[presel][i] = new TH1D( Form( "h_%s_diph_vtxz_%d", kin_regions[presel].c_str(), i ), "Diphoton vertex z@@Events@@cm", 40, -20., 20. );
      h_diph_vtxz[presel_noneleveto][i] = (TH1D*)h_diph_vtxz[presel][i]->Clone();
      h_diph_vtxz[incl][i] = (TH1D*)h_diph_vtxz[presel][i]->Clone();
      h_diph_vtxz[singlevtx][i] = (TH1D*)h_diph_vtxz[presel][i]->Clone();
      h_diph_vtxz[elastic][i] = new TH1D( Form( "h_%s_diph_vtxz_%d", kin_regions[elastic].c_str(), i ), "Diphoton vertex z@@Events@@cm", 20, -20., 20. );
      h_diph_vtxz[xicomp][i] = new TH1D( Form( "h_%s_diph_vtxz_%d", kin_regions[xicomp].c_str(), i ), "Diphoton vertex z@@Events@@cm", 20, -20., 20. );
      h_diph_vtxz[xitight][i] = new TH1D( Form( "h_%s_diph_vtxz_%d", kin_regions[xitight].c_str(), i ), "Diphoton vertex z@@Events@@cm", 20, -20., 20. );
    }

    auto tree = s.tree();
    TreeEvent ev;
    ev.attach( tree, s.type() == Sample::kData );
    const unsigned long long num_entries = tree->GetEntriesFast();

    string sname = s.name();
    rescaled[i] = false;
    if ( sample_type == mc_signal ) {
      if ( scaling_signal > 1. && sample_weight < 1.e-5 ) {
      //if ( scaling_signal > 1. ) {
        sname += Form( " (#times%.0f)", scaling_signal );
        sample_weight *= scaling_signal;
        rescaled[i] = true;
      }
    }
    cout << "(" << (int)sample_type << ")" << s.name() << ": " << sample_weight << endl;
    const double unrescaled_sample_weight = sample_weight;
    if ( s.type() == Sample::kBackground && s.name().find( "Inclusive #gamma#gamma" ) != string::npos )
      if ( !rescale_background_parametric )
        sample_weight *= rescaling_background;
    // loop on events
    for ( unsigned long long j = 0; j < num_entries; ++j ) {
      if ( fmod( j*1./num_entries, 0.1 ) == 0 ) cerr << "   event " << j << " / " << num_entries << endl;

      tree->GetEntry( j );

      if ( ev.hlt_accept[0] == 0 )
        continue;

      bool has_diphoton_cand = false, has_elastic_diphoton_cand = false;
      bool has_elastic_xicomp_cand = false, has_elastic_xitight_cand = false;
      bool has_incl_diphoton_event = false;

      //if ( s.type() == Sample::kData && !ev_selector.isSelected( ev.run_id, ev.lumisection, ev.event_number ) ) continue;

      //const float s_weight = ( s.type() == Sample::kData ) ? 1. : ( sample_weight * ev.pileup_weight );
      const double sp_weight = ev.pileup_weight;

      unsigned short num_preselected_diphotons = 0;

      // loop on diphotons
      for ( unsigned short k = 0; k < ev.num_diphoton; ++k ) {

        // EB: 0 < |eta| < 1.4442
        // EE: |eta| > 1.566
        unsigned short ev_class = TreeEvent::invalid;
        if ( fabs( ev.diphoton_eta1[k] ) <= eta_cut && fabs( ev.diphoton_eta2[k] ) <= eta_cut ) {
          if ( fabs( ev.diphoton_eta1[k] ) < min_etaveto && fabs( ev.diphoton_eta2[k] ) > max_etaveto ) ev_class = TreeEvent::ebee;
          else if ( fabs( ev.diphoton_eta2[k] ) < min_etaveto && fabs( ev.diphoton_eta1[k] ) > max_etaveto ) ev_class = TreeEvent::ebee;
          else if ( fabs( ev.diphoton_eta1[k] ) < min_etaveto && fabs( ev.diphoton_eta2[k] ) < min_etaveto ) ev_class = TreeEvent::ebeb;
          else if ( fabs( ev.diphoton_eta1[k] ) > max_etaveto && fabs( ev.diphoton_eta2[k] ) > max_etaveto ) ev_class = TreeEvent::eeee;
        }

        //----- only keep EBEE and EBEB diphoton events

        if ( ev_class == TreeEvent::invalid ) continue;
        if ( ev_class == TreeEvent::eeee ) continue; //FIXME FIXME
        //if ( ev_class == TreeEvent::ebee ) continue; //FIXME FIXME
        //cout << ">> evt class: " << TreeEvent::classes[ev_class] << "  " << ev.diphoton_eta1[k] << "\t" << ev.diphoton_eta2[k] << endl;

//        if ( ev.diphoton_r91[k] < 0.8 || ev.diphoton_r92[k] < 0.8 ) continue; //FIXME for safety (avoids a bug in the TTree producer...)

        const float acop = 1.-fabs( ev.diphoton_dphi[k]/M_PI );
        const float lacop = log10( acop );
        const float xip = ( ev.diphoton_pt1[k]*exp( +ev.diphoton_eta1[k] )+ev.diphoton_pt2[k]*exp( +ev.diphoton_eta2[k] ) ) / sqrt_s;
        const float xim = ( ev.diphoton_pt1[k]*exp( -ev.diphoton_eta1[k] )+ev.diphoton_pt2[k]*exp( -ev.diphoton_eta2[k] ) ) / sqrt_s;

        const TVector3 diph_vtx( ev.vertex_x[ev.diphoton_vertex_id[k]], ev.vertex_y[ev.diphoton_vertex_id[k]], ev.vertex_z[ev.diphoton_vertex_id[k]] );
        const TVector3 sc_pho1( ev.diphoton_supercluster_x1[k], ev.diphoton_supercluster_y1[k], ev.diphoton_supercluster_z1[k] );
        const TVector3 sc_pho2( ev.diphoton_supercluster_x2[k], ev.diphoton_supercluster_y2[k], ev.diphoton_supercluster_z2[k] );

        float s_weight = 1.;
        if ( s.type() != Sample::kData ) {
          const float eff_pho1 = loosemvaeff( ev.diphoton_r91[k], sc_pho1.Eta() );
          const float eff_pho2 = loosemvaeff( ev.diphoton_r92[k], sc_pho2.Eta() );
          s_weight = sp_weight * ( eff_pho1*eff_pho2 );
        }
        if ( rescale_background_parametric && s.type() == Sample::kBackground && s.name().find( "Inclusive #gamma#gamma" ) != string::npos )
          s_weight *= incl_fit_func->Eval( lacop );

        //----- look at surrounding objects

        const float min_lep_vtx_dist = 2.0; // in cm
        float closest_lep_vtx_dist = 999.; // in cm

        unsigned short num_matched_ele = 0, num_matched_mu = 0, num_matched_jets = 0;
        for ( unsigned short l = 0; l < ev.num_electron; ++l ) {
          TLorentzVector ele; ele.SetPtEtaPhiE( ev.electron_pt[l], ev.electron_eta[l], ev.electron_phi[l], ev.electron_energy[l] );
          const TVector3 ele_vtx( ev.electron_vtx_x[l], ev.electron_vtx_y[l], ev.electron_vtx_z[l] );
          const float ele_dist = ( ele_vtx-diph_vtx ).Mag();
          if ( ele_dist < closest_lep_vtx_dist ) closest_lep_vtx_dist = ele_dist;
          if ( ele_dist < min_lep_vtx_dist ) num_matched_ele++;
        }
        for ( unsigned short l = 0; l < ev.num_muon; ++l ) {
          TLorentzVector mu; mu.SetPtEtaPhiE( ev.muon_pt[l], ev.muon_eta[l], ev.muon_phi[l], ev.muon_energy[l] );
          const TVector3 mu_vtx( ev.muon_vtx_x[l], ev.muon_vtx_y[l], ev.muon_vtx_z[l] );
          const float mu_dist = ( mu_vtx-diph_vtx ).Mag();
          if ( mu_dist < closest_lep_vtx_dist ) closest_lep_vtx_dist = mu_dist;
          if ( mu_dist < min_lep_vtx_dist ) num_matched_mu++;
        }
        for ( unsigned short l = 0; l < ev.num_jet; ++l ) {
          if ( ev.jet_dipho_match[l] != k ) continue;
          if ( ev.jet_pt[l] < 500. ) continue;
          TLorentzVector jet; jet.SetPtEtaPhiE( ev.jet_pt[l], ev.jet_eta[l], ev.jet_phi[l], ev.jet_energy[l] );
          num_matched_jets++;
        }

        //--- no selection
        h_ptlead[nosel][i]->Fill( ev.diphoton_pt1[k], s_weight );
        h_ptsublead[nosel][i]->Fill( ev.diphoton_pt2[k], s_weight );
        h_pt[nosel][i]->Fill( ev.diphoton_pt1[k], s_weight );
        h_pt[nosel][i]->Fill( ev.diphoton_pt2[k], s_weight );
        h_pt_varbins[nosel][i]->Fill( ev.diphoton_pt1[k], s_weight );
        h_pt_varbins[nosel][i]->Fill( ev.diphoton_pt2[k], s_weight );
        h_etalead[nosel][i]->Fill( ev.diphoton_eta1[k], s_weight );
        h_etasublead[nosel][i]->Fill( ev.diphoton_eta2[k], s_weight );
        h_eta[nosel][i]->Fill( ev.diphoton_eta1[k], s_weight );
        h_eta[nosel][i]->Fill( ev.diphoton_eta2[k], s_weight );
        h_r9lead[nosel][i]->Fill( ev.diphoton_r91[k], s_weight );
        h_r9sublead[nosel][i]->Fill( ev.diphoton_r92[k], s_weight );
        h_r9[nosel][i]->Fill( ev.diphoton_r91[k], s_weight );
        h_r9[nosel][i]->Fill( ev.diphoton_r92[k], s_weight );
        h_mass[nosel][i]->Fill( ev.diphoton_mass[k], s_weight );
        h_mass_rebin[nosel][i]->Fill( ev.diphoton_mass[k], s_weight );
        h_ptpair[nosel][i]->Fill( ev.diphoton_pt[k], s_weight );
        h_dphi[nosel][i]->Fill( acop, s_weight );
        h_dphi_varbins[nosel][i]->Fill( acop, s_weight );
        h_dphi_logbins[nosel][i]->Fill( lacop, s_weight );
        h_dphi_tightbins[nosel][i]->Fill( acop*1.e3, s_weight );
        h_met[nosel][i]->Fill( ev.met, s_weight );
        h_xip[nosel][i]->Fill( xip, s_weight );
        h_xim[nosel][i]->Fill( xim, s_weight );
        h_diph_vtxz[nosel][i]->Fill( diph_vtx.z(), s_weight );
        h_diph_numjets[nosel][i]->Fill( num_matched_jets, s_weight );
        h_diph_numleptons[nosel][i]->Fill( num_matched_ele+num_matched_mu, s_weight );
        if ( ev.num_vertex == 1 ) {
          h_ptlead[singlevtx][i]->Fill( ev.diphoton_pt1[k], s_weight );
          h_ptsublead[singlevtx][i]->Fill( ev.diphoton_pt2[k], s_weight );
          h_pt[singlevtx][i]->Fill( ev.diphoton_pt1[k], s_weight );
          h_pt[singlevtx][i]->Fill( ev.diphoton_pt2[k], s_weight );
          h_pt_varbins[singlevtx][i]->Fill( ev.diphoton_pt1[k], s_weight );
          h_pt_varbins[singlevtx][i]->Fill( ev.diphoton_pt2[k], s_weight );
          h_etalead[singlevtx][i]->Fill( ev.diphoton_eta1[k], s_weight );
          h_etasublead[singlevtx][i]->Fill( ev.diphoton_eta2[k], s_weight );
          h_eta[singlevtx][i]->Fill( ev.diphoton_eta1[k], s_weight );
          h_eta[singlevtx][i]->Fill( ev.diphoton_eta2[k], s_weight );
          h_r9lead[singlevtx][i]->Fill( ev.diphoton_r91[k], s_weight );
          h_r9sublead[singlevtx][i]->Fill( ev.diphoton_r92[k], s_weight );
          h_r9[singlevtx][i]->Fill( ev.diphoton_r91[k], s_weight );
          h_r9[singlevtx][i]->Fill( ev.diphoton_r92[k], s_weight );
          h_mass[singlevtx][i]->Fill( ev.diphoton_mass[k], s_weight );
          h_mass_rebin[singlevtx][i]->Fill( ev.diphoton_mass[k], s_weight );
          h_ptpair[singlevtx][i]->Fill( ev.diphoton_pt[k], s_weight );
          h_dphi[singlevtx][i]->Fill( acop, s_weight );
          h_dphi_varbins[singlevtx][i]->Fill( acop, s_weight );
          h_dphi_logbins[singlevtx][i]->Fill( lacop, s_weight );
          h_dphi_tightbins[singlevtx][i]->Fill( acop*1.e3, s_weight );
          h_met[singlevtx][i]->Fill( ev.met, s_weight );
          h_xip[singlevtx][i]->Fill( xip, s_weight );
          h_xim[singlevtx][i]->Fill( xim, s_weight );
          h_diph_vtxz[singlevtx][i]->Fill( diph_vtx.z(), s_weight );
          h_diph_numjets[singlevtx][i]->Fill( num_matched_jets, s_weight );
          h_diph_numleptons[singlevtx][i]->Fill( num_matched_ele+num_matched_mu, s_weight );
        }
//        h2_eleveto[i]->Fill( ev.diphoton_ele_veto1[k], ev.diphoton_ele_veto2[k] );

        evts_sel[i][nosel] += s_weight * sample_weight;

        //--- QCD control region

        bool is_preselected = true;
        if ( fabs( ev.diphoton_eta1[k] ) > min_etaveto && fabs( ev.diphoton_eta1[k] ) < max_etaveto ) is_preselected = false;
        if ( fabs( ev.diphoton_eta2[k] ) > min_etaveto && fabs( ev.diphoton_eta2[k] ) < max_etaveto ) is_preselected = false;
        is_preselected &= ( fabs( ev.diphoton_eta1[k] ) <= eta_cut );
        is_preselected &= ( fabs( ev.diphoton_eta2[k] ) <= eta_cut );

        //if ( ( ev.diphoton_r91[k] < 0.85 || ev.diphoton_r92[k] < 0.85 ) && ev.diphoton_mass[k] < 800. ) {
        if ( is_preselected && ( ev.diphoton_r91[k] < 0.85 || ev.diphoton_r92[k] < 0.85 ) && ev.diphoton_pt1[k] < 150. && ev.diphoton_pt2[k] < 150. ) {
        //if ( is_preselected && ( ev.diphoton_r91[k] < 0.85 && ev.diphoton_r92[k] < 0.85 ) && ev.diphoton_pt1[k] < 150. && ev.diphoton_pt2[k] < 150. ) {
        //if ( ( ev.diphoton_r91[k] < 0.85 || ev.diphoton_r92[k] < 0.85 ) && ev.diphoton_pt1[k] < 150. && ev.diphoton_pt2[k] < 150. && ev.diphoton_pt[k] < 50. ) {
        //if ( is_preselected && ( ev.diphoton_r91[k] < 0.85 || ev.diphoton_r92[k] < 0.85 ) && ev.diphoton_pt[k] < 50. ) {
          h_ptlead[qcd][i]->Fill( ev.diphoton_pt1[k], s_weight );
          h_ptsublead[qcd][i]->Fill( ev.diphoton_pt2[k], s_weight );
          h_pt[qcd][i]->Fill( ev.diphoton_pt1[k], s_weight );
          h_pt[qcd][i]->Fill( ev.diphoton_pt2[k], s_weight );
          h_pt_varbins[qcd][i]->Fill( ev.diphoton_pt1[k], s_weight );
          h_pt_varbins[qcd][i]->Fill( ev.diphoton_pt2[k], s_weight );
          h_etalead[qcd][i]->Fill( ev.diphoton_eta1[k], s_weight );
          h_etasublead[qcd][i]->Fill( ev.diphoton_eta2[k], s_weight );
          h_eta[qcd][i]->Fill( ev.diphoton_eta1[k], s_weight );
          h_eta[qcd][i]->Fill( ev.diphoton_eta2[k], s_weight );
          h_r9lead[qcd][i]->Fill( ev.diphoton_r91[k], s_weight );
          h_r9sublead[qcd][i]->Fill( ev.diphoton_r92[k], s_weight );
          h_r9[qcd][i]->Fill( ev.diphoton_r91[k], s_weight );
          h_r9[qcd][i]->Fill( ev.diphoton_r92[k], s_weight );
          h_mass[qcd][i]->Fill( ev.diphoton_mass[k], s_weight );
          h_mass_rebin[qcd][i]->Fill( ev.diphoton_mass[k], s_weight );
          h_ptpair[qcd][i]->Fill( ev.diphoton_pt[k], s_weight );
          h_dphi[qcd][i]->Fill( acop, s_weight );
          h_dphi_varbins[qcd][i]->Fill( acop, s_weight );
          h_dphi_logbins[qcd][i]->Fill( lacop, s_weight );
          h_dphi_tightbins[qcd][i]->Fill( acop*1.e3, s_weight );
          h_met[qcd][i]->Fill( ev.met, s_weight );
          h_xip[qcd][i]->Fill( xip, s_weight );
          h_xim[qcd][i]->Fill( xim, s_weight );
          h_diph_vtxz[qcd][i]->Fill( diph_vtx.z(), s_weight );
          h_diph_numjets[qcd][i]->Fill( num_matched_jets, s_weight );
          h_diph_numleptons[qcd][i]->Fill( num_matched_ele+num_matched_mu, s_weight );
          evts_sel[i][qcd] += s_weight * sample_weight;
        }

        //--- preselection definition

//        is_preselected &= ( ev.diphoton_id1[k] >= -0.2 && ev.diphoton_id2[k] >= -0.2 );

        //if ( is_preselected && ev.diphoton_pt1[k] >= 300. && ev.diphoton_pt2[k] >= 300. ) {
        //if ( is_preselected && ev.diphoton_pt1[k] >= 250. && ev.diphoton_pt2[k] >= 250. ) {
        //if ( is_preselected && ev.diphoton_pt1[k] >= 200. && ev.diphoton_pt2[k] >= 200. ) {
        if ( is_preselected && ev.diphoton_pt1[k] >= 200. && ev.diphoton_pt2[k] >= 200. ) {
          //FIXME h_dphi[incl][i]->Fill( acop, s_weight );
          h_dphi_varbins[incl][i]->Fill( acop, s_weight );
          h_dphi_logbins[incl][i]->Fill( lacop, s_weight );
          h_dphi_tightbins[incl][i]->Fill( acop*1.e3, s_weight );
        }
        if ( is_preselected && ev.diphoton_pt1[k] >= 200. && ev.diphoton_pt2[k] >= 200. && acop > 0.025 ) {
        //if ( is_preselected && ev.diphoton_pt1[k] >= 200. && ev.diphoton_pt2[k] >= 200. && acop > 0.005 ) {
          h_ptlead[incl][i]->Fill( ev.diphoton_pt1[k], s_weight );
          h_ptsublead[incl][i]->Fill( ev.diphoton_pt2[k], s_weight );
          h_pt[incl][i]->Fill( ev.diphoton_pt1[k], s_weight );
          h_pt[incl][i]->Fill( ev.diphoton_pt2[k], s_weight );
          h_pt_varbins[incl][i]->Fill( ev.diphoton_pt1[k], s_weight );
          h_pt_varbins[incl][i]->Fill( ev.diphoton_pt2[k], s_weight );
          h_etalead[incl][i]->Fill( ev.diphoton_eta1[k], s_weight );
          h_etasublead[incl][i]->Fill( ev.diphoton_eta2[k], s_weight );
          h_eta[incl][i]->Fill( ev.diphoton_eta1[k], s_weight );
          h_eta[incl][i]->Fill( ev.diphoton_eta2[k], s_weight );
          h_r9lead[incl][i]->Fill( ev.diphoton_r91[k], s_weight );
          h_r9sublead[incl][i]->Fill( ev.diphoton_r92[k], s_weight );
          h_r9[incl][i]->Fill( ev.diphoton_r91[k], s_weight );
          h_r9[incl][i]->Fill( ev.diphoton_r92[k], s_weight );
          h_mass_rebin[incl][i]->Fill( ev.diphoton_mass[k], s_weight );
          h_ptpair[incl][i]->Fill( ev.diphoton_pt[k], s_weight );
          h_dphi[incl][i]->Fill( acop, s_weight );
          h_met[incl][i]->Fill( ev.met, s_weight );
          h_xip[incl][i]->Fill( xip, s_weight );
          h_xim[incl][i]->Fill( xim, s_weight );
          h_diph_vtxz[incl][i]->Fill( diph_vtx.z(), s_weight );
          h_diph_numjets[incl][i]->Fill( num_matched_jets, s_weight );
          h_diph_numleptons[incl][i]->Fill( num_matched_ele+num_matched_mu, s_weight );
          evts_sel[i][incl] += s_weight * sample_weight;

          has_incl_diphoton_event = true;
        }

        is_preselected &= ( ev.diphoton_pt1[k] >= pt_cut && ev.diphoton_pt2[k] >= pt_cut );
        is_preselected &= ( ev.diphoton_r91[k] >= r9_cut );
        is_preselected &= ( ev.diphoton_r92[k] >= r9_cut );
        is_preselected &= ( ev.diphoton_mass[k] >= mass_cut );

        //--- preselection

        if ( !is_preselected ) continue;

        h_ptlead[presel_noneleveto][i]->Fill( ev.diphoton_pt1[k], s_weight );
        h_ptsublead[presel_noneleveto][i]->Fill( ev.diphoton_pt2[k], s_weight );
        h_pt[presel_noneleveto][i]->Fill( ev.diphoton_pt1[k], s_weight );
        h_pt[presel_noneleveto][i]->Fill( ev.diphoton_pt2[k], s_weight );
        h_pt_varbins[presel_noneleveto][i]->Fill( ev.diphoton_pt1[k], s_weight );
        h_pt_varbins[presel_noneleveto][i]->Fill( ev.diphoton_pt2[k], s_weight );
        h_etalead[presel_noneleveto][i]->Fill( ev.diphoton_eta1[k], s_weight );
        h_etasublead[presel_noneleveto][i]->Fill( ev.diphoton_eta2[k], s_weight );
        h_eta[presel_noneleveto][i]->Fill( ev.diphoton_eta1[k], s_weight );
        h_eta[presel_noneleveto][i]->Fill( ev.diphoton_eta2[k], s_weight );
        h_r9lead[presel_noneleveto][i]->Fill( ev.diphoton_r91[k], s_weight );
        h_r9sublead[presel_noneleveto][i]->Fill( ev.diphoton_r92[k], s_weight );
        h_r9[presel_noneleveto][i]->Fill( ev.diphoton_r91[k], s_weight );
        h_r9[presel_noneleveto][i]->Fill( ev.diphoton_r92[k], s_weight );
        h_mass[presel_noneleveto][i]->Fill( ev.diphoton_mass[k], s_weight );
        h_mass_rebin[presel_noneleveto][i]->Fill( ev.diphoton_mass[k], s_weight );
        h_ptpair[presel_noneleveto][i]->Fill( ev.diphoton_pt[k], s_weight );
        h_dphi[presel_noneleveto][i]->Fill( acop, s_weight );
        h_dphi_varbins[presel_noneleveto][i]->Fill( acop, s_weight );
        h_dphi_logbins[presel_noneleveto][i]->Fill( lacop, s_weight );
        h_dphi_tightbins[presel_noneleveto][i]->Fill( acop*1.e3, s_weight );
        h_met[presel_noneleveto][i]->Fill( ev.met, s_weight );
        h_xip[presel_noneleveto][i]->Fill( xip, s_weight );
        h_xim[presel_noneleveto][i]->Fill( xim, s_weight );
        h_diph_vtxz[presel_noneleveto][i]->Fill( diph_vtx.z(), s_weight );
        h_diph_numjets[presel_noneleveto][i]->Fill( num_matched_jets, s_weight );
        h_diph_numleptons[presel_noneleveto][i]->Fill( num_matched_ele+num_matched_mu, s_weight );
        evts_sel[i][presel_noneleveto] += s_weight * sample_weight;

        //if ( ev.diphoton_id1[k] > -.2 && ev.diphoton_id2[k] > -.2 ) //FIXME FIXME FIXME
        h2_eleveto[i]->Fill( ev.diphoton_ele_veto1[k], ev.diphoton_ele_veto2[k] );
        if ( ev.diphoton_ele_veto1[k] != 1 ) continue;
        if ( ev.diphoton_ele_veto2[k] != 1 ) continue;

        num_preselected_diphotons++;

        h_ptlead[presel][i]->Fill( ev.diphoton_pt1[k], s_weight );
        h_ptsublead[presel][i]->Fill( ev.diphoton_pt2[k], s_weight );
        h_pt[presel][i]->Fill( ev.diphoton_pt1[k], s_weight );
        h_pt[presel][i]->Fill( ev.diphoton_pt2[k], s_weight );
        h_pt_varbins[presel][i]->Fill( ev.diphoton_pt1[k], s_weight );
        h_pt_varbins[presel][i]->Fill( ev.diphoton_pt2[k], s_weight );
        h_etalead[presel][i]->Fill( ev.diphoton_eta1[k], s_weight );
        h_etasublead[presel][i]->Fill( ev.diphoton_eta2[k], s_weight );
        h_eta[presel][i]->Fill( ev.diphoton_eta1[k], s_weight );
        h_eta[presel][i]->Fill( ev.diphoton_eta2[k], s_weight );
        h_r9lead[presel][i]->Fill( ev.diphoton_r91[k], s_weight );
        h_r9sublead[presel][i]->Fill( ev.diphoton_r92[k], s_weight );
        h_r9[presel][i]->Fill( ev.diphoton_r91[k], s_weight );
        h_r9[presel][i]->Fill( ev.diphoton_r92[k], s_weight );
        h_mass[presel][i]->Fill( ev.diphoton_mass[k], s_weight );
        h_mass_rebin[presel][i]->Fill( ev.diphoton_mass[k], s_weight );
        h_ptpair[presel][i]->Fill( ev.diphoton_pt[k], s_weight );
        h_dphi[presel][i]->Fill( acop, s_weight );
        h_dphi_varbins[presel][i]->Fill( acop, s_weight );
        h_dphi_logbins[presel][i]->Fill( lacop, s_weight );
        h_dphi_tightbins[presel][i]->Fill( acop*1.e3, s_weight );
        h_met[presel][i]->Fill( ev.met, s_weight );
        h_xip[presel][i]->Fill( xip, s_weight );
        h_xim[presel][i]->Fill( xim, s_weight );
        h_diph_vtxz[presel][i]->Fill( diph_vtx.z(), s_weight );
        h_diph_numjets[presel][i]->Fill( num_matched_jets, s_weight );
        h_diph_numleptons[presel][i]->Fill( num_matched_ele+num_matched_mu, s_weight );
        evts_sel[i][presel] += s_weight * sample_weight;

        //----- from that point on, we have a diphoton candidate

        if ( s.name().find( "QCD" ) != string::npos ) continue; //FIXME skip QCD from now on, as MC weight is incredible
        has_diphoton_cand = true;

        if ( acop < 0.005 ) {
          //--- elastic region
          h_ptlead[elastic][i]->Fill( ev.diphoton_pt1[k], s_weight );
          h_ptsublead[elastic][i]->Fill( ev.diphoton_pt2[k], s_weight );
          h_pt[elastic][i]->Fill( ev.diphoton_pt1[k], s_weight );
          h_pt[elastic][i]->Fill( ev.diphoton_pt2[k], s_weight );
          h_pt_varbins[elastic][i]->Fill( ev.diphoton_pt1[k], s_weight );
          h_pt_varbins[elastic][i]->Fill( ev.diphoton_pt2[k], s_weight );
          h_etalead[elastic][i]->Fill( ev.diphoton_eta1[k], s_weight );
          h_etasublead[elastic][i]->Fill( ev.diphoton_eta2[k], s_weight );
          h_eta[elastic][i]->Fill( ev.diphoton_eta1[k], s_weight );
          h_eta[elastic][i]->Fill( ev.diphoton_eta2[k], s_weight );
          h_r9lead[elastic][i]->Fill( ev.diphoton_r91[k], s_weight );
          h_r9sublead[elastic][i]->Fill( ev.diphoton_r92[k], s_weight );
          h_r9[elastic][i]->Fill( ev.diphoton_r91[k], s_weight );
          h_r9[elastic][i]->Fill( ev.diphoton_r92[k], s_weight );
          h_mass[elastic][i]->Fill( ev.diphoton_mass[k], s_weight );
          h_mass_rebin[elastic][i]->Fill( ev.diphoton_mass[k], s_weight );
          h_ptpair[elastic][i]->Fill( ev.diphoton_pt[k], s_weight );
          h_dphi[elastic][i]->Fill( acop*1.e3, s_weight );
          h_dphi_varbins[elastic][i]->Fill( acop, s_weight );
          h_dphi_logbins[elastic][i]->Fill( lacop, s_weight );
          h_dphi_tightbins[elastic][i]->Fill( acop*1.e3, s_weight );
          h_met[elastic][i]->Fill( ev.met, s_weight );
          h_xip[elastic][i]->Fill( xip, s_weight );
          h_xim[elastic][i]->Fill( xim, s_weight );
          h_diph_vtxz[elastic][i]->Fill( diph_vtx.z(), s_weight );
          h_diph_numjets[elastic][i]->Fill( num_matched_jets, s_weight );
          h_diph_numleptons[elastic][i]->Fill( num_matched_ele+num_matched_mu, s_weight );
          //h2_eleveto[i]->Fill( ev.diphoton_ele_veto1[k], ev.diphoton_ele_veto2[k] );
          evts_sel[i][elastic] += s_weight * sample_weight;

          has_elastic_diphoton_cand = true;

          //--- count the vertices around diphoton vertex
          if ( s.type() == Sample::kData ) {
            h_diph_nvtxaround[0]->Fill( ev.diphoton_vertex_vtx1mmdist[k] );
            h_diph_nvtxaround[1]->Fill( ev.diphoton_vertex_vtx2mmdist[k] );
            h_diph_nvtxaround[2]->Fill( ev.diphoton_vertex_vtx5mmdist[k] );
            h_diph_nvtxaround[3]->Fill( ev.diphoton_vertex_vtx1cmdist[k] );
          }

          if ( ( ( xim > pots_accept["56N"] || xim > pots_accept["56F"] ) && xim < 0.15 )
            && ( ( xip > pots_accept["45N"] || xip > pots_accept["45F"] ) && xip < 0.15 ) ) {
            // xi(gg) within the PPS pots acceptance
            h_ptlead[xicomp][i]->Fill( ev.diphoton_pt1[k], s_weight );
            h_ptsublead[xicomp][i]->Fill( ev.diphoton_pt2[k], s_weight );
            h_pt[xicomp][i]->Fill( ev.diphoton_pt1[k], s_weight );
            h_pt[xicomp][i]->Fill( ev.diphoton_pt2[k], s_weight );
            h_pt_varbins[xicomp][i]->Fill( ev.diphoton_pt1[k], s_weight );
            h_pt_varbins[xicomp][i]->Fill( ev.diphoton_pt2[k], s_weight );
            h_etalead[xicomp][i]->Fill( ev.diphoton_eta1[k], s_weight );
            h_etasublead[xicomp][i]->Fill( ev.diphoton_eta2[k], s_weight );
            h_eta[xicomp][i]->Fill( ev.diphoton_eta1[k], s_weight );
            h_eta[xicomp][i]->Fill( ev.diphoton_eta2[k], s_weight );
            h_r9lead[xicomp][i]->Fill( ev.diphoton_r91[k], s_weight );
            h_r9sublead[xicomp][i]->Fill( ev.diphoton_r92[k], s_weight );
            h_r9[xicomp][i]->Fill( ev.diphoton_r91[k], s_weight );
            h_r9[xicomp][i]->Fill( ev.diphoton_r92[k], s_weight );
            h_mass[xicomp][i]->Fill( ev.diphoton_mass[k], s_weight );
            h_mass_rebin[xicomp][i]->Fill( ev.diphoton_mass[k], s_weight );
            h_ptpair[xicomp][i]->Fill( ev.diphoton_pt[k], s_weight );
            h_dphi[xicomp][i]->Fill( acop*1.e3, s_weight );
            h_dphi_varbins[xicomp][i]->Fill( acop, s_weight );
            h_dphi_logbins[xicomp][i]->Fill( lacop, s_weight );
            h_dphi_tightbins[xicomp][i]->Fill( acop*1.e3, s_weight );
            h_met[xicomp][i]->Fill( ev.met, s_weight );
            h_xip[xicomp][i]->Fill( xip, s_weight );
            h_xim[xicomp][i]->Fill( xim, s_weight );
            h_diph_vtxz[xicomp][i]->Fill( diph_vtx.z(), s_weight );
            h_diph_numjets[xicomp][i]->Fill( num_matched_jets, s_weight );
            h_diph_numleptons[xicomp][i]->Fill( num_matched_ele+num_matched_mu, s_weight );
            evts_sel[i][xicomp] += s_weight * sample_weight;

            has_elastic_xicomp_cand = true;
            if ( ( ( xim > pots_accept_tight["56N"] || xim > pots_accept_tight["56F"] ) && xim < 0.15 )
              && ( ( xip > pots_accept_tight["45N"] || xip > pots_accept_tight["45F"] ) && xip < 0.15 ) ) {
              // xi(gg) within tight PPS pots acceptance
              h_ptlead[xitight][i]->Fill( ev.diphoton_pt1[k], s_weight );
              h_ptsublead[xitight][i]->Fill( ev.diphoton_pt2[k], s_weight );
              h_pt[xitight][i]->Fill( ev.diphoton_pt1[k], s_weight );
              h_pt[xitight][i]->Fill( ev.diphoton_pt2[k], s_weight );
              h_pt_varbins[xitight][i]->Fill( ev.diphoton_pt1[k], s_weight );
              h_pt_varbins[xitight][i]->Fill( ev.diphoton_pt2[k], s_weight );
              h_etalead[xitight][i]->Fill( ev.diphoton_eta1[k], s_weight );
              h_etasublead[xitight][i]->Fill( ev.diphoton_eta2[k], s_weight );
              h_eta[xitight][i]->Fill( ev.diphoton_eta1[k], s_weight );
              h_eta[xitight][i]->Fill( ev.diphoton_eta2[k], s_weight );
              h_r9lead[xitight][i]->Fill( ev.diphoton_r91[k], s_weight );
              h_r9sublead[xitight][i]->Fill( ev.diphoton_r92[k], s_weight );
              h_r9[xitight][i]->Fill( ev.diphoton_r91[k], s_weight );
              h_r9[xitight][i]->Fill( ev.diphoton_r92[k], s_weight );
              h_mass[xitight][i]->Fill( ev.diphoton_mass[k], s_weight );
              h_mass_rebin[xitight][i]->Fill( ev.diphoton_mass[k], s_weight );
              h_ptpair[xitight][i]->Fill( ev.diphoton_pt[k], s_weight );
              h_dphi[xitight][i]->Fill( acop*1.e3, s_weight );
              h_dphi_varbins[xitight][i]->Fill( acop, s_weight );
              h_dphi_logbins[xitight][i]->Fill( lacop, s_weight );
              h_dphi_tightbins[xitight][i]->Fill( acop*1.e3, s_weight );
              h_met[xitight][i]->Fill( ev.met, s_weight );
              h_xip[xitight][i]->Fill( xip, s_weight );
              h_xim[xitight][i]->Fill( xim, s_weight );
              h_diph_vtxz[xitight][i]->Fill( diph_vtx.z(), s_weight );
              h_diph_numjets[xitight][i]->Fill( num_matched_jets, s_weight );
              h_diph_numleptons[xitight][i]->Fill( num_matched_ele+num_matched_mu, s_weight );
              has_elastic_xitight_cand = true;
              evts_sel[i][xitight] += s_weight * sample_weight;
            }
          }
        }

        //h2_excl_sel[i]->Fill( num_matched_ele+num_matched_mu+num_matched_jets, acop, 1./num_entries );
        h2_excl_sel[sample_type]->Fill( closest_lep_vtx_dist, acop, s_weight );
        h2_excl_acop_dpt[sample_type]->Fill( fabs( ev.diphoton_pt1[k]-ev.diphoton_pt2[k] ), acop, s_weight );
        h2_excl_jet[sample_type]->Fill( num_matched_jets, acop, s_weight );

        if ( sample_type == mc_inclusive && has_elastic_diphoton_cand ) {
          for ( const auto& pot : pots_accept ) {
            if ( pot.first.find( "45" ) != string::npos && xip > pot.second )
              in_pot_accept[pot.first] += s_weight;
            if ( pot.first.find( "56" ) != string::npos && xim > pot.second )
              in_pot_accept[pot.first] += s_weight;
          }
        }

      } // loop on diphotons

      h_ndiph[nosel][i]->Fill( ev.num_diphoton, sp_weight );
      h_nvtx[nosel][i]->Fill( ev.num_vertex, sp_weight );
      h_nvtx_unc[nosel][i]->Fill( ev.num_vertex, 1. );
      num_after[nosel][i] += 1./s.num_events();
      /*if ( num_preselected_diphotons == 0 ) {
        h_ndiph[qcd][i]->Fill( ev.num_diphoton, sp_weight );
        h_nvtx[qcd][i]->Fill( ev.num_vertex, sp_weight );
        h_nvtx_unc[qcd][i]->Fill( ev.num_vertex, 1. );
        num_after[qcd][i]++;
      }*/
        //num_after[presel_noneleveto][i]++;
      if ( has_incl_diphoton_cand )
        num_after[incl][i] += 1./s.num_events();
      if ( has_diphoton_cand ) {
        h_ndiph[presel][i]->Fill( ev.num_diphoton, sp_weight );
        h_nvtx[presel][i]->Fill( ev.num_vertex, sp_weight );
        h_nvtx_unc[presel][i]->Fill( ev.num_vertex, 1. );
        num_after[presel][i] += 1./s.num_events();
        if ( has_elastic_diphoton_cand ) {
          h_ndiph[elastic][i]->Fill( ev.num_diphoton, sp_weight );
          h_nvtx[elastic][i]->Fill( ev.num_vertex, sp_weight );
          h_nvtx_unc[elastic][i]->Fill( ev.num_vertex, 1. );
          num_after[elastic][i] += 1./s.num_events();
          if ( has_elastic_xicomp_cand ) {
            num_after[xicomp][i] += 1./s.num_events();
            if ( has_elastic_xitight_cand ) {
              num_after[xitight][i] += 1./s.num_events();
            }
          }
        }
      }
    } // loop on events

    //for ( unsigned short j = 0; j < num_regions-3; ++j ) {
    for ( size_t j = 0; j < regions_in_cutflow.size(); ++j ) {
      const int reg = (int)regions_in_cutflow.at(j);
      for ( unsigned short k = 0; k < TreeEvent::num_classes-2; ++k ) {
        h_cutflow[i]->SetBinContent( j+1, evts_sel[i][reg] );
        h_cutflow[i]->GetXaxis()->SetBinLabel( j+1, kin_regions_label_small.at((int)regions_in_cutflow.at(j)).c_str() );
        //h_cutflow[i]->GetXaxis()->SetBinLabel( j+1, kin_regions_label.at(reg).c_str() );
      }
    }
    for ( unsigned short j = 0; j < num_regions; ++j ) {
      h_nvtx[j][i]->Scale( sample_weight );
      h_nvtx_unc[j][i]->Scale( sample_weight );
      h_ndiph[j][i]->Scale( sample_weight );
      h_mass[j][i]->Scale( sample_weight );
      h_mass_rebin[j][i]->Scale( sample_weight );
      h_ptpair[j][i]->Scale( sample_weight );
      h_met[j][i]->Scale( sample_weight );
      h_ptlead[j][i]->Scale( sample_weight );
      h_ptsublead[j][i]->Scale( sample_weight );
      h_pt[j][i]->Scale( sample_weight );
      h_pt_varbins[j][i]->Scale( sample_weight );
      h_etalead[j][i]->Scale( sample_weight );
      h_etasublead[j][i]->Scale( sample_weight );
      h_eta[j][i]->Scale( sample_weight );
      h_r9lead[j][i]->Scale( sample_weight );
      h_r9sublead[j][i]->Scale( sample_weight );
      h_r9[j][i]->Scale( sample_weight );
      h_dphi[j][i]->Scale( sample_weight );
      h_dphi_logbins[j][i]->Scale( sample_weight );
      h_dphi_varbins[j][i]->Scale( sample_weight );
      h_dphi_tightbins[j][i]->Scale( sample_weight );
      h_xip[j][i]->Scale( sample_weight );
      h_xim[j][i]->Scale( sample_weight );
      h_diph_vtxz[j][i]->Scale( sample_weight );
      h_diph_numjets[j][i]->Scale( sample_weight );
      h_diph_numleptons[j][i]->Scale( sample_weight );

      if ( s.type() == Sample::kBackground ) {
        h_nvtx[j][i]->SetFillColorAlpha( s.colour(), alpha );
        h_nvtx[j][i]->SetLineColor( s.colour() );
        h_nvtx_unc[j][i]->SetFillColorAlpha( s.colour(), alpha );
        h_nvtx_unc[j][i]->SetLineColor( s.colour() );
        h_ndiph[j][i]->SetFillColorAlpha( s.colour(), alpha );
        h_mass[j][i]->SetFillColorAlpha( s.colour(), alpha );
        h_mass_rebin[j][i]->SetFillColorAlpha( s.colour(), alpha );
        h_ptpair[j][i]->SetFillColorAlpha( s.colour(), alpha );
        h_met[j][i]->SetFillColorAlpha( s.colour(), alpha );
        h_ptlead[j][i]->SetFillColorAlpha( s.colour(), alpha );
        h_ptsublead[j][i]->SetFillColorAlpha( s.colour(), alpha );
        h_pt[j][i]->SetFillColorAlpha( s.colour(), alpha );
        h_pt_varbins[j][i]->SetFillColorAlpha( s.colour(), alpha );
        h_etalead[j][i]->SetFillColorAlpha( s.colour(), alpha );
        h_etasublead[j][i]->SetFillColorAlpha( s.colour(), alpha );
        h_eta[j][i]->SetFillColorAlpha( s.colour(), alpha );
        h_r9lead[j][i]->SetFillColorAlpha( s.colour(), alpha );
        h_r9sublead[j][i]->SetFillColorAlpha( s.colour(), alpha );
        h_r9[j][i]->SetFillColorAlpha( s.colour(), alpha );
        h_dphi[j][i]->SetFillColorAlpha( s.colour(), alpha );
        h_dphi_varbins[j][i]->SetFillColorAlpha( s.colour(), alpha );
        h_dphi_logbins[j][i]->SetFillColorAlpha( s.colour(), alpha );
        h_dphi_tightbins[j][i]->SetFillColorAlpha( s.colour(), alpha );
        h_xip[j][i]->SetFillColorAlpha( s.colour(), alpha );
        h_xim[j][i]->SetFillColorAlpha( s.colour(), alpha );
        h_diph_vtxz[j][i]->SetFillColorAlpha( s.colour(), alpha );
        h_diph_numjets[j][i]->SetFillColorAlpha( s.colour(), alpha );
        h_diph_numleptons[j][i]->SetFillColorAlpha( s.colour(), alpha );
      }

      hm_nvtx[j][sample_type].emplace_back( sname, h_nvtx[j][i] );
      hm_nvtx_unc[j][sample_type].emplace_back( sname, h_nvtx_unc[j][i] );
      hm_ndiph[j][sample_type].emplace_back( sname, h_ndiph[j][i] );
      hm_mass[j][sample_type].emplace_back( sname, h_mass[j][i] );
      hm_mass_rebin[j][sample_type].emplace_back( sname, h_mass_rebin[j][i] );
      hm_ptpair[j][sample_type].emplace_back( sname, h_ptpair[j][i] );
      hm_met[j][sample_type].emplace_back( sname, h_met[j][i] );
      hm_ptlead[j][sample_type].emplace_back( sname, h_ptlead[j][i] );
      hm_ptsublead[j][sample_type].emplace_back( sname, h_ptsublead[j][i] )  ;
      hm_pt[j][sample_type].emplace_back( sname, h_pt[j][i] )  ;
      hm_pt_varbins[j][sample_type].emplace_back( sname, h_pt_varbins[j][i] )  ;
      hm_etalead[j][sample_type].emplace_back( sname, h_etalead[j][i] );
      hm_etasublead[j][sample_type].emplace_back( sname, h_etasublead[j][i] )  ;
      hm_eta[j][sample_type].emplace_back( sname, h_eta[j][i] )  ;
      hm_r9lead[j][sample_type].emplace_back( sname, h_r9lead[j][i] );
      hm_r9sublead[j][sample_type].emplace_back( sname, h_r9sublead[j][i] );
      hm_r9[j][sample_type].emplace_back( sname, h_r9[j][i] );
      hm_dphi[j][sample_type].emplace_back( sname, h_dphi[j][i] );
      hm_dphi_varbins[j][sample_type].emplace_back( sname, h_dphi_varbins[j][i] );
      hm_dphi_logbins[j][sample_type].emplace_back( sname, h_dphi_logbins[j][i] );
      hm_dphi_tightbins[j][sample_type].emplace_back( sname, h_dphi_tightbins[j][i] );
      hm_xip[j][sample_type].emplace_back( sname, h_xip[j][i] );
      hm_xim[j][sample_type].emplace_back( sname, h_xim[j][i] );
      hm_diph_vtxz[j][sample_type].emplace_back( sname, h_diph_vtxz[j][i] );
      hm_diph_numjets[j][sample_type].emplace_back( sname, h_diph_numjets[j][i] );
      hm_diph_numleptons[j][sample_type].emplace_back( sname, h_diph_numleptons[j][i] );
    }
    h_cutflow[i]->Sumw2();
    //h_cutflow[i]->Scale( sample_weight );
    if ( s.type() == Sample::kBackground )
      h_cutflow[i]->SetFillColorAlpha( s.colour(), alpha );
    hm_cutflow[sample_type].emplace_back( sname, h_cutflow[i] );
  } // loop on samples

  //-------------------------------------------------------------------------------------
  // HERE BEGINS THE PLOTTING PART
  //-------------------------------------------------------------------------------------

  {
    const auto region = xitight;
    //auto& hists_mass = hm_mass[region][0][1]->GetHists()
    auto hm_dat = hm_mass_rebin[region][0];
    if ( hm_dat.size() != 1 )
      throw length_error( "Invalid size for data mass distribution" );
    auto h_data = (TH1D*)hm_dat.begin()->second->Clone();

    auto hm_sig = hm_mass_rebin[region][2];
    //if ( hm_sig.size() != 1 )
    //  throw length_error( "Invalid size for signal mass distribution" );
    //FIXME!!!
    if ( !hm_sig.empty() ) {
      auto h_signal = (TH1D*)hm_sig.begin()->second->Clone();
      h_signal->Scale( 1./scaling_signal );

      auto hm_bck = hm_mass_rebin[region][1];
      TH1D* h_backgr = nullptr;
      unsigned short i = 0;
      for ( const auto& obj : hm_bck ) {
        if ( i++ == 0 ) h_backgr = (TH1D*)hm_bck.begin()->second->Clone();
        else h_backgr->Add( obj.second );
      }
      //auto h_mass_sum = (TH1D*)( hm_mass_rebin[xitight][0][1]->GetStack().Last() );
      TFile f( "sum_bckg_mass.root", "recreate" );
      //h_data->Write( "data_obs" );
      h_signal->Write( "aaaa" );
      h_backgr->Write( "incl" );
    }
  }

  gStyle->SetOptStat( 0 );
  //Plotter plt( "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp", Form( "CMS Preliminary 2016, #sqrt{s} = 13 TeV, L = %.1f fb^{-1}", the_lumi*1.e-3 ) );
  //Plotter plt( "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp", Form( "#sqrt{s} = 13 TeV, L = %g fb^{-1}", the_lumi*1.e-3 ) );
  Plotter plt( "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp", Form( "%.1f fb^{-1} (13 TeV)", the_lumi*1.e-3 ) );

  plt.draw_multiplot( "cutflow", hm_cutflow[0], hm_cutflow[1], hm_cutflow[2], "", false, true /*logy*/, true /*legend*/, false/*, 0.1, 1.9*/ );

  for ( unsigned short i = 0; i < num_regions; ++i ) {
    const string class_name = Form( "#font[52]{%s}", kin_regions_label[i].c_str() );
    // filename, data, mc, sig, label="", colours=1, logy=0, legend=1, logx=0, min_y=-0.4, max_y=2.4
    plt.draw_multiplot( kin_regions[i]+"_diphoton_mass", hm_mass[i][0], hm_mass[i][1], hm_mass[i][2], class_name, false, false );
    plt.draw_multiplot( kin_regions[i]+"_diphoton_mass_logscale", hm_mass[i][0], hm_mass[i][1], hm_mass[i][2], class_name, false, true );
    plt.draw_multiplot( kin_regions[i]+"_diphoton_ptpair", hm_ptpair[i][0], hm_ptpair[i][1], hm_ptpair[i][2], class_name, false, false );
    plt.draw_multiplot( kin_regions[i]+"_diphoton_ptpair_logscale", hm_ptpair[i][0], hm_ptpair[i][1], hm_ptpair[i][2], class_name, false, true );
    plt.draw_multiplot( kin_regions[i]+"_diphoton_met", hm_met[i][0], hm_met[i][1], hm_met[i][2], class_name, false );
    plt.draw_multiplot( kin_regions[i]+"_diphoton_pt_singlepho", hm_pt[i][0], hm_pt[i][1], hm_pt[i][2], class_name, false );
    plt.draw_multiplot( kin_regions[i]+"_diphoton_pt_singlepho_varbins", hm_pt_varbins[i][0], hm_pt_varbins[i][1], hm_pt_varbins[i][2], class_name, false );
    plt.draw_multiplot( kin_regions[i]+"_diphoton_pt_leadpho", hm_ptlead[i][0], hm_ptlead[i][1], hm_ptlead[i][2], class_name, false );
    plt.draw_multiplot( kin_regions[i]+"_diphoton_pt_subleadpho", hm_ptsublead[i][0], hm_ptsublead[i][1], hm_ptsublead[i][2], class_name, false );
    plt.draw_multiplot( kin_regions[i]+"_diphoton_eta_singlepho", hm_eta[i][0], hm_eta[i][1], hm_eta[i][2], class_name, false );
    plt.draw_multiplot( kin_regions[i]+"_diphoton_eta_leadpho", hm_etalead[i][0], hm_etalead[i][1], hm_etalead[i][2], class_name, false );
    plt.draw_multiplot( kin_regions[i]+"_diphoton_eta_subleadpho", hm_etasublead[i][0], hm_etasublead[i][1], hm_etasublead[i][2], class_name, false );
    plt.draw_multiplot( kin_regions[i]+"_diphoton_r9_singlepho", hm_r9[i][0], hm_r9[i][1], hm_r9[i][2], class_name, false, false );
    plt.draw_multiplot( kin_regions[i]+"_diphoton_r9_singlepho_logscale", hm_r9[i][0], hm_r9[i][1], hm_r9[i][2], class_name, false, true );
    plt.draw_multiplot( kin_regions[i]+"_diphoton_r9_leadpho", hm_r9lead[i][0], hm_r9lead[i][1], hm_r9lead[i][2], class_name, false );
    plt.draw_multiplot( kin_regions[i]+"_diphoton_r9_subleadpho", hm_r9sublead[i][0], hm_r9sublead[i][1], hm_r9sublead[i][2], class_name, false );
    plt.draw_multiplot( kin_regions[i]+"_diphoton_dphi", hm_dphi[i][0], hm_dphi[i][1], hm_dphi[i][2], class_name, false );
    //plt.draw_multiplot( kin_regions[i]+"_diphoton_dphi_varbins", hm_dphi_varbins[i][0], hm_dphi_varbins[i][1], hm_dphi_varbins[i][2], class_name, false, true, false, true, -2.1, 4.1 );
    plt.draw_multiplot( kin_regions[i]+"_diphoton_dphi_varbins", hm_dphi_varbins[i][0], hm_dphi_varbins[i][1], hm_dphi_varbins[i][2], class_name+" (no acop.sel.)", false, true, false, true );
    plt.draw_multiplot( kin_regions[i]+"_diphoton_dphi_logbins", hm_dphi_logbins[i][0], hm_dphi_logbins[i][1], hm_dphi_logbins[i][2], class_name+" (no acop.sel.)", false, true, false, false );
    plt.draw_multiplot( kin_regions[i]+"_diphoton_dphi_tightbins", hm_dphi_tightbins[i][0], hm_dphi_tightbins[i][1], hm_dphi_tightbins[i][2], class_name+" (no acop.sel.)", false, true, false, false, 0.24, 1.76 );
    plt.draw_multiplot( kin_regions[i]+"_diphoton_xip", hm_xip[i][0], hm_xip[i][1], hm_xip[i][2], class_name, false );
    plt.draw_multiplot( kin_regions[i]+"_diphoton_xim", hm_xim[i][0], hm_xim[i][1], hm_xim[i][2], class_name, false );
    plt.draw_multiplot( kin_regions[i]+"_diphoton_vtx_z", hm_diph_vtxz[i][0], hm_diph_vtxz[i][1], hm_diph_vtxz[i][2], class_name, false );
    plt.draw_multiplot( kin_regions[i]+"_diphoton_numjets", hm_diph_numjets[i][0], hm_diph_numjets[i][1], hm_diph_numjets[i][2], class_name, false, true );
    plt.draw_multiplot( kin_regions[i]+"_diphoton_numleptons", hm_diph_numleptons[i][0], hm_diph_numleptons[i][1], hm_diph_numleptons[i][2], class_name, false );
    plt.draw_multiplot( kin_regions[i]+"_diphoton_mult", hm_ndiph[i][0], hm_ndiph[i][1], hm_ndiph[i][2], class_name, false, true );
    plt.draw_multiplot( kin_regions[i]+"_diphoton_nvtx", hm_nvtx[i][0], hm_nvtx[i][1], hm_nvtx[i][2], class_name, false, false );
    plt.draw_multiplot( kin_regions[i]+"_diphoton_nvtx_logscale", hm_nvtx[i][0], hm_nvtx[i][1], hm_nvtx[i][2], class_name, false, true, false );
    plt.draw_multiplot( kin_regions[i]+"_diphoton_nvtx_uncorr_logscale", hm_nvtx_unc[i][0], hm_nvtx_unc[i][1], hm_nvtx_unc[i][2], class_name, false, true, false );
  }
  if ( !rescale_background_parametric ) {
    gStyle->SetOptFit( 1 );
    Canvas c( "data_mc_dphi_presel", Form( "%.1f fb^{-1} (13 TeV)", the_lumi*1.e-3 ), "Preliminary" );
    const auto reg = incl;
//    c.SetLogx();
    THStack hs_mc;
    TH1D* h_mc = nullptr;
    unsigned short i = 0;
    for ( const auto& plt : hm_dphi_logbins[reg][1] ) {
      if ( i == 0 )
        h_mc = (TH1D*)plt.second->Clone();
      else
        h_mc->Add( plt.second );
      ++i;
    }
    auto h_data = hm_dphi_logbins[reg][0].begin()->second;
    auto h_ratio = (TH1D*)h_data->Clone();
    h_ratio->Divide( h_mc );
    h_ratio->Draw( "p" );
    auto fit_func = new TF1( "fit_func", "pol1", -3., 0. );
    h_ratio->Fit( fit_func );
    cout << "integral [0.,0.005] = " << fit_func->Integral( -4., log10( 0.005 ) ) << endl;
    cout << "integral [0.25,1.0] = " << fit_func->Integral( log10(0.25), log10( 1. ) ) << endl;
    fit_func->SaveAs( "incl_fit_result.root" );
    c.Prettify( h_ratio );
    h_ratio->GetYaxis()->SetTitle( "Data/Pred." );
    PaveText::topLabel( ( kin_regions_label[reg]+" (no acop.sel.)" ).c_str() );
    c.Save( "png,pdf", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
  }
  {
    Canvas c( "num_diph_vtxaround", Form( "%.1f fb^{-1} (13 TeV)", the_lumi*1.e-3 ), "Preliminary" );
    THStack hs;
    c.SetLegendY1( 0.75 );
    TString label[4] = { "at 1 mm distance", "at 2 mm distance", "at 5 mm distance", "at 1 cm distance" };
    for ( unsigned short i = 0; i < 4; ++i ) {
      h_diph_nvtxaround[i]->SetMarkerStyle( 20+i );
      h_diph_nvtxaround[i]->SetMarkerColor( Canvas::colour_pool[i] );
      h_diph_nvtxaround[i]->SetLineColor( kBlack );
      h_diph_nvtxaround[i]->SetLineWidth( 2 );
      c.AddLegendEntry( h_diph_nvtxaround[i], label[i], "elp" );
      hs.Add( h_diph_nvtxaround[i] );
    }
    hs.Draw( "e0,nostack" );
    hs.GetHistogram()->SetTitle( h_diph_nvtxaround[0]->GetTitle() );
    //hs.GetXaxis()->SetTitle( h_diph_nvtxaround[0]->GetXaxis()->GetTitle() );
    //hs.GetYaxis()->SetTitle( h_diph_nvtxaround[0]->GetYaxis()->GetTitle() );
    c.Prettify( hs.GetHistogram() );
    hs.GetHistogram()->SetTitle( "" );
    hs.GetHistogram()->SetMaximum( 1.e3 );
    c.SetLogy();
    PaveText::topLabel( kin_regions_label[elastic].c_str() );
    c.Save( "png,pdf", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
    gStyle->SetOptFit( 0 );
  }

  for ( unsigned short i = 0; i < num_regions; ++i ) {
    cout << ">>>>>>>>>>> REGION: " << kin_regions[i] << " <<<<<<<<<<<" << endl;
    double n_data = 0., n_data_err = 0.;
    double n_mc = 0., n_mc_err = 0.;
    double n_sublead_mc = 0., n_sublead_mc_err = 0.;
    double n_sig = 0., n_sig_err = 0.;
    for ( const auto& h : hm_xim[i][0] ) { // data
      double n1, n2;
      n1 = h.second->IntegralAndError( 0, h.second->GetNbinsX()+1, n2 );
      n_data += n1;
      //n_data_err += n2*n2;
    }
    for ( const auto& h : hm_xim[i][1] ) { // MC
      double n1, n2;
      n1 = h.second->IntegralAndError( 0, h.second->GetNbinsX()+1, n2 );
      if ( h.first.find( "#gamma#gamma" ) == string::npos ) {
        n_sublead_mc += n1;
        n_sublead_mc_err += n2*n2;
      }
      n_mc += n1;
      n_mc_err += n2*n2;
      cout << setw( 30 ) << h.first << " >> " << n1 << " +/- " << n2 << endl;
    }
    //for ( const auto& h : hm_xim[i][2] ) { // signal
    if ( !hm_xim[i][2].empty() ) {
      const auto& h = *hm_xim[i][2].begin();
      double n1, n2;
      n1 = h.second->IntegralAndError( 0, h.second->GetNbinsX()+1, n2 );
      n_sig += n1;
      n_sig_err += n2*n2;
    }
    cout << ">>> TOTAL::EVENTS <<<\n"
         << "signal: " << n_sig << " +/- " << sqrt( n_sig_err ) << "\n"
         << "  data: " << n_data << " +/- " << sqrt( n_data_err ) << "\n"
         << "    MC: " << n_mc << " +/- " << sqrt( n_mc_err ) << "\n"
         << "sub MC: " << n_sublead_mc << " +/- " << sqrt( n_sublead_mc_err ) << "\n"
         << " ratio: " << n_data/n_mc << " +/- " << ( n_data/n_mc*sqrt( n_data_err/n_data/n_data+n_mc_err/n_mc/n_mc ) ) << endl;
    cout << "events fractions:\n";
    size_t j = 0;
    for ( const auto& num : num_after[i] ) {
      const auto& smp = dsh.sample( j++ );
      cout << smp.name() << ":: "<< num << endl;
    }
    cout << endl;
  }
  TF1* f_xim = nullptr, *f_xip = nullptr;
  TGraph* gr_1sig_xim = nullptr, *gr_1sig_xip = nullptr;
  TH1D* h_mc_xim = nullptr, *h_mc_xip = nullptr;
  TH1D* h_data_xim = nullptr, *h_data_xip = nullptr;
  TFitResult* fit_xim = nullptr, *fit_xip = nullptr;
  TVirtualFitter::SetDefaultFitter( "Minuit" );
  for ( auto& hm : map<string,Plotter::HistsMap*>{ { "xim_incl", hm_xim[incl] }, { "xip_incl", hm_xip[incl] } } ) {
    gStyle->SetOptFit( 1 );
    Canvas c( hm.first.c_str(), Form( "%.1f fb^{-1} (13 TeV)", the_lumi*1.e-3 ), "Preliminary" );
    auto h_sum_data = (TH1D*)hm.second[the_data].begin()->second->Clone();
    h_sum_data->SetMarkerSize( 1.0 );
    h_sum_data->SetTitle( Form( "#xi_{#gamma#gamma}^{%s}@@Events@@?.2f", hm.first == "xim_incl" ? "-" : "+" ) );
    auto h_sum_mc = (TH1D*)h_sum_data->Clone();
    h_sum_mc->Clear();
    for ( const auto& h : hm.second[mc_inclusive] ) {
      if ( h.first.find( "QCD" ) != string::npos )
        continue;
      h_sum_mc->Add( h.second );
    }
    for ( const auto& h : hm.second[the_data] )
      h_sum_data->Add( h.second );
    h_sum_data->SetLineColor( kBlack );
    h_sum_data->SetMarkerStyle( 20 );
    //h_sum_data->SetMarkerColor( kGray+1 );
    h_sum_data->Draw( "p" );
    //h_sum_mc->SetTitle( Form( ";#xi_{#gamma#gamma}^{%s};Events", hm.first == "xim_incl" ? "-" : "+" ) );
//    c.AddLegendEntry( h_sum_mc, "#Sigma MC", "lp" );
    /*h_sum_data->SetLineColor( kBlack );
    h_sum_data->SetMarkerStyle( 20 );
    h_sum_data->SetMarkerColor( kBlack );
    h_sum_data->Draw( "p" );
    c.AddLegendEntry( h_sum_data, "Data", "lp" );*/
    //TFitResultPtr fit = h_sum_mc->Fit( "[0]+exp([1]*x)", "sr", "", 0.025, 0.5 );
    //TFitResultPtr fit_mc = h_sum_mc->Fit( "expo", "sr", "", 0.012, 0.14 );
    //TFitResultPtr fit_mc = h_sum_mc->Fit( &f_expo, "sr", "", 350./13000., 0.2 );
    //TFitResultPtr fit_data = h_sum_data->Fit( &f_expo, "sr", "", 350./13000., 0.2 );
    TFitResultPtr fit_data = h_sum_data->Fit( &f_expo, "sr", "", 350./13000., 0.15 );
    //TFitResultPtr fit_data = h_sum_data->Fit( &f_expo, "sr", "", 0.045, 0.175 );
    /////TFitResultPtr fit_data = h_sum_data->Fit( &f_expo, "sr", "", 0.065, 0.3 );
    cout << hm.first << ", prob = " << fit_data->Prob() << endl;
    TString label;
    if ( hm.first == "xim_incl" ) {
      f_xim = (TF1*)f_expo.Clone( "fit_xim" );
      h_mc_xim = (TH1D*)h_sum_mc->Clone( "sum_mc_xim" );
      h_data_xim = (TH1D*)h_sum_data->Clone( "sum_data_xim" );
      fit_xim = (TFitResult*)fit_data->Clone( "res_fit_xim" );
      label = "#xi_{#gamma#gamma}^{-}";
    }
    else {
      f_xip = (TF1*)f_expo.Clone( "fit_xip" );
      h_mc_xip = (TH1D*)h_sum_mc->Clone( "sum_mc_xip" );
      h_data_xip = (TH1D*)h_sum_mc->Clone( "sum_data_xip" );
      fit_xip = (TFitResult*)fit_data->Clone( "res_fit_xip" );
      label = "#xi_{#gamma#gamma}^{+}";
    }
    fit_data->GetCorrelationMatrix().Print();
    //TFitResultPtr fit_mc = h_sum_mc->Fit( "expo", "s" );
    //TFitResultPtr fit_data = h_sum_data->Fit( "expo", "s+" );
    /*if ( fit.Get() ) {
      const double* params_fit = fit->GetParams();
      c.AddLegendEntry( h_sum, Form( "%g+Sect.45, Mean = %.2e, RMS = %.3f", params_fit[0], params_fit[1] ), "f" );
    }*/
    auto top = PaveText::topLabel( ( kin_regions_label[incl]+" (data)" ).c_str() );
    c.Prettify( h_sum_data );
    c.Save( "png,pdf", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
    //--- parameters correlation
    Canvas c2( ( hm.first+"_param_corr" ).c_str(), Form( "%.1f fb^{-1} (13 TeV)", the_lumi*1.e-3 ), "Preliminary" );
    c2.SetLegendX1( 0.6 );
    c2.SetLegendY1( 0.2 );
    gMinuit->SetErrorDef( 4 );
    auto gr2 = dynamic_cast<TGraph*>( gMinuit->Contour( 100, 0, 1 ) );
    //gr2->SetFillColor( 42 );
    //gr2->SetFillColor( 5 );
    gr2->SetFillColor( kYellow-4 );
    gr2->Draw( "alf" );
    gMinuit->SetErrorDef( 1 );
    auto gr1 = dynamic_cast<TGraph*>( gMinuit->Contour( 100, 0, 1 ) );
    //gr1->SetFillColor( 38 );
    //gr1->SetFillColor( 3 );
    if ( hm.first == "xim_incl" )
      gr_1sig_xim = (TGraph*)gr1->Clone( "fit_xim_1sig" );
    else
      gr_1sig_xip = (TGraph*)gr1->Clone( "fit_xip_1sig" );
    gr1->SetFillColor( kGreen+1 );
    gr1->Draw( "lf" );
    auto params = fit_data->GetParams(), err_params = fit_data->GetErrors();
    TGraphErrors gr_central;
    gr_central.SetPoint( 0, params[0], params[1] );
    gr_central.SetPointError( 0, err_params[0], err_params[1] );
    gr_central.SetLineWidth( 3 );
    //gr_central.SetMarkerStyle( 20 );
    //gr_central.SetMarkerColor( kBlack );
    gr_central.Draw( "e,same" );
    gr2->GetHistogram()->SetTitle( ";Intercept;Slope" );
    c2.Prettify( gr2->GetHistogram() );
    top->Draw( "same" );
    c2.AddLegendEntry( &gr_central, "1D uncert.", "l" );
    c2.AddLegendEntry( gr1, "#pm 1 #sigma", "f" );
    c2.AddLegendEntry( gr2, "#pm 2 #sigma", "f" );
    PaveText pt( 0.18, 0.7 );
    pt.AddText( Form( "%s (#rho = %.3f)", label.Data(), fit_data->GetCorrelationMatrix()[0][1] ) );
    pt.SetTextAlign( kHAlignLeft+kVAlignBottom );
    pt.SetTextSize( 0.045 );
    pt.Draw( "same" );
    c2.Save( "png,pdf", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
  }
  auto f_fits = TFile::Open( "fits_results.root", "recreate" );
  f_xim->Write();
  f_xip->Write();
  gr_1sig_xim->Write();
  gr_1sig_xip->Write();
  fit_xim->Write();
  fit_xip->Write();
  h_mc_xim->Write();
  h_mc_xip->Write();
  h_data_xim->Write();
  h_data_xip->Write();
  //delete f_fits;

  /*FIXME FIXME FIXME FIXME FIXME FIXME FIXME
  plot_2ddiscrim( "presel_excl_selection", h2_excl_sel, true );
  plot_2ddiscrim( "presel_excl_dpt_vs_acop", h2_excl_acop_dpt, true );
  plot_2ddiscrim( "presel_excl_matched_selection", h2_excl_jet, false );
  */

  for ( unsigned short i = 0; i < num_samples; ++i ) {
    Canvas c( Form( "nosel_eleveto_%d", i ), Form( "%.1f fb^{-1} (13 TeV)", the_lumi*1.e-3 ), "Preliminary" );
    c.Pad()->SetRightMargin( 0.12 );
    gStyle->SetPalette( kLightTemperature );
    h2_eleveto[i]->Scale( 1./h2_eleveto[i]->Integral() );
    h2_eleveto[i]->Draw( "colz text" );
    h2_eleveto[i]->SetMarkerSize( 2. );
    PaveText::topLabel( dsh.sample( i ).name() );
    c.Prettify( h2_eleveto[i] );
    c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
  }

  {
    Canvas c( "presel_diphoton_pt_sigovbckg", Form( "%.1f fb^{-1} (13 TeV)", the_lumi*1.e-3 ), "Preliminary" );
    TH1D* hist_bck = (TH1D*)hm_ptpair[presel][1][0].second->Clone( "bck" ),
         *hist_sig = (TH1D*)hm_ptpair[presel][2][0].second->Clone( "sig" );
    hist_bck->Scale( 0. ); hist_sig->Scale( 0. );
    for ( const auto& hm : hm_ptpair[presel][1] ) { hist_bck->Add( hm.second ); }
    for ( const auto& hm : hm_ptpair[presel][2] ) { hist_sig->Add( hm.second ); }

    TH1D* hist_signif = (TH1D*)hist_bck->Clone( "signif" ); hist_signif->Scale( 0. );
    for ( int j = 1; j < hist_signif->GetNbinsX(); ++j ) {
      const double bck_value = hist_bck->GetBinContent( j ), sig_value = hist_sig->GetBinContent( j );
      const double signif = ( sig_value+bck_value != 0. ) ? sig_value / sqrt( sig_value+bck_value ) : 0.;
      //cout << " bin " << j << ": background: " << bck_value << ", signal: " << sig_value << " --> significance: " << signif << endl;
      hist_signif->SetBinContent( j, signif );
    }
    hist_signif->Draw( "p" );
    hist_signif->SetMarkerStyle( 20 );
    c.Prettify( hist_signif );
    c.SetLogy();
    c.Save( "png,pdf", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/mc_comparison" );
  }
  cout << "aaaaaaa" << endl;
}


void logarithmicBins( TAxis* axis )
{
  int bins = axis->GetNbins();

  Axis_t from = axis->GetXmin();
  Axis_t to = axis->GetXmax();
  Axis_t width = ( to-from )/bins;
  Axis_t *new_bins = new Axis_t[bins+1];
  for ( int i = 0; i <= bins; ++i )
    new_bins[i] = TMath::Power( 10, from+i*width );
  axis->Set(bins, new_bins);
  delete [] new_bins;
}

void plot_2ddiscrim( const char* name, TH2D* h2[], bool logx )
{
  Canvas c( name, "13 TeV", "Simulation" );
  THStack st;
  for ( unsigned short i = 1; i < num_types; ++i ) {
    st.Add( h2[i] );
    st.SetTitle( h2[i]->GetTitle() );
    h2[i]->SetLineColor( 1+i );
    h2[i]->SetMarkerStyle( 19+i );
    h2[i]->SetMarkerColor( 1+i );
    h2[i]->SetMarkerSize( 0.5 );
  }
  st.Draw( "p,nostack" );
  //st.Draw( "box,nostack" );
  //c.AddLegendEntry( h2[the_data], "Data" );
  c.AddLegendEntry( h2[mc_inclusive], "Incl. backgrounds" );
  c.AddLegendEntry( h2[mc_signal], "Exclusive signal" );
  if ( logx ) c.SetLogx();
  c.SetLogy();
  TLine cut( st.GetHistogram()->GetXaxis()->GetXmin(), 0.005, st.GetHistogram()->GetXaxis()->GetXmax(), 0.005 );
  cut.SetLineColor( kBlack );
  cut.SetLineWidth( 3 );
  cut.SetLineStyle( 2 );
  cut.Draw();
  c.Prettify( st.GetHistogram() );
  c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
}

