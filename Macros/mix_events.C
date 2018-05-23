#include "Canvas.h"
#include "pot_alignment.h"
#include "xi_reconstruction.h"
#include "diproton_candidate.h"
#include "tree_reader.h"

#include "TH2.h"
#include "TF1.h"
#include "TStyle.h"

//#define OUTPUT_DIR "/afs/cern.ch/user/l/lforthom/www/private/twophoton/test"
#define OUTPUT_DIR "/afs/cern.ch/user/l/lforthom/www/private/twophoton/test_15may"

void mix_events()
{
//  const double exp_yield = 586.877;
//  const double exp_yield = 525;
  const double exp_yield = 566.994;
  const double rel_err_mass = 0.02, rel_err_rap = 0.074;

  //----- inclusive selection
  // 45
  auto mix_45 = new TF1( "mix_45", "expo", 0., 0.2 );
  /*mix_45->SetParameter( 0, 4.90128 ); //mix_45->SetParError( 0, 0.0136147 );
  mix_45->SetParameter( 1, -10.4410 ); //mix_45->SetParError( 1, 0.165792 );
  mix_45->SetParameter( 0, 5.69505e+00 );
  mix_45->SetParameter( 1, -1.05041e+01 );*/
  mix_45->SetParameter( 0, 6.12193e+00 ); mix_45->SetParError( 0, 3.05054e-02 );
  mix_45->SetParameter( 1, -1.64478e+01 ); mix_45->SetParError( 1, 4.91350e-01 );
  // 56
  auto mix_56 = new TF1( "mix_56", "expo", 0., 0.2 );
  /*mix_56->SetParameter( 0, 5.10791 ); //mix_56->SetParError( 0, 0.0174085 );
  mix_56->SetParameter( 1, -12.4708 ); //mix_56->SetParError( 1, 0.176941 );
  mix_56->SetParameter( 0, 5.87960e+00 );
  mix_56->SetParameter( 1, -1.17574e+01 );*/
  mix_56->SetParameter( 0, 6.16033e+00 ); mix_56->SetParError( 0, 3.12646e-02 );
  mix_56->SetParameter( 1, -1.72669e+01 ); mix_56->SetParError( 1, 4.99068e-01 );

  TFile f( "Samples/output_Run2016BCG_looseCuts_28jun.root" );
  //TFile f( "Samples/output_alignmentrun_mbtree.root" );
  treeinfo ev;
  ev.read( dynamic_cast<TTree*>( f.Get( "ntp" ) ) );
  const unsigned long long num_toys = 100000;
  const unsigned long long num_events = ev.tree->GetEntriesFast(); // will correspond to the toys number
  //const unsigned long long num_toys = ev.tree->GetEntriesFast(); // will correspond to the toys number
  xi_reco::load_file( "TreeProducer/data/optics_jun22.root" );
  pot_align::load_file( "TreeProducer/data/alignment_collection_v2.out" );

  map<unsigned short,TH2D*> m_h2_xicorr;
  map<unsigned short,TH1D*> m_h_xireco;
  map<unsigned short,const char*> m_pot_names = { { 2, "45N" }, { 3, "45F" }, { 102, "56N" }, { 103, "56F" } };
  //map<unsigned short,double> m_pot_limits = { { 2, 0.067 }, { 3, 0.066 }, { 102, 0.070 }, { 103, 0.061 } }; // < 10% radiation damage
  map<unsigned short,double> m_pot_limits = { { 2, 0.034 }, { 3, 0.023 }, { 102, 0.042 }, { 103, 0.032 } };
  const double min_xigg = 0.02;
  for ( const auto& p : m_pot_names ) {
    m_h2_xicorr[p.first] = new TH2D( Form( "xi_corr_%d", p.first ), Form( ";#xi_{%s};#xi_{#gamma#gamma}", p.second ), 50, 0., 0.2, 50, 0., 0.2 );
    m_h_xireco[p.first] = new TH1D( Form( "xi_reco_%d", p.first ), Form( ";$xi_{%s};Events", p.second ), 50, 0., 0.2 );
  }
  //auto h2_mcorr = new TH2D( "mcorr", ";m_{pp} (GeV);m_{#gamma#gamma} (GeV)", 40, 250., 1850., 40, 250., 1850. ),
  auto h2_mcorr = new TH2D( "mcorr", ";m_{pp} (GeV);m_{#gamma#gamma} (GeV)", 40, 450., 2050., 40, 450., 2050. ),
       h2_mcorr_2sigma = (TH2D*)h2_mcorr->Clone( "mcorr_2sigma" ),
       h2_mcorr_3sigma = (TH2D*)h2_mcorr->Clone( "mcorr_3sigma" );
  //auto h2_ycorr = new TH2D( "ycorr", ";y_{pp};y_{#gamma#gamma}", 50, -2., 2., 50, -2., 2. ),
  auto h2_ycorr = new TH2D( "ycorr", ";y_{pp};y_{#gamma#gamma}", 50, -1., 1., 50, -1., 1. ),
       h2_ycorr_2sigma = (TH2D*)h2_ycorr->Clone( "ycorr_2sigma" ),
       h2_ycorr_3sigma = (TH2D*)h2_ycorr->Clone( "ycorr_3sigma" );

  double num_match2d_2sigma = 0., num2_match2d_2sigma = 0.;
  double num_match2d_3sigma = 0., num2_match2d_3sigma = 0.;

  cout << "number of toys: " << num_toys << endl;
  const double int_45 = mix_45->Integral( min_xigg, 0.15 );
  const double int_56 = mix_56->Integral( min_xigg, 0.15 );
  const double mean_int = hypot( int_45, int_56 );
  /*const double int_45 = mix_45->Integral( m_pot_limits[103], 0.15 );
  const double int_56 = mix_56->Integral( m_pot_limits[3], 0.15 );*/

  const double weight = exp_yield/num_toys/mean_int; //FIXME

  for ( unsigned long long i = 0; i < num_toys; ++i ) {
    ev.tree->GetEntry( rand()*1./RAND_MAX*num_events );
    //ev.tree->GetEntry( i );
    if ( fmod( i*1., num_toys*0.05 ) == 0 )
      cout << "event " << i << endl;

    auto align = pot_align::get_alignments( ev.fill_number );

    const double xi_45_rnd = mix_45->GetRandom( min_xigg, 0.15 );
    const double xi_56_rnd = mix_56->GetRandom( min_xigg, 0.15 );
    const double m_diph = 13.e3*sqrt( xi_45_rnd*xi_56_rnd );
    const double y_diph = 0.5*log( xi_45_rnd/xi_56_rnd );

    map<unsigned short,pair<double,double> > xi_meas45, xi_meas56;
    for ( unsigned short j = 0; j < ev.num_proton_track; ++j ) {
      const unsigned short pot_id = 100*ev.proton_track_side[j]+ev.proton_track_pot[j];
      auto al = align[pot_id];
//      pot_align::align_t al;
      double xi, xi_err;
      xi_reco::reconstruct( ev.proton_track_x[j]+al.x, ev.proton_track_side[j], ev.proton_track_pot[j], xi, xi_err );
      if ( xi < m_pot_limits[pot_id] ) continue;
      /*if ( ev.proton_track_side[j] == 0 && xi_45_rnd < m_pot_limits[pot_id] ) continue;
      if ( ev.proton_track_side[j] == 1 && xi_56_rnd < m_pot_limits[pot_id] ) continue;*/
      m_h_xireco[pot_id]->Fill( xi );
      if ( ev.proton_track_side[j] == 0 ) {
        m_h2_xicorr[pot_id]->Fill( xi, xi_45_rnd );
        xi_meas45[pot_id] = { xi, xi_err };
      }
      if ( ev.proton_track_side[j] == 1 ) {
        m_h2_xicorr[pot_id]->Fill( xi, xi_56_rnd );
        xi_meas56[pot_id] = { xi, xi_err };
      }
    }
    /*double m_diph = -999., y_diph = -999.;
    for ( unsigned short j = 0; j < ev.num_diphoton; ++j ) {
      if ( ev.diphoton_mass[j] > m_diph ) {
        m_diph = ev.diphoton_mass[j];
        y_diph = ev.diphoton_rapidity[j];
      }
    }*/
    for ( const auto& m45 : xi_meas45 ) {
      for ( const auto& m56 : xi_meas56 ) {
        diproton_candidate_t pp( m45.second.first, m45.second.second, m56.second.first, m56.second.second );
        const double m_pp = pp.mass(), err_m_pp = pp.mass_error();
        const double y_pp = pp.rapidity(), err_y_pp = pp.rapidity_error();
        h2_mcorr->Fill( m_pp, m_diph, weight );
        h2_ycorr->Fill( y_pp, y_diph, weight );
        if ( is_matched( 2.0, m_diph, m_pp, m_diph*rel_err_mass, err_m_pp )
          && is_matched( 2.0, y_diph, y_pp, y_diph*rel_err_rap, err_y_pp ) ) {
          num_match2d_2sigma += weight;
          num2_match2d_2sigma += weight*weight;
          h2_mcorr_2sigma->Fill( m_pp, m_diph, weight );
          h2_ycorr_2sigma->Fill( y_pp, y_diph, weight );
        }
        if ( is_matched( 3.0, m_diph, m_pp, m_diph*rel_err_mass, err_m_pp )
          && is_matched( 3.0, y_diph, y_pp, y_diph*rel_err_rap, err_y_pp ) ) {
          num_match2d_3sigma += weight;
          num2_match2d_3sigma += weight*weight;
          h2_mcorr_3sigma->Fill( m_pp, m_diph, weight );
          h2_ycorr_3sigma->Fill( y_pp, y_diph, weight );
        }
      }
    }
  }
  const double err_num_match2d_2sigma = sqrt( num2_match2d_2sigma );
  const double err_num_match2d_3sigma = sqrt( num2_match2d_3sigma );

  cout << "num matched at 2 sigma: " << num_match2d_2sigma << " +/- " << err_num_match2d_2sigma << endl;
  cout << "num matched at 3 sigma: " << num_match2d_3sigma << " +/- " << err_num_match2d_3sigma << endl;

  gStyle->SetOptStat( 0 );

  const string title = "CMS-TOTEM Pseudo-experiments";
  for ( const auto& p : m_pot_names ) {
    {
      Canvas c( Form( "toy_xi_corr_%s", p.second ), title.c_str() );
      m_h2_xicorr[p.first]->Draw( "colz" );
      c.Prettify( m_h2_xicorr[p.first] );
      c.Save( "pdf,png", OUTPUT_DIR );
    }
    {
      Canvas c( Form( "xi_reco_%s", p.second ), "CMS-TOTEM Preliminary 2016, #sqrt{s} = 13 TeV" );
      m_h_xireco[p.first]->Draw( "p" );
      m_h_xireco[p.first]->SetMarkerStyle( 24 );
      c.Prettify( m_h_xireco[p.first] );
      c.Save( "pdf,png", OUTPUT_DIR );
    }
  }
  auto diag = new TF1( "diag", "x" );
  unsigned short i = 0;
  for ( auto& h_mcorr : { h2_mcorr, h2_mcorr_2sigma, h2_mcorr_3sigma } ) {
    Canvas c( ( i == 0 ) ? "toy_m_corr" : Form( "toy_m_corr_%dsigma", i+1 ), title.c_str() );
    h_mcorr->Draw( "col" );
    c.Prettify( h_mcorr );
    //diag->SetLineColor( kGray );
    const double min_x = h_mcorr->GetXaxis()->GetXmin(), max_x = h_mcorr->GetXaxis()->GetXmax();
    diag->DrawF1( min_x, max_x, "same" );
    const double err = hypot( 0.055, rel_err_mass );
    auto diag_m1 = new TF1( "diag", "x*(1-[0])" ); diag_m1->SetParameter( 0, err );
    auto diag_p1 = new TF1( "diag", "x*(1+[0])" ); diag_p1->SetParameter( 0, err );
    auto diag_m2 = new TF1( "diag", "x*(1-2*[0])" ); diag_m2->SetParameter( 0, err );
    auto diag_p2 = new TF1( "diag", "x*(1+2*[0])" ); diag_p2->SetParameter( 0, err );
    auto diag_m3 = new TF1( "diag", "x*(1-3*[0])" ); diag_m3->SetParameter( 0, err );
    auto diag_p3 = new TF1( "diag", "x*(1+3*[0])" ); diag_p3->SetParameter( 0, err );
    diag_m1->DrawF1( min_x, max_x, "same" );
    diag_p1->DrawF1( min_x, max_x, "same" );
    diag_m2->DrawF1( min_x, max_x, "same" );
    diag_p2->DrawF1( min_x, max_x, "same" );
    diag_m3->DrawF1( min_x, max_x, "same" );
    diag_p3->DrawF1( min_x, max_x, "same" );
    diag_m1->SetLineColor( kGray+2 );
    diag_p1->SetLineColor( kGray+2 );
    diag_m2->SetLineColor( kGray+1 );
    diag_p2->SetLineColor( kGray+1 );
    diag_m3->SetLineColor( kGray+0 );
    diag_p3->SetLineColor( kGray+0 );
    diag_m2->SetLineStyle( 2 );
    diag_p2->SetLineStyle( 2 );
    diag_m3->SetLineStyle( 3 );
    diag_p3->SetLineStyle( 3 );
    c.SetLegendX1( 0.17 );
    c.SetLegendY1( 0.78 );
    c.AddLegendEntry( diag_m1->GetHistogram(), "#pm 1 #sigma", "l" );
    c.AddLegendEntry( diag_m2->GetHistogram(), "#pm 2 #sigma", "l" );
    c.AddLegendEntry( diag_m3->GetHistogram(), "#pm 3 #sigma", "l" );
    //diag_m1->SetLineColor( kGreen );
    //diag_p1->SetLineColor( kGreen );
    //diag_m2->SetLineColor( kYellow );
    //diag_p2->SetLineColor( kYellow );
    c.Save( "pdf,png", OUTPUT_DIR );
    ++i;
  }
  i = 0;
  for ( auto& h_ycorr : { h2_ycorr, h2_ycorr_2sigma, h2_ycorr_3sigma } ) {
    Canvas c( ( i == 0 ) ? "toy_y_corr" : "toy_y_corr_2sigma", title.c_str() );
    h_ycorr->Draw( "col" );
    c.Prettify( h_ycorr );
    diag->SetLineColor( kRed );
    double err = hypot( 0.055, rel_err_rap );
    const double min_x = h_ycorr->GetXaxis()->GetXmin(), max_x = h_ycorr->GetXaxis()->GetXmax();
    diag->DrawF1( min_x, max_x, "same" );
    auto diag_m1 = new TF1( "diag", "x+[0]" ); diag_m1->SetParameter( 0, err );
    auto diag_p1 = new TF1( "diag", "x-[0]" ); diag_p1->SetParameter( 0, err );
    auto diag_m2 = new TF1( "diag", "x+2*[0]" ); diag_m2->SetParameter( 0, err );
    auto diag_p2 = new TF1( "diag", "x-2*[0]" ); diag_p2->SetParameter( 0, err );
    auto diag_m3 = new TF1( "diag", "x+3*[0]" ); diag_m3->SetParameter( 0, err );
    auto diag_p3 = new TF1( "diag", "x-3*[0]" ); diag_p3->SetParameter( 0, err );
    diag_m1->DrawF1( min_x, max_x, "same" );
    diag_p1->DrawF1( min_x, max_x, "same" );
    diag_m2->DrawF1( min_x, max_x, "same" );
    diag_p2->DrawF1( min_x, max_x, "same" );
    diag_m3->DrawF1( min_x, max_x, "same" );
    diag_p3->DrawF1( min_x, max_x, "same" );
    diag_m1->SetLineColor( kGray+2 );
    diag_p1->SetLineColor( kGray+2 );
    diag_m2->SetLineColor( kGray+1 );
    diag_p2->SetLineColor( kGray+1 );
    diag_m3->SetLineColor( kGray );
    diag_p3->SetLineColor( kGray );
    diag_m2->SetLineStyle( 2 );
    diag_p2->SetLineStyle( 2 );
    diag_m3->SetLineStyle( 3 );
    diag_p3->SetLineStyle( 3 );
    c.SetLegendX1( 0.17 );
    c.SetLegendY1( 0.78 );
    c.AddLegendEntry( diag_m1->GetHistogram(), "#pm 1 #sigma", "l" );
    c.AddLegendEntry( diag_m2->GetHistogram(), "#pm 2 #sigma", "l" );
    c.AddLegendEntry( diag_m3->GetHistogram(), "#pm 3 #sigma", "l" );
    c.Save( "pdf,png", OUTPUT_DIR );
    ++i;
  }
}
