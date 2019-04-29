#include "Canvas.h"
#include "pot_alignment.h"
#include "xi_reconstruction.h"
#include "diproton_candidate.h"
#include "DiphotonAnalyzer/TreeProducer/interface/ProtonInfoEvent.h"

#include "TH2.h"
#include "TF1.h"
#include "TStyle.h"

//#define OUTPUT_DIR "/afs/cern.ch/user/l/lforthom/www/private/twophoton/test"
#define OUTPUT_DIR "/afs/cern.ch/user/l/lforthom/www/private/twophoton/test_15may"

void mix_events( const char* input_filename = "fits_results.root", bool scan = false, unsigned long long num_toys = 5000000 )
{
//  const double exp_yield = 586.877;
//  const double exp_yield = 525;
  const double rel_err_mass = 0.02, rel_err_rap = 0.074;
//  const double exp_yield = 566.994;
//  const double exp_yield = 112.653; // xicomp
//  const double exp_yield = 3.26289; // xitight (MC)
  double exp_yield = 2.62498; // xitight (MC)
  exp_yield /= (1.-0.5277); // fraction of events with at least 1 proton track

  auto in_file = TFile::Open( input_filename );
  auto mix_45 = (TF1*)in_file->Get( "fit_xip" )->Clone(), mix_56 = (TF1*)in_file->Get( "fit_xim" )->Clone();
  auto fit_1sig_45 = (TGraph*)in_file->Get( "fit_xip_1sig" )->Clone(), fit_1sig_56 = (TGraph*)in_file->Get( "fit_xim_1sig" )->Clone();
  delete in_file;

  //TFile f( "/eos/cms/store/user/lforthom/ProtonTree/DoubleEG/proton_ntuple-Run2016BCG_94Xrereco_v1.root" );
  TFile f( "proton_ntuple-Run2016BCG_94Xrereco_v1.root" );
  auto tree = dynamic_cast<TTree*>( f.Get( "protonTreeProducer/ntp" ) );
  ProtonInfoEvent ev;
  ev.attach( tree, {
    "fill_number",
    "num_fwd_track", "fwd_track_arm", "fwd_track_pot", "fwd_track_x"
  } );
  const unsigned long long num_events = tree->GetEntriesFast();
  //xi_reco::load_file( "TreeProducer/data/optics_jun22.root" );
  xi_reco::load_optics_file( "TreeProducer/data/optics_17may22.root" );
  pot_align::load_file( "TreeProducer/data/alignment_collection_v2.out" );

  map<unsigned short,TH2D*> m_h2_xicorr;
  map<unsigned short,TH1D*> m_h_xireco;
  map<unsigned short,const char*> m_pot_names = { { 2, "45N" }, { 3, "45F" }, { 102, "56N" }, { 103, "56F" } };
  map<unsigned short,double> m_pot_limits = { { 2, 0.067 }, { 3, 0.066 }, { 102, 0.070 }, { 103, 0.061 } }; // < 10% radiation damage
  //map<unsigned short,double> m_pot_limits = { { 2, 0.034 }, { 3, 0.023 }, { 102, 0.042 }, { 103, 0.032 } };
  const double max_xigg = 0.15;
  for ( const auto& p : m_pot_names ) {
    m_h2_xicorr[p.first] = new TH2D( Form( "xi_corr_%d", p.first ), Form( ";#xi_{%s};#xi_{#gamma#gamma}", p.second ), 50, 0., 0.2, 50, 0., 0.2 );
    m_h_xireco[p.first] = new TH1D( Form( "xi_reco_%d", p.first ), Form( ";#xi_{%s};Events", p.second ), 50, 0., 0.2 );
  }
  //auto h2_mcorr = new TH2D( "mcorr", ";m_{pp} (GeV);m_{#gamma#gamma} (GeV)", 40, 250., 1850., 40, 250., 1850. ),
  auto h2_mcorr = new TH2D( "mcorr", ";m_{pp} (GeV);m_{#gamma#gamma} (GeV)", 80, 450., 2050., 80, 450., 2050. ),
       h2_mcorr_2sigma = (TH2D*)h2_mcorr->Clone( "mcorr_2sigma" ),
       h2_mcorr_3sigma = (TH2D*)h2_mcorr->Clone( "mcorr_3sigma" );
  //auto h2_ycorr = new TH2D( "ycorr", ";y_{pp};y_{#gamma#gamma}", 50, -2., 2., 50, -2., 2. ),
  auto h2_ycorr = new TH2D( "ycorr", ";y_{pp};y_{#gamma#gamma}", 100, -1., 1., 100, -1., 1. ),
       h2_ycorr_2sigma = (TH2D*)h2_ycorr->Clone( "ycorr_2sigma" ),
       h2_ycorr_3sigma = (TH2D*)h2_ycorr->Clone( "ycorr_3sigma" );

  cout << "number of toys: " << num_toys << endl;

  //const double weight = exp_yield/num_toys/mean_int; //FIXME
  const double weight = exp_yield/num_toys; //FIXME
  cout << "weight=" << weight << endl;

  unsigned long long idx = TMath::Max( 0ll, (long long)( rand()*1./RAND_MAX*num_events-num_toys ) );

  unsigned short num_try = 1;
  if ( scan )
    num_try = fit_1sig_45->GetN();

  double max_nomatch = 0., err_max_nomatch = 0., max_match2d_2sigma = 0., err_max_match2d_2sigma = 0., max_match2d_3sigma = 0., err_max_match2d_3sigma = 0.;
  double min_nomatch = 9999., err_min_nomatch = 9999., min_match2d_2sigma = 9999., err_min_match2d_2sigma = 9999., min_match2d_3sigma = 9999., err_min_match2d_3sigma = 9999.;
  for ( unsigned short l = 0; l < num_try; ++l ) {
    double num_nomatch = 0., num2_nomatch = 0.;
    double num_match2d_2sigma = 0., num2_match2d_2sigma = 0.;
    double num_match2d_3sigma = 0., num2_match2d_3sigma = 0.;

    if ( scan ) {
      cout << "scanning point #" << l << endl;
      mix_45->SetParameter( 0, fit_1sig_45->GetX()[l] );
      mix_45->SetParameter( 1, fit_1sig_45->GetY()[l] );
      mix_56->SetParameter( 0, fit_1sig_56->GetX()[l] );
      mix_56->SetParameter( 1, fit_1sig_56->GetY()[l] );
    }
    for ( unsigned long long i = 0; i < num_toys; ++i ) {
      //tree->GetEntry( rand()*1./RAND_MAX*num_events );
      tree->GetEntry( idx++ );
      if ( !scan && fmod( i*1., num_toys*0.1 ) == 0 )
        cout << "event " << i << endl;

      const double xi_45_rnd = mix_45->GetRandom( min( m_pot_limits[  2], m_pot_limits[  3] ), max_xigg );
      const double xi_56_rnd = mix_56->GetRandom( min( m_pot_limits[102], m_pot_limits[103] ), max_xigg );

      const double m_diph = 13.e3*sqrt( xi_45_rnd*xi_56_rnd );
      const double y_diph = 0.5*log( xi_45_rnd/xi_56_rnd );

      auto align = pot_align::get_alignments( ev.fill_number );

      map<unsigned short,pair<double,double> > xi_meas45, xi_meas56;
      //cout << ev.num_fwd_track << endl;
      for ( unsigned short j = 0; j < ev.num_fwd_track; ++j ) {
        const unsigned short pot_id = 100*ev.fwd_track_arm[j]+ev.fwd_track_pot[j];
        auto al = align[pot_id];
//        pot_align::align_t al;
        double xi, xi_err;
        xi_reco::reconstruct( ev.fwd_track_x[j]+al.x, ev.fwd_track_arm[j], ev.fwd_track_pot[j], xi, xi_err );
        if ( xi < m_pot_limits[pot_id] || xi > 0.15 ) continue; //FIXME
        m_h_xireco[pot_id]->Fill( xi );
        if ( ev.fwd_track_arm[j] == 0 ) {
          m_h2_xicorr[pot_id]->Fill( xi, xi_45_rnd );
          xi_meas45[pot_id] = { xi, xi_err };
        }
        if ( ev.fwd_track_arm[j] == 1 ) {
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
          const bool m_matched_2sig = is_matched( 2., m_diph, m_pp, m_diph*rel_err_mass, err_m_pp );
          const bool m_matched_3sig = is_matched( 3., m_diph, m_pp, m_diph*rel_err_mass, err_m_pp );
          const bool y_matched_2sig = is_matched( 2., y_diph, y_pp, rel_err_rap, err_y_pp );
          const bool y_matched_3sig = is_matched( 3., y_diph, y_pp, rel_err_rap, err_y_pp );
          num_nomatch += weight;
          num2_nomatch += weight*weight;
          h2_mcorr->Fill( m_pp, m_diph, weight );
          h2_ycorr->Fill( y_pp, y_diph, weight );
          if ( m_matched_2sig && y_matched_2sig ) {
            num_match2d_2sigma += weight;
            num2_match2d_2sigma += weight*weight;
            h2_mcorr_2sigma->Fill( m_pp, m_diph, weight );
            h2_ycorr_2sigma->Fill( y_pp, y_diph, weight );
          }
          if ( m_matched_3sig && y_matched_3sig ) {
            num_match2d_3sigma += weight;
            num2_match2d_3sigma += weight*weight;
            h2_mcorr_3sigma->Fill( m_pp, m_diph, weight );
            h2_ycorr_3sigma->Fill( y_pp, y_diph, weight );
          }
        }
      }
    }
    const double err_num_nomatch = sqrt( num2_nomatch );
    const double err_num_match2d_2sigma = sqrt( num2_match2d_2sigma );
    const double err_num_match2d_3sigma = sqrt( num2_match2d_3sigma );

    cout << "num unmatched: " << num_nomatch << " +/- " << err_num_nomatch << endl;
    cout << "num matched at 2 sigma: " << num_match2d_2sigma << " +/- " << err_num_match2d_2sigma << endl;
    cout << "num matched at 3 sigma: " << num_match2d_3sigma << " +/- " << err_num_match2d_3sigma << endl;
    if ( max_nomatch < num_nomatch ) {
      max_nomatch = num_nomatch;
      err_max_nomatch = err_num_nomatch;
    }
    if ( min_nomatch > num_nomatch ) {
      min_nomatch = num_nomatch;
      err_min_nomatch = err_num_nomatch;
    }
    if ( max_match2d_2sigma < num_match2d_2sigma ) {
      max_match2d_2sigma = num_match2d_2sigma;
      err_max_match2d_2sigma = err_num_match2d_2sigma;
    }
    if ( min_match2d_2sigma > num_match2d_2sigma ) {
      min_match2d_2sigma = num_match2d_2sigma;
      err_min_match2d_2sigma = err_num_match2d_2sigma;
    }
    if ( max_match2d_3sigma < num_match2d_3sigma ) {
      max_match2d_3sigma = num_match2d_3sigma;
      err_max_match2d_3sigma = err_num_match2d_3sigma;
    }
    if ( min_match2d_3sigma > num_match2d_3sigma ) {
      min_match2d_3sigma = num_match2d_3sigma;
      err_min_match2d_3sigma = err_num_match2d_3sigma;
    }
  }
  cout << "min num unmatched: " << min_nomatch << " +/- " << err_min_nomatch << endl;
  cout << "min num matched at 2 sigma: " << min_match2d_2sigma << " +/- " << err_min_match2d_2sigma << endl;
  cout << "min num matched at 3 sigma: " << min_match2d_3sigma << " +/- " << err_min_match2d_3sigma << endl;
  cout << "max num unmatched: " << max_nomatch << " +/- " << err_max_nomatch << endl;
  cout << "max num matched at 2 sigma: " << max_match2d_2sigma << " +/- " << err_max_match2d_2sigma << endl;
  cout << "max num matched at 3 sigma: " << max_match2d_3sigma << " +/- " << err_max_match2d_3sigma << endl;
  gStyle->SetOptStat( 0 );

  const string title = "2016 conditions (13 TeV)";
  for ( const auto& p : m_pot_names ) {
    {
      Canvas c( Form( "toy_xi_corr_%s", p.second ), title.c_str() );
      m_h2_xicorr[p.first]->Draw( "colz" );
      c.Prettify( m_h2_xicorr[p.first] );
      c.Save( "pdf,png", OUTPUT_DIR );
    }
    {
      Canvas c( Form( "xi_reco_%s", p.second ), "9.4 fb^{-1} (13 TeV)", "Preliminary" );
      m_h_xireco[p.first]->Draw( "p" );
      m_h_xireco[p.first]->SetMarkerStyle( 24 );
      c.Prettify( m_h_xireco[p.first] );
      c.Save( "pdf,png", OUTPUT_DIR );
    }
  }
  auto diag = new TF1( "diag", "x" );
  unsigned short i = 0;
  for ( auto& h_mcorr : { h2_mcorr, h2_mcorr_2sigma, h2_mcorr_3sigma } ) {
    Canvas c( ( i == 0 ) ? "toy_m_corr" : Form( "toy_m_corr_%dsigma", i+1 ), title.c_str(), "Pseudo-experiments", false, Canvas::Align::right );
    c.Divide( 2, 2 );
    c.cd( 3 );
    h_mcorr->Draw( "col" );
    c.Prettify( h_mcorr );
    h_mcorr->GetXaxis()->SetLabelSize( 16 );
    h_mcorr->GetYaxis()->SetLabelSize( 16 );
    //diag->SetLineColor( kGray );
    const double min_x = h_mcorr->GetXaxis()->GetXmin(), max_x = h_mcorr->GetXaxis()->GetXmax();
    diag->DrawF1( min_x, max_x, "same" );
    const double err = hypot( 0.055, rel_err_mass );
    auto diag_m1 = new TF1( "diag", "x*(1-  [0])" ); diag_m1->SetParameter( 0, err ); diag_m1->DrawF1( min_x, max_x, "same" ); diag_m1->SetLineColor( kBlack );
    auto diag_p1 = new TF1( "diag", "x*(1+  [0])" ); diag_p1->SetParameter( 0, err ); diag_p1->DrawF1( min_x, max_x, "same" ); diag_p1->SetLineColor( kBlack );
    auto diag_m2 = new TF1( "diag", "x*(1-2*[0])" ); diag_m2->SetParameter( 0, err ); diag_m2->DrawF1( min_x, max_x, "same" ); diag_m2->SetLineColor( kBlack ); diag_m2->SetLineStyle( 2 );
    auto diag_p2 = new TF1( "diag", "x*(1+2*[0])" ); diag_p2->SetParameter( 0, err ); diag_p2->DrawF1( min_x, max_x, "same" ); diag_p2->SetLineColor( kBlack ); diag_p2->SetLineStyle( 2 );
    auto diag_m3 = new TF1( "diag", "x*(1-3*[0])" ); diag_m3->SetParameter( 0, err ); diag_m3->DrawF1( min_x, max_x, "same" ); diag_m3->SetLineColor( kBlack ); diag_m3->SetLineStyle( 3 );
    auto diag_p3 = new TF1( "diag", "x*(1+3*[0])" ); diag_p3->SetParameter( 0, err ); diag_p3->DrawF1( min_x, max_x, "same" ); diag_p3->SetLineColor( kBlack ); diag_p3->SetLineStyle( 3 );
    c.cd( 4 );
    auto py = h_mcorr->ProjectionY();
    py->Draw( "hist,hbar" );
    py->SetFillColorAlpha( kBlack, 0.3 );
    py->GetYaxis()->SetLabelFont( h_mcorr->GetXaxis()->GetLabelFont() );
    py->GetYaxis()->SetLabelSize( h_mcorr->GetXaxis()->GetLabelSize() );
    c.cd( 1 );
    auto px = h_mcorr->ProjectionX();
    px->Draw( "hist,bar" );
    px->SetFillColorAlpha( kBlack, 0.3 );
    px->GetYaxis()->SetLabelFont( h_mcorr->GetYaxis()->GetLabelFont() );
    px->GetYaxis()->SetLabelSize( h_mcorr->GetYaxis()->GetLabelSize() );
    c.cd();
    c.SetLegendX1( 0.67 );
    c.SetLegendY1( 0.63 );
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
    Canvas c( ( i == 0 ) ? "toy_y_corr" : Form( "toy_y_corr_%dsigma", i+1 ), title.c_str(), "Pseudo-experiments", false, Canvas::Align::right );
    c.Divide( 2, 2 );
    c.cd( 3 );
    h_ycorr->Draw( "col" );
    c.Prettify( h_ycorr );
    h_ycorr->GetXaxis()->SetLabelSize( 16 );
    h_ycorr->GetYaxis()->SetLabelSize( 16 );
    diag->SetLineColor( kRed );
    double err = hypot( 0.055, rel_err_rap );
    const double min_x = h_ycorr->GetXaxis()->GetXmin(), max_x = h_ycorr->GetXaxis()->GetXmax();
    diag->DrawF1( min_x, max_x, "same" );
    auto diag_m1 = new TF1( "diag", "x+  [0]" ); diag_m1->SetParameter( 0, err ); diag_m1->DrawF1( min_x, max_x, "same" ); diag_m1->SetLineColor( kBlack );
    auto diag_p1 = new TF1( "diag", "x-  [0]" ); diag_p1->SetParameter( 0, err ); diag_p1->DrawF1( min_x, max_x, "same" ); diag_p1->SetLineColor( kBlack );
    auto diag_m2 = new TF1( "diag", "x+2*[0]" ); diag_m2->SetParameter( 0, err ); diag_m2->DrawF1( min_x, max_x, "same" ); diag_m2->SetLineColor( kBlack ); diag_m2->SetLineStyle( 2 );
    auto diag_p2 = new TF1( "diag", "x-2*[0]" ); diag_p2->SetParameter( 0, err ); diag_p2->DrawF1( min_x, max_x, "same" ); diag_p2->SetLineColor( kBlack ); diag_p2->SetLineStyle( 2 );
    auto diag_m3 = new TF1( "diag", "x+3*[0]" ); diag_m3->SetParameter( 0, err ); diag_m3->DrawF1( min_x, max_x, "same" ); diag_m3->SetLineColor( kBlack ); diag_m3->SetLineStyle( 3 );
    auto diag_p3 = new TF1( "diag", "x-3*[0]" ); diag_p3->SetParameter( 0, err ); diag_p3->DrawF1( min_x, max_x, "same" ); diag_p3->SetLineColor( kBlack ); diag_p3->SetLineStyle( 3 );
    c.cd( 4 );
    auto py = h_ycorr->ProjectionY( "_py", 0, -1, "e" );
    py->Draw( "hist,hbar" );
    py->SetFillColorAlpha( kBlack, 0.3 );
    py->GetYaxis()->SetLabelFont( h_ycorr->GetXaxis()->GetLabelFont() );
    py->GetYaxis()->SetLabelSize( h_ycorr->GetXaxis()->GetLabelSize() );
    /*py->GetYaxis()->SetTitle( "Pseudo-events" );
    py->GetYaxis()->SetTitleFont( h_ycorr->GetXaxis()->GetTitleFont() );
    py->GetYaxis()->SetTitleSize( h_ycorr->GetXaxis()->GetTitleSize() );*/
    c.cd( 1 );
    auto px = h_ycorr->ProjectionX( "_px", 0, -1, "e" );
    px->Draw( "hist,bar" );
    px->SetFillColorAlpha( kBlack, 0.3 );
    px->GetYaxis()->SetLabelFont( h_ycorr->GetYaxis()->GetLabelFont() );
    px->GetYaxis()->SetLabelSize( h_ycorr->GetYaxis()->GetLabelSize() );
    //px->GetYaxis()->SetTitle( "Pseudo-events" );
    c.cd();
    /*c.SetLegendX1( 0.17 );
    c.SetLegendY1( 0.78 );*/
    c.SetLegendX1( 0.67 );
    c.SetLegendY1( 0.63 );
    c.AddLegendEntry( diag_m1->GetHistogram(), "#pm 1 #sigma", "l" );
    c.AddLegendEntry( diag_m2->GetHistogram(), "#pm 2 #sigma", "l" );
    c.AddLegendEntry( diag_m3->GetHistogram(), "#pm 3 #sigma", "l" );
    c.Save( "pdf,png", OUTPUT_DIR );
    ++i;
  }
}
