#include "Canvas.h"
#include "pot_alignment.h"
#include "xi_reconstruction.h"
#include "diproton_candidate.h"

#include "TFile.h"
#include "TTree.h"
#include "TF1.h"

#define MAX_PROTONS 30

void xi_accept()
{
  TFile f( "Samples/output_Run2016BCG_looseCuts_28jun.root" );
  TTree* tr = dynamic_cast<TTree*>( f.Get( "ntp" ) );

  xi_reco::load_file( "TreeProducer/data/optics_jun22.root" );
  pot_align::load_file( "TreeProducer/data/alignment_collection_v2.out" );

  const float num_sigma = 2.0;

  map<unsigned short,string> pot_names = { { 2, "45N" }, { 3, "45F" }, { 102, "56N" }, { 103, "56F" } };
  map<unsigned short,float> xi_accept = { { 2, 0.033 }, { 3, 0.024 }, { 102, 0.050 }, { 103, 0.037 } };
  map<unsigned short,float> xi_accept_new = { { 2, 0.033 }, { 3, 0.024 }, { 102, 0.040 }, { 103, 0.032 } };

//  const vector<unsigned short> fills_to_study = { 4947, 4985, 5017, 5030 };

  unsigned int fill_number;
  tr->SetBranchAddress( "fill_number", &fill_number );

  unsigned int num_proton_track;
  float proton_track_x[MAX_PROTONS], proton_track_y[MAX_PROTONS];
  unsigned int proton_track_side[MAX_PROTONS], proton_track_pot[MAX_PROTONS];
  float proton_track_normchi2[MAX_PROTONS];
  tr->SetBranchAddress( "num_proton_track", &num_proton_track );
  tr->SetBranchAddress( "proton_track_x", proton_track_x );
  tr->SetBranchAddress( "proton_track_y", proton_track_y );
  tr->SetBranchAddress( "proton_track_side", proton_track_side );
  tr->SetBranchAddress( "proton_track_pot", proton_track_pot );
  tr->SetBranchAddress( "proton_track_normchi2", proton_track_normchi2 );

  typedef map<unsigned short,TH1D*> plots_t;
  map<unsigned short,plots_t> plots_vs_fill;
  plots_t combined_plots;
  for ( const auto& pot : pot_names ) {
    combined_plots[pot.first] = new TH1D( Form( "xi_spectrum_%s_combined", pot.second.c_str() ), Form( "#xi(%s)@@Events", pot.second.c_str() ), 800, 0., 0.4 );
  }

  for ( unsigned long long i=0; i<tr->GetEntriesFast(); i++ ) {
    tr->GetEntry( i );

    auto align = pot_align::get_alignments( fill_number );
    if ( plots_vs_fill.find( fill_number )==plots_vs_fill.end() ) { // fill info not yet created
      // create the plots for this fill
      plots_t plots;
      for ( const auto& pot : pot_names ) {
        plots[pot.first] = new TH1D( Form( "xi_spectrum_%s_fill%d", pot.second.c_str(), fill_number ), Form( "#xi(%s)@@Events", pot.second.c_str() ), 100, 0., 0.25 );
      }
      plots_vs_fill.insert( make_pair( fill_number, plots ) );
    }
    // retrieve the plots for this fill
    auto& fill_plots = plots_vs_fill[fill_number];

    for ( unsigned short j=0; j<num_proton_track; j++ ) {
      double xi, xi_err;
      const unsigned short pot_id = proton_track_side[j]*100+proton_track_pot[j];
      const auto& al = align[pot_id];
      xi_reco::reconstruct( proton_track_x[j]+al.x, proton_track_side[j], proton_track_pot[j], xi, xi_err );

      //float xi_cut = xi_accept[pot_id];
      fill_plots[pot_id]->Fill( xi );
      combined_plots[pot_id]->Fill( xi );
      /*if ( xi>=xi_cut ) {
      }
      else { // outside pot acceptance
      }*/
    }
  }

  //----- plotting part

  unsigned short pot_id = 0;
  Canvas c_combined( "xiaccept_limitstudy_allpots", "CMS+TOTEM Preliminary 2016, #sqrt{s} = 13 TeV, L = 9.4 fb^{-1}" );
  c_combined.Divide( 2, 2 );
  c_combined.SetLegendX1( 0.18 );
  c_combined.SetLegendY1( 0.82 );

  const float min_frac = 0.01;

  for ( const auto& pot : xi_accept ) {
    const char* pot_name = pot_names[pot.first].c_str();
    /*{
      Canvas c( Form( "xiaccept_%s", pot_name ), "CMS+TOTEM Preliminary 2016, #sqrt{s} = 13 TeV, L = 9.4 fb^{-1}" );
      THStack st;
      unsigned short i = 0;
      //for ( auto& fill : plots_vs_fill ) {
      for ( auto& fill : fills_to_study ) {
        if ( i>10 ) break; //FIXME
        //auto plot = fill.second[pot.first];
        auto plot = plots_vs_fill[fill][pot.first];
        plot->SetLineWidth( 2 );
        plot->SetLineColor( i+1 );
        st.Add( plot );
        //c.AddLegendEntry( plot, Form( "Fill %d", fill.first ) );
        c.AddLegendEntry( plot, Form( "Fill %d", fill ) );
        if ( i==0 ) st.SetTitle( plot->GetTitle() );
        i++;
      }
      st.Draw( "nostack" );
      c.Prettify( st.GetHistogram() );
      c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
    }*/
    {
      c_combined.cd( 1+pot_id );
      TPad* pad = (TPad*)c_combined.GetPad( 1+pot_id );
      cout << "new plot for " << pot_name << endl;
      double y_max = 1.3;
      TGraph* gr_accept = new TGraph; // brilliant ROOT, brilliant...
      unsigned short i = 0;
      for ( double xi_cut = 0.; xi_cut<=0.16 ; xi_cut+=0.01 ) {
        TH1D* plot = combined_plots[pot.first];
        double num_inside = plot->Integral( plot->GetXaxis()->FindBin( xi_cut ), plot->GetXaxis()->FindBin( plot->GetXaxis()->GetXmax() )+1 );
        double num_total = plot->Integral();
        gr_accept->SetPoint( i, xi_cut, 1.-num_inside/num_total );
        i++;
      }
      float xi_cut = 0.0, dist_goal = 99.;
      for ( float xi=0.01; xi<0.1; xi+=0.001 ) {
        double frac = gr_accept->Eval( xi );
        if ( fabs( frac-min_frac )<dist_goal ) {
          dist_goal = fabs( frac-min_frac );
          xi_cut = xi;
        }
      }
      gr_accept->SetPoint( i, 0.1, 0. );
      gr_accept->Draw( "acf" );
      gr_accept->GetXaxis()->SetTitle( Form( "Lower cut on #xi(%s)", pot_name ) );
      gr_accept->GetYaxis()->SetTitle( "Ev.fraction outside pot accept." );
      gr_accept->GetXaxis()->SetRangeUser( 0., 0.15 );
      gr_accept->GetYaxis()->SetRangeUser( 0., y_max );
      gr_accept->SetFillColor( kBlack );
      gr_accept->SetFillStyle( 3004 );
      gr_accept->Draw( "c" );
      c_combined.Prettify( gr_accept->GetHistogram() );
      auto cut_new = new TLine( xi_cut, 0., xi_cut, y_max ); // damn ROOT, damn...
      cut_new->SetLineWidth( 3 );
      cut_new->SetLineColor( kBlue-2 );
      cut_new->Draw( "same" );
      /*TLine* cut_new_horiz = new TLine( 0., min_frac, xi_cut, min_frac ); // damn ROOT, damn...
      cut_new_horiz->SetLineWidth( 2 );
      cut_new_horiz->SetLineColor( kBlue-2 );
      cut_new_horiz->Draw( "same" );*/
      auto cut_orig = new TLine( pot.second, 0., pot.second, y_max ); // damn ROOT, damn...
      cut_orig->SetLineWidth( 3 );
      cut_orig->SetLineStyle( 2 );
      cut_orig->SetLineColor( kRed+1 );
      cut_orig->Draw( "same" );
      if ( pot_id==0 ) {
        c_combined.AddLegendEntry( cut_orig, "Predicted from optics", "l" );
        c_combined.AddLegendEntry( cut_new, "Adapted cuts", "l" );
      }
      pad->SetLogy();
    }
    pot_id++;
  }
  TLegend* leg = c_combined.GetLegend();
  if ( leg ) {
    leg->SetFillStyle( 0 );
    leg->SetLineWidth( 0 );
    leg->SetTextSize( 0.065 );
  }
  c_combined.Save( "png,pdf", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
}

