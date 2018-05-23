#include "xi_reconstruction.h"
#include "Canvas.h"

void validate_new_reco()
{
  xi_reco::load_file( "TreeProducer/data/optics_jun22.root" );
  xi_reco::load_optics_file( "TreeProducer/data/optics_17may22.root" );

  const map<unsigned short,const char*> pot_names = {
    { 2, "45 near" }, { 3, "45 far" }, { 102, "56 near" }, { 103, "56 far" }
  };
  map<unsigned short,TGraphErrors*> gr_xi_old, gr_xi_new;
  map<unsigned short,TGraph*> gr_rel_xi_old, gr_rel_xi_new;
  double xi, err_xi;

  for ( const auto& p : pot_names ) {
    gr_xi_old[p.first] = new TGraphErrors;
    gr_xi_new[p.first] = new TGraphErrors;
    gr_rel_xi_old[p.first] = new TGraph;
    gr_rel_xi_new[p.first] = new TGraph;
    for ( double x = 0.; x <= 0.025; x += 5.e-6 ) {
      xi_reco::reconstruct_old( x, p.first/100, p.first%100, xi, err_xi );
      const unsigned short id = gr_xi_old[p.first]->GetN();
      gr_xi_old[p.first]->SetPoint( id, x, xi );
      gr_xi_old[p.first]->SetPointError( id, 0., err_xi );
      if ( xi > 0.02 && xi < 0.15 )
        gr_rel_xi_old[p.first]->SetPoint( gr_rel_xi_old[p.first]->GetN(), xi, err_xi/xi*100. );

      xi_reco::reconstruct( x, p.first/100, p.first%100, xi, err_xi );
      gr_xi_new[p.first]->SetPoint( id, x, xi );
      gr_xi_new[p.first]->SetPointError( id, 0., err_xi );
      if ( xi > 0.02 && xi < 0.15 )
        gr_rel_xi_new[p.first]->SetPoint( gr_rel_xi_new[p.first]->GetN(), xi, err_xi/xi*100. );
    }
    {
      Canvas c( Form( "xi_reco_valid_%d", p.first ), Form( "%s strips - TOTEM parameterisation 2016", p.second ) );
      TMultiGraph mg;
      gr_xi_new[p.first]->SetLineWidth( 4 );
      gr_xi_new[p.first]->SetLineColor( kRed+1 );
      gr_xi_new[p.first]->SetFillColorAlpha( kRed+1, 0.3 );
      //gr_xi_new[p.first]->SetFillStyle( 3004 );
      mg.Add( gr_xi_new[p.first] );
      gr_xi_old[p.first]->SetLineWidth( 4 );
      gr_xi_old[p.first]->SetFillColor( kBlack );
      gr_xi_old[p.first]->SetFillStyle( 3004 );
      mg.Add( gr_xi_old[p.first] );
      mg.SetTitle( ";x (m);#xi" );
      mg.Draw( "al3" );
      c.SetLegendX1( 0.18 );
      c.AddLegendEntry( gr_xi_old[p.first], "Previous approach", "lf" );
      c.AddLegendEntry( gr_xi_new[p.first], "All effects accounted for", "lf" );
      c.Prettify( mg.GetHistogram() );
      c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/xi_study/" );
    }
    {
      Canvas c( Form( "xi_reco_reldxi_valid_%d", p.first ), Form( "%s strips - TOTEM parameterisation 2016", p.second ) );
      TMultiGraph mg;
      gr_rel_xi_old[p.first]->SetLineWidth( 4 );
      mg.Add( gr_rel_xi_old[p.first] );
      gr_rel_xi_new[p.first]->SetLineWidth( 4 );
      gr_rel_xi_new[p.first]->SetLineColor( kRed+1 );
      mg.Add( gr_rel_xi_new[p.first] );
      mg.SetTitle( ";#xi;#delta#xi/#xi (%)" );
      mg.Draw( "al" );
      mg.GetXaxis()->SetRangeUser( 0.02, 0.15 );
      mg.GetYaxis()->SetRangeUser( 0., 15. );
      c.SetLegendX1( 0.4 );
      c.AddLegendEntry( gr_rel_xi_old[p.first], "Previous approach", "l" );
      c.AddLegendEntry( gr_rel_xi_new[p.first], "All effects accounted for", "l" );
      c.Prettify( mg.GetHistogram() );
      c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/xi_study/" );
    }
  }

}

