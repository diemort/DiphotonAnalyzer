#include "TFile.h"
#include "TSpline.h"
#include <map>
#include "Canvas.h"

void xi_resol()
{
  std::map<unsigned short,TGraph*> g_x_vs_xi, g_Lx_vs_xi, g_vx_vs_xi, g_xi_vs_x;
  std::map<unsigned short,double> disp_rel_err;
  std::map<unsigned short,const char*> m_pots = {
    { 2, "45-near" },
    { 3, "45-far" },
    { 102, "56-near" },
    { 103, "56-far" }
  };

  TFile file( "TreeProducer/data/optics_17may22.root" );
  std::map<unsigned short,const char*> paths = {
    { 2, "ctppsPlotOpticalFunctions_45/ip5_to_station_150_h_1_lhcb2" }, // 45N
    { 3, "ctppsPlotOpticalFunctions_45/ip5_to_station_150_h_2_lhcb2" }, // 45F
    { 102, "ctppsPlotOpticalFunctions_56/ip5_to_station_150_h_1_lhcb1" }, // 56N
    { 103, "ctppsPlotOpticalFunctions_56/ip5_to_station_150_h_2_lhcb1" } // 56F
  };
  for ( const auto& p : paths ) {
    g_x_vs_xi[p.first] = dynamic_cast<TGraph*>( file.Get( Form( "%s/g_x0_vs_xi", p.second ) )->Clone() );
    g_Lx_vs_xi[p.first] = dynamic_cast<TGraph*>( file.Get( Form( "%s/g_L_x_vs_xi", p.second ) )->Clone() );
    g_vx_vs_xi[p.first] = dynamic_cast<TGraph*>( file.Get( Form( "%s/g_v_x_vs_xi", p.second ) )->Clone() );
    g_xi_vs_x[p.first] = dynamic_cast<TGraph*>( file.Get( Form( "%s/g_xi_vs_xso", p.second ) )->Clone() );
  };
  disp_rel_err = std::map<unsigned short,double>{
    { 2, 0.0387 }, { 3, 0.0435 }, { 102, 0.0526 }, { 103, 0.0614 }
  };

  // a few parameters
  const double si_x_alignment = 150.0e-6; // in m, alignment uncertainty
  const double si_vertex = 10.e-6;
  const double si_beam_div = 20.e-6;
  const double si_sensor_resol = 12.e-6;
  const double si_beam_mom = 1.e-3;
  const double de_xi = 1.e-3;

  const double min_xi = 0.02, max_xi = 0.15;

  for ( const auto& pot : m_pots ) {
    unsigned short i = 0;
    TGraph g_dx, g_vx, g_thx, g_align, g_sensor, g_beammom, g_tot;

    for ( double xi = min_xi; xi <= max_xi; xi += 1.e-4 ) {
      const double de_x = g_x_vs_xi.at( pot.first )->Eval( xi+de_xi/2. ) - g_x_vs_xi.at( pot.first )->Eval( xi-de_xi/2. );
      const double D = de_x / de_xi;  // m
      const double Lx = g_Lx_vs_xi.at( pot.first )->Eval( xi ); // m
      const double vx = g_vx_vs_xi.at( pot.first )->Eval( xi );
      const double si_de_disp = disp_rel_err.at( pot.first );

      // determine uncertainty of xi
      const double si_D_x = si_de_disp * xi;
      const double si_sensor = si_sensor_resol / fabs( D );
      const double si_alig = si_x_alignment / fabs( D );
      const double si_term_th_x = fabs( Lx/D ) * si_beam_div;
      const double si_term_vtx_x = fabs( vx/D ) * si_vertex;

      g_dx.SetPoint( i, xi, si_D_x / xi * 100 );
      g_sensor.SetPoint( i, xi, si_sensor / xi * 100 );
      g_align.SetPoint( i, xi, si_alig / xi * 100 );
      g_thx.SetPoint( i, xi, si_term_th_x / xi * 100 );
      g_vx.SetPoint( i, xi, si_term_vtx_x / xi * 100 );
      g_beammom.SetPoint( i, xi, si_beam_mom / xi * 100 );

      double xi_err = 0.;
      for ( const auto &c : { si_beam_mom, si_D_x, si_sensor, si_alig, si_term_th_x, si_term_vtx_x } )
        xi_err += c*c;
      xi_err = sqrt( xi_err );

      g_tot.SetPoint( i, xi, xi_err / xi * 100 );

      ++i;
    }
    Canvas c( Form( "rel_dxi_%d", pot.first ), "PPS optics parameterisation 2016" );
    TMultiGraph mg;

    g_dx.SetLineColor( kRed+1 );
    g_dx.SetLineWidth( 3 );
    g_dx.SetLineStyle( 2 );
    mg.Add( &g_dx );

    g_sensor.SetLineColor( kBlue+1 );
    g_sensor.SetLineWidth( 3 );
    g_sensor.SetLineStyle( 3 );
    mg.Add( &g_sensor );

    g_align.SetLineColor( kGreen-2 );
    g_align.SetLineWidth( 3 );
    g_align.SetLineStyle( 4 );
    mg.Add( &g_align );

    g_vx.SetLineColor( kOrange+1 );
    g_vx.SetLineWidth( 3 );
    g_vx.SetLineStyle( 5 );
    mg.Add( &g_vx );

    g_thx.SetLineColor( kCyan+2 );
    g_thx.SetLineWidth( 3 );
    g_thx.SetLineStyle( 6 );
    mg.Add( &g_thx );

    g_beammom.SetLineColor( kAzure+1 );
    g_beammom.SetLineWidth( 3 );
    g_beammom.SetLineStyle( 7 );
    mg.Add( &g_beammom );

    g_tot.SetLineWidth( 5 );
    mg.Add( &g_tot );

    if ( pot.first == 2 ) {
      c.AddLegendEntry( &g_dx, "D_{x}", "l" );
      c.AddLegendEntry( &g_sensor, "Sensor resol.", "l" );
      c.AddLegendEntry( &g_align, "Alignment", "l" );
      c.AddLegendEntry( &g_vx, "Vertex smearing", "l" );
      c.AddLegendEntry( &g_thx, "Angular term", "l" );
      c.AddLegendEntry( &g_beammom, "Beam momentum", "l" );
      c.AddLegendEntry( &g_tot, "Combined (in quadrature)", "l" );
    }

    mg.Draw( "al" );
    mg.GetHistogram()->SetTitle( ";#xi;#delta#xi / #xi (%)" );
    c.Prettify( mg.GetHistogram() );
    TLegend* leg = c.GetLegend();
    if ( leg ) {
      leg->SetNColumns( 2 );
      leg->SetX1( 0.18 );
      leg->SetX2( 0.78 );
      leg->SetY1( 0.65 );
      leg->SetY2( 0.85 );
    }
    auto lab = new PaveText( 0.18, 0.89, 0.2, 0.92 );
    lab->SetTextSize( 0.05 );
    lab->SetTextAlign( kVAlignTop+kHAlignLeft );
    lab->AddText( pot.second );
    lab->Draw();

    mg.GetXaxis()->SetLimits( min_xi, max_xi );
    mg.GetYaxis()->SetRangeUser( 0., 15. );
    c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/xi_study/" );
  }
}

