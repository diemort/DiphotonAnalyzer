#ifndef DiphotonAnalyzer_Macros_xi_reconstruction_h
#define DiphotonAnalyzer_Macros_xi_reconstruction_h

#include "TFile.h"
#include "TSpline.h"
#include <map>

namespace xi_reco
{
  std::map<unsigned short,TSpline*> splines, r_splines;
  std::map<unsigned short,TGraph*> g_x0_vs_xi, g_x_vs_xi, g_Lx_vs_xi, g_vx_vs_xi, g_xi_vs_x;
  std::map<unsigned short,double> disp_rel_err;
  std::map<unsigned short,std::pair<double,double> > dispersions; // RP -> ( disp, err_disp )

  void set_dispersions( const std::map<unsigned short,std::pair<double,double> >& disp ) { dispersions = disp; }

  void
  load_file( const char* filename )
  {
    auto file = TFile::Open( filename );
    splines = std::map<unsigned short,TSpline*>( {
      { 2, dynamic_cast<TSpline*>( file->Get( "s_x_to_xi_L_1_N" )->Clone() ) }, //45N
      { 3, dynamic_cast<TSpline*>( file->Get( "s_x_to_xi_L_1_F" )->Clone() ) }, //45F
      { 102, dynamic_cast<TSpline*>( file->Get( "s_x_to_xi_R_1_N" )->Clone() ) }, //56N
      { 103, dynamic_cast<TSpline*>( file->Get( "s_x_to_xi_R_1_F" )->Clone() ) } //56F
    } );
    delete file;
  }

  void
  load_optics_file( const char* filename )
  {
    TFile file( filename );
    std::map<unsigned short,const char*> paths = {
      { 2, "ctppsPlotOpticalFunctions_45/ip5_to_station_150_h_1_lhcb2" }, // 45N
      { 3, "ctppsPlotOpticalFunctions_45/ip5_to_station_150_h_2_lhcb2" }, // 45F
      { 102, "ctppsPlotOpticalFunctions_56/ip5_to_station_150_h_1_lhcb1" }, // 56N
      { 103, "ctppsPlotOpticalFunctions_56/ip5_to_station_150_h_2_lhcb1" } // 56F
    };
    for ( const auto& p : paths ) {
      g_x0_vs_xi[p.first] = dynamic_cast<TGraph*>( file.Get( Form( "%s/g_x0_vs_xi", p.second ) )->Clone() );
      g_x_vs_xi[p.first] = new TGraph;
      for ( int i = 0; i < g_x0_vs_xi[p.first]->GetN(); ++i )
        g_x_vs_xi[p.first]->SetPoint( i, g_x0_vs_xi[p.first]->GetX()[i], g_x0_vs_xi[p.first]->GetY()[i]-g_x0_vs_xi[p.first]->GetY()[0] );
      g_Lx_vs_xi[p.first] = dynamic_cast<TGraph*>( file.Get( Form( "%s/g_L_x_vs_xi", p.second ) )->Clone() );
      g_vx_vs_xi[p.first] = dynamic_cast<TGraph*>( file.Get( Form( "%s/g_v_x_vs_xi", p.second ) )->Clone() );
      g_xi_vs_x[p.first] = dynamic_cast<TGraph*>( file.Get( Form( "%s/g_xi_vs_xso", p.second ) )->Clone() );
      // build spline
      r_splines[p.first] = new TSpline3( "", g_xi_vs_x[p.first]->GetX(), g_xi_vs_x[p.first]->GetY(), g_xi_vs_x[p.first]->GetN() );
    };
    disp_rel_err = std::map<unsigned short,double>{
      { 2, 0.0387 }, { 3, 0.0435 }, { 102, 0.0526 }, { 103, 0.0614 }
    };
  }

  unsigned short
  get_id( unsigned short arm, unsigned short pot )
  {
    return 100*arm+pot;
  }

  TGraph*
  get_graph( unsigned int arm, unsigned int pot )
  {
    auto it_gr = g_x_vs_xi.find( get_id( arm, pot ) );
    if ( it_gr == g_x_vs_xi.end() ) return nullptr;
    return it_gr->second;
  }

  TSpline*
  get_spline( unsigned int arm, unsigned int pot )
  {
    auto it_sp = splines.find( get_id( arm, pot ) );
    if ( it_sp == splines.end() ) return nullptr;
    return it_sp->second;
  }

  void
  reconstruct_old( double x, unsigned int arm, unsigned int pot, double& xi, double& xi_err )
  {
    xi = xi_err = 0.;
    TSpline* sp = get_spline( arm, pot );
    if ( !sp ) return;
    xi = sp->Eval( x );
    // old
    /*double de_xi = sp->Eval( x+0.4e-3 )-xi, de_dx = 0.1*xi;
    xi_err = sqrt( de_xi*de_xi+de_dx*de_dx );*/
    // determine uncertainty of xi
    const double si_x_alignment = 150.0e-6; // in m, alignment uncertainty
    const double si_x_neglected_angle = 150.0e-6; // in m, to (approximately) account for the neglected angular term in proton transport
    const double si_rel_D = 0.055; // 1, relative uncertainty of dispersion
    const double si_x = hypot( si_x_alignment, si_x_neglected_angle );
    const double si_xi_from_x = sp->Eval( x+si_x )-xi;
    const double si_xi_from_D_x = si_rel_D * xi;
    xi_err = hypot( si_xi_from_x, si_xi_from_D_x );
  }

  void
  reconstruct( double x, unsigned int arm, unsigned int pot, double& xi, double& xi_err )
  {
    xi = xi_err = 0.;
    const unsigned short sid = get_id( arm, pot );
    if ( r_splines.count( sid ) == 0 )
      return;
    const auto sp = r_splines.at( sid );
    xi = sp->Eval( x );
    // prepare the run for xi uncertainty
    const double de_xi = 1.e-3;
    const double de_x = g_x0_vs_xi.at( sid )->Eval( xi+de_xi/2. ) - g_x0_vs_xi.at( sid )->Eval( xi-de_xi/2. );
    const double D = de_x / de_xi;  // m
    const double Lx = g_Lx_vs_xi.at( sid )->Eval( xi ); // m
    const double vx = g_vx_vs_xi.at( sid )->Eval( xi );
    // a few parameters
    const double si_x_alignment = 150.0e-6; // in m, alignment uncertainty
    const double si_vertex = 10.e-6;
    const double si_beam_div = 20.e-6;
    const double si_sensor_resol = 12.e-6;
    const double si_de_disp = disp_rel_err.at( get_id( arm, pot ) );
    const double si_beam_mom = 1.e-3;
    // determine uncertainty of xi
    const double si_D_x = si_de_disp * xi;
    const double si_sensor = si_sensor_resol / fabs( D );
    const double si_alig = si_x_alignment / fabs( D );
    const double si_term_th_x = fabs( Lx/D ) * si_beam_div;
    const double si_term_vtx_x = fabs( vx/D ) * si_vertex;
    for ( const auto &c : { si_beam_mom, si_D_x, si_sensor, si_alig, si_term_th_x, si_term_vtx_x } )
      xi_err += c*c;
    xi_err = sqrt( xi_err );
  }

  void
  reconstruct_lin( double x, unsigned int arm, unsigned int pot, double& xi, double& xi_err )
  {
    xi = xi_err = 0;
    if ( dispersions.count( 100*arm+pot ) == 0 ) return;

    const auto disp = dispersions.at( 100*arm+pot ); // ( disp, err_disp )
    const double de_x = 150.e-6; // alignment uncertainty

    xi = x/disp.first;
    xi_err = hypot( de_x/disp.first, disp.second*xi );
    //xi_err = xi * hypot( de_x/x, disp.second );
  }
}

#endif
