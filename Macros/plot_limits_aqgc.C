#include "Canvas.h"

#define OUT_PATH "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp/limits"

void plot_limits_aqgc( const char* path = "/afs/cern.ch/work/l/lforthom/private/twophoton/mcprod/fpmc/build/scan_result_aqgc.dat" )
{
  gStyle->SetPalette( kBeach );
  //const double xsec_limit = 12.64e-3; //FIXME
  const double xsec_limit = 13.352e-3; //FIXME
  const double sel_eff = 0.415;
  ifstream file( path );
  //TGraph2D gr_xs( path, "%lg\t%lg\t%lg" );
  TGraph2D gr_xs;
  TGraph gr_xs_z1, gr_xs_z2;
  double z1, z2, xs;
  double min_z1 = 100., min_z2 = 100., max_z1 = -100., max_z2 = -100.;
  while ( file >> z1 >> z2 >> xs ) {
    z1 *= 1.e12;
    z2 *= 1.e12;
    xs *= sel_eff;
    min_z1 = min( min_z1, z1 );
    max_z1 = max( max_z1, z1 );
    min_z2 = min( min_z2, z2 );
    max_z2 = max( max_z2, z2 );
    //if ( fabs( z1 ) < 3.5e-12 && fabs( z2 ) < 3.5e-12 )
      gr_xs.SetPoint( gr_xs.GetN(), z1, z2, xs );
    if ( fabs( z1 ) < 1.e-6 )
      gr_xs_z2.SetPoint( gr_xs_z2.GetN(), z2, xs );
    if ( fabs( z2 ) < 1.e-6 )
      gr_xs_z1.SetPoint( gr_xs_z1.GetN(), z1, xs );
  }

  gStyle->SetOptFit( 1111 );
  const char* top_label = "FPMC #gamma#gamma #rightarrow #gamma#gamma + AQGC (13 TeV)";
  double lim1_zeta1, lim2_zeta1, lim1_zeta2, lim2_zeta2;
  {
    Canvas c( "limits_aqgc_1d", top_label, "Simulation preliminary" );
    TMultiGraph mg;
    gr_xs_z1.SetLineWidth( 3 );
    gr_xs_z1.SetLineStyle( 1 );
    gr_xs_z1.SetLineColor( kBlack );
    gr_xs_z1.SetMarkerStyle( 24 );
    gr_xs_z2.SetLineWidth( 3 );
    gr_xs_z2.SetLineStyle( 2 );
    gr_xs_z2.SetLineColor( kRed+1 );
    gr_xs_z2.SetMarkerColor( kRed+1 );
    gr_xs_z2.SetMarkerStyle( 25 );
    mg.Add( &gr_xs_z1 );
    mg.Add( &gr_xs_z2 );
    mg.Draw( "ap" );
    mg.GetHistogram()->SetTitle( ";#zeta_{1,2} (#times 10^{-12});#sigma_{obs} (pb)" );
    c.Prettify( mg.GetHistogram() );
    gStyle->SetFuncColor( kBlack );
    gr_xs_z1.Fit( "pol2", "0" );
    gr_xs_z2.Fit( "pol2", "0" );
    gPad->Update();
    auto fit1 = (TF1*)gr_xs_z1.GetListOfFunctions()->FindObject( "pol2" );
    fit1->SetLineColor( kBlack );
    fit1->SetLineStyle( 1 );
    fit1->Draw( "same" );
    lim1_zeta1 = fit1->GetX( xsec_limit, min_z1, 0. );
    lim2_zeta1 = fit1->GetX( xsec_limit, 0., max_z1 );
    cout << "zeta1 --> " << lim1_zeta1 << "|" << lim2_zeta1 << endl;
    auto fit2 = (TF1*)gr_xs_z2.GetListOfFunctions()->FindObject( "pol2" );
    fit2->SetLineColor( kRed+1 );
    fit2->SetLineStyle( 2 );
    fit2->Draw( "same" );
    lim1_zeta2 = fit2->GetX( xsec_limit, min_z2, 0. );
    lim2_zeta2 = fit2->GetX( xsec_limit, 0., max_z2 );
    cout << "zeta2 --> " << lim1_zeta2 << "|" << lim2_zeta2 << endl;
    auto stat1 = (TPaveStats*)gr_xs_z1.GetListOfFunctions()->FindObject( "stats" );
    stat1->SetTextSize( 0.025 );
    stat1->SetY1NDC( 0.73 );
    stat1->SetY2NDC( 0.93 );
    stat1->Draw();
    auto stat2 = (TPaveStats*)gr_xs_z2.GetListOfFunctions()->FindObject( "stats" );
    //stat2->SetX1NDC( 0.6 );
    //stat2->SetX2NDC( 0.875 );
    stat2->SetY1NDC( 0.52 );
    stat2->SetY2NDC( 0.72 );
    stat2->SetTextColor( kRed+1 );
    stat2->SetTextSize( 0.025 );
    stat2->Draw();

    /*TF1 f_pol2_1( "pol2_1", "pol2", min_z1, max_z1 );
    TF1 f_pol2_2( "pol2_2", "pol2", min_z2, max_z2 );
    gr_xs_z1.Fit( &f_pol2_1 );
    gr_xs_z2.Fit( &f_pol2_2 );
    gPad->Update();*/
    /*unsigned short i = 0;
    for ( auto& ps : vector<pair<TGraph,int> >{ { gr_xs_z1, kBlack }, { gr_xs_z2, kRed+1 } } ) {
      //ps.first.GetListOfFunctions()->Dump();
      auto fit = (TF1*)ps.first.GetListOfFunctions()->FindObject(Form( "pol2_%d", i+1 ));
      fit->SetLineColor( ps.second );
      auto stat = (TPaveStats*)ps.first.GetListOfFunctions()->FindObject("stats");
      stat->SetX1NDC( 0.6 );
      stat->SetX2NDC( 0.875 );
      stat->SetY1NDC( 0.55+(1-i)*0.18 );
      stat->SetY2NDC( 0.72+(1-i)*0.18 );
      stat->SetTextColor( ps.second );
      stat->SetTextSize( 0.03 );
      ++i;
    }*/
    auto xaxis = mg.GetHistogram()->GetXaxis();
    TF1 f_limit( "limit", "[0]", xaxis->GetXmin(), xaxis->GetXmax() );
    f_limit.SetParameter( 0, xsec_limit );
    f_limit.SetLineColor( kRed );
    f_limit.SetLineWidth( 3 );
    f_limit.Draw( "same" );
    c.SetLegendX1( 0.18 );
    c.SetLegendY1( 0.16 );
    c.AddLegendEntry( &f_limit, "#sigma_{95%}^{lim}", "l" );
    c.AddLegendEntry( &gr_xs_z2, "#sigma_{obs}(#zeta_{1} = 0, #zeta_{2})", "lp" );
    c.AddLegendEntry( &gr_xs_z1, "#sigma_{obs}(#zeta_{1}, #zeta_{2} = 0)", "lp" );
    mg.SetMaximum( 2. );
    mg.SetMinimum( 1.4e-4 );
    c.SetLogy();
    c.SetGrid();
    c.Save( "pdf,png", OUT_PATH );
  }
  auto cont = gr_xs.GetContourList( xsec_limit );
  if ( !cont || cont->GetSize() > 1 )
    throw runtime_error( "invalid contour retrieved!" );
  auto gr_cont = static_cast<TGraph*>( cont->At( 0 ) );
  if ( !gr_cont )
    throw runtime_error( "failed to build contour graph!" );
  {
    Canvas c( "limits_aqgc_2d", "9.4 fb^{-1} (13 TeV)", "Preliminary", false, Canvas::Align::right );
    gr_cont->SetTitle( ";#zeta_{1} (#times 10^{-12});#zeta_{2} (#times 10^{-12})" );
    gr_cont->SetLineWidth( 3 );
    gr_cont->SetLineStyle( 2 );
    gr_cont->SetLineColor( 1 );
    gr_cont->Draw( "al" );
    TGraphAsymmErrors gr_1d;
    gr_1d.SetPoint( 0, 0., 0. );
    gr_1d.SetPointError( 0, -lim1_zeta1, lim2_zeta1, -lim1_zeta2, lim2_zeta2 );
    gr_1d.SetLineWidth( 5 );
    gr_1d.Draw( "e same" );
    TGraph gr_sm;
    gr_sm.SetPoint( 0, 0., 0. );
    gr_sm.SetMarkerStyle( 20 );
    gr_sm.Draw( "p same" );
    c.SetLegendX1( 0.16 );
    c.SetLegendY1( 0.14 );
    c.AddLegendEntry( &gr_sm, "Standard Model", "p" );
    c.AddLegendEntry( gr_cont, "CMS 95% confidence region", "l" );
    c.AddLegendEntry( &gr_1d, "CMS 1-D limit, 95% confidence region", "l" );
    c.Prettify( (TH1D*)gr_cont->GetHistogram() );
    c.SetLogz();
    //gr_xs.GetXaxis()->SetRangeUser( -2.5e-12, 2.5e-12 );
    c.Save( "pdf,png", OUT_PATH );
  }
  {
    Canvas c( "scan_xsec_aqgc", top_label, "Simulation preliminary" );
    c.SetRightMargin( 0.15 );
    //gr_xs.Draw( "cont4z" );
    gr_xs.Draw( "cont4z" );
    //gr_xs.Draw( "colz" );
    gr_xs.SetTitle( ";#zeta_{1} (#times 10^{-12});#zeta_{2} (#times 10^{-12});#sigma_{obs} (pb)" );
    /*gr_cont->SetLineWidth( 3 );
    gr_cont->SetLineStyle( 2 );
    gr_cont->SetLineColor( 1 );
    gr_cont->Draw( "l same" );
    TGraphAsymmErrors gr_1d;
    gr_1d.SetPoint( 0, 0., 0. );
    gr_1d.SetPointError( 0, -lim1_zeta1, lim2_zeta1, -lim1_zeta2, lim2_zeta2 );
    gr_1d.SetLineWidth( 5 );
    gr_1d.Draw( "e same" );
    TGraph gr_sm;
    gr_sm.SetPoint( 0, 0., 0. );
    gr_sm.SetMarkerStyle( 20 );
    gr_sm.Draw( "p same" );
    c.SetLegendX1( 0.16 );
    c.SetLegendY1( 0.14 );
    c.AddLegendEntry( &gr_sm, "Standard Model", "p" );
    c.AddLegendEntry( gr_cont, "CMS 95% confidence region", "l" );
    c.AddLegendEntry( &gr_1d, "CMS 1-D limit, 95% confidence region", "l" );*/
    c.Prettify( (TH1D*)gr_xs.GetHistogram() );
    c.SetLogz();
    c.SetGrid();
    //gr_xs.GetXaxis()->SetRangeUser( -2.5e-12, 2.5e-12 );
    c.Save( "pdf,png", OUT_PATH );
  }
}
