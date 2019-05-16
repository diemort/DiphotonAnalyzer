#include "Canvas.h"

void plot_limits()
{
  //const string combine_path = "/afs/cern.ch/work/l/lforthom/private/twophoton/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit";
  const string combine_path = "/afs/cern.ch/work/l/lforthom/private/twophoton/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/results";
  const string method = "HybridNew";
  vector<int> masses = { 500, 750, 1000, 1250, 1500, 1750, 2000 };
  map<string,string> m_sig_lab = {
    { "obs", "" },
    { "exp", ".quant0.500" },
    { "m1s", ".quant0.160" },
    { "p1s", ".quant0.840" },
    { "m2s", ".quant0.025" },
    { "p2s", ".quant0.975" },
  };
  TGraph g_obs, g_exp;
  TGraphAsymmErrors g_exp_1sig, g_exp_2sig;
  unsigned short i = 0;
  for ( const auto& m : masses ) {
    map<string,double> m_vals;
    for ( const auto& s : m_sig_lab ) {
      const string out_combine_path = Form( "%s/higgsCombineTest.%s.mH%d%s.root", combine_path.c_str(), method.c_str(), m, s.second.c_str() );
      //cout << s.first << "|" << out_combine_path << endl;
      auto output_combine = TFile::Open( out_combine_path.c_str() );
      if ( !output_combine )
        continue;
      auto tree_combine = (TTree*)output_combine->Get( "limit" );
      if ( !tree_combine )
        continue;
      double limit = 0.;
      tree_combine->SetBranchAddress( "limit", &limit );
      tree_combine->GetEntry( 0 );
      m_vals[s.first] = limit;
      cout << s.first << "|" << limit << endl;
    }
    if ( m_vals.count( "obs" ) > 0 && m_vals.at( "obs" ) != 0. )
      g_obs.SetPoint( i, m, m_vals.at( "obs" ) );
    if ( m_vals.count( "exp" ) > 0 && m_vals.at( "exp" ) != 0. ) {
      g_exp.SetPoint( i, m, m_vals.at( "exp" ) );
      g_exp_1sig.SetPoint( i, m, m_vals.at( "exp" ) );
      if ( m_vals.count( "p1s" ) > 0 && m_vals.count( "m1s" ) > 0 )
        g_exp_1sig.SetPointError( i, 0., 0., fabs( m_vals.at( "exp" )-m_vals.at( "m1s" ) ), fabs( m_vals.at( "exp" )-m_vals.at( "p1s" ) ) );
      g_exp_2sig.SetPoint( i, m, m_vals.at( "exp" ) );
      if ( m_vals.count( "p2s" ) > 0 && m_vals.count( "m2s" ) > 0 )
        g_exp_2sig.SetPointError( i, 0., 0., fabs( m_vals.at( "exp" )-m_vals.at( "m2s" ) ), fabs( m_vals.at( "exp" )-m_vals.at( "p2s" ) ) );
    }
    ++i;
  }
  {
    Canvas c( "limits_alp", "9.4 fb^{-1} (13 TeV)" );
    g_exp_2sig.SetFillColor( kYellow-4 );
    g_exp_1sig.SetFillColor( kGreen+1 );
    g_exp_2sig.Draw( "a3" );
    g_exp_1sig.Draw( "3" );
    /*g_exp_2sig.Draw( "a4" );
    g_exp_1sig.Draw( "4" );*/
    g_obs.SetLineStyle( 1 );
    g_obs.SetLineWidth( 3 );
    g_exp.SetLineStyle( 2 );
    g_exp.SetLineWidth( 3 );
    g_exp.Draw( "l" );
    g_obs.Draw( "l" );
    c.AddLegendEntry( &g_obs, "Observed", "l" );
    c.AddLegendEntry( &g_exp, "Expected", "l" );
    c.AddLegendEntry( &g_exp_1sig, "68% expected", "f" );
    c.AddLegendEntry( &g_exp_2sig, "95% expected", "f" );
    g_exp_2sig.SetTitle( ";m_{a} (GeV);Upper limit on #sigma(#gamma#gamma #rightarrow a #rightarrow #gamma#gamma)" );
    c.Prettify( g_exp_2sig.GetHistogram() );
    g_exp_2sig.SetMaximum( 1.e3 );
    c.SetLogy();
    c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
  }
}
