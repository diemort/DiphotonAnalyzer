#include "Canvas.h"

void plot_limits_asym()
{
  string k_nm, out_file;
  string combine_path = "/afs/cern.ch/work/l/lforthom/private/twophoton/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit";
  const string method = "AsymptoticLimits";
  map<float,string> m_sig_lab = {
    { -1., "obs" },
    { 0.5, "exp" },
    { 0.16, "m1s" },
    { 0.84, "p1s" },
    { 0.025, "m2s" },
    { 0.975, "p2s" },
  };
  struct value_t {
    string combine_path, k_nm, out_file;
    vector<int> masses;
  };
  vector<value_t> samples = {
    { "results_f_5e-1", "5 #times 10^{-1}", "limits_alp_5e-1", { 500, 750, 1000, 1250, 1500, 1750, 2000 } },
    { "results_f_1e-1", "10^{-1}", "limits_alp_1e-1", { 500, 750, 1000, 1250, 1500, 1750 } },
    { "results_f_1e-1p5", "10^{-1.5}", "limits_alp_1e-1p5", { 500, 750, 1000, 1250, 1500, 1750 } }
  };
  for ( const auto& smp : samples ) {
    TGraph g_obs, g_exp, g_one;
    TGraphAsymmErrors g_exp_1sig, g_exp_2sig;
    unsigned short i = 0;
    for ( const auto& m : smp.masses ) {
      map<string,double> m_vals;
      const string out_combine_path = Form( "%s/%s/higgsCombineTest.%s.mH%d.root", combine_path.c_str(), smp.combine_path.c_str(), method.c_str(), m );
      auto output_combine = TFile::Open( out_combine_path.c_str() );
      if ( !output_combine )
        continue;
      auto tree_combine = (TTree*)output_combine->Get( "limit" );
      if ( !tree_combine )
        continue;
      double limit = 0.;
      float quantileExpected = 0.;
      tree_combine->SetBranchAddress( "limit", &limit );
      tree_combine->SetBranchAddress( "quantileExpected", &quantileExpected );
      for ( unsigned long long j = 0; j < tree_combine->GetEntriesFast(); ++j ) {
        tree_combine->GetEntry( j );
        m_vals[m_sig_lab.at( quantileExpected )] = limit;
      }
      g_one.SetPoint( i, m, 1. );
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
      Canvas c( smp.out_file.c_str(), "9.4 fb^{-1} (13 TeV)" );
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
      g_one.SetLineColor( kRed+1 );
      g_one.SetLineStyle( 2 );
      g_one.SetLineWidth( 3 );
      g_one.Draw( "l" );
      g_exp.Draw( "l" );
      g_obs.Draw( "l" );
      c.AddLegendEntry( &g_obs, "Observed", "l" );
      c.AddLegendEntry( &g_exp, "Expected", "l" );
      c.AddLegendEntry( &g_exp_1sig, "68% expected", "f" );
      c.AddLegendEntry( &g_exp_2sig, "95% expected", "f" );
      //g_exp_2sig.SetTitle( ";m_{a} (GeV);Upper limit on #sigma/#sigma(#gamma#gamma#rightarrowa#rightarrow#gamma#gamma) (pb)" );
      g_exp_2sig.SetTitle( ";m_{a} (GeV);Upper limit on #sigma/#sigma_{ALP}" );
      c.Prettify( g_exp_2sig.GetHistogram() );
      g_exp_2sig.GetXaxis()->SetRangeUser( *smp.masses.begin(), *smp.masses.rbegin() );
      g_exp_2sig.SetMaximum( 1.e3 );
      c.SetLogy();
      PaveText::topLabel( Form( "f^{-1} = %s GeV^{-1}", smp.k_nm.c_str() ) );
      c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
    }
  }
}
