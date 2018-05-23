#include "Canvas.h"
#include "pot_alignment.h"
#include "xi_reconstruction.h"

#include "DiphotonAnalyzer/TreeProducer/interface/TreeEvent.h"

#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH2.h"

#include <fstream>

#define MAX_PROTONS 30

map<unsigned short,const char*> pots_names = { { 2, "Sector 45 - near pot" }, { 3, "Sector 45 - far pot" }, { 102, "Sector 56 - near pot" }, { 103, "Sector 56 - far pot" } };

void test_xireco()
{
  TFile f( "Samples/output_Run2016BCG_looseCuts_28jun.root" );
  TTree* tr = dynamic_cast<TTree*>( f.Get( "ntp" ) );

  xi_reco::load_file( "TreeProducer/data/optics_jun22.root" );
  pot_align::load_file( "TreeProducer/data/alignment_collection_v2.out" );

  TreeEvent ev;
  ev.attach( tr, true );

  map<unsigned short,TProfile*> m_xi;
  for ( const auto& pot : pots_names )
    m_xi[pot.first] = new TProfile( Form( "dxi_vs_xi_%d", pot.first ), Form( "%s;#xi(%s);#Delta#xi/#xi", pot.second, pot.second ), 32, 0.02, 0.18 );

  const unsigned long long num_events = tr->GetEntriesFast();
  for ( unsigned long long i = 0; i < num_events; ++i ) {
    tr->GetEntry( i );
    //cout << "event " << i << ": " << ev.num_proton_track << " proton tracks, " << ev.num_diphoton << " diphoton candidates" << endl;

    // first loop to identify the tracks and their respective pot

    auto align = pot_align::get_alignments( ev.fill_number );

    for ( unsigned short j = 0; j < ev.num_proton_track; ++j ) {
      const unsigned short pot_id = 100*ev.proton_track_side[j]+ev.proton_track_pot[j];
      const auto& al = align[pot_id];

      //----- reconstruct the kinematics
      double xi = 0., xi_err = 0.;
      xi_reco::reconstruct( ev.proton_track_x[j]+al.x, ev.proton_track_side[j], ev.proton_track_pot[j], xi, xi_err );

      //----- associate each track to a RP
      m_xi[pot_id]->Fill( xi, xi_err/xi );
    }

  } // loop on events
  Canvas c( "dxi_vs_xi", "CMS+TOTEM Preliminary 2016, #sqrt{s} = 13 TeV" );
  THStack hs;
  int colours[] = { kBlack, kRed+1, kBlue+1, kGreen-2 };
  unsigned short i = 0;
  const unsigned short num_points = 100;
  const double min_x = 0.02, max_x = 0.2;
  const double si_x_alignment = 150.0e-6; // in m, alignment uncertainty
  const double si_x_neglected_angle = 150.0e-6; // in m, to (approximately) account for the neglected angular term in proton transport
  const double si_rel_D = 0.055; // 1, relative uncertainty of dispersion
  for ( const auto& pot : pots_names ) {
    m_xi[pot.first]->SetLineWidth( 3 );
    const unsigned short arm = pot.first/100, pot = pot%100;
    auto sp = xi_reco::get_spline( arm, pot );
    for ( unsigned short j = 0; j < num_points; ++j ) {
      const double xi = min_x+( max_x-min_x )*j/( num_points-1 );
      const double si_x = hypot( si_x_alignment, si_x_neglected_angle );
      const double si_xi_from_x = sp->Eval( x+si_x )-xi;
      const double si_xi_from_D_x = si_rel_D * xi;
      xi_err = hypot( si_xi_from_x, si_xi_from_D_x );
    }

    m_xi[pot.first]->SetLineColor( colours[i] );
    //hs.Add( m_xi[pot.first], "hist c" );
    hs.Add( m_xi[pot.first] );
    c.AddLegendEntry( m_xi[pot.first], pot.second );
    ++i;
  }
  hs.Draw( "nostack" );
  //hs.GetXaxis()->SetTitle( m_xi.begin()->second->GetXaxis()->GetTitle() );
  hs.GetXaxis()->SetTitle( "#xi" );
  hs.GetYaxis()->SetTitle( m_xi.begin()->second->GetYaxis()->GetTitle() );
  c.Prettify( hs.GetHistogram() );
  c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
}
