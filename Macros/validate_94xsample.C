#include "Canvas.h"
#include "../TreeProducer/interface/TreeEvent.h"
//#include "/afs/cern.ch/work/l/lforthom/private/twophoton/CMSSW_8_0_26_patch1/src/DiphotonAnalyzer/TreeProducer/interface/TreeEvent.h"
#include "../../../../CMSSW_8_0_26_patch1/src/DiphotonAnalyzer/TreeProducer/interface/TreeEvent.h"

void validate_94xsample()
{
  //auto file_94x = TFile::Open( "samples/ntuple-Run2016B_94Xrereco_v2.root" );
  auto file_94x = TFile::Open( "/eos/cms/store/group/phys_pps/diphoton/DoubleEG/ntuple-Run2016BCG_94Xrereco_v1.root" );
  auto file_80x = TFile::Open( "/eos/cms/store/user/lforthom/twophoton/samples_80x/output_Run2016BCG_looseCuts_28jun.root" );
  auto tree_94x = dynamic_cast<TTree*>( file_94x->Get( "ntp" ) );
  auto tree_80x = dynamic_cast<TTree*>( file_80x->Get( "ntp" ) );

  const float eta_cut = 2.5, min_etaveto = 1.4442, max_etaveto = 1.566;
  bool is_elastic = true;

  vector<const char*> v_type = { "8_0_26_patch1", "9_4_9", "9_4_9 (+ele.veto)" };
  vector<TH1D*> v_h_ptg, v_h_etag, v_h_mgg, v_h_ptgg, v_h_dphigg, v_h_dphigg_zoom;
  vector<double> pt_bins{ 70., 80., 100., 120., 150., 180., 220., 250., 300., 400., 500., 600., 700., 800., 1000., 1200., 1500. };
  vector<double> acop_bins{ 0., 0.01, 0.02, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 1. };
  for ( unsigned short i = 0; i < 3; ++i ) {
    if ( is_elastic )
      v_h_ptg.emplace_back( new TH1D( Form( "ptg_%s", v_type[i] ), "p_{T}^{#gamma}@@Events@@GeV", 50, 50., 550. ) );
    else
      //v_h_ptg.emplace_back( new TH1D( Form( "ptg_%s", v_type[i] ), "p_{T}^{#gamma}@@Events@@GeV", 45, 75., 975. ) );
      v_h_ptg.emplace_back( new TH1D( Form( "ptg_%s", v_type[i] ), "p_{T}^{#gamma}@@Events@@GeV", pt_bins.size()-1, pt_bins.data() ) );
    v_h_etag.emplace_back( new TH1D( Form( "etag_%s", v_type[i] ), "#eta^{#gamma}@@Events@@?.g", 25, -2.5, 2.5 ) );
    v_h_mgg.emplace_back( new TH1D( Form( "mgg_%s", v_type[i] ), "m_{#gamma#gamma}@@Events@@GeV", 50, 300., 2300. ) );
    v_h_ptgg.emplace_back( new TH1D( Form( "ptgg_%s", v_type[i] ), "p_{T}^{#gamma#gamma}@@Events@@GeV", 50, 0., 250. ) );
    //v_h_ptgg.emplace_back( new TH1D( Form( "ptgg_%s", v_type[i] ), "p_{T}^{#gamma#gamma}@@Events@@GeV", 100, 0., 500. ) );
    v_h_dphigg.emplace_back( new TH1D( Form( "dphigg_%s", v_type[i] ), "1-|#Delta#phi(#gamma#gamma)/#pi|@@Events@@?.g", 50, 0., 1. ) );
    //v_h_dphigg.emplace_back( new TH1D( Form( "dphigg_%s", v_type[i] ), "1-|#Delta#phi(#gamma#gamma)/#pi|@@Events@@?.g", acop_bins.size()-1, acop_bins.data() ) );
    v_h_dphigg_zoom.emplace_back( new TH1D( Form( "dphigg_zoom_%s", v_type[i] ), "1-|#Delta#phi(#gamma#gamma)/#pi| (#times 10^{-3})@@Events@@?.g", 10, 0., 5. ) );
  }

  cout << "opening 80x version" << endl;
  TreeEvent evt_80x;
  evt_80x.attach( tree_80x, true );
  for ( unsigned long long i = 0; i < tree_80x->GetEntriesFast(); ++i ) {
    tree_80x->GetEntry( i );
    //if ( evt_80x.run_id > 275400 ) continue;
    for ( unsigned short j = 0; j < evt_80x.num_diphoton; ++j ) {
      unsigned short ev_class = gggg::TreeEvent::invalid;
      if ( fabs( evt_80x.diphoton_eta1[j] ) <= eta_cut && fabs( evt_80x.diphoton_eta2[j] ) <= eta_cut ) {
        if ( fabs( evt_80x.diphoton_eta1[j] ) < min_etaveto && fabs( evt_80x.diphoton_eta2[j] ) > max_etaveto ) ev_class = gggg::TreeEvent::ebee;
        else if ( fabs( evt_80x.diphoton_eta2[j] ) < min_etaveto && fabs( evt_80x.diphoton_eta1[j] ) > max_etaveto ) ev_class = gggg::TreeEvent::ebee;
        else if ( fabs( evt_80x.diphoton_eta1[j] ) < min_etaveto && fabs( evt_80x.diphoton_eta2[j] ) < min_etaveto ) ev_class = gggg::TreeEvent::ebeb;
        else if ( fabs( evt_80x.diphoton_eta1[j] ) > max_etaveto && fabs( evt_80x.diphoton_eta2[j] ) > max_etaveto ) ev_class = gggg::TreeEvent::eeee;
      }
      //----- only keep EBEE and EBEB diphoton events
      if ( ev_class == gggg::TreeEvent::invalid ) continue;
      if ( ev_class == gggg::TreeEvent::eeee ) continue; //FIXME FIXME

      if ( evt_80x.diphoton_pt1[j] < 75. ) continue;
      if ( evt_80x.diphoton_pt2[j] < 75. ) continue;
      if ( fabs( evt_80x.diphoton_eta1[j] ) > eta_cut || ( fabs( evt_80x.diphoton_eta1[j] ) > min_etaveto && fabs( evt_80x.diphoton_eta1[j] ) < max_etaveto ) ) continue;
      if ( fabs( evt_80x.diphoton_eta2[j] ) > eta_cut || ( fabs( evt_80x.diphoton_eta2[j] ) > min_etaveto && fabs( evt_80x.diphoton_eta2[j] ) < max_etaveto ) ) continue;
      if ( evt_80x.diphoton_r91[j] < 0.94 ) continue;
      if ( evt_80x.diphoton_r92[j] < 0.94 ) continue;
      if ( evt_80x.diphoton_mass[j] < 350. ) continue;
      //----- back-to-back photons
      if ( is_elastic && 1.-fabs( evt_80x.diphoton_dphi[j] )/M_PI > 0.005 ) continue;

      //----- ELASTIC REGION
      v_h_ptg[0]->Fill( evt_80x.diphoton_pt1[j] );
      v_h_ptg[0]->Fill( evt_80x.diphoton_pt2[j] );
      v_h_etag[0]->Fill( evt_80x.diphoton_eta1[j] );
      v_h_etag[0]->Fill( evt_80x.diphoton_eta2[j] );
      v_h_mgg[0]->Fill( evt_80x.diphoton_mass[j] );
      v_h_ptgg[0]->Fill( evt_80x.diphoton_pt[j] );
      v_h_dphigg[0]->Fill( 1.-fabs( evt_80x.diphoton_dphi[j]/M_PI ) );
      v_h_dphigg_zoom[0]->Fill( ( 1.-fabs( evt_80x.diphoton_dphi[j]/M_PI ) )*1.e3 );
    }
  }

  cout << "opening 94x version" << endl;
  gggg::TreeEvent evt_94x;
  evt_94x.attach( tree_94x, true );
  for ( unsigned long long i = 0; i < tree_94x->GetEntriesFast(); ++i ) {
    tree_94x->GetEntry( i );
    for ( unsigned short j = 0; j < evt_94x.num_diphoton; ++j ) {
      unsigned short ev_class = gggg::TreeEvent::invalid;
      if ( fabs( evt_94x.diphoton_eta1[j] ) <= eta_cut && fabs( evt_94x.diphoton_eta2[j] ) <= eta_cut ) {
        if ( fabs( evt_94x.diphoton_eta1[j] ) < min_etaveto && fabs( evt_94x.diphoton_eta2[j] ) > max_etaveto ) ev_class = gggg::TreeEvent::ebee;
        else if ( fabs( evt_94x.diphoton_eta2[j] ) < min_etaveto && fabs( evt_94x.diphoton_eta1[j] ) > max_etaveto ) ev_class = gggg::TreeEvent::ebee;
        else if ( fabs( evt_94x.diphoton_eta1[j] ) < min_etaveto && fabs( evt_94x.diphoton_eta2[j] ) < min_etaveto ) ev_class = gggg::TreeEvent::ebeb;
        else if ( fabs( evt_94x.diphoton_eta1[j] ) > max_etaveto && fabs( evt_94x.diphoton_eta2[j] ) > max_etaveto ) ev_class = gggg::TreeEvent::eeee;
      }
      //----- only keep EBEE and EBEB diphoton events
      if ( ev_class == gggg::TreeEvent::invalid ) continue;
      if ( ev_class == gggg::TreeEvent::eeee ) continue; //FIXME FIXME

      if ( evt_94x.diphoton_pt1[j] < 75. ) continue;
      if ( evt_94x.diphoton_pt2[j] < 75. ) continue;
      if ( fabs( evt_94x.diphoton_eta1[j] ) > eta_cut || ( fabs( evt_94x.diphoton_eta1[j] ) > min_etaveto && fabs( evt_94x.diphoton_eta1[j] ) < max_etaveto ) ) continue;
      if ( fabs( evt_94x.diphoton_eta2[j] ) > eta_cut || ( fabs( evt_94x.diphoton_eta2[j] ) > min_etaveto && fabs( evt_94x.diphoton_eta2[j] ) < max_etaveto ) ) continue;
      if ( evt_94x.diphoton_r91[j] < 0.94 ) continue;
      if ( evt_94x.diphoton_r92[j] < 0.94 ) continue;
      if ( evt_94x.diphoton_mass[j] < 350. ) continue;
      //----- back-to-back photons
      if ( is_elastic && 1.-fabs( evt_94x.diphoton_dphi[j] )/M_PI > 0.005 ) continue;

      //----- ELASTIC REGION
      v_h_ptg[1]->Fill( evt_94x.diphoton_pt1[j] );
      v_h_ptg[1]->Fill( evt_94x.diphoton_pt2[j] );
      v_h_etag[1]->Fill( evt_94x.diphoton_eta1[j] );
      v_h_etag[1]->Fill( evt_94x.diphoton_eta2[j] );
      v_h_mgg[1]->Fill( evt_94x.diphoton_mass[j] );
      v_h_ptgg[1]->Fill( evt_94x.diphoton_pt[j] );
      v_h_dphigg[1]->Fill( 1.-fabs( evt_94x.diphoton_dphi[j]/M_PI ) );
      v_h_dphigg_zoom[1]->Fill( ( 1.-fabs( evt_94x.diphoton_dphi[j]/M_PI ) )*1.e3 );

      if ( evt_94x.diphoton_ele_veto1[j] == 1 && evt_94x.diphoton_ele_veto2[j] == 1 ) {
        v_h_ptg[2]->Fill( evt_94x.diphoton_pt1[j] );
        v_h_ptg[2]->Fill( evt_94x.diphoton_pt2[j] );
        v_h_etag[2]->Fill( evt_94x.diphoton_eta1[j] );
        v_h_etag[2]->Fill( evt_94x.diphoton_eta2[j] );
        v_h_mgg[2]->Fill( evt_94x.diphoton_mass[j] );
        v_h_ptgg[2]->Fill( evt_94x.diphoton_pt[j] );
        v_h_dphigg[2]->Fill( 1.-fabs( evt_94x.diphoton_dphi[j]/M_PI ) );
        v_h_dphigg_zoom[2]->Fill( ( 1.-fabs( evt_94x.diphoton_dphi[j]/M_PI ) )*1.e3 );
      }
    }
  }

  string prepend = "";
  if ( !is_elastic )
    prepend = "nosel_";
  unordered_map<string,vector<TH1D*> > plots = {
    { "comp_ptg", v_h_ptg },
    { "comp_etag", v_h_etag },
    { "comp_mgg", v_h_mgg },
    { "comp_ptgg", v_h_ptgg },
    { "comp_dphigg_full", v_h_dphigg },
    { "comp_dphigg", v_h_dphigg_zoom }
  };
  for ( auto& sv : plots ) {
    gStyle->SetOptStat( 111111 );
    //Canvas c( sv.first.c_str(), "2016B - 4.3 fb^{-1} (13 TeV)", "Preliminary", true );
    Canvas c( sv.first.c_str(), "9.4 fb^{-1} (13 TeV)", "Preliminary", true );
    THStack hs;
    hs.SetTitle( sv.second[0]->GetTitle() );
    //for ( size_t i = 0; i < sv.second.size(); ++i )
    for ( size_t i = 0; i < 2; ++i )
      hs.Add( sv.second[i], i < 1 ? "hist sames" : "p sames" );
    sv.second[0]->SetLineColor( kRed+1 );
    sv.second[0]->SetLineStyle( 2 );
    sv.second[0]->SetLineWidth( 3 );
/*    sv.second[1]->SetLineColor( kBlack );
    sv.second[1]->SetLineStyle( 1 );
    sv.second[1]->SetLineWidth( 2 );
    sv.second[2]->SetLineColor( kBlack );
    sv.second[2]->SetMarkerStyle( 24 );
    sv.second[2]->SetLineWidth( 2 );
    sv.second[2]->Sumw2();*/
    sv.second[1]->SetLineColor( kBlack );
    sv.second[1]->SetMarkerStyle( 24 );
    sv.second[1]->SetLineWidth( 2 );
    sv.second[1]->Sumw2();
    hs.Draw( "nostack" );
    gPad->Update();
    //for ( size_t i = 0; i < sv.second.size(); ++i ) {
    for ( size_t i = 0; i < 2; ++i ) {
      auto stat = (TPaveStats*)sv.second[i]->FindObject("stats")->Clone();
      stat->SetX1NDC( 0.6 );
      stat->SetX2NDC( 0.875 );
      stat->SetY1NDC( 0.55+(1-i)*0.18 );
      stat->SetY2NDC( 0.72+(1-i)*0.18 );
      stat->SetTextColor( sv.second[i]->GetLineColor() );
      stat->SetTextSize( 0.03 );
      stat->Draw();
    }
    c.Modified();
    c.Prettify( hs.GetHistogram() );
    if ( !is_elastic && sv.first.find( "ptg_" ) != 0 ) {
      c.GetPad( 1 )->SetLogy();
      hs.SetMaximum( hs.GetMaximum()*10. ); // should double it...
    }
    else
      hs.SetMaximum( hs.GetMaximum()*0.85 ); // should double it...
    hs.SetTitle( "" );
    c.cd( 2 );
    gStyle->SetOptStat( 0 );
    Canvas::HistsMap hm;
    hm.emplace_back( "80x", sv.second[0] );
    sv.second[0]->Sumw2();
    hm.emplace_back( "94x", sv.second[1] );
    //c.RatioPlot( hm, -0.4, 2.4, "94X/80X", 1.0 );
    c.RatioPlot( hm, -1.6, 3.6, "94X/80X", 1.0 );
    c.cd();
    auto lab = new PaveText( 0.135, 0.95, 0.2, 0.96 );
    lab->SetTextAlign( kVAlignBottom+kHAlignLeft );
    if ( is_elastic )
      lab->AddText( "#font[52]{Elastic selection}" );
    else
      lab->AddText( "#font[52]{Preselection}" );
    lab->Draw( "same" );
    c.Pad()->SetFillStyle( 0 );
    c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
  }
}
