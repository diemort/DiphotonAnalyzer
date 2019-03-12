#include "Macros/Canvas.h"
#include "TreeProducer/interface/TreeEvent.h"
//#include "/afs/cern.ch/work/l/lforthom/private/twophoton/CMSSW_8_0_26_patch1/src/DiphotonAnalyzer/TreeProducer/interface/TreeEvent.h"
#include "../../../CMSSW_8_0_26_patch1/src/DiphotonAnalyzer/TreeProducer/interface/TreeEvent.h"

void validate_94xsample()
{
  //auto file_94x = TFile::Open( "samples/ntuple-Run2016B_94Xrereco_v2.root" );
  auto file_94x = TFile::Open( "/eos/cms/store/group/phys_pps/diphoton/DoubleEG/ntuple-Run2016BCG_94Xrereco_v1.root" );
  auto file_80x = TFile::Open( "../../../CMSSW_8_0_26_patch1/src/DiphotonAnalyzer/Samples/output_Run2016BCG_looseCuts_28jun.root" );
  auto tree_94x = dynamic_cast<TTree*>( file_94x->Get( "ntp" ) );
  auto tree_80x = dynamic_cast<TTree*>( file_80x->Get( "ntp" ) );
  gggg::TreeEvent evt_94x;
  evt_94x.attach( tree_94x, true );
  TreeEvent evt_80x;
  evt_80x.attach( tree_80x, true );

  const float eta_cut = 2.5, min_etaveto = 1.4442, max_etaveto = 1.566;

  vector<const char*> v_type = { "8_0_26_patch1", "9_4_9" };
  vector<TH1D*> v_h_ptg, v_h_etag, v_h_mgg, v_h_ptgg, v_h_dphigg;
  for ( unsigned short i = 0; i < 2; ++i ) {
    v_h_ptg.emplace_back( new TH1D( Form( "ptg_%s", v_type[i] ), "p_{T}^{#gamma}@@Events@@GeV", 50, 50., 550. ) );
    v_h_etag.emplace_back( new TH1D( Form( "etag_%s", v_type[i] ), "#eta^{#gamma}@@Events@@?.g", 25, -2.5, 2.5 ) );
    v_h_mgg.emplace_back( new TH1D( Form( "mgg_%s", v_type[i] ), "m_{#gamma#gamma}@@Events@@GeV", 50, 300., 2300. ) );
    v_h_ptgg.emplace_back( new TH1D( Form( "ptgg_%s", v_type[i] ), "p_{T}^{#gamma#gamma}@@Events@@GeV", 50, 0., 250. ) );
    v_h_dphigg.emplace_back( new TH1D( Form( "dphigg_%s", v_type[i] ), "1-|#Delta#phi(#gamma#gamma)/#pi| (#times 10^{-3})@@Events@@?.g", 10, 0., 5. ) );
  }

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
      if ( 1.-fabs( evt_80x.diphoton_dphi[j] )/M_PI > 0.005 ) continue;

      //----- ELASTIC REGION
      v_h_ptg[0]->Fill( evt_80x.diphoton_pt1[j] );
      v_h_ptg[0]->Fill( evt_80x.diphoton_pt2[j] );
      v_h_etag[0]->Fill( evt_80x.diphoton_eta1[j] );
      v_h_etag[0]->Fill( evt_80x.diphoton_eta2[j] );
      v_h_mgg[0]->Fill( evt_80x.diphoton_mass[j] );
      v_h_ptgg[0]->Fill( evt_80x.diphoton_pt[j] );
      v_h_dphigg[0]->Fill( ( 1.-fabs( evt_80x.diphoton_dphi[j]/M_PI ) )*1.e3 );
    }
  }

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
      if ( 1.-fabs( evt_94x.diphoton_dphi[j] )/M_PI > 0.005 ) continue;

      //----- ELASTIC REGION
      v_h_ptg[1]->Fill( evt_94x.diphoton_pt1[j] );
      v_h_ptg[1]->Fill( evt_94x.diphoton_pt2[j] );
      v_h_etag[1]->Fill( evt_94x.diphoton_eta1[j] );
      v_h_etag[1]->Fill( evt_94x.diphoton_eta2[j] );
      v_h_mgg[1]->Fill( evt_94x.diphoton_mass[j] );
      v_h_ptgg[1]->Fill( evt_94x.diphoton_pt[j] );
      v_h_dphigg[1]->Fill( ( 1.-fabs( evt_94x.diphoton_dphi[j]/M_PI ) )*1.e3 );
    }
  }

  unordered_map<const char*,vector<TH1D*> > plots = {
    { "comp_ptg", v_h_ptg },
    { "comp_etag", v_h_etag },
    { "comp_mgg", v_h_mgg },
    { "comp_ptgg", v_h_ptgg },
    { "comp_dphigg", v_h_dphigg }
  };
  for ( auto& sv : plots ) {
    gStyle->SetOptStat( 111111 );
    //Canvas c( sv.first, "2016B - 4.3 fb^{-1} (13 TeV)", "Preliminary", true );
    Canvas c( sv.first, "9.41 fb^{-1} (13 TeV)", "Preliminary", true );
    THStack hs;
    hs.SetTitle( sv.second[0]->GetTitle() );
    for ( size_t i = 0; i < sv.second.size(); ++i )
      hs.Add( sv.second[i], i > 0 ? "p sames" : "hist sames" );
    sv.second[0]->SetLineColor( kRed+1 );
    sv.second[0]->SetLineStyle( 2 );
    sv.second[0]->SetLineWidth( 3 );
    sv.second[1]->SetLineColor( kBlack );
    sv.second[1]->SetMarkerStyle( 24 );
    sv.second[1]->SetLineWidth( 2 );
    sv.second[1]->Sumw2();
    hs.Draw( "nostack" );
    gPad->Update();
    for ( size_t i = 0; i < sv.second.size(); ++i ) {
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
    hs.SetMaximum( hs.GetMaximum()*0.85 ); // should double it...
    hs.SetTitle( "" );
    c.cd( 2 );
    gStyle->SetOptStat( 0 );
    Canvas::HistsMap hm;
    hm.emplace_back( "80x", sv.second[0] );
    sv.second[0]->Sumw2();
    hm.emplace_back( "94x", sv.second[1] );
    c.RatioPlot( hm, -0.4, 2.4, "94X/80X", 1.0 );
    c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
  }
}
