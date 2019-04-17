#include "Canvas.h"
#include "Plotter.h"

void plot_resol( const char* title, TH1D* h_resol, const char* top_title="" );
const char* output_dir = "/afs/cern.ch/user/l/lforthom/www/private/twophoton/mc_study";
//const char* file = "Samples/output_0000_GammaGammaToGammaGamma_13TeV_fpmc.root";
//const char* file = "Samples/output_GammaGammaToGammaGamma_fpmc_v2.root";
//const char* file = "Samples/output_GammaGammaToGammaGamma_fpmc_5jul_pufix.root";
//const char* file = "Samples/output_GammaGammaToGammaGamma_fpmc_zeta1em14.root";
//const char* file = "Samples/output_DiPhotonJetsBox_MGG-80toInf_Sherpa_v3.root";
/*const char* file = "/eos/cms/store/user/lforthom/twophoton/samples_80x/output_GammaGammaToGammaGamma_fpmc_justin_sm.root";
const string typ_lab = "FPMC, SM pred.";*/
const char* file = "/eos/cms/store/user/lforthom/twophoton/samples_80x/output_GammaGammaToGammaGamma_13TeV_superchic2_27jun.root";
const string typ_lab = "SuperChic 2, SM pred.";
unique_ptr<PaveText> lab;

void mc_study()
{
  const double pt_bins[] = { 50., 150., 500., 1000. };
  const unsigned short num_pt_bins = sizeof( pt_bins )/sizeof( pt_bins[0] )-1;
  double mean_pt_bins[num_pt_bins], err_pt_bins[num_pt_bins];
  TH1D* h_ptdist[num_pt_bins];
  for ( unsigned short i=0; i<num_pt_bins; i++ ) {
    h_ptdist[i] = new TH1D( Form( "ptdist_%d", i ), "", 100, -5., 5. );
    mean_pt_bins[i] = ( pt_bins[i+1]+pt_bins[i] )/2.;
    err_pt_bins[i] = ( pt_bins[i+1]-pt_bins[i] )/2.;
  }

  const double diph_pt_bins[] = { 0., 2., 10., 20., 50. };
  const unsigned short num_diph_pt_bins = sizeof( diph_pt_bins )/sizeof( diph_pt_bins[0] )-1;
  double mean_diph_pt_bins[num_diph_pt_bins], err_diph_pt_bins[num_diph_pt_bins];
  for ( unsigned short i=0; i<num_diph_pt_bins; i++ ) {
    mean_diph_pt_bins[i] = ( diph_pt_bins[i+1]+diph_pt_bins[i] )/2.;
    err_diph_pt_bins[i] = ( diph_pt_bins[i+1]-diph_pt_bins[i] )/2.;
  }
  TH1D* h_diph_ptdist_numer = new TH1D( "diph_ptdist_numer", "", num_diph_pt_bins, diph_pt_bins ),
       *h_diph_ptdist_denom = (TH1D*)h_diph_ptdist_numer->Clone( "diph_ptdist_denom" );

  const float dr_cut = 0.1;

  TFile f( file );
  TTree* t = dynamic_cast<TTree*>( f.Get( "ntp" ) );

  unsigned int num_dipho;
  const unsigned short max_dipho = 50;
  float photon1_pt[max_dipho], photon1_eta[max_dipho], photon1_phi[max_dipho], photon1_energy[max_dipho], photon1_r9[max_dipho],
        photon2_pt[max_dipho], photon2_eta[max_dipho], photon2_phi[max_dipho], photon2_energy[max_dipho], photon2_r9[max_dipho];
  float gen_photon1_pt[max_dipho], gen_photon1_eta[max_dipho], gen_photon1_phi[max_dipho], gen_photon1_energy[max_dipho],
        gen_photon2_pt[max_dipho], gen_photon2_eta[max_dipho], gen_photon2_phi[max_dipho], gen_photon2_energy[max_dipho];
  float dipho_vx[max_dipho], dipho_vy[max_dipho], dipho_vz[max_dipho];
  float dipho_m[max_dipho], dipho_dphi[max_dipho];
  t->SetBranchAddress( "num_diphoton", &num_dipho );
  t->SetBranchAddress( "diphoton_mass", dipho_m );
  t->SetBranchAddress( "diphoton_dphi", dipho_dphi );
  t->SetBranchAddress( "diphoton_pt1", photon1_pt );
  t->SetBranchAddress( "diphoton_eta1", photon1_eta );
  t->SetBranchAddress( "diphoton_phi1", photon1_phi );
  t->SetBranchAddress( "diphoton_energy1", photon1_energy );
  t->SetBranchAddress( "diphoton_pt2", photon2_pt );
  t->SetBranchAddress( "diphoton_eta2", photon2_eta );
  t->SetBranchAddress( "diphoton_phi2", photon2_phi );
  t->SetBranchAddress( "diphoton_energy2", photon2_energy );
  t->SetBranchAddress( "diphoton_r91", photon1_r9 );
  t->SetBranchAddress( "diphoton_r92", photon2_r9 );
  t->SetBranchAddress( "diphoton_genpt1", gen_photon1_pt );
  t->SetBranchAddress( "diphoton_geneta1", gen_photon1_eta );
  t->SetBranchAddress( "diphoton_genphi1", gen_photon1_phi );
  t->SetBranchAddress( "diphoton_genenergy1", gen_photon1_energy );
  t->SetBranchAddress( "diphoton_genpt2", gen_photon2_pt );
  t->SetBranchAddress( "diphoton_geneta2", gen_photon2_eta );
  t->SetBranchAddress( "diphoton_genphi2", gen_photon2_phi );
  t->SetBranchAddress( "diphoton_genenergy2", gen_photon2_energy );
  t->SetBranchAddress( "diphoton_vertex_x", dipho_vx );
  t->SetBranchAddress( "diphoton_vertex_y", dipho_vy );
  t->SetBranchAddress( "diphoton_vertex_z", dipho_vz );

  unsigned int num_genpho;
  const unsigned short max_genpho = 10;
  float genpho_pt[max_genpho], genpho_eta[max_genpho], genpho_phi[max_genpho], genpho_energy[max_genpho];
  t->SetBranchAddress( "num_gen_photon", &num_genpho );
  t->SetBranchAddress( "gen_photon_pt", genpho_pt );
  t->SetBranchAddress( "gen_photon_eta", genpho_eta );
  t->SetBranchAddress( "gen_photon_phi", genpho_phi );
  t->SetBranchAddress( "gen_photon_energy", genpho_energy );

  const unsigned short max_genpart = 10;
  unsigned int num_genpart;
  int genpart_pdgid[max_genpart];
  float genpart_pt[max_genpart], genpart_eta[max_genpart], genpart_phi[max_genpart], genpart_energy[max_genpart];
  t->SetBranchAddress( "num_gen_part", &num_genpart );
  t->SetBranchAddress( "gen_part_pdgid", genpart_pdgid );
  t->SetBranchAddress( "gen_part_pt", genpart_pt );
  t->SetBranchAddress( "gen_part_eta", genpart_eta );
  t->SetBranchAddress( "gen_part_phi", genpart_phi );
  t->SetBranchAddress( "gen_part_energy", genpart_energy );

  float dipho_vtx_smear_x, dipho_vtx_smear_y, dipho_vtx_smear_z;
  t->SetBranchAddress( "diphoton_genvertex_smeared_x", &dipho_vtx_smear_x );
  t->SetBranchAddress( "diphoton_genvertex_smeared_y", &dipho_vtx_smear_y );
  t->SetBranchAddress( "diphoton_genvertex_smeared_z", &dipho_vtx_smear_z );

  float pho1_sc_x[max_dipho], pho1_sc_y[max_dipho], pho1_sc_z[max_dipho];
  float pho2_sc_x[max_dipho], pho2_sc_y[max_dipho], pho2_sc_z[max_dipho];
  t->SetBranchAddress( "diphoton_supercluster_x1", pho1_sc_x );
  t->SetBranchAddress( "diphoton_supercluster_y1", pho1_sc_y );
  t->SetBranchAddress( "diphoton_supercluster_z1", pho1_sc_z );
  t->SetBranchAddress( "diphoton_supercluster_x2", pho2_sc_x );
  t->SetBranchAddress( "diphoton_supercluster_y2", pho2_sc_y );
  t->SetBranchAddress( "diphoton_supercluster_z2", pho2_sc_z );

  float pu_weight;
  t->SetBranchAddress( "pileup_weight", &pu_weight );

  TH1D* h_ptsingle_gen = new TH1D( "pt_gen", "Single photon p_{T}@@Events fraction@@GeV", 100, 0., 1000. ),
       *h_ptsingle_reco = (TH1D*) h_ptsingle_gen->Clone( "pt_reco" );
  TH1D* h_ptpair_gen = new TH1D( "ptpair_gen", "Diphoton p_{T}@@Events fraction@@GeV?.2f", 100, 0., 25. ),
       *h_ptpair_reco = (TH1D*) h_ptpair_gen->Clone( "ptpair_reco" );
  TH1D* h_etasingle_gen = new TH1D( "etasingle_gen", "Leading photon #eta@@Events fraction@@?.2f", 50, -2.5, 2.5 ),
       *h_etasingle_reco = (TH1D*) h_etasingle_gen->Clone( "etasingle_reco" );
  TH1D* h_phisingle_gen = new TH1D( "phisingle_gen", "Leading photon #phi@@Events fraction@@?.2f", 20, -TMath::Pi(), TMath::Pi() ),
       *h_phisingle_reco = (TH1D*) h_phisingle_gen->Clone( "phisingle_reco" );
  TH1D* h_energysingle_gen = new TH1D( "energysingle_gen", "Leading photon energy@@Events fraction@@GeV?.2f", 50, 0., 500. ),
       *h_energysingle_reco = (TH1D*) h_energysingle_gen->Clone( "energysingle_reco" );
  TH1D* h_ptdiff = new TH1D( "ptdiff", "p_{T}^{reco} - p_{T}^{gen}@@Events fraction@@GeV?.2f", 100, -25., 25. );
  TH1D* h_etadiff = new TH1D( "etadiff", "#eta_{T}^{reco} - #eta_{T}^{gen}@@Events fraction@@GeV?.2f", 100, -0.5, 0.5 );
  TH2D* h_ptevol = new TH2D( "ptevol", "Generated leading photon p_{T}@@Reconstructed - generated leading photon p_{T}", 250, 0., 1000., 100, -25., 25. );
  TH1D* h_vtxdist_gen = new TH1D( "vtxdist_gen", "d(generated/reconstructed diphoton vertex)@@Events fraction@@cm?.2f", 200, 0., 20. ),
       *h_vtxdist_gen_purw = (TH1D*) h_vtxdist_gen->Clone( "vtxdist_gen_purw" );
  TH1D* h_match_mass = new TH1D( "match_mass", "(m_{reco}(#gamma#gamma)-m_{gen}(#gamma#gamma))/m_{gen}(#gamma#gamma)@@Events fraction@@?.3f", 200, -0.5, 0.5 ),
       *h_match_ptpair = new TH1D( "match_pt", "(p_{T}^{reco}(#gamma#gamma)-p_{T}^{gen}(#gamma#gamma))/p_{T}^{gen}(#gamma#gamma)@@Events fraction", 100, -15., 15. );
  TH1D* h_diphvtxz_gen = new TH1D( "diphvtx_gen", "Diphoton vertex z@@Events fraction@@cm", 60, -30., 30. ),
       *h_diphvtxz_reco = (TH1D*)h_diphvtxz_gen->Clone( "diphvtx_reco" );
  TH1D* h_mresol = new TH1D( "mass_resol", "(m_{reco}-m_{gen})/m_{gen}@@Events fraction@@?.3f", 100, -0.25, 0.25 );
  TH1D* h_yresol = new TH1D( "rapidity_resol", "(y_{reco}-y_{gen})/y_{gen}@@Events fraction@@?.3f", 100, -0.25, 0.25 );
  TH1D* h_phiresol = new TH1D( "phi_resol", "(#phi_{reco}-#phi_{gen})/#phi_{gen}@@Events fraction@@?.4f", 150, -0.015, 0.015 );
  TH1D* h_xiresol1 = new TH1D( "xi_resol_p1", "(#xi_{reco}-#xi_{gen})/#xi_{gen}@@Events fraction@@?.3f", 100, -0.25, 0.25 ),
       *h_xiresol2 = (TH1D*)h_xiresol1->Clone( "xi_resol_p2" );

  const unsigned short num_vtx_pos = 11;
  const float min_vtx_dist = -10., max_vtx_dist = 10.;
  float vtx_shifts[num_vtx_pos];
  TH1D* h_mdiff_shift[num_vtx_pos];
  for ( unsigned short i=0; i<num_vtx_pos; i++ ) {
    vtx_shifts[i] = min_vtx_dist+( max_vtx_dist-min_vtx_dist )*i/( num_vtx_pos-1 );
    h_mdiff_shift[i] = new TH1D( Form( "mdiff_shift_%d", i ), "#Deltam(#gamma#gamma)/m_{gen}(#gamma#gamma)@@Events@@GeV?.2f", 2000, -0.5, 0.5 );
    //cout << "vtx shift " << i << " = " << vtx_shifts[i] << endl;
  }

  const unsigned int num_entries = t->GetEntries();
  for ( unsigned long long i=0; i<num_entries; ++i ) {
    t->GetEntry( i );
    bool has_pho1_match = false,
         has_pho2_match = false;

    TLorentzVector inpro1, inpro2, inpho1, inpho2;
    const double p_pro = sqrt( 6500.*6500.-0.938272013*0.938272013 );
    inpro1.SetPxPyPzE( 0., 0., +p_pro, 6500. );
    inpro2.SetPxPyPzE( 0., 0., -p_pro, 6500. );
    inpho1.SetPtEtaPhiE( genpart_pt[0], genpart_eta[0], genpart_phi[0], genpart_energy[0] );
    inpho2.SetPtEtaPhiE( genpart_pt[1], genpart_eta[1], genpart_phi[1], genpart_energy[1] );
    TLorentzVector outpro1 = inpro1-inpho1, outpro2 = inpro2-inpho2;
    const double xi_pro1 = 1.-outpro1.Pz()/inpro1.Pz(), xi_pro2 = 1.-outpro2.Pz()/inpro2.Pz();

    for ( unsigned int j=0; j<num_dipho; ++j ) {
      TLorentzVector pho1, pho2,
                     pho1_gen, pho2_gen;
      TVector3 dipho_vtx( dipho_vx[j], dipho_vy[j], dipho_vz[j] ),
               dipho_gen_vtx( dipho_vtx_smear_x, dipho_vtx_smear_y, dipho_vtx_smear_z );

      if ( dipho_m[j] < 350. ) continue;
      if ( photon1_pt[j] < 75. || photon2_pt[j] < 75. ) continue;
      if ( fabs( photon1_eta[j] ) > 2.5 || fabs( photon2_eta[j] ) > 2.5 ) continue;
      if ( fabs( photon1_eta[j] ) > 1.4442 && fabs( photon1_eta[j] ) < 1.566 ) continue;
      if ( fabs( photon2_eta[j] ) > 1.4442 && fabs( photon2_eta[j] ) < 1.566 ) continue;
      if ( photon1_r9[j] < 0.94 || photon2_r9[j] < 0.94 ) continue;
      if ( 1.-fabs( dipho_dphi[j]/M_PI ) > 0.005 ) continue;

      const TVector3 pho1_sc( pho1_sc_x[j], pho1_sc_y[j], pho1_sc_z[j] );
      const TVector3 pho2_sc( pho2_sc_x[j], pho2_sc_y[j], pho2_sc_z[j] );

      /*pho1.SetPtEtaPhiM( photon1_pt[j], photon1_eta[j], photon1_phi[j], 0. );
      pho2.SetPtEtaPhiM( photon2_pt[j], photon2_eta[j], photon2_phi[j], 0. );*///FIXME
      pho1.SetPtEtaPhiE( photon1_pt[j], photon1_eta[j], photon1_phi[j], photon1_energy[j] );
      pho2.SetPtEtaPhiE( photon2_pt[j], photon2_eta[j], photon2_phi[j], photon2_energy[j] );
      pho1_gen.SetPtEtaPhiE( gen_photon1_pt[j], gen_photon1_eta[j], gen_photon1_phi[j], gen_photon1_energy[j] );
      pho2_gen.SetPtEtaPhiE( gen_photon2_pt[j], gen_photon2_eta[j], gen_photon2_phi[j], gen_photon2_energy[j] );

      //cout << pho1.M() << " | " << pho1_gen.M() << "\t" << pho2.M() << " | " << pho2_gen.M() << endl;

      const float xip = ( photon1_pt[j]*exp( +photon1_eta[j] ) + photon2_pt[j]*exp( +photon2_eta[j] ) ) / 13000.,
                  xim = ( photon1_pt[j]*exp( -photon1_eta[j] ) + photon2_pt[j]*exp( -photon2_eta[j] ) ) / 13000.;
      h_xiresol1->Fill( ( xip-xi_pro1 )/xip, pu_weight/num_entries );
      h_xiresol2->Fill( ( xim-xi_pro2 )/xim, pu_weight/num_entries );

      h_mresol->Fill( ( ( pho1+pho2 ).M()-( pho1_gen+pho2_gen ).M() )/( pho1_gen+pho2_gen ).M(), pu_weight/num_entries );
      h_yresol->Fill( ( ( pho1+pho2 ).Rapidity()-( pho1_gen+pho2_gen ).Rapidity() )/( pho1_gen+pho2_gen ).Rapidity(), pu_weight/num_entries );
      h_phiresol->Fill( ( pho1.Phi()-pho1_gen.Phi() )/pho1_gen.Phi(), pu_weight/num_entries );
      h_phiresol->Fill( ( pho2.Phi()-pho2_gen.Phi() )/pho2_gen.Phi(), pu_weight/num_entries );

      h_ptsingle_reco->Fill( photon1_pt[j], pu_weight/num_entries );
      h_ptsingle_gen->Fill( gen_photon1_pt[j], pu_weight/num_entries );
      h_etasingle_reco->Fill( photon1_eta[j], pu_weight/num_entries );
      h_etasingle_gen->Fill( gen_photon1_eta[j], pu_weight/num_entries );
      h_phisingle_reco->Fill( photon1_phi[j], pu_weight/num_entries );
      h_phisingle_gen->Fill( gen_photon1_phi[j], pu_weight/num_entries );
      h_energysingle_reco->Fill( photon1_energy[j], pu_weight/num_entries );
      h_energysingle_gen->Fill( gen_photon1_energy[j], pu_weight/num_entries );

      h_ptdiff->Fill( photon1_pt[j]-gen_photon1_pt[j], pu_weight/num_entries );
      h_ptdiff->Fill( photon2_pt[j]-gen_photon2_pt[j], pu_weight/num_entries );
      h_etadiff->Fill( photon1_eta[j]-gen_photon1_eta[j], pu_weight/num_entries );
      h_etadiff->Fill( photon2_eta[j]-gen_photon2_eta[j], pu_weight/num_entries );

      h_ptevol->Fill( gen_photon1_pt[j], photon1_pt[j]-gen_photon1_pt[j], pu_weight/num_entries );

      int bin_id = -1;
      for ( unsigned int l=0; l<num_pt_bins; l++ ) {
        if ( gen_photon1_pt[j]<pt_bins[l] || gen_photon1_pt[j]>pt_bins[l+1] ) continue;
        bin_id = l;
      }
      if ( bin_id>=0 ) {
        h_ptdist[bin_id]->Fill( ( photon1_pt[j]-gen_photon1_pt[j] )/gen_photon1_pt[j], pu_weight/num_entries );
      }

      h_diphvtxz_gen->Fill( dipho_gen_vtx.z(), pu_weight/num_entries );
      h_diphvtxz_reco->Fill( dipho_vtx.z(), pu_weight/num_entries );

      const double reco_gen_dist = ( dipho_vtx-dipho_gen_vtx ).Mag();
      //const double reco_gen_dist = fabs( dipho_vz[j]-dipho_vtx_smear_z ); //FIXME
//cout << reco_gen_dist << "\t" << dipho_vz[j] << "\t" << dipho_vtx_smear_z << endl;
      h_vtxdist_gen->Fill( reco_gen_dist, 1./num_entries );
      h_vtxdist_gen_purw->Fill( reco_gen_dist, pu_weight/num_entries );

      const TLorentzVector pair_gen( pho1_gen+pho2_gen ),
                           pair_reco( pho1+pho2 );

      h_diph_ptdist_numer->Fill( pair_reco.Pt(), pu_weight/num_entries );
      h_diph_ptdist_denom->Fill( pair_gen.Pt(), pu_weight/num_entries );

      h_ptpair_gen->Fill( pair_gen.Pt(), pu_weight/num_entries );
      h_ptpair_reco->Fill( pair_reco.Pt(), pu_weight/num_entries );
      h_match_mass->Fill( ( pair_reco.M()-pair_gen.M() )/pair_gen.M(), pu_weight/num_entries );
      h_match_ptpair->Fill( ( pair_reco.Pt()-pair_gen.Pt() )/pair_gen.Pt(), pu_weight/num_entries );

      const double m_gen = ( pho1_gen+pho2_gen ).M();
      for ( unsigned short k=0; k<num_vtx_pos; k++ ) {
        const TVector3 shift_pos( 0., 0., vtx_shifts[k] );
	//dipho_gen_vtx.Print(); ( dipho_gen_vtx+shift_pos ).Print();
        TVector3 p1_shift = ( pho1_sc-( dipho_gen_vtx+shift_pos ) ).Unit() * photon1_energy[j];
        TVector3 p2_shift = ( pho2_sc-( dipho_gen_vtx+shift_pos ) ).Unit() * photon2_energy[j];
        TLorentzVector corr_p1_shift( p1_shift.x(), p1_shift.y(), p1_shift.z(), photon1_energy[j] );
        TLorentzVector corr_p2_shift( p2_shift.x(), p2_shift.y(), p2_shift.z(), photon2_energy[j] );
	const double m_reco = ( corr_p1_shift+corr_p2_shift ).M();
        h_mdiff_shift[k]->Fill( ( m_reco-m_gen ) / m_gen );
	//cout << ( pho1_gen+pho2_gen ).M() << "\t" << ( corr_p1_shift+corr_p2_shift ).M() << endl;
      }
    }
  }
  gStyle->SetOptStat( 0 );
  //lab.reset( new PaveText( 0.8, 0.6, 0.87, 0.7 ) );
  lab.reset( new PaveText( 0.8, 0.8, 0.87, 0.9 ) );
  //PaveText lab( 0.8, 0.6, 0.85, 0.7 );
  lab->SetTextSize( 0.035 );
  lab->SetFillStyle( 0 );
  lab->SetLineWidth( 0 );
  lab->AddText( "Elastic #gamma#gamma#rightarrow#gamma#gamma" );
  lab->AddText( typ_lab.c_str() );
  {
    Canvas c( "mc_study_mresol_vs_vtxpos", "2016 PU cond. (13 TeV)", "Simulation" );
    TGraphErrors gr;
    for ( unsigned short i=0; i<num_vtx_pos; i++ ) {
      const double mean = h_mdiff_shift[i]->GetMean();
      const double rms = h_mdiff_shift[i]->GetRMS();
      const unsigned short n = gr.GetN();
      cout << vtx_shifts[i] << "\t" << mean << "\t" << rms << endl;
      gr.SetPoint( n, vtx_shifts[i], mean );
      gr.SetPointError( n, 0., rms );
    }
    gr.Draw( "ap" );
    gr.SetMarkerStyle( 20 );
    gr.GetHistogram()->SetTitle( "z shift to the generated vertex position (cm)@@m_{gen}(#gamma#gamma)-m_{reco}(#gamma#gamma) (GeV)" );
    gr.GetYaxis()->SetRangeUser( -0.1, 0.1 );
    c.Prettify( gr.GetHistogram() );
    c.Save( "pdf,png", output_dir );
  }
  for ( unsigned short i=0; i<num_vtx_pos; i++ ) {
    Canvas c( Form( "xxx_%d", i ), "2016 PU cond. (13 TeV)", "Simulation" );
    gStyle->SetOptStat( 1111 );
    h_mdiff_shift[i]->Draw();
    c.Prettify( h_mdiff_shift[i] );
    c.Save( "pdf,png", output_dir );
  }
  gStyle->SetOptStat( 0 );
  {
    Canvas c( "mc_study_ptevol", "2016 PU cond. (13 TeV)", "Simulation" );
    h_ptevol->Draw( "colz" );
    c.Prettify( h_ptevol );
    c.Save( "pdf,png", output_dir );
  }
  {
    Canvas c( "mc_study_ptresol", "2016 PU cond. (13 TeV)", "Simulation" );
    double pt_error_mean[num_pt_bins], pt_error_rms[num_pt_bins];
    for ( unsigned short i=0; i<num_pt_bins; i++ ) {
      pt_error_mean[i] = h_ptdist[i]->GetMean();
      pt_error_rms[i] = h_ptdist[i]->GetRMS();
    }
    TGraphErrors gr( num_pt_bins, mean_pt_bins, pt_error_mean, err_pt_bins, pt_error_rms );
    gr.Draw( "ap" );
    gr.SetMarkerStyle( 20 );
    gr.SetTitle( "Single photon p_{T} region@@(p_{T}^{reco}-p_{T}^{gen})/p_{T}^{gen}" );
    gr.SetLineWidth( 2 );
    gr.GetXaxis()->SetRangeUser( pt_bins[0], pt_bins[num_pt_bins] );
    c.Prettify( gr.GetHistogram() );
    c.SetGrid();
    c.Save( "pdf,png", output_dir );
  }
  {
    gStyle->SetOptStat( 1111 );
    Canvas c( "mc_study_match_mass", "2016 PU cond. (13 TeV)", "Simulation" );
    h_match_mass->Draw( "hist" );
    h_match_mass->SetMinimum( 1.e-4 );
    c.Prettify( h_match_mass );
    c.SetLogy();
    c.Save( "pdf,png", output_dir );
    gStyle->SetOptStat( 0 );
  }
  {
    Canvas c( "mc_study_ptpair_ratio", "2016 PU cond. (13 TeV)", "Simulation" );
    h_diph_ptdist_numer->Divide( h_diph_ptdist_denom );
    h_diph_ptdist_numer->Sumw2();
    h_diph_ptdist_numer->SetMarkerStyle( 20 );
    h_diph_ptdist_numer->Draw( "p" );
    c.Prettify( h_diph_ptdist_numer );
    c.Save( "pdf,png", output_dir );
  }
  {
    Canvas c( "mc_study_match_ptpair", "2016 PU cond. (13 TeV)", "Simulation" );
    h_match_ptpair->Draw( "hist" );
    c.Prettify( h_match_ptpair );
    c.SetLogy();
    c.Save( "pdf,png", output_dir );
  }
  gStyle->SetOptStat( 11111 );
  {
    Canvas c( "mc_study_ptdiff", "2016 PU cond. (13 TeV)", "Simulation" );
    h_ptdiff->Draw();
    c.Prettify( h_ptdiff );
    c.Save( "pdf,png", output_dir );
  }
  {
    Canvas c( "mc_study_etadiff", "2016 PU cond. (13 TeV)", "Simulation" );
    h_etadiff->Draw();
    c.Prettify( h_etadiff );
    c.Save( "pdf,png", output_dir );
  }
  gStyle->SetOptStat( 0 );
  {
    Canvas c( "mc_study_genvtx_dist", "2016 PU cond. (13 TeV)", "Simulation" );
    h_vtxdist_gen->Draw( "hist" );
    c.Prettify( h_vtxdist_gen );
    c.SetLogy();
    c.Save( "png,pdf", output_dir );
  }
  {
    //Plotter plt( output_dir, "CMS Simulation - Excl. #gamma#gamma #rightarrow #gamma#gamma, #sqrt{s} = 13 TeV, 2016 PU cond." );
    bool ltor = true;
    TH1D* cumul = dynamic_cast<TH1D*>( h_vtxdist_gen->GetCumulative( ltor ) ),
         *cumul_purw = dynamic_cast<TH1D*>( h_vtxdist_gen_purw->GetCumulative( ltor ) );
    Canvas c( "mc_study_genvtx_dist_cumul", "2016 PU cond. (13 TeV)", "Simulation" );
    THStack st;
    cumul->SetLineColor( kBlack );
    cumul->SetLineWidth( 2 );
    cumul->SetFillColor( kGray );
    cumul->SetFillStyle( 0 );
    st.Add( cumul, "hist" );
    c.AddLegendEntry( cumul, "No PU reweighting", "f" );
    //cumul_purw->SetMarkerStyle( 24 );
    cumul_purw->SetLineColor( kRed+1 );
    cumul_purw->SetLineWidth( 2 );
    st.Add( cumul_purw, "hist" );
    c.AddLegendEntry( cumul_purw, "PU reweighting", "l" );
    st.Draw( "nostack" );
    st.SetTitle( "d(diphoton vtx^{gen}/diphoton vtx^{reco})@@Events fraction (cumulative)@@cm?.2f" );
    //c.SetLegendY1( 0.1 );
    st.SetMaximum( 1.35 );
    c.Prettify( st.GetHistogram() );
    st.SetTitle("");
    c.SetGrid();
    c.Save( "pdf,png", output_dir );
    //Plotter::HistsMap hm;
    //plt.BaseCanvas()->SetLegendY1( 0.1 );
    //hm.push_back( make_pair( "No PU reweighting", cumul ) );
    //hm.push_back( make_pair( "PU reweighting", cumul_purw ) );
    //plt.draw_multiplot( "mc_study_genvtx_dist_cumul", hm );
  }
  const char* top_title = "2016 PU cond. (13 TeV)";
  Plotter plt( output_dir, top_title );
  {
    Plotter::HistsMap hm;
    hm.push_back( make_pair( "Gen level", h_ptpair_gen ) );
    hm.push_back( make_pair( "Reco level", h_ptpair_reco ) );
    plt.draw_multiplot( "mc_study_ptpair_comp", hm, false, true, 0.00005 );
  }
  {
    Plotter::HistsMap hm;
    hm.push_back( make_pair( "Gen level", h_diphvtxz_gen ) );
    hm.push_back( make_pair( "Reco level", h_diphvtxz_reco ) );
    plt.draw_multiplot( "mc_study_diphoton_vtxz", hm );
  }
  {
    Plotter::HistsMap hm;
    hm.push_back( make_pair( "Gen level", h_ptsingle_gen ) );
    hm.push_back( make_pair( "Reco level", h_ptsingle_reco ) );
    plt.draw_multiplot( "mc_study_ptcomp", hm );
  }
  {
    Plotter::HistsMap hm;
    hm.push_back( make_pair( "Gen level", h_etasingle_gen ) );
    hm.push_back( make_pair( "Reco level", h_etasingle_reco ) );
    plt.draw_multiplot( "mc_study_etacomp", hm );
  }
  {
    Plotter::HistsMap hm;
    hm.push_back( make_pair( "Gen level", h_phisingle_gen ) );
    hm.push_back( make_pair( "Reco level", h_phisingle_reco ) );
    plt.draw_multiplot( "mc_study_phicomp", hm );
  }
  {
    Plotter::HistsMap hm;
    hm.push_back( make_pair( "Gen level", h_energysingle_gen ) );
    hm.push_back( make_pair( "Reco level", h_energysingle_reco ) );
    plt.draw_multiplot( "mc_study_energycomp", hm );
  }
  plot_resol( "mc_study_massresolution", h_mresol, top_title );
  plot_resol( "mc_study_rapidityresolution", h_yresol, top_title );
  plot_resol( "mc_study_singlephiresolution", h_phiresol, top_title );

  {
    Canvas c( "mc_study_xiresolution", "2016 PU cond. (13 TeV)", "Simulation" );
    c.SetLegendX1( 0.16 );
    c.SetLegendY1( 0.58 );
    c.SetLegendSizeY( 0.25 );
    h_xiresol1->Draw( "hist" );
    h_xiresol1->SetLineColor( kBlack );
    h_xiresol1->SetLineWidth( 2 );
    h_xiresol1->SetFillColor( kBlack );
    h_xiresol1->SetFillStyle( 3004 );
    //h_xiresol1->SetMaximum( 1.15e3 );
    //FIXME h_xiresol1->SetMaximum( 0.19 );
    h_xiresol1->SetMaximum( h_xiresol1->GetMaximum()*1.2 );
    TFitResultPtr fit1 = h_xiresol1->Fit( "gaus", "NS" );
    if ( fit1.Get() ) {
      const double* params_fit1 = fit1->GetParams();
      //c.AddLegendEntry( h_xiresol1, Form( "Sect.45, #sigma = %.4f", params_fit1[2] ), "f" );
      //c.AddLegendEntry( h_xiresol1, Form( "Sect.45, Mean = %.2e, RMS = %.3f", h_xiresol1->GetMean(), h_xiresol1->GetRMS() ), "f" );
      c.AddLegendEntry( h_xiresol1, "#font[72]{Sector 45}", "f" );
      c.AddLegendEntry( nullptr, Form( "Mean = %.2e", h_xiresol1->GetMean() ), "" );
      c.AddLegendEntry( nullptr, Form( "RMS = %.g", h_xiresol1->GetRMS() ), "" );
    }
    h_xiresol2->Draw( "hist same" );
    h_xiresol2->SetLineColor( kRed+1 );
    h_xiresol2->SetFillStyle( 3005 );
    h_xiresol2->SetFillColor( kRed+1 );
    h_xiresol2->SetLineWidth( 2 );
    //h_xiresol2->SetLineStyle( 2 );
    TFitResultPtr fit2 = h_xiresol2->Fit( "gaus", "NS" );
    if ( fit2.Get() ) {
      const double* params_fit2 = fit2->GetParams();
      //c.AddLegendEntry( h_xiresol2, Form( "Sect.56, #sigma = %.4f", params_fit2[2] ), "f" );
      //c.AddLegendEntry( h_xiresol2, Form( "Sect.56, Mean = %.2e, RMS = %.3f", h_xiresol2->GetMean(), h_xiresol2->GetRMS() ), "f" );
      c.AddLegendEntry( h_xiresol2, "#font[72]{Sector 56}", "f" );
      c.AddLegendEntry( nullptr, Form( "Mean = %.2e", h_xiresol2->GetMean() ), "" );
      c.AddLegendEntry( nullptr, Form( "RMS = %.g", h_xiresol2->GetRMS() ), "" );
    }
    c.Prettify( h_xiresol1 );

    if ( lab ) lab->Draw( "same" );

    c.GetLegend()->SetLineColor( kWhite );
    c.Save( "pdf,png", output_dir );
  }
}

void plot_resol( const char* title, TH1D* h_resol, const char* top_title="" ) {
  Canvas c( title, top_title, "Simulation" );
  c.SetLegendX1( 0.16 );
  //c.SetLegendY1( 0.6 );
  //c.SetLegendSizeY( 0.3 );
  h_resol->Draw( "hist" );
  h_resol->SetLineColor( kBlack );
  h_resol->SetLineWidth( 2 );
  h_resol->SetFillColor( kBlack );
  h_resol->SetFillStyle( 3004 );
  //TFitResultPtr fit1 = h_resol->Fit( "gaus", "NS" );
  TFitResultPtr fit1 = h_resol->Fit( "gaus", "S+" );
  if ( fit1.Get() ) {
    const double* params_fit1 = fit1->GetParams();
    c.AddLegendEntry( h_resol, Form( "Mean = %.1e", h_resol->GetMean() ), "f" );
    c.AddLegendEntry( nullptr, Form( "RMS = %.g", h_resol->GetRMS() ), "" );
    //c.AddLegendEntry( h_resol, Form( "Mean = %.1e, RMS = %.3f", h_resol->GetMean(), h_resol->GetRMS() ), "f" );
    //c.AddLegendEntry( h_resol, Form( "Mean = %.2e, RMS = %.3f", params_fit1[1], params_fit1[2] ), "f" );
  }
  c.Prettify( h_resol );

  /*PaveText lab( 0.8, 0.6, 0.87, 0.7 );
  lab.SetTextSize( 0.035 );
  lab.SetFillStyle( 0 );
  lab.SetLineWidth( 0 );
  lab.AddText( "Elastic #gamma#gamma#rightarrow#gamma#gamma" );
  lab.AddText( typ_lab.c_str() );*/
  if ( lab )
    lab->Draw( "same" );

  c.GetLegend()->SetY1( 0.725 );
  c.GetLegend()->SetY2( 0.825 );
  c.Save( "pdf,png", output_dir );
}
