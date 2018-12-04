#ifndef DiphotonAnalyzer_Macros_tree_reader_h
#define DiphotonAnalyzer_Macros_tree_reader_h

#include "TTree.h"
#include <memory>

#define MAX_PROTONS 30
#define MAX_DIPH 20
#define MAX_ELE 50
#define MAX_MUON 50
#define MAX_JET 200

struct treeinfo
{
  unsigned int fill_number;
  unsigned int num_proton_track;
  float proton_track_x[MAX_PROTONS], proton_track_y[MAX_PROTONS];
  float proton_track_normchi2[MAX_PROTONS];
  unsigned int proton_track_arm[MAX_PROTONS], proton_track_pot[MAX_PROTONS];
  unsigned int num_diphoton;
  float diphoton_pt1[MAX_DIPH], diphoton_pt2[MAX_DIPH];
  float diphoton_eta1[MAX_DIPH], diphoton_eta2[MAX_DIPH];
  float diphoton_r91[MAX_DIPH], diphoton_r92[MAX_DIPH];
  float diphoton_mass[MAX_DIPH], diphoton_pt[MAX_DIPH], diphoton_rapidity[MAX_DIPH], diphoton_dphi[MAX_DIPH];
  float diphoton_vertex_x[MAX_DIPH], diphoton_vertex_y[MAX_DIPH], diphoton_vertex_z[MAX_DIPH];
  unsigned int num_electron;
  float electron_vtx_x[MAX_ELE], electron_vtx_y[MAX_ELE], electron_vtx_z[MAX_ELE];
  unsigned int num_muon;
  float muon_vtx_x[MAX_ELE], muon_vtx_y[MAX_ELE], muon_vtx_z[MAX_ELE];
  unsigned int num_jet;
  float jet_pt[MAX_JET];
  int jet_dipho_match[MAX_JET];

  void read( TTree* tr ) {
    tree = unique_ptr<TTree>( tr );
    tree->SetBranchAddress( "fill_number", &fill_number );
    tree->SetBranchAddress( "num_proton_track", &num_proton_track );
    tree->SetBranchAddress( "proton_track_x", proton_track_x );
    tree->SetBranchAddress( "proton_track_y", proton_track_y );
    tree->SetBranchAddress( "proton_track_arm", proton_track_arm );
    tree->SetBranchAddress( "proton_track_pot", proton_track_pot );
    tree->SetBranchAddress( "proton_track_normchi2", proton_track_normchi2 );
    tree->SetBranchAddress( "num_diphoton", &num_diphoton );
    tree->SetBranchAddress( "diphoton_pt1", diphoton_pt1 );
    tree->SetBranchAddress( "diphoton_pt2", diphoton_pt2 );
    tree->SetBranchAddress( "diphoton_eta1", diphoton_eta1 );
    tree->SetBranchAddress( "diphoton_eta2", diphoton_eta2 );
    tree->SetBranchAddress( "diphoton_r91", diphoton_r91 );
    tree->SetBranchAddress( "diphoton_r92", diphoton_r92 );
    tree->SetBranchAddress( "diphoton_mass", diphoton_mass );
    tree->SetBranchAddress( "diphoton_pt", diphoton_pt );
    tree->SetBranchAddress( "diphoton_rapidity", diphoton_rapidity );
    tree->SetBranchAddress( "diphoton_dphi", diphoton_dphi );
    tree->SetBranchAddress( "diphoton_vertex_x", diphoton_vertex_x );
    tree->SetBranchAddress( "diphoton_vertex_y", diphoton_vertex_y );
    tree->SetBranchAddress( "diphoton_vertex_z", diphoton_vertex_z );
    tree->SetBranchAddress( "num_electron", &num_electron );
    tree->SetBranchAddress( "electron_vtx_x", electron_vtx_x );
    tree->SetBranchAddress( "electron_vtx_y", electron_vtx_y );
    tree->SetBranchAddress( "electron_vtx_z", electron_vtx_z );
    tree->SetBranchAddress( "num_muon", &num_muon );
    tree->SetBranchAddress( "muon_vtx_x", muon_vtx_x );
    tree->SetBranchAddress( "muon_vtx_y", muon_vtx_y );
    tree->SetBranchAddress( "muon_vtx_z", muon_vtx_z );
    tree->SetBranchAddress( "num_jet", &num_jet );
    tree->SetBranchAddress( "jet_pt", jet_pt );
    tree->SetBranchAddress( "jet_dipho_match", jet_dipho_match );
  }
  unique_ptr<TTree> tree;
};

#endif
