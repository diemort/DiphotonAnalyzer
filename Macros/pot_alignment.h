#ifndef DiphotonAnalyzer_Macros_pot_alignment_h
#define DiphotonAnalyzer_Macros_pot_alignment_h

#include <iostream>
#include <fstream>
#include <regex>
#include <map>

namespace pot_align
{
  regex rgx_header( "\\[\\w+\\/fill_(\\d+)\\/(.+)" );
  regex rgx_align( "id=(\\d+),sh_x=([0-9.+-]+),sh_x_unc=([0-9.+-]+),sh_y=([0-9.+-]+),sh_y_unc=([0-9.+-]+)" );

  typedef struct {
    float x, x_unc, y, y_unc;
  } align_t;
  map<unsigned int,map<unsigned int,align_t> > alignments;

  void load_file( const char* file ) {
    ifstream f( file );
    string buf;
    smatch match;
    unsigned int fill_num = 0;
    while ( getline( f, buf ) ) {
      if ( regex_match( buf, match, rgx_header ) ) {
        fill_num = atoi( match[1].str().c_str() );
        alignments.insert( make_pair( fill_num, map<unsigned int, align_t>() ) );
        continue;
      }
      if ( regex_match( buf, match, rgx_align ) && fill_num>0 ) {
        map<unsigned int,align_t>& last_align = alignments[fill_num];
        align_t align;
        align.x = atof( match[2].str().c_str() )*1.e-3;
        align.x_unc = atof( match[3].str().c_str() )*1.e-3;
        align.y = atof( match[4].str().c_str() )*1.e-3;
        align.y_unc = atof( match[5].str().c_str() )*1.e-3;
        last_align.insert( make_pair( atoi( match[1].str().c_str() ), align ) );
      }
    }
    f.close();
  }

  map<unsigned int, align_t> get_alignments( unsigned int fill ) {
    const auto& al = alignments.find( fill );
    if ( al!=alignments.end() ) return al->second;
    if ( fill<alignments.begin()->first ) return alignments.begin()->second;
    if ( fill>alignments.rbegin()->first ) return alignments.rbegin()->second;
    return map<unsigned int,align_t>();
  }

  void print() {
    for ( const auto& fill_info : alignments ) {
      cout << ">>> fill #" << fill_info.first << endl;
      for ( const auto& pot_info : fill_info.second ) {
        cout << "    pot " << pot_info.first << " -> alignment: x=" << pot_info.second.x << " +/- " << pot_info.second.x_unc << ", y=" << pot_info.second.y << " +/- " << pot_info.second.y_unc << endl;
      }
    }
  }
}

#endif
