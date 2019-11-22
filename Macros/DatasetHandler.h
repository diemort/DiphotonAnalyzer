#ifndef DatasetHandler_h
#define DatasetHandler_h

#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"
#include "TColor.h"
#include "TROOT.h"

#include <fstream>
#include <vector>
#include <streambuf>

#include "EventFilter/Utilities/interface/reader.h"
#include "FWCore/Utilities/interface/Exception.h"

class Sample
{
  public:
    typedef enum {
      kInvalidType = 0,
      kData,
      kSignal,
      kBackground
    } SampleType;

  public:
    Sample() : name_( "[no name]" ), type_( kInvalidType ), tree_( 0 ), xsection_( -1. ), nevts_( 0 ), colour_( kWhite ) {}
    Sample( const char* name, SampleType type, TTree* tree, float xsection, unsigned int nevts, const Color_t& col ) :
      name_( name ), type_( type ), tree_( tree ), xsection_( xsection ), nevts_( nevts ), colour_( col ) {}
    ~Sample() {}

    const std::string& name() const { return name_; }
    SampleType type() const { return type_; }
    TTree* tree() { return tree_; }
    float cross_section() const { return xsection_; }
    unsigned int num_events() const { return nevts_; }
    Color_t colour() const { return colour_; }

  private:
    std::string name_;
    SampleType type_;
    TTree* tree_;
    float xsection_;
    unsigned int nevts_;
    Color_t colour_;
};

class DatasetHandler
{
  public:
    DatasetHandler() {}
    DatasetHandler( const char* json_file ) {
      std::ifstream file( json_file );
      std::string str( ( std::istreambuf_iterator<char>( file ) ), std::istreambuf_iterator<char>() );
      Json::Reader reader;
      Json::Value tree;

      reader.parse( str, tree );
      if ( !tree.isMember( "datasets" ) ) throw cms::Exception( "InvalidJSON" ) << "The datasets collection was not found!";
      Json::Value ds_list = tree["datasets"];
      for ( unsigned short i=0; i<ds_list.size(); i++ ) {
        Json::Value ds = ds_list[i];
        Sample::SampleType ds_type = Sample::kInvalidType;
        if ( !ds.isMember( "path" ) || !ds["path"].isString() ) {
          continue;
        }

        // retrieve the TTree
        const char* ds_path = ds["path"].asCString();
        if ( !strcmp( ds_path, "" ) ) {
          throw cms::Exception( "InvalidFile" ) << "Invalid file name provided: " << ds_path;
        }
        TFile* file = new TFile( ds_path );
        files_.push_back( file );

        // retrieve the sample name
        std::string name = "[no name]";
        if ( ds.isMember( "name" ) && ds["name"].isString() ) {
          name = ds["name"].asString();
        }

        // retrieve the sample cross section
        float xsec = -1.;
        if ( ds.isMember( "cross_section" ) ) {
          xsec = ds["cross_section"].asDouble();
        }

        // retrieve the sample size
        unsigned int nevts = 0;
        if ( ds.isMember( "num_events" ) ) {
          nevts = ds["num_events"].asUInt();
        }

        // retrieve the sample plotting colour
        //Color_t col( static_cast<TColor>( *gROOT->GetColor( kWhite ) ) );
        Color_t col = kWhite;
        if ( ds.isMember( "colour" ) ) {
          col = static_cast<Color_t>( ds["colour"].asUInt() );
          //col = static_cast<TColor>( *gROOT->GetColor( ds["colour"].asUInt() ) );
        }

        // retrieve the sample type
        if ( ds.isMember( "type" ) && ds["type"].isString() ) {
          const std::string type = ds["type"].asString();
          if ( type=="data" ) { ds_type = Sample::kData; }
          else if ( type=="signal" || type=="sig" ) { ds_type = Sample::kSignal; }
          else if ( type=="background" || type=="bck" ) { ds_type = Sample::kBackground; }
        }
        samples_.emplace_back( name.c_str(), ds_type, dynamic_cast<TTree*>( file->Get( "ntp" ) ), xsec, nevts, col );
      }
    }
    ~DatasetHandler() {
      for ( unsigned short i=0; i<files_.size(); i++ ) {
        files_[i]->Close();
        delete files_[i];
      }
      files_.clear();
    }

    unsigned short size() const { return samples_.size(); }
    Sample& sample( const unsigned short i ) {
      if ( i>samples_.size() ) throw cms::Exception( "AboveRange" ) << "Trying to access a sample above range";
      return samples_[i];
    }

  private:
    std::vector<TFile*> files_;
    std::vector<Sample> samples_;

};

#endif
