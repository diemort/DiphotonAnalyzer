#include "TFile.h"
#include "TH2.h"
#include <memory>

class PhotonScalesParser
{
  public:
    PhotonScalesParser( const char* filename, const char* hist_name = "EGamma_SF2D" )
    {
      TFile file( filename );
      hist_effmc_.reset( (TH2D*)file.Get( hist_name )->Clone() );
    }

    float efficiency( float pt, float eta_sc ) const {
      if ( !hist_effmc_ ) return 0.;
      int bin_x = min( max( hist_effmc_->GetXaxis()->FindBin( eta_sc ), 1 ), hist_effmc_->GetXaxis()->GetNbins() ),
          bin_y = min( max( hist_effmc_->GetYaxis()->FindBin( pt ), 1 ), hist_effmc_->GetYaxis()->GetNbins() );
      printf( "%d:%d\n", bin_x, bin_y );
      return hist_effmc_->GetBinContent( bin_x, bin_y );
    }

  private:
    std::unique_ptr<TH2D> hist_effmc_;
};
