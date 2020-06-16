#ifndef DiphotonAnalyzer_Macros_xi_accept_h
#define DiphotonAnalyzer_Macros_xi_accept_h

namespace xi_accept
{
  TF1* acc_near_arm0, *acc_near_arm1;
  void load_file( const char* filename ) {
    auto file = TFile::Open( filename );
    /*for ( size_t i = 0; i < 2; ++i ) {
      // loop on sectors
      auto func_near = file->Get<TF1>( Form( "erf_highxi_%lu", i*100+2 ) );
      auto func_far = file->Get<TF1>( Form( "erf_highxi_%lu", i*100+3 ) );
    }*/
    acc_near_arm0 = (TF1*)file->Get( "erf_highxi_2" );
    acc_near_arm1 = (TF1*)file->Get( "erf_highxi_102" );
  }

  void vary_uncertainty( double nsig ) {
    for ( auto& acc : { acc_near_arm0, acc_near_arm1 } )
      for ( int i = 0; i < acc->GetNpar(); ++i )
        acc->SetParameter( i, acc->GetParameter( i )+nsig*acc->GetParError( i ) );
  }

  double weight( double xi1, double xi2 ) {
    return ( xi1 > 0.12 ? acc_near_arm0->Eval( xi1 ) : 1. )
         * ( xi2 > 0.12 ? acc_near_arm1->Eval( xi2 ) : 1. );
  }
}

#endif
