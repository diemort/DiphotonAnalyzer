#ifndef DiphotonAnalyzer_Macro_diproton_candidate_h
#define DiphotonAnalyzer_Macro_diproton_candidate_h

#define sqrt_s 13.0e3 //FIXME

class diproton_candidate_t
{
  public:
    diproton_candidate_t( float xi1, float xi1_err, float xi2, float xi2_err ) :
      xi1_( xi1 ), xi2_( xi2 ), xi1_err_( xi1_err ), xi2_err_( xi2_err ) {}
    double mass() const { return sqrt_s * sqrt( xi1_*xi2_ ); }
    double mass_error() const { return mass()*0.5 * sqrt( pow( xi1_err_/xi1_, 2 ) + pow( xi2_err_/xi2_, 2 ) ); }
    double rapidity() const { return 0.5*log( xi1_/xi2_ ); }
    double rapidity_error() const { return sqrt( pow( xi1_err_/xi1_, 2 ) + pow( xi2_err_/xi2_, 2 ) ) * 0.5; }
  private:
    float xi1_, xi2_, xi1_err_, xi2_err_;
};

bool is_matched( int n_sigma, float x1, float x2, float err_x1, float err_x2 )
{
  const double combined_error = std::hypot( err_x1, err_x2 );
  const double delta = fabs( x1-x2 );
  return ( delta/combined_error <= n_sigma );
}

#endif
