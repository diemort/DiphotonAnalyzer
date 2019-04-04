#ifndef Canvas_h
#define Canvas_h

#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TH1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TObjArray.h"
#include "TObjString.h"

#include <memory>
#include <string.h>

#define EXPERIMENT "CMS"

class PaveText : public TPaveText
{
 public:
  inline PaveText( float x1, float y1, float x2=-999., float y2=-999. ) :
    TPaveText( x1, y1, ( x2 != -999. ) ? x2 : x1+0.25, ( y2 != -999. ) ? y2 : y1+0.15, "NB NDC" ) {
    TPaveText::SetFillStyle( 0 );
    TPaveText::SetFillColor( 0 );
    TPaveText::SetLineWidth( 0 );
    TPaveText::SetLineStyle( 0 );
    TPaveText::SetTextFont( 42 );
    TPaveText::SetTextSize( 0.033 );
    TPaveText::SetTextAlign( kHAlignRight+kVAlignBottom );
  }
};

class Canvas : public TCanvas
{
 public:
  enum struct Align { left = 0, right = 1 };
  inline Canvas( const char* name, const char* title = "", const char* type = "Preliminary", bool ratio = false, Align label_align = Align::left ) :
    //TCanvas(name, "", 450, 450),
    TCanvas( name, "", 600, 600 ),
    fTitle( title ), fPlotType( type ),
    fLegXpos( 0.5 ), fLegYpos( 0.775 ), fLegXsize( 0.35 ), fLegYsize( 0.15 ),
    fRatio( ratio ), divided_( false ), align_( label_align )
  {
    Build();
  }

  inline void Prettify( TH1* obj ) {
    TAxis* x = dynamic_cast<TAxis*>( obj->GetXaxis() ),
          *y = dynamic_cast<TAxis*>( obj->GetYaxis() );
    x->SetLabelFont( 43 );
    x->SetTitleFont( 43 );
    y->SetLabelFont( 43 );
    y->SetTitleFont( 43 );
    if ( !divided_ ) {
      x->SetLabelSize( 20 );
      x->SetTitleSize( 26 );
      y->SetLabelSize( 20 );
      y->SetTitleSize( 26 );
      y->SetTitleOffset( 1.4 );
    }
    else {
      x->SetLabelSize( 16 );
      x->SetTitleSize( 20 );
      y->SetLabelSize( 15 );
      y->SetTitleSize( 20 );
      x->SetTitleOffset( 1.8 );
      x->SetLabelOffset( 0.02 );
      y->SetTitleOffset( 2.1 );
    }
    x->SetTitleColor( kBlack );
    if ( fRatio ) {
      x->SetTitleOffset( 3.2 );
      x->SetLabelOffset( 0.02 );
    }
    bool ss_bins = true;
    double first_bin_size = x->GetBinWidth( 1 );
    for ( unsigned short i = 2; i < x->GetNbins(); ++i ) {
      if ( x->GetBinWidth( i ) != first_bin_size ) {
        ss_bins = false;
        break;
      }
    }

    // axis titles
    TString ttle = obj->GetTitle();
    if ( ttle.Contains( "@@" ) ) {
      TObjArray* tok = ttle.Tokenize( "@@" );
      TString x_title = "", y_title = "", unit = "", form_spec = "";
      if ( tok->GetEntries()>0 ) x_title = dynamic_cast<TObjString*>( tok->At( 0 ) )->String();
      if ( tok->GetEntries()>1 ) y_title = dynamic_cast<TObjString*>( tok->At( 1 ) )->String();
      if ( tok->GetEntries()>2 ) {
        unit = dynamic_cast<TObjString*>( tok->At( 2 ) )->String();
        if ( unit.Contains( "?" ) ) { // extract format specifier
          TObjArray* tok2 = unit.Tokenize( "?" );
          if ( tok2->GetEntries() > 1 ) {
            unit = dynamic_cast<TObjString*>( tok2->At( 0 ) )->String();
            form_spec = dynamic_cast<TObjString*>( tok2->At( 1 ) )->String();
          }
          else {
            unit = "";
            form_spec = dynamic_cast<TObjString*>( tok2->At( 0 ) )->String();
          }
        }
        if ( !unit.IsNull() || !form_spec.IsNull() ) {
          if ( !unit.IsNull() ) x_title = Form( "%s (%s)", x_title.Data(), unit.Data() );
          if ( !ss_bins ) {
            y_title = Form( "%s / bin", y_title.Data() );
          }
          else {
            if ( !form_spec.IsNull() ) {
              TString format = Form( "%%s / %%%s %%s", form_spec.Data() );
              y_title = Form( format.Data(), y_title.Data(), GetBinning( obj ), unit.Data() );
            }
            else y_title = Form( "%s / %d %s", y_title.Data(), static_cast<unsigned int>( GetBinning( obj ) ), unit.Data() );
          }
        }
      }
      obj->GetXaxis()->SetTitle( x_title );
      obj->GetYaxis()->SetTitle( y_title );
      obj->SetTitle( "" );
    }
  }

  inline void Divide( int num_cols, int num_lines = 1, float xmargin = 0.01, float ymargin = 0.01, int color = 0 ) override {
    if ( fRatio ) return;
    TCanvas::Divide( num_cols, num_lines, xmargin, ymargin, color );
    double top_margin = 0.055;
    double pad_x = 1./num_cols, pad_y = ( 1.-top_margin )/num_lines;
    for ( unsigned short l = 0; l < num_lines; ++l ) {
      for ( unsigned short c = 0; c < num_cols; ++c ) { // fetch one line by one line
        TPad* p = dynamic_cast<TPad*>( TCanvas::GetPad( 1+l*num_cols+c ) );
        p->SetPad( 0.+pad_x*c, 1.-top_margin-pad_y*( l+1 ), 0.+pad_x*( c+1 ), 1.-top_margin-pad_y*l );
        //p->SetFillStyle( 4000 );
        p->SetFillStyle( 0 );
        p->SetLeftMargin( 0.15 );
        if ( c == num_cols-1 ) p->SetRightMargin( TCanvas::GetRightMargin() );
        p->SetTopMargin( 0. );
        p->SetBottomMargin( 0.15 );
        p->SetTicks( 1, 1 );
      }
    }
    divided_ = true;
  }

  typedef std::vector< std::pair<std::string,TH1*> > HistsMap;
  inline TH1* RatioPlot( HistsMap hm, float ymin = -999., float ymax = -999., const char* yaxis = "", float xline = -999. ) {
    if ( !fRatio )
      return nullptr;
    TH1* denom = hm.begin()->second,
        *numer = 0;
    denom->GetXaxis()->SetTitle( "" );
    TCanvas::cd( 2 );

    TH1D* denom_err = (TH1D*)denom->Clone(),
         *denom_err2 = (TH1D*)denom->Clone();
    denom_err2->Sumw2( false );
    denom_err->Divide( denom_err2 );

    unsigned short i = 0;
    for ( HistsMap::const_iterator it = hm.begin()+1; it != hm.end(); ++it ) {
      numer = dynamic_cast<TH1*>( it->second->Clone() );
      //ratio1->Sumw2(); ratio2->Sumw2();
      numer->Divide( denom );
      numer->Draw( ( i == 0 ) ? "p" : "p same" );
      //numer->Draw( "p same" );
      if ( ymin != ymax ) {
        numer->GetYaxis()->SetRangeUser( ymin, ymax );
      }
      Prettify( numer );
      if ( strcmp( yaxis, "" ) == 0 ) numer->GetYaxis()->SetTitle( Form( "Ratio%s", ( hm.size() > 2 ) ? "s" : "" ) );
      else numer->GetYaxis()->SetTitle( yaxis );
      i++;
    }
    denom_err->Draw( "e2same" );
    denom_err->SetFillColor( kBlack );
    denom_err->SetFillStyle( 3004 );

    if ( xline != -999. ) {
      TLine* l = new TLine( denom->GetXaxis()->GetXmin(), xline, denom->GetXaxis()->GetXmax(), xline );
      l->SetLineColor( kBlack );
      l->SetLineWidth( 1 );
      l->Draw();
    }
    Prettify( denom );
    TCanvas::cd();
    return denom;
  }

  inline void SetTopLabel( const char* lab = "" ) {
    TCanvas::cd();
    if ( strcmp( lab, "" ) != 0 ) fTitle = lab;
    if ( !fTopLabel ) BuildTopLabel();
    else fTopLabel->Clear();
    fTopLabel->AddText( fTitle );
    //fTopLabel->Draw();
  }

  inline TLegend* GetLegend() { return fLeg.get(); }
  inline void SetLegendX1( double x ) { fLegXpos = x; }
  inline void SetLegendY1( double y ) { fLegYpos = y; }
  inline void SetLegendSizeX( double x ) { fLegXsize = x; }
  inline void SetLegendSizeY( double y ) { fLegYsize = y; }
  inline void AddLegendEntry( const TObject* obj, const char* title, Option_t* option = "lpf" ) {
    if ( !fLeg ) CreateLegend();
    fLeg->AddEntry( obj, title, option );
    const unsigned int num_entries = fLeg->GetNRows();
    if ( num_entries > 3 ) {
      fLegYsize += ( num_entries-3 )*0.015;
      SetLegendY1( fLeg->GetY1NDC() );
    }
  }

  inline void Save( const char* ext, const char* out_dir = "." ) {
    TCanvas::cd();
    if ( fLeg && TCanvas::FindObject( fLeg.get() ) == 0 )
      fLeg->Draw();
    if ( fBanner && TCanvas::FindObject( fBanner.get() ) == 0 )
      fBanner->Draw();
    if ( fTopLabel && TCanvas::FindObject( fTopLabel.get() ) == 0 )
      fTopLabel->Draw();

    TCanvas::SetFillStyle( 0 );
    const TString ext_str( ext );
    TObjArray* tok = TString( ext ).Tokenize( "," );
    for ( unsigned short i = 0; i < tok->GetEntries(); ++i ) {
      const TString ext_str = dynamic_cast<TObjString*>( tok->At( i ) )->String();
      bool good_format = false;
      good_format |= ( ext_str == "pdf" );
      good_format |= ( ext_str == "png" );
      if ( !good_format ) continue;
      TCanvas::SaveAs( Form( "%s/%s.%s", out_dir, TCanvas::GetName(), ext_str.Data() ) );
    }
  }

 private:
  inline void Build() {
    TCanvas::SetLeftMargin( 0.14 );
    TCanvas::SetTopMargin( 0.06 );
    TCanvas::SetRightMargin( 0.1 );
    TCanvas::SetBottomMargin( 0.12 );
    TCanvas::SetTicks( 1, 1 );
    //TCanvas::SetFillStyle( 4000 );
    TCanvas::SetFillStyle( 0 );

    if ( align_ == Align::left ) {
      fBanner.reset( new PaveText( 0.17, 0.855, 0.185, 0.932 ) );
      fBanner->SetTextAlign( kHAlignLeft+kVAlignTop );
    }
    else if ( align_ == Align::right ) {
      fBanner.reset( new PaveText( 0.82, 0.855, 0.835, 0.932 ) );
      fBanner->SetTextAlign( kHAlignRight+kVAlignTop );
    }
    fBanner->AddText( Form( "#font[62]{%s}", EXPERIMENT ) );
    fBanner->AddText( Form( "#scale[0.75]{#font[52]{%s}}", fPlotType.Data() ) );
    fBanner->SetTextSize( 0.04 );

    SetTopLabel();
    if ( fRatio ) DivideCanvas();
  }

  inline void DivideCanvas() {
    TCanvas::Divide( 1, 2 );
    TPad* p1 = dynamic_cast<TPad*>( TCanvas::GetPad( 1 ) ),
         *p2 = dynamic_cast<TPad*>( TCanvas::GetPad( 2 ) );
    p1->SetPad( 0., 0.3, 1., 1. );
    p2->SetPad( 0., 0.0, 1., 0.3 );
    p1->SetLeftMargin( TCanvas::GetLeftMargin() );
    p1->SetRightMargin( TCanvas::GetRightMargin() );
    p2->SetLeftMargin( TCanvas::GetLeftMargin() );
    p2->SetRightMargin( TCanvas::GetRightMargin() );
    p1->SetTopMargin( TCanvas::GetTopMargin()+0.025 );
    p1->SetBottomMargin( 0.02 );
    p2->SetTopMargin( 0.02 );
    p2->SetBottomMargin( TCanvas::GetBottomMargin()+0.225 );
    p1->SetTicks( 1, 1 );
    p2->SetTicks( 1, 1 );
    p2->SetGrid( 0, 1 );
    TCanvas::cd( 1 );
  }

  inline void BuildTopLabel() {
    TCanvas::cd();
    fTopLabel.reset( new PaveText( 0.5, 0.95, 0.925, 0.96 ) );
  }

  inline void CreateLegend() {
    if ( fLeg ) return;
    if ( fRatio ) TCanvas::cd( 1 );
    fLeg.reset( new TLegend( fLegXpos, fLegYpos, fLegXpos+fLegXsize, fLegYpos+fLegYsize ) );
    fLeg->SetFillStyle( 0 );
    //fLeg->SetLineColor( kWhite );
    //fLeg->SetLineColor( kGray );
    //fLeg->SetLineWidth( 1 );
    fLeg->SetLineWidth( 0 );
    fLeg->SetTextSize( 0.035 );
    fLeg->Draw();
  }
  inline float GetBinning( const TH1* h ) {
    return ( h->GetXaxis()->GetXmax() - h->GetXaxis()->GetXmin() ) / h->GetXaxis()->GetNbins();
  }

  TString fTitle, fPlotType;
  std::unique_ptr<PaveText> fBanner, fTopLabel;
  std::unique_ptr<TLegend> fLeg;
  double fLegXpos, fLegYpos, fLegXsize, fLegYsize;
  bool fRatio;
  bool divided_;
  Align align_;
};

TH1*
WithOverflow( TH1* h )
{
  //function to paint the histogram h with an extra bin for overflows
  unsigned short nx = h->GetNbinsX()+1;
  double* xbins = new double[nx+1];
  for ( unsigned short i = 0; i < nx; ++i ) { xbins[i] = h->GetBinLowEdge( i+1 ); }
  xbins[nx] = xbins[nx-1] + h->GetBinWidth( nx );
  //book a temporary histogram having extra bins for overflows
  TH1D* htmp = new TH1D( Form( "%s_withoverflow", h->GetName() ), h->GetTitle(), nx, xbins );
  htmp->Sumw2();
  //fill the new histogram including the overflows
  for ( unsigned short i = 1; i <= nx; ++i ) {
    htmp->SetBinContent( htmp->FindBin( htmp->GetBinCenter( i ) ), h->GetBinContent( i ) );
    htmp->SetBinError( htmp->FindBin( htmp->GetBinCenter( i ) ), h->GetBinError( i ) );
  }
  htmp->SetBinContent( htmp->FindBin( h->GetBinLowEdge( 1 )-1 ), h->GetBinContent( 0 ) );
  htmp->SetBinError( htmp->FindBin( h->GetBinLowEdge( 1 )-1 ), h->GetBinError( 0 ) );
  // Restore the number of entries
  htmp->SetEntries( h->GetEffectiveEntries() );
  return htmp;
}

#endif
