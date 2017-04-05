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

#include <string.h>
#include <iostream> //FIXME

class Canvas : public TCanvas
{
  public:
    typedef std::vector< std::pair<const char*,TH1*> > HistsMap;

  public:
    // Default constructor
    inline Canvas() :
      TCanvas( "invalid", "", 600, 600 ),
      fTitle( "" ), fTopLabel( 0 ),
      fLeg( 0 ), fLegXSize( 0.35 ), fLegYSize( 0.15 ),
      fRatio( false )
    {}

    inline Canvas( const char* name, const char* title="", bool ratio=false ) :
      //TCanvas(name, "", 450, 450),
      TCanvas( name, "", 600, 600 ),
      fTitle( title ), fTopLabel( 0 ),
      fLeg( 0 ), fLegXSize( 0.35 ), fLegYSize( 0.15 ),
      fRatio( ratio )
    {
      Build();
    }
    /// Copy constructor
    inline Canvas( const Canvas& c, const char* name="" ) :
      //TCanvas( name, c.fTitle, c.GetWw(), c.GetWh() ),
      TCanvas( ( strcmp( name, "" )!=0 ) ? name : c.GetName(), c.fTitle, c.GetWw(), c.GetWh() ),
      fTitle( c.fTitle ), fTopLabel( c.fTopLabel ),
      fLeg( c.fLeg ), fLegXSize( c.fLegXSize ), fLegYSize( c.fLegYSize ),
      fRatio( c.fRatio )
    {
      Build();
    }

    inline ~Canvas() {
cout << __PRETTY_FUNCTION__ << " for " << GetName() << endl;
      if ( fLeg ) delete fLeg;
      if ( fTopLabel ) delete fTopLabel;
    }

    inline Canvas* Clone( const char* name ) const { return new Canvas( *this, name ); }

    void SetRatioPlot()
    {
      if ( fRatio ) return;
      fRatio = true;
      DivideCanvas();
    }
    inline void SetTitle( const char* title ) { fTitle = title; }

    inline void Prettify( TH1& obj ) { Prettify( &obj ); }
    inline void Prettify( TH1* obj )
    {
      TAxis* x = dynamic_cast<TAxis*>( obj->GetXaxis() ),
            *y = dynamic_cast<TAxis*>( obj->GetYaxis() );
      x->SetLabelFont( 43 ); x->SetLabelSize( 20 );
      x->SetTitleFont( 43 ); x->SetTitleSize( 26 );
      y->SetLabelFont( 43 ); y->SetLabelSize( 20 );
      y->SetTitleFont( 43 ); y->SetTitleSize( 26 );
      x->SetTitleColor( kBlack );
      if ( fRatio ) {
        x->SetTitleOffset( 3.2 );
        x->SetLabelOffset( 0.02 );
      }
      y->SetTitleOffset( 1.4 );

      // axis titles
      TString ttle = obj->GetTitle();
      if ( ttle.Contains( "\\" ) ) {
        TObjArray* tok = ttle.Tokenize( "\\" );
        TString x_title = "", y_title = "", unit = "", form_spec = "";
        if ( tok->GetEntries()>0 ) x_title = dynamic_cast<TObjString*>( tok->At( 0 ) )->String();
        if ( tok->GetEntries()>1 ) y_title = dynamic_cast<TObjString*>( tok->At( 1 ) )->String();
        if ( tok->GetEntries()>2 ) {
          unit = dynamic_cast<TObjString*>( tok->At( 2 ) )->String();
          if ( unit.Contains( "?" ) ) { // extract format specifier
            TObjArray* tok2 = unit.Tokenize( "?" );
            if ( tok2->GetEntries()>1 ) {
              unit = dynamic_cast<TObjString*>( tok2->At( 0 ) )->String();
              form_spec = dynamic_cast<TObjString*>( tok2->At( 1 ) )->String();
            }
            else {
              unit = "";
              form_spec = dynamic_cast<TObjString*>( tok2->At( 0 ) )->String();
            }
          }
          if ( !unit.IsNull() or !form_spec.IsNull() ) {
            if ( !unit.IsNull() ) x_title = Form( "%s (%s)", x_title.Data(), unit.Data() );
            if ( !form_spec.IsNull() ) {
              TString format = Form( "%%s / %%%s %%s", form_spec.Data() );
              y_title = Form( format.Data(), y_title.Data(), GetBinning( obj ), unit.Data() );
            }
            else y_title = Form( "%s / %d %s", y_title.Data(), static_cast<unsigned int>( GetBinning( obj ) ), unit.Data() );
          }
        }
        obj->GetXaxis()->SetTitle( x_title );
        obj->GetYaxis()->SetTitle( y_title );
        obj->SetTitle( "" );
      }
    }

    inline void RatioPlot( HistsMap hm, float ymin=-999., float ymax=-999., float xline=-999. )
    {
      if ( !fRatio ) return;
      TH1* denom = hm.begin()->second,
          *numer = 0;
      denom->GetXaxis()->SetTitle( "" );
      TCanvas::cd( 2 );
      unsigned short i = 0;
      for ( HistsMap::const_iterator it=hm.begin()+1; it!=hm.end(); it++ ) {
        numer = dynamic_cast<TH1*>( it->second->Clone() );
        //ratio1->Sumw2(); ratio2->Sumw2();
        numer->Divide( denom );
        numer->Draw( ( i==0 ) ? "p" : "p same" );
        //numer->Draw( "p same" );
        if ( ymin!=ymax ) {
          numer->GetYaxis()->SetRangeUser( ymin, ymax );
        }
        Prettify( numer );
        numer->GetYaxis()->SetTitle( Form( "Ratio%s", ( hm.size()>2 ) ? "s" : "" ) );
        i++;
      }
      if ( xline!=-999. ) {
        TLine* l = new TLine( denom->GetXaxis()->GetXmin(), xline, denom->GetXaxis()->GetXmax(), xline );
        l->SetLineColor( kRed );
        l->SetLineWidth( 1 );
        l->Draw();
      }
      Prettify( denom );
      TCanvas::cd();
    }

    inline void SetTopLabel( const char* lab="" )
    {
      TCanvas::cd();
      if ( strcmp( lab, "" )!=0 ) fTitle = lab;
      if ( !fTopLabel ) BuildTopLabel();
      else fTopLabel->Clear();
      fTopLabel->AddText( fTitle );
      //fTopLabel->Draw();
    }

    inline TLegend* GetLegend() { return fLeg; }
    inline void SetLegendX1( double x )
    {
      if ( !fLeg ) CreateLegend();
      fLeg->SetX1NDC( x );
      fLeg->SetX2NDC( x+fLegXSize );
    }
    inline void SetLegendY1( double y )
    {
      if ( !fLeg ) CreateLegend();
      fLeg->SetY1NDC( y );
      fLeg->SetY2NDC( y+fLegYSize );
    }
    inline void AddLegendEntry( const TObject* obj, const char* title, Option_t* option="lpf" )
    {
      if ( !fLeg ) CreateLegend();
std::cout << __PRETTY_FUNCTION__ << "::" << fLeg << std::endl;
      fLeg->AddEntry( obj, title, option );
      if ( fLeg->GetNRows()>4 ) {
        fLegYSize += 0.08;
        SetLegendY1( fLeg->GetY1NDC() );
      }
    }

    inline void Save( const char* ext, const char* out_dir="." )
    {
      TCanvas::cd();
      if ( fLeg and TCanvas::FindObject( fLeg )==0 ) {
        fLeg->Draw();
      }
      if ( fTopLabel and TCanvas::FindObject(fTopLabel)==0 ) {
        fTopLabel->Draw();
      }

      const TString ext_str( ext );
      TObjArray* tok = TString( ext ).Tokenize( "," );
      for ( unsigned short i=0; i<tok->GetEntries(); i++ ) {
        const TString ext_str = dynamic_cast<TObjString*>( tok->At( i ) )->String();
        if ( ext_str!="pdf" && ext_str!="png" ) continue;
        TCanvas::SaveAs( Form( "%s/%s.%s", out_dir, TCanvas::GetName(), ext_str.Data() ) );
      }
    }

  private:
    inline void Build()
    {
      TCanvas::SetLeftMargin( 0.14 );
      TCanvas::SetTopMargin( 0.06 );
      TCanvas::SetRightMargin( 0.1 );
      TCanvas::SetBottomMargin( 0.12 );
      TCanvas::SetTicks( 1, 1 );

      SetTopLabel();
      if ( fRatio ) DivideCanvas();
    }

    inline void DivideCanvas()
    {
      TCanvas::Divide(1,2);
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

    inline void BuildTopLabel()
    {
      TCanvas::cd();
      fTopLabel = new TPaveText( 0.5, 0.95, 0.925, 0.96, "NB NDC" );
      fTopLabel->SetFillStyle( 0 );
      fTopLabel->SetFillColor( 0 );
      fTopLabel->SetLineColor( 0 );
      fTopLabel->SetLineStyle( 0 );
      fTopLabel->SetTextFont( 42 );
      fTopLabel->SetTextSize( 0.033 );
      fTopLabel->SetTextAlign( kHAlignRight+kVAlignBottom );
    }

    inline void CreateLegend()
    {
      if ( fLeg ) return;
      if ( fRatio ) TCanvas::cd(1);
      fLeg = new TLegend( 0.5, 0.75, 0.5+fLegXSize, 0.75+fLegYSize );
      //fLeg->SetFillStyle( 0 );
      //fLeg->SetLineColor( kWhite );
      fLeg->SetLineColor( kGray );
      fLeg->SetLineWidth( 1 );
      fLeg->SetTextSize( 0.032 );
    }
    inline float GetBinning( const TH1* h ) { return ( h->GetXaxis()->GetXmax() - h->GetXaxis()->GetXmin() ) / h->GetXaxis()->GetNbins(); }

    TString fTitle;
    TPaveText* fTopLabel;
    TLegend* fLeg;
    double fLegXSize, fLegYSize;
    bool fRatio;
};

TH1*
WithOverflow( TH1* h )
{
  //function to paint the histogram h with an extra bin for overflows
  unsigned short nx = h->GetNbinsX()+1;
  double* xbins = new double[nx+1];
  for ( unsigned short i=0; i<nx; i++ ) { xbins[i] = h->GetBinLowEdge( i+1 ); }
  xbins[nx] = xbins[nx-1] + h->GetBinWidth( nx );
  //book a temporary histogram having extra bins for overflows
  TH1D* htmp = new TH1D( Form( "%s_withoverflow", h->GetName() ), h->GetTitle(), nx, xbins );
  htmp->Sumw2();
  //fill the new histogram including the overflows
  for ( unsigned short i=1; i<=nx; i++ ) {
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
