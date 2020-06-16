#include "Canvas.h"

void test_canvas()
{
  Canvas c( "test_canvas", "Testing the canvas", "Experiment", "Preliminary", false );
  TH1D h_test( "h_test", ";x;y", 100, 0., 1. );
  h_test.Draw();
  c.AddLegendEntry( &h_test, "test" );
  c.Prettify( &h_test );
  c.Save( "pdf,png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
}
