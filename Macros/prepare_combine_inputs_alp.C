#include "TH1.h"

void prepare_combine_inputs_alp()
{
  struct value_t {
    string path;
    double coupling, mass, xsec;
  };
  vector<value_t> samples = {
    { "samples/ntuple-alp-m500_f1e-1.root",  1.e-1,  500., 2.3101621130799566E-003 },
    { "samples/ntuple-alp-m750_f1e-1.root",  1.e-1,  750., 1.9776221405498604E-003 },
    { "samples/ntuple-alp-m1000_f1e-1.root", 1.e-1, 1000., 9.5337063012939575E-004 },
    { "samples/ntuple-alp-m1250_f1e-1.root", 1.e-1, 1250., 4.6240501744627136E-004 },
    { "samples/ntuple-alp-m1500_f1e-1.root", 1.e-1, 1500., 1.9682574485681756E-004 },
    { "samples/ntuple-alp-m1750_f1e-1.root", 1.e-1, 1750., 5.8186796297957804E-005 },
    { "samples/ntuple-alp-m2000_f1e-1.root", 1.e-1, 2000., 3.7767336768070358E-007 },
    { "samples/ntuple-alp-m500_f5e-1.root",  5.e-1,  500., 5.8835732216095324E-002 },
    { "samples/ntuple-alp-m750_f5e-1.root",  5.e-1,  750., 5.0151589355963207E-002 },
    { "samples/ntuple-alp-m1000_f5e-1.root", 5.e-1, 1000., 2.3823877047682900E-002 },
    { "samples/ntuple-alp-m1250_f5e-1.root", 5.e-1, 1250., 1.0890909682316622E-002 },
    { "samples/ntuple-alp-m1500_f5e-1.root", 5.e-1, 1500., 4.8234860094796576E-003 },
    { "samples/ntuple-alp-m1750_f5e-1.root", 5.e-1, 1750., 1.5037045828166359E-003 },
    { "samples/ntuple-alp-m2000_f5e-1.root", 5.e-1, 2000., 2.1620443123727138E-004 },
    /*{ "samples/ntuple-alp-m500_f1e-1p5.root", 1.e-1.5,  500., 2.2969922909241779E-004 },
    { "samples/ntuple-alp-m750_f1e-1p5.root",  1.e-1.5,  750., 2.0008675939864182E-004 },
    { "samples/ntuple-alp-m1000_f1e-1p5.root", 1.e-1.5, 1000., 9.5952100756794203E-005 },
    { "samples/ntuple-alp-m1250_f1e-1p5.root", 1.e-1.5, 1250., 4.6541893851199770E-005 },
    { "samples/ntuple-alp-m1500_f1e-1p5.root", 1.e-1.5, 1500., 1.9678083746387744E-005 },
    { "samples/ntuple-alp-m1750_f1e-1p5.root", 1.e-1.5, 1750., 5.7993168158861949E-006 },
    { "samples/ntuple-alp-m2000_f1e-1p5.root", 1.e-1.5, 2000., 3.5524088265405244E-009 },*/
  };
}
