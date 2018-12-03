#ifndef RecoCTPPS_ProtonProducer_XiInterpolator_h
#define RecoCTPPS_ProtonProducer_XiInterpolator_h

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"
#include "DataFormats/Provenance/interface/RunID.h"

#include "DiphotonAnalyzer/TreeProducer/interface/ProtonKinematicsUtils.h"

#include "DiphotonAnalyzer/TreeProducer/interface/RPAlignmentConstants.h"
#include "DiphotonAnalyzer/TreeProducer/interface/RPCalibrationConstants.h"

#include "TFile.h"
#include "TGraph.h"
#include "TSpline.h"

#include <iostream>

namespace ProtonUtils
{
  class XiInterpolator
  {
    public:
      XiInterpolator();
      XiInterpolator( const char* filename );
      ~XiInterpolator();

      void loadInterpolationGraphs( const char* filename, bool conversion=false );
      void setCalibrationConstants( const edm::RunNumber_t& run_id );
      void setAlignmentConstants( const CTPPSAlCa::RPAlignmentConstants& ac ) { align_ = ac; }

      void computeXiLinear( const CTPPSDetId& detid, const CTPPSLocalTrackLite& trk, float& xi, float& err_xi );
      void computeXiSpline( const CTPPSDetId& detid, const CTPPSLocalTrackLite& trk, float& xi, float& err_xi );

    private:
      void extractSpline( const TGraph* gr_in, const char* name_out, TGraph* gr_out, TSpline3* sp_out );

      TSpline3 *isLF_, *isLN_, *isRF_, *isRN_;
      CTPPSAlCa::RPAlignmentConstants align_;
      CTPPSAlCa::RPCalibrationConstants calib_;
  };
}

#endif
