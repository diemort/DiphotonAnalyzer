#!/bin/sh
#JSON_FILE=/afs/cern.ch/user/l/lforthom/work/yyanalysis/CMSSW_8_1_0_pre8/src/Cert_279760-280385_13TeV_PromptReco_Collisions16_JSON_NoL1T_PPSruns.txt
#JSON_FILE=utils/processed_runs_BCG_PPSruns.json
JSON_FILE=/afs/cern.ch/user/l/lforthom/public/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_PPSruns_preTS2.txt
PILEUP_FILE=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt
#pileupCalc.py -i ${JSON_FILE} --inputLumiJSON ${PILEUP_FILE} --calcMode true --minBiasXsec 69200 --maxPileupBin 75 --numPileupBins 75 pileup_data16BCG_PPSruns.root
pileupCalc.py -i ${JSON_FILE} --inputLumiJSON ${PILEUP_FILE} --calcMode true --minBiasXsec 69200 --maxPileupBin 99 --numPileupBins 99 pileup_data16BCG_PPSruns.root
