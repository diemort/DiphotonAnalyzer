#!/usr/bin/python

from ROOT import TH1D
#from SimGeneral.MixingModule.mix_CSA14_50ns_PoissonOOTPU_cfi import mix
# QCD
#from SimGeneral.MixingModule.mix_2016_25ns_Moriond17MC_PoissonOOTPU_cfi import mix
from SimGeneral.MixingModule.mix_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU_cfi import mix

nbins = mix.input.nbPileupEvents.probFunctionVariable[-1]+1

prob = TH1D('pileup', '', nbins, 0., nbins*1.)

i = 0
for value in mix.input.nbPileupEvents.probValue:
    prob.Fill(i, value)
    i += 1

prob.SaveAs('pileup_mc.root')

