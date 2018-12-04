import ROOT
from ROOT import *
from array import array

ROOT.gROOT.SetBatch(1)
'''
gr = ROOT.TGraph2D('scan_alps.dat')
gr.SetTitle('Generator-level #sigma(m_{A},f^{-1}) (pb), #sqrt{s} = 13 TeV, p_{T}^{#gamma} > 75 GeV;log_{10}(m_{a}/GeV);log_{10}(f^{-1}/GeV^{-1})')
gr.SetMaximum(10.)
'''
gr = ROOT.TGraph2D('scan_aqgc.dat')
gr.SetTitle('Generator-level #sigma(#zeta_{1},#zeta_{2}) (pb), #sqrt{s} = 13 TeV, p_{T}^{#gamma} > 75 GeV;#zeta_{1};#zeta_{2}')

#for l in open('scan.dat'):
#    z1, z2, xsec = [float(a) for a in l.split()]
#    print xsec
#    gr.SetPoint(gr.GetN(), z1, z2, xsec)
#for i in range(gr.GetN()):
#    print gr.GetX()[i], gr.GetY()[i], gr.GetZ()[i]

ROOT.gStyle.SetPalette(ROOT.kLightTemperature)
#ROOT.gStyle.SetPalette(ROOT.kGreyScale)
#ROOT.gStyle.SetPalette(ROOT.kBeach)

#c = ROOT.TCanvas('', '', 180, 180)
c = ROOT.TCanvas()
c.SetLeftMargin(0.135)
c.SetRightMargin(0.125)
c.SetTopMargin(0.125)
c.SetBottomMargin(0.15)
c.SetLogz()
#gr.Draw('arr')
gr.Draw('colz')
cont = array('d', [2.96e-3])
h2 = gr.GetHistogram().Clone()
h2.SetContour(len(cont), cont)
h2.SetLineWidth(4)
h2.SetLineColor(1)
h2.SetLineStyle(7)
h2.Draw('same cont2 list')
#ROOT.gPad.Update()
#contours = ROOT.TObjArray(ROOT.gROOT.GetListOfSpecials().FindObject('contours'))
contours = ROOT.gROOT.GetListOfSpecials().FindObject('contours')
#cnt = contours.At(0).First()
#cnt.SetLineStyle(7)
#cnt.Draw('same')
#gr.Draw('surf3')
#gr.Draw('pcol')
#gr.GetXaxis().SetRangeUser(-1.e-10, 1.e-10)
#gr.GetYaxis().SetRangeUser(-1.e-10, 1.e-10)
gr.GetXaxis().SetTitleSize(0.065)
gr.GetYaxis().SetTitleSize(0.065)
c.SaveAs('scan.png')
c.SaveAs('scan.pdf')

