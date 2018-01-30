import ROOT 
import sys 
import numpy as np
import scipy.stats as stats
from scipy import mean as MEAN
import array
ROOT.gROOT.SetBatch(1)
import math

# read in a TTree and calculate moments to be spat out in the SL pythom model file - :)
# Run with -> python toys2ModelFile.py in.root > out.py 
makePlots = False

def getCoefficients(m1,m2,m3):
  pi  = math.pi;

  if (8*m2*m2*m2 >= m3*m3) :C = -2*((2*m2)**0.5)*math.cos(4*pi/3. + (1./3.)*math.atan(((8*m2*m2*m2-m3*m3)/(m3*m3))**0.5) );
  else: C = -2*((2*m2)**0.5)*math.cosh((-1./3)*math.atanh(((-8*m2*m2*m2+m3*m3)/(m3*m3))**0.5)) ;
  if (m2 < (C)*(C)/2.) : print "Oh No! ?"#*B = TMath::Sqrt(m2-(*C)*(*C)/2.);
  B = (m2-C*C/2)**0.5
  A = m1-C/2

  return A,B,C

def mean(tree,x):

  nevt = tree.GetEntries()
  mean = 0
  ncount = 0 
  for i in range(nevt):
    tree.GetEntry(i)
    y = getattr(tree,x)
    #if y<=0: continue 
    mean += y
    ncount+=1
  return mean/ncount

def covariance(tree,x,mx,y,my):
  
  nevt = tree.GetEntries()
  cvar = 0 
  ncount = 0 
  for i in range(nevt):
    tree.GetEntry(i)
    xv = getattr(tree,x)
    #if xv<=0: continue 
    yv = getattr(tree,y)
    #if yv<=0: continue 

    dx = xv - mx
    dy = yv - my
    cvar+=dx*dy
    ncount+=1
  return cvar/(ncount-1)

def skew(tree,x,mx):
  
  nevt = tree.GetEntries()
  skew = 0 
  ncount = 0 
  for i in range(nevt):
    tree.GetEntry(i)
    y = getattr(tree,x)
    #if y<=0: continue 
    dx = y - mx
    skew+=dx*dx*dx
    ncount +=1
  return skew/(ncount)

def plotCompare(tree,b,mean,var,skew):
  
  tree.GetEntry(0)
  mle = getattr(tree,"mle_b_%d"%(b))

  ROOT.gStyle.SetOptStat(0)
  tree.Draw("b_%d>>h_%d"%(b,b));
  h = ROOT.gROOT.FindObject("h_%d"%b);  h.SetLineColor(1)
  h.Rebin(4)
  hg = h.Clone(); hg.SetName("hg_%d"%(b)); hg.SetLineColor(2)
  hq = h.Clone(); hq.SetName("hq_%d"%(b)); hq.SetLineColor(ROOT.kGreen+2)
  hg.SetLineWidth(2)
  hq.SetLineWidth(2)
  for ib in range(h.GetNbinsX()):
    hg.SetBinContent(ib,0)
    hq.SetBinContent(ib,0)
    hg.SetBinError(ib,0)
    hq.SetBinError(ib,0)
  r = ROOT.TRandom3()
  
  A,B,C = getCoefficients(mean,var,skew)

  mynewcalc = []
  for i in range(100000):
    rx = r.Gaus(0,1)
    hg.Fill(mean*(1+rx*(var**0.5)/mean))
    mvq = A+B*rx+(C/2)*rx*rx
    mynewcalc.append(mvq)
    hq.Fill(mvq)

  #mean2 = MEAN(mynewcalc)
  #var2  = stats.moment(mynewcalc,moment=2)
  #skew2  = stats.moment(mynewcalc,moment=3)
  #print "Bin = ",b, "means = ", mean,mean2, "var = ", var, var2, "skew = ", skew,skew2

  h.SetTitle("")
  h.GetXaxis().SetTitle("#hat{#it{n}}_{%d}"%(b))
  h.GetYaxis().SetTitleOffset(1.4)
  h.GetYaxis().SetTitle("Arbitrary units")
  h.GetXaxis().SetTitleSize(0.05)
  h.GetXaxis().SetTitleOffset(0.82)
  h.GetYaxis().SetTitleOffset(1.6)
  h.Scale(1./h.Integral())
  hq.Scale(1./hq.Integral())
  hg.Scale(1./hg.Integral())
  c = ROOT.TCanvas("x%d"%b,"x",900,800)
  c.SetLeftMargin(0.12)
  c.SetRightMargin(0.05)
  h.SetMaximum(h.GetMaximum()*1.2)
  h.SetMinimum(0.0)
  h.SetLineWidth(2)
  #for bi in range(h.GetNbinsX()): h.SetBinError(bi+1,0);
  h.SetMarkerStyle(20)
  #h.Draw("Psame")
  h.Draw("P")
  hg.Draw("samehist")
  hq.Draw("samehist")

  tlat = ROOT.TLatex()
  tlat.SetTextSize(0.032)
  tlat.SetNDC()
  tlat.SetTextFont(42)
  tlat.DrawLatex(0.12,0.92,"#it{m}_{1}=%.2f, #it{m}_{2}=%.2f, #it{m}_{3}=%.2f"%(mean,var,skew))

  lmle = ROOT.TLine(mle,h.GetMinimum(),mle,h.GetMaximum())
  lmle.SetLineStyle(2)
  lmle.SetLineWidth(2)
  #lmle.Draw()
  
  lgle = ROOT.TLine(mean,h.GetMinimum(),mean,h.GetMaximum())
  lgle.SetLineStyle(2)
  lgle.SetLineWidth(2)
  lgle.SetLineColor(2)
  #lgle.Draw()
  
  lple = ROOT.TLine(A,h.GetMinimum(),A,h.GetMaximum())
  lple.SetLineStyle(2)
  lple.SetLineWidth(2)
  lple.SetLineColor(ROOT.kGreen+2)
  #lple.Draw()
  
  pad = ROOT.TPad("p1%d"%b,"p1",0.59,0.44,0.99,0.99)
  pad.SetRightMargin(0.05)
  pad.SetTopMargin(0.16)
  pad.SetLeftMargin(0.2)
  pad.SetBottomMargin(0.2)
  pad.Draw()
  pad.cd()
  # draw the polynomial
  f1 = ROOT.TF1("myf_%d"%b,"%g+%g*x+%g*x*x"%(A,B,C/2),-5,5)
  f1.GetXaxis().SetNdivisions(511)
  f1.GetXaxis().SetLabelSize(0.07)
  f1.GetYaxis().SetNdivisions(511)
  f1.GetYaxis().SetLabelSize(0.06)
  f1.GetYaxis().SetTitleSize(0.08)
  f1.GetXaxis().SetTitleSize(0.08)
  f1.GetYaxis().SetTitleOffset(1.3)
  f1.SetLineColor(ROOT.kGreen+2)

  f1.SetTitle("");
  f1.GetXaxis().SetTitle("#it{#theta}")
  f1.GetYaxis().SetTitle("#it{n}(#it{#theta})")
  
  f1.Draw("L")
  f2 = ROOT.TF1("myf2_%d"%b,"%g+%g*x"%(mean,var**0.5),-5,5)
  f2.SetLineColor(2)
  
  tlat.SetTextSize(0.055)
  tlat.SetTextColor(2); tlat.DrawLatex(0.2,0.93,"#it{n}(#it{#theta}) = %.2f + %.2f#it{#theta}"%(mean,var**0.5))
  tlat.SetTextColor(ROOT.kGreen+2); tlat.DrawLatex(0.2,0.88,"#it{n}(#it{#theta}) = %.2f + %.2f#it{#theta} + %.2f#it{#theta}^{2}"%(A,B,C/2))
  #f1.GetYaxis().SetTitle("#it{n}(#it{#theta}) = %.2f + %.2f#it{#theta} + %.2f#it{#theta}^{2}"%(A,B,C/2))

  f2.Draw("Lsame")

  c.cd()
  c.SaveAs("distribution_%d.png"%(b))
  c.SaveAs("distribution_%d.pdf"%(b))

if len(sys.argv) != 3:
    print "Usage python toys2ModelFile.py inputFile.root outputFile.py"
    exit()

fi = ROOT.TFile.Open(sys.argv[1])
tree = fi.Get("toys")
fo = sys.argv[2]
nbins = 90
allData = []
for iEntry in range(tree.GetEntries()):
    tree.GetEntry(iEntry)
    allDataTemp = []
    for branch in ["b_%d"%(i+1) for i in range(nbins)]:
        allDataTemp.append(getattr(tree,branch))
    allData.append(allDataTemp)
allData = np.array(allData)
allData = allData
means = allData.mean(0)
cov = np.cov(allData.T)
allDataMinusMean = allData - means
allDataMinusMean3 = allDataMinusMean**3
skews = allDataMinusMean3.mean(0)

means        = list(means)
covariance2D = [ [cov[i,j] for j in range(nbins)] for i in range(nbins) ]
skews        = list(skews)
covariance   = [ covariance2D[i][j] for i in range(nbins) for j in range(nbins) ]

fother = ROOT.TFile.Open("histos.root")
dataH = fother.Get("data")
signalH = fother.Get("signal")

data   = [dataH.GetBinContent(b+1) for b in range(nbins)]
signal = [signalH.GetBinContent(b+1) for b in range(nbins)]

# now, for each bin, lets make a plot comparing the toys, a gaussian and the quadratic
if makePlots:
    for b in range(nbins): plotCompare(tree,b+1,means[b],covariance2D[b][b],skews[b]) 

outString = []
outString.append( "import numpy as np")
outString.append( "import array")
outString.append( "name = 'Generated Model' ")
outString.append( "nbins = %d"%nbins)
outString.append( "data = array.array('d',\t"+str(data)+"\t)")
outString.append( "background = array.array('d',\t"+str(means)+"\t)")
outString.append( "covariance   = array.array('d',\t"+str(covariance)+"\t)")
outString.append( "third_moment = array.array('d',\t"+str(skews)+"\t)")
outString.append( "signal = array.array('d',\t"+str(signal)+"\t)")

outString = "\n".join(outString)
with open(fo,'w') as f:
    f.write(outString)

