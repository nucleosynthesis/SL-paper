import ROOT 
import sys 
import numpy as np
from scipy import stats
import array
ROOT.gROOT.SetBatch(1)

# read in a TTree and calculate moments to be spat out in the SL pythom model file - :)
# Run with -> python toys2ModelFile.py in.root > out.py 
makePlots = False

def getCoefficients(m1,m2,m3):

  inside = np.complex(m3*m3 - 8*m2*m2*m2,0)
  root = inside**0.5;
  c_m3 = np.complex(-m3,0);
  k3   = c_m3+root;
  k    = k3**(1./3);

  j2 = np.complex(1,(3)**0.5);
  j  = np.complex(1,-(3)**0.5);
      
  c = -j*m2/k-0.5*j2*k;

  C = np.real(c);
  if (m2 < (C*C)/2): 
  	B = m2**0.5
	C = 0
  else:
        B = (m2 - (C*C)/2)**0.5;
  A = m1 - C/2;

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
  hg = h.Clone(); hg.SetName("hg_%d"%(b)); hg.SetLineColor(2)
  hq = h.Clone(); hq.SetName("hq_%d"%(b)); hq.SetLineColor(ROOT.kGreen+2)
  for ib in range(h.GetNbinsX()):
    hg.SetBinContent(ib,0)
    hq.SetBinContent(ib,0)
    hg.SetBinError(ib,0)
    hq.SetBinError(ib,0)
  r = ROOT.TRandom3()
  
  A,B,C = getCoefficients(mean,var,skew)

  for i in range(10000):
    rx = r.Gaus(0,1)
    hg.Fill(mean*(1+rx*(var**0.5)/mean))
    hq.Fill(A+B*rx+(C/2)*rx*rx)

  h.SetTitle("")
  h.GetXaxis().SetTitle("Bin - %d"%(b))
  h.GetYaxis().SetTitleOffset(1.4)
  h.GetYaxis().SetTitle("Arbitrary units")
  h.Scale(1./h.Integral())
  hq.Scale(1./hq.Integral())
  hg.Scale(1./hg.Integral())
  c = ROOT.TCanvas("x%d"%b,"x",900,800)
  h.SetMaximum(h.GetMaximum()*1.2)
  h.Draw()
  hg.Draw("samehist")
  hq.Draw("samehist")

  tlat = ROOT.TLatex()
  tlat.SetTextSize(0.032)
  tlat.SetNDC()
  tlat.SetTextFont(42)
  tlat.DrawLatex(0.1,0.94,"#mu_{1}=%g, #mu_{2}=%g, #mu_{3}=%g"%(mean,var,skew))

  lmle = ROOT.TLine(mle,h.GetMinimum(),mle,h.GetMaximum())
  lmle.SetLineStyle(2)
  lmle.SetLineWidth(2)
  lmle.Draw()
  
  lgle = ROOT.TLine(mean,h.GetMinimum(),mean,h.GetMaximum())
  lgle.SetLineStyle(2)
  lgle.SetLineWidth(2)
  lgle.SetLineColor(2)
  lgle.Draw()
  
  lple = ROOT.TLine(A,h.GetMinimum(),A,h.GetMaximum())
  lple.SetLineStyle(2)
  lple.SetLineWidth(2)
  lple.SetLineColor(ROOT.kGreen+2)
  lple.Draw()
  
  pad = ROOT.TPad("p1%d"%b,"p1",0.72,0.72,0.99,0.99)
  pad.Draw()
  pad.cd()
  # draw the polynomial
  f1 = ROOT.TF1("myf_%d"%b,"%g+%g*x+%g*x*x"%(A,B,C/2),-5,5)
  f1.GetXaxis().SetNdivisions(511)
  f1.GetXaxis().SetLabelSize(0.06)
  f1.GetYaxis().SetNdivisions(511)
  f1.GetYaxis().SetLabelSize(0.06)
  f1.Draw("L")

  f1.GetXaxis().SetTitleSize(0.06)
  f1.GetXaxis().SetTitle("x = #theta_{%d}"%(b))
  c.cd()
  c.SaveAs("distribution_%d.png"%(b))

if len(sys.argv) != 3:
    print "Usage python toys2ModelFileCombine.py inputFile.root outputFile.py"
    exit()
fi = ROOT.TFile.Open(sys.argv[1])
fo = sys.argv[2]
tree = fi.Get("tree_prefit")
nbins = 90
allData = []
for iEntry in range(tree.GetEntries()):
    tree.GetEntry(iEntry)
    allDataTemp = []
    for ibin in range(nbins):
        icat = int(ibin/30) + 1
        iBinInCat = (ibin % 30) + 1
        branch = "n_exp_final_binch{0}_proc_background_{1}".format(icat,iBinInCat)
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

