import ROOT 
import sys 

# read in a TTree and calculate moments to be spat out in the SL pythom model file - :)
# Run with -> python toys2ModelFile.py in.root > out.py 

def mean(tree,x):

  nevt = tree.GetEntries()
  mean = 0 
  for i in range(nevt):
    tree.GetEntry(i)
    mean += getattr(tree,x)
  return mean/nevt

def covariance(tree,x,mx,y,my):
  
  nevt = tree.GetEntries()
  cvar = 0 
  for i in range(nevt):
    tree.GetEntry(i)
    dx = getattr(tree,x) - mx
    dy = getattr(tree,y) - my
    cvar+=dx*dy
  return cvar/(nevt-1)

def skew(tree,x,mx):
  
  nevt = tree.GetEntries()
  skew = 0 
  for i in range(nevt):
    tree.GetEntry(i)
    dx = getattr(tree,x) - mx
    skew+=dx*dx*dx

  return skew/(nevt)


fi = ROOT.TFile.Open(sys.argv[1])
tree = fi.Get("toys")
nbins = 90

means      = [ mean(tree,"b_%d"%i) for i in range(nbins) ]
covariance = [ [covariance(tree,"b_%d"%i,means[i],"b_%d"%j,means[j]) for j in range(nbins)] for i in range(nbins) ]
skews      = [ skew(tree,"b_%d"%i,means[i]) for i in range(nbins) ]
covariance = [ covariance[i][j] for i in range(nbins) for j in range(nbins) ]

fother = ROOT.TFile.Open("histos.root")
dataH = fother.Get("data")
signalH = fother.Get("signal")

data   = [dataH.GetBinContent(b+1) for b in range(nbins)]
signal = [signalH.GetBinContent(b+1) for b in range(nbins)]

import numpy
import array
import numpy
import array
print "import numpy"
print "import array"
print "name = 'Generated Model' "
print "nbins = %d"%nbins
print "data = array.array('d',",data,")"
print "background = array.array('d',",means,")"
print "covariance   = array.array('d',",covariance,")"
print "third_moment = array.array('d',",skews,")"
print "signal = array.array('d',",signal,")"

