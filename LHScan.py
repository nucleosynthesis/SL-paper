import ROOT
import array,numpy,sys
from matplotlib import pyplot as plt

noSys = False

fi = ROOT.TFile.Open("workspace.root")
ws = fi.Get("w")

nbins = 90

if noSys: 
  allN = ws.genobj("nuisances")
  iterN = allN.createIterator()
  nNuis = allN.getSize()
  for i in range(nNuis):
    obj = iterN.Next()
    ws.var(obj.GetName()).setConstant(True)

mu   = ws.var("mu"); mu.setConstant()
data = ws.data("obsdata")
pdf  = ws.pdf("combined_pdf")
con  = ws.pdf("nuisance_pdf")

data.Print()
con.Print()
pdf.Print()

nll = pdf.createNLL(data,ROOT.RooFit.ExternalConstraints(ROOT.RooArgSet(con)))
mu.setConstant(False)
minimF = ROOT.RooMinimizer(nll)
minimF.minimize("Minuit","minimize")
nll0 = nll.getVal()

mu.setConstant(True)
minim = ROOT.RooMinimizer(nll)
# do a scan
rmin = -1.
rmax = 2
np = 20

R = numpy.linspace(rmin, rmax, np)
C = []
B = []
for r in R:
  mu.setVal(r)
  minim.minimize("Minuit","minimize")
  C.append(nll.getVal()-nll0)
  thisB = []
  for b in range(nbins): thisB.append(ws.function("expected_background_bin%d"%(b+1)).getVal())
  B.append(thisB)

# also make the usual tree 
fout = ROOT.TFile("full-LH.root","RECREATE")
tree = ROOT.TTree("limit","limit")

dnll = array.array('d',[0])
r    = array.array('d',[0])
BS = []
for i in range(nbins): BS.append(array.array('d',[0]))

tree.Branch("r",r,"r/D")
tree.Branch("deltaNLL",dnll,"deltaNLL/D")
for i in range(nbins) : tree.Branch("bkg_bin_%d"%(i+1),BS[i],"bkg_bin_%d/D"%(i+1))

#nll0 = min(C)
for i in range(len(R)) :
  dnll[0]= C[i]
  r[0]=R[i]
  for b in range(nbins): BS[b][0]=B[i][b]
  tree.Fill()

fout.cd()
tree.Write()
  

plt.plot(R,C)
plt.ylabel("-Log(L)")
plt.xlabel("$\mu$")
plt.show()
