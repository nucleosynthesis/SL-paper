import ROOT
import array,numpy,sys
from matplotlib import pyplot as plt

fi = ROOT.TFile.Open("workspace.root")
ws = fi.Get("w")

mu   = ws.var("mu"); mu.setConstant()
data = ws.data("obsdata")
pdf  = ws.pdf("combined_pdf")
con  = ws.pdf("nuisance_pdf")

data.Print()
con.Print()
pdf.Print()

nll = pdf.createNLL(data,ROOT.RooFit.ExternalConstraints(ROOT.RooArgSet(con)))
minim = ROOT.RooMinimizer(nll)
# do a scan
rmin = -1.
rmax = 2
np = 20

R = numpy.linspace(rmin, rmax, np)
C = []
for r in R:
  mu.setVal(r)
  minim.minimize("Minuit","minimize")
  C.append(2*nll.getVal())

# also make the usual tree 
fout = ROOT.TFile("full-LH.root","RECREATE")
tree = ROOT.TTree("limit","limit")

dnll = array.array('d',[0])
r    = array.array('d',[0])

tree.Branch("r",r,"r/D")
tree.Branch("deltaNLL",dnll,"deltaNLL/D")

nll0 = min(C)
for i in range(len(R)) :
  dnll[0]= (C[i]-nll0)/2
  r[0]=R[i]
  tree.Fill()

fout.cd()
tree.Write()
  

plt.plot(R,C)
plt.ylabel("-2 Log(L)")
plt.xlabel("$\mu$")
plt.show()
