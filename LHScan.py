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
rmin = -2.
rmax = 6
np = 20

R = numpy.linspace(rmin, rmax, np)
C = []
for r in R:
  mu.setVal(r)
  minim.minimize("Minuit","minimize")
  C.append(2*nll.getVal())

plt.plot(R,C)
plt.ylabel("-2 Log(L)")
plt.xlabel("$\mu$")
plt.show()
