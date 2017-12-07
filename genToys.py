import ROOT
import array

fi = ROOT.TFile.Open("workspace.root")
ws = fi.Get("w")

mu   = ws.var("mu"); mu.setVal(0); mu.setConstant()
#data = ws.data("obsdata")
#pdf  = ws.pdf("combined_pdf")
con  = ws.pdf("nuisance_pdf")
nuis = ws.genobj("nuisances")

fout = ROOT.TFile("toy_trees.root","RECREATE");
tree = ROOT.TTree("toys","toys");

ntoys = 1000
nbins = 90

bV    = [array.array('d',[0]) for it in range(nbins)] 

for b in range(nbins):
  tree.Branch("b_%d"%b,bV[b],"b_%d/D"%b)

for i in range(ntoys): 
 
 genData = con.generate(nuis,ROOT.RooFit.NumEvents(1))
 itern = nuis.createIterator()
 for ni in range(nuis.getSize()):
   n = itern.Next()
   ws.var(n.GetName()).setVal(genData.get(0).getRealValue(n.GetName()))
   #iprint "Setting ", n.GetName(), ws.var(n.GetName()).getVal()
   for b in range(nbins):
     bV[b][0] = ws.function("expected_background_bin%d"%b).getVal()
   
 tree.Fill()
fout.cd()
tree.Write()
