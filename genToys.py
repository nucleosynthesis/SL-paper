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

ntoys = 5000
nbins = 90

bV       = [array.array('d',[0]) for it in range(nbins)] 
bmle     = [array.array('d',[ws.function("expected_background_bin%d"%(it+1)).getVal()]) for it in range(nbins)] 

for b in range(nbins):
  tree.Branch("b_%d"%(b+1),bV[b],"b_%d/D"%(b+1))
  tree.Branch("mle_b_%d"%(b+1),bmle[b],"mle_b_%d/D"%(b+1))

itern = nuis.createIterator()
nuisV = [array.array('d',[0]) for it in range(nuis.getSize())] 

for ni in range(nuis.getSize()):
   n = itern.Next()
   nme = n.GetName()
   tree.Branch("%s"%(nme),nuisV[ni],"%s/D"%(nme))

genData = con.generate(nuis,ROOT.RooFit.NumEvents(ntoys))

for i in range(ntoys): 

 # set nuisances to the generated values
 itern = nuis.createIterator()
 for ni in range(nuis.getSize()):
   n = itern.Next()
   nv = genData.get(i).getRealValue(n.GetName())
   ws.var(n.GetName()).setVal(nv)
   nuisV[ni][0]=nv

 # now find the bin values in this toy
 for b in range(nbins):
   bV[b][0] = ws.function("expected_background_bin%d"%(b+1)).getVal()

 # fill the tree  
 tree.Fill()
 #end loop

fout.cd()
tree.Write()
