import ROOT
ROOT.gROOT.SetBatch(1)
# open the file which contains all the info 
fi = ROOT.TFile.Open("histos.root")

"""
TFile**		histos.root
 TFile*		histos.root
  KEY: TH1F	data;1
  KEY: TH1F	nominal_bkg;1
  KEY: TH1F	signal;1
  KEY: TH1F	jes_dn_bkg;1
  KEY: TH1F	jes_up_bkg;1
  KEY: TH1F	isr_dn_bkg;1
  KEY: TH1F	isr_up_bkg;1
  KEY: TH1F	mc_unc;1
"""

def integrate(h,x,y):
 I = 0
 for i in range(x,y+1):
   I+=h.GetBinContent(i+1)
 return I

mu = ROOT.RooRealVar("mu","#mu",0,-1,5)
ONE = ROOT.RooConstVar("One","",1)
allNuisancePdfs = []
allNuisanceParameters = ROOT.RooArgSet();

data = fi.Get("data")
nbins = data.GetNbinsX()
# next we play the usual game of making RooDataHist for *all* the bins. 
sampleType  = ROOT.RooCategory("bin_number","Bin Number");
observation = ROOT.RooRealVar("observed","Observed Events bin",1);

for b in range(nbins):
  sampleType.defineType("%d"%b,b)
  sampleType.setIndex(b)

obsargset = ROOT.RooArgSet(observation,sampleType)
obsdata = ROOT.RooDataSet("obsdata","Data in all Bins",obsargset)

# Interpolation !
def interpString(B,C): 

 return  "TMath::Max(0, (TMath::Abs(@0)<=1)*(1+%g*@0+%g*@0*@0) + (@0<-1)*(%g+%g*@0) + (@0>1)*(%g+%g*@0)  )"%(B,C,1-C,B-2*C,1-C,B+2*C)

# Fill the dataset
for b in range(1,nbins+1):
  sampleType.setIndex(b-1)
  observation.setVal(data.GetBinContent(b));
  obsdata.add(ROOT.RooArgSet(observation,sampleType));

signal = fi.Get("signal")
bkg = fi.Get("nominal_bkg")
# Ok, that was easy, now we want to make an "expectation" for each of the bins 0->89

# 1, for the JES and ISR systematics, lets make the template variations 
# JES SYSTEMATIC ####################################################################
hjesu = fi.Get("jes_up_bkg")
hjesd = fi.Get("jes_dn_bkg")

nuis_JES    = ROOT.RooRealVar("nuis_JES","nuis_JES",0,-5,5)
nuis_JES_In = ROOT.RooRealVar("nuis_JES_In","nuis_JES_In",0,-5,5); nuis_JES_In.setConstant()
nuis_JES_pdf = ROOT.RooGaussian("nuis_JES_PDF","",nuis_JES,nuis_JES_In,ONE)
allNuisanceParameters.add(nuis_JES)
allNuisancePdfs.append(nuis_JES_pdf)

allSysdF_JES = []
allSysdN_JES = []
for c in range(3):
  
  nco = integrate(bkg,c*30,c*30+29)
  ncu = integrate(hjesu,c*30,c*30+29)
  ncd = integrate(hjesd,c*30,c*30+29)

  ku = (ncu/nco)-1
  kd = (ncd/nco)-1
  
  dN = ROOT.RooFormulaVar("dN_JES_c%d"%(c)," ","TMath::Power(1+%g,@0)*(@0>=0) + TMath::Power(1+%g,-1*@0)*(@0<0)"%(ku,kd),ROOT.RooArgList(nuis_JES))
  allSysdN_JES.append(dN)

  for x in range(30):
     b = c*30+x
     fo = bkg.GetBinContent(b+1)/nco
     fu = hjesu.GetBinContent(b+1)/ncu
     fd = hjesd.GetBinContent(b+1)/ncd

     B = (fu-fd)/(2*fo)
     C = 1 - fu/(2*fo) - fd/(2*fo)

     #df = ROOT.RooFormulaVar("df_JES_c%d_b%d"%(c,b),"TMath::Max(0,1+%g*@0+%g*@0*@0)"%(B,C),ROOT.RooArgList(nuis_JES))
     df = ROOT.RooFormulaVar("df_JES_c%d_b%d"%(c,b),interpString(B,C),ROOT.RooArgList(nuis_JES))
     plot = nuis_JES.frame()
     df.plotOn(plot)
     cv = ROOT.TCanvas();
     plot.Draw()
     cv.SaveAs("plot_df_JES_c%d_b%d.png"%(c,b))
     allSysdF_JES.append(df) 
   
#####################################################################################
# ISR SYSTEMATIC ####################################################################
hisru = fi.Get("isr_up_bkg")
hisrd = fi.Get("isr_dn_bkg")


nuis_ISR = ROOT.RooRealVar("nuis_ISR","nuis_ISR",0,-5,5)
nuis_ISR_In = ROOT.RooRealVar("nuis_ISR_In","nuis_ISR_In",0,-5,5); nuis_ISR_In.setConstant()
nuis_ISR_pdf = ROOT.RooGaussian("nuis_ISR_PDF","",nuis_ISR,nuis_ISR_In,ONE)
allNuisanceParameters.add(nuis_ISR)
allNuisancePdfs.append(nuis_ISR_pdf)
allSysdF_ISR = []
allSysdN_ISR = []
for c in range(3):
  
  nco = integrate(bkg,c*30,c*30+29)
  ncu = integrate(hisru,c*30,c*30+29)
  ncd = integrate(hisrd,c*30,c*30+29)

  ku = (ncu/nco)-1
  kd = (ncd/nco)-1
  
  dN = ROOT.RooFormulaVar("dN_ISR_c%d"%(c),"TMath::Power(1+%g,@0)*(@0>=0) + TMath::Power(1+%g,-1*@0)*(@0<0)"%(ku,kd),ROOT.RooArgList(nuis_ISR))
  allSysdN_ISR.append(dN)

  for x in range(30):
     b = c*30+x
     fo = bkg.GetBinContent(b+1)/nco
     fu = hisru.GetBinContent(b+1)/ncu
     fd = hisrd.GetBinContent(b+1)/ncd

     B = (fu-fd)/(2*fo)
     C = 1 - fu/(2*fo) - fd/(2*fo)

     df = ROOT.RooFormulaVar("df_ISR_c%d_b%d"%(c,b),interpString(B,C),ROOT.RooArgList(nuis_ISR))
     # ok why not make a plot of each of these guys too?
     plot = nuis_ISR.frame()
     df.plotOn(plot)
     cv = ROOT.TCanvas();
     plot.Draw()
     cv.SaveAs("plot_df_ISR_c%d_b%d.png"%(c,b))

     allSysdF_ISR.append(df) 
   
#####################################################################################
# Next we make the "fraction in each category" 
fractions = []
sumFrac   = []
for c in range(3):
  
  nco = integrate(bkg,c*30,c*30+29)
  sums = ROOT.RooArgList()

  for x in range(30):
    b = c*30+x
    fo = bkg.GetBinContent(b+1)/nco

    f = ROOT.RooFormulaVar("f_c%d_b%d"%(c,b),"%g*@0*@1"%fo,ROOT.RooArgList(allSysdF_JES[b],allSysdF_ISR[b]))
    fractions.append(f)
    sums.add(f)

  total = ROOT.RooAddition("total_FRACNORMTERM_c%d"%c,"total for fraction in cat %d"%c,sums)
  sumFrac.append(total)
#####################################################################################
# MC - stat uncerts
# This one is a little easier since we just do a log-normal multiploed for each of the bins
allSysdB_MC = []
allNuis_MC   = []
allNuisIn_MC = []

for b in range(nbins):
    
    nuis_MC = ROOT.RooRealVar("nuis_MC_b%d"%b,"nuis_MC_b%d"%b,0,-5,5)
    nuis_MC_In = ROOT.RooRealVar("nuis_MC_In_b%d"%b,"nuis_MC_In_b%d"%b,0,-5,5); nuis_MC_In.setConstant()

    allNuisIn_MC.append(nuis_MC_In)
    allNuis_MC.append(nuis_MC)
    allNuisanceParameters.add(nuis_MC)
    nuis_MC_pdf = ROOT.RooGaussian("nuis_MC_PDF_b%d"%b,"",allNuis_MC[-1],allNuisIn_MC[-1],ONE)
    allNuisancePdfs.append(nuis_MC_pdf)
    dB = ROOT.RooFormulaVar("dB_mcstat_B%d"%(b),"TMath::Power((1+%g),@0)"%fo,ROOT.RooArgList(allNuis_MC[-1]))
    allSysdB_MC.append(dB)
#####################################################################################
# make a prodPdf of the nuisances 
listPdfs = ROOT.RooArgList(); 
for pdf in allNuisancePdfs: listPdfs.add(pdf)

listPdfs.setName("constraintPdfs")
nuisancePdf = ROOT.RooProdPdf("nuisance_pdf","",listPdfs)

# Ok finally, we are ready to make the expectation per bin and make a Poisson PDF  
expected_Backgrounds = []
expected_Signals     = []
expected_Totals      = []
all_poissons 	     = []
combined_pdf = ROOT.RooSimultaneous("combined_pdf","combined_pdf",sampleType);
for c in range(3):
  
  nco = integrate(bkg,c*30,c*30+29)
  dNISR = allSysdN_ISR[c]
  dNJES = allSysdN_JES[c]
  sumFracNorm = sumFrac[c]

  for x in range(30):
    
    b = c*30+x

    dfMC  = allSysdB_MC[b]
    dfJES = allSysdF_JES[b]
    dfISR = allSysdF_ISR[b]
    fo = bkg.GetBinContent(b+1)/nco
    
    ns = signal.GetBinContent(b+1)

    expectation_b = ROOT.RooFormulaVar("expected_background_bin%d"%b,"(%g*@0*@1*@2)*(%g/@3)*@4*@5"%(nco,fo),ROOT.RooArgList(dNISR,dNJES,dfMC,sumFracNorm,dfJES,dfISR))
    expectation_s = ROOT.RooFormulaVar("expectation_signal_bin%d"%b,"%g*@0"%ns,ROOT.RooArgList(mu)) 
    expected_Backgrounds.append(expectation_b)
    expected_Signals.append(expectation_s)

    expectation   = ROOT.RooAddition("total_expected_bin_%d"%b,"",ROOT.RooArgList(expectation_s,expectation_b))
    expected_Totals.append(expectation)
    
    pois = ROOT.RooPoisson("pdf_bin_%d"%b,"Poisson in bin %d"%b,observation,expectation); 
    all_poissons.append(pois)
    combined_pdf.addPdf(pois,"%d"%b)
  
combined_pdf.Print("v")

output = ROOT.TFile("workspace.root","RECREATE")
wks = ROOT.RooWorkspace("w","w")
getattr(wks,"import")(combined_pdf)
getattr(wks,"import")(nuisancePdf)
getattr(wks,"import")(listPdfs,listPdfs.GetName())
getattr(wks,"import")(allNuisanceParameters,"nuisances")
getattr(wks,"import")(obsdata)
output.WriteTObject(wks)
