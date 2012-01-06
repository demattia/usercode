import ROOT as r
from GetParameters import *

# Range of search [10,500] and maximal number of Signal Events [30] HardCoded Here!

# return RooWorkspace given the mass hypothesis

def GetWorkspace(M,LeptonType):

	# GetParameters for mass hypothesis
	pars = GetParameters(M,"dataFile"+LeptonType)
	
	# Create RooRealVars for input parameters (hardcoded for now)
	Bkg = r.RooRealVar("Bkg","Bkg",pars['Bkg'])
	#BkgE = r.RooRealVar("BkgE","BkgE",pars['BkgE']) # for Gaussian
	BkgE = r.RooRealVar("BkgE","BkgE",1 + pars['BkgE']/pars['Bkg']) # for Lognormal
	SigWidth = r.RooRealVar("SigWidth","SigWidth",pars['SigWidth'])
	#SigWidthE = r.RooRealVar("SigWidthE","SigWidthE",pars['SigWidthE']) # for Gaussian 
	SigWidthE = r.RooRealVar("SigWidthE","SigWidthE",1 + pars['SigWidthE']/pars['SigWidth']) # for Lognormal

	# not actually used but needed for dummies
	SigEff = r.RooRealVar("SigEff","SigEff",1)
	#SigEffE = r.RooRealVar("SigEffE","SigEffE",0) # for Gaussian
	SigEffE = r.RooRealVar("SigEffE","SigEffE",1) # for Lognormal

	# observable 
	mass = r.RooRealVar("mass","mass",15,500)

	# some variables to work with, s, b, mean, sigWidth, sigEff
	s = r.RooRealVar("S","S",0,0,50)
	b = r.RooRealVar("B","B",Bkg.getVal(),max(0,Bkg.getVal()-10*BkgE.getVal()),Bkg.getVal()+10*BkgE.getVal())
	b.Print()
	mean = r.RooRealVar("mean","mean",M)
	sigWidth = r.RooRealVar("sigWidth","sigWidth",SigWidth.getVal(),max(0,SigWidth.getVal()-10*pars['SigWidthE']),SigWidth.getVal()+10*pars['SigWidthE'])
	sigEff = r.RooRealVar("sigEff","sigEff",SigEff.getVal(),0,1)
	seff = r.RooFormulaVar("seff","@0*@1",r.RooArgList(s,sigEff))

	# prior PDFs on nuisance parameters Gaussian or Lognormal -> Lognormal preffered by Stat Committee ( no jokes with them )
	#prior_sigWidth = r.RooGaussian("prior_sigWidth","prior_sigWidth",sigWidth,SigWidth,SigWidthE)
	prior_sigWidth = r.RooLognormal("prior_sigWidth","prior_sigWidth",sigWidth,SigWidth,SigWidthE)
	#prior_b = r.RooGaussian("prior_b","prior_b",b,Bkg,BkgE)
	prior_b = r.RooLognormal("prior_b","prior_b",b,Bkg,BkgE)
	#prior_sigEff = r.RooGaussian("prior_sigEff","prior_sigEFf",sigEff,SigEff,SigEffE)
	prior_sigEff = r.RooLognormal("prior_sigEff","prior_sigEFf",sigEff,SigEff,SigEffE)
	#prior_nuisance = r.RooProdPdf("prior_nuisance","prior_nuisance",r.RooArgList(prior_sigWidth,prior_b,prior_sigEff))
	prior_nuisance = r.RooProdPdf("prior_nuisance","prior_nuisance",r.RooArgList(prior_sigWidth,prior_b))

	# shapes of signal and background
	signalPdf = r.RooGaussian("signalPdf","signalPdf",mass,mean,sigWidth)

	weight = r.RooRealVar("weight","weight",1,0,50)
	bkgdataTmp = r.RooDataSet.read("data/masses_backgroundMC_"+LeptonType+".txt",r.RooArgList(mass,weight))
	bkgdata = r.RooDataSet(bkgdataTmp.GetName(),bkgdataTmp.GetTitle(),bkgdataTmp,bkgdataTmp.get(),"",weight.GetName())

	if LeptonType=='Electrons':

		E_fraction = r.RooRealVar("fraction","fraction",0.9,0,1)
		E_mean = r.RooRealVar("m","m",90,80,100)
		E_sigma = r.RooRealVar("s","s",2,0.1,10)
		E_uni = r.RooUniform("uni","uni",r.RooArgSet(mass))
		E_gaus = r.RooGaussian("gaus","gaus",mass,E_mean,E_sigma)
		bkgPdf = r.RooAddPdf("pdf","pdf",E_gaus,E_uni,E_fraction)
		bkgPdf.fitTo(bkgdata,r.RooFit.SumW2Error(r.kTRUE))

	if LeptonType=='Muons':
		bkgPdf = r.RooUniform("bkgPdf","bkgPdf",r.RooArgSet(mass))

	# model 
	#model_no_nuis = r.RooAddPdf("model_no_nuis","model_no_nuis",r.RooArgList(signalPdf,bkgPdf),r.RooArgList(seff,b))
	model_no_nuis = r.RooAddPdf("model_no_nuis","model_no_nuis",r.RooArgList(signalPdf,bkgPdf),r.RooArgList(s,b))
	model = r.RooProdPdf("model","model",model_no_nuis,prior_nuisance)

	# define parameter sets
	observables = r.RooArgSet(mass)
	#nuisanceParameters = r.RooArgSet(b,sigWidth,sigEff)
	nuisanceParameters = r.RooArgSet(b,sigWidth)
	POI = r.RooArgSet(s)

	# modelConfig
	w = r.RooWorkspace('WS'+str(int(M)))
	modelCfg = r.RooStats.ModelConfig('modelCfg',w)
	modelCfg.SetPdf(model)
	modelCfg.SetObservables(observables)
	modelCfg.SetNuisanceParameters(nuisanceParameters)
	modelCfg.SetParametersOfInterest(POI)

	w.Print()

	# Some debugging code -- uncomment to see plots
	if (0):
		# Plot PDFs and plot background data points
		c = r.TCanvas("c","c",600,1200)
		signalPdf.createHistogram("mass").Draw();
		c.Update()
		print "signal pdf"
		r.cin.get()
		bkgPdf.createHistogram("mass").Draw();
		c.Update()
		print "background pdf"
		r.cin.get()
		
		frame = mass.frame()
		
		frame.SetTitle("Background")
		bkgdata.plotOn(frame)
		frame.Draw()
		c.Update()
		print "background points"
		r.cin.get()

	return w
