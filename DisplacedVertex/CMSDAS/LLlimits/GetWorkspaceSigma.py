import ROOT as r
from GetParameters import *

# Be aware of the following hardcoded parameters:
#  * range of mass fit [15,500]
#  * details of electron background PDF (see below)
#
# The upper bound of the sigma search space is now
# dynamically determined -- given that we expect a limit
# of approximately 3/(lumi*eff) I've set the upper bound at
# 10x this.

# return RooWorkspace given the mass hypothesis

def GetWorkspaceSigma(M,LeptonType,eff,efferr):

	# GetParameters for mass hypothesis
	pars = GetParameters(M,"dataFile"+LeptonType)
	
	# Create RooRealVars for input parameters (hardcoded for now)
	Bkg = r.RooRealVar("Bkg","Bkg",pars['Bkg'])
	#BkgE = r.RooRealVar("BkgE","BkgE",pars['BkgE']) # for Gaussian
	BkgE = r.RooRealVar("BkgE","BkgE",1 + pars['BkgE']/pars['Bkg']) # for lognormal
	SigWidth = r.RooRealVar("SigWidth","SigWidth",pars['SigWidth'])
	#SigWidthE = r.RooRealVar("SigWidthE","SigWidthE",pars['SigWidthE']) # for Gaussian 
	SigWidthE = r.RooRealVar("SigWidthE","SigWidthE",1 + pars['SigWidthE']/pars['SigWidth']) # for lognormal
	SigEff = r.RooRealVar("SigEff","SigEff",eff)

	# now we read the sigeff directly from the command line
	#SigEffE = r.RooRealVar("SigEffE","SigEffE",pars['SigEffE']*eff) # for Gaussian
	#SigEffE = r.RooRealVar("SigEffE","SigEffE",1 + pars['SigEffE']) # for lognormal
	#SigEffE = r.RooRealVar("SigEffE","SigEffE",efferr) # for Gaussian
	SigEffE = r.RooRealVar("SigEffE","SigEffE",1 + efferr/eff) # for lognormal

	Lumi = r.RooRealVar("Lumi", "Lumi", pars['Lumi'])
	#LumiE = r.RooRealVar("LumiE", "LumiE", pars['LumiErr']*Lumi.getVal()) # for Gaussian
	LumiE = r.RooRealVar("LumiE", "LumiE", 1+pars['LumiErr']) # for lognormal

	# observable 
	mass = r.RooRealVar("mass","mass",15,500)

	# some variables to work with, s, b, mean, sigWidth, sigEff
	sUpperBound = 30/(pars['Lumi']*eff)
	print "upper limit of search space is: " + str(sUpperBound)
	s = r.RooRealVar("S","S",0,0,sUpperBound)
	b = r.RooRealVar("B","B",Bkg.getVal(),max(0,Bkg.getVal()-5*pars['BkgE'],Bkg.getVal()+5*pars['BkgE']))
	b.Print()
	mean = r.RooRealVar("mean","mean",M)
	sigWidth = r.RooRealVar("sigWidth","sigWidth",SigWidth.getVal(),SigWidth.getVal()-5*pars['SigWidthE'],SigWidth.getVal()+5*pars['SigWidthE'])
	sigEff = r.RooRealVar("sigEff","sigEff",SigEff.getVal(),SigEff.getVal()*(1-5*pars['SigEffE']),SigEff.getVal()*(1+5*pars['SigEffE']))
	lumi = r.RooRealVar("lumi","lumi",Lumi.getVal(),Lumi.getVal()*(1-5*pars['LumiErr']),Lumi.getVal()*(1+5*pars['LumiErr']))
	seff = r.RooFormulaVar("seff","@0*@1*@2",r.RooArgList(s,sigEff,lumi))

	# prior PDFs on nuisance parameters; now use lognormal as preferred by statistics committee
	#prior_sigWidth = r.RooGaussian("prior_sigWidth","prior_sigWidth",sigWidth,SigWidth,SigWidthE)
	prior_sigWidth = r.RooLognormal("prior_sigWidth","prior_sigWidth",sigWidth,SigWidth,SigWidthE)
	#prior_b = r.RooGaussian("prior_b","prior_b",b,Bkg,BkgE)
	prior_b = r.RooLognormal("prior_b","prior_b",b,Bkg,BkgE)
	#prior_sigEff = r.RooGaussian("prior_sigEff","prior_sigEFf",sigEff,SigEff,SigEffE)
	prior_sigEff = r.RooLognormal("prior_sigEff","prior_sigEFf",sigEff,SigEff,SigEffE)
	#prior_lumi = r.RooGaussian("prior_lumi", "prior_lumi", lumi, Lumi, LumiE)
	prior_lumi = r.RooLognormal("prior_lumi", "prior_lumi", lumi, Lumi, LumiE)

	prior_nuisance = r.RooProdPdf("prior_nuisance","prior_nuisance",r.RooArgList(prior_sigWidth,prior_b,prior_sigEff,prior_lumi))

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
	model_no_nuis = r.RooAddPdf("model_no_nuis","model_no_nuis",r.RooArgList(signalPdf,bkgPdf),r.RooArgList(seff,b))
	model = r.RooProdPdf("model","model",model_no_nuis,prior_nuisance)

	# define parameter sets
	observables = r.RooArgSet(mass)
	nuisanceParameters = r.RooArgSet(b,sigWidth,sigEff,lumi)
	POI = r.RooArgSet(s)

	# modelConfig
	w = r.RooWorkspace('WS'+str(int(M)))
	modelCfg = r.RooStats.ModelConfig('modelCfg',w)
	modelCfg.SetPdf(model)
	modelCfg.SetObservables(observables)
	modelCfg.SetNuisanceParameters(nuisanceParameters)
	modelCfg.SetParametersOfInterest(POI)

	w.Print()

	return w
