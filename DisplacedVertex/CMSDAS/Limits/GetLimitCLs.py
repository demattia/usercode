import ROOT as r
import sys,os.path

from GetParameters import *

def GetLimit(M, eff, efferr, LeptonType, fileName):



	# GetParameters for mass hypothesis
	pars = GetParameters(M,"dataFile"+LeptonType)
	
	# Create RooRealVars for input parameters (hardcoded for now)
	Bkg = r.RooRealVar("Bkg","Bkg",pars['Bkg'])
	BkgE = r.RooRealVar("BkgE","BkgE",1 + pars['BkgE']/pars['Bkg']) # for lognormal
	SigWidth = r.RooRealVar("SigWidth","SigWidth",pars['SigWidth'])
	SigWidthE = r.RooRealVar("SigWidthE","SigWidthE",1 + pars['SigWidthE']/pars['SigWidth']) # for lognormal
	SigEff = r.RooRealVar("SigEff","SigEff",eff)

	# now we read the sigeff directly from the command line
	SigEffE = r.RooRealVar("SigEffE","SigEffE",1 + efferr/eff) # for lognormal

	Lumi = r.RooRealVar("Lumi", "Lumi", pars['Lumi'])
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

	prior_sigWidth = r.RooLognormal("prior_sigWidth","prior_sigWidth",sigWidth,SigWidth,SigWidthE)
	prior_b = r.RooLognormal("prior_b","prior_b",b,Bkg,BkgE)
	prior_sigEff = r.RooLognormal("prior_sigEff","prior_sigEFf",sigEff,SigEff,SigEffE)
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





	# print("Output will be written to " +fileName)
	# M = w.var("mean").getVal()
	# nuisanceParameters = w.set('modelCfg_NuisParams')
	# observables = w.set('modelCfg_Observables')
	# POI = w.set('modelCfg_POI')
	#w.Print()

	# mass = w.arg("mass")
	data = r.RooDataSet.read("data/masses_data_"+LeptonType+".txt",r.RooArgList(mass)) 
	
	'''
	frame = mass.frame()
	data.plotOn(frame)

	bkgpdf = w.pdf("model")
	bkgpdf.plotOn(frame)
	frame.Draw()
	'''

	# modelCfg = r.RooStats.ModelConfig('modelCfg',w)
	# modelCfg.SetPdf(w.pdf("model"))
	# modelCfg.SetNuisanceParameters(nuisanceParameters)
	# modelCfg.SetObservables(observables)
	# modelCfg.SetParametersOfInterest(POI)
	
	
	# Bayesian limits
	prior_S = r.RooUniform("prior_S","prior_S",r.RooArgSet(w.var("S")))
	w.var("S").Print()
	modelCfg.SetPriorPdf(prior_S)
	bc = r.RooStats.BayesianCalculator(data,modelCfg)
	bc.SetConfidenceLevel(0.95)
	bc.SetLeftSideTailFraction(0.0)
	IntBC = bc.GetInterval()
	
	# Get the result
	outf = open(fileName+'BC.txt','write')
	print "Bayes" + str(IntBC.LowerLimit()) + " " + str(IntBC.UpperLimit())
	outf.write(str(M) + " " + str(IntBC.LowerLimit()) + " " + str(IntBC.UpperLimit()) +'\n')
	

	# CLs limits
	data.SetName("modelData")
        getattr(w,'import')(data)
        getattr(w,'import')(modelCfg)

        # CLs calculation
        r.gROOT.LoadMacro("roostats_cl95.C+")
	Int = r.RunInverter(w, "modelCfg","","modelData",0,0,50,0,IntBC.UpperLimit()*5,2000,True)

	q = [0]*8
        q[0] = M
        q[1] = Int.UpperLimit()
        q[2] = Int.UpperLimitEstimatedError()
        q[3] = Int.GetExpectedUpperLimit(0)
        q[4] = Int.GetExpectedUpperLimit(-1)
	q[5] = Int.GetExpectedUpperLimit(1)
        q[6] = Int.GetExpectedUpperLimit(-2)
        q[7] = Int.GetExpectedUpperLimit(2)

	outf = open(fileName+'.txt','write')
	import pickle
	pickle.dump(q,outf)

