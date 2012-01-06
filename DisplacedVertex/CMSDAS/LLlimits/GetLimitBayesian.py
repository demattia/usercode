import ROOT as r
import sys,os.path

def GetLimit(w,LeptonType,fileName):
	print("Output will be written to " +fileName)

	M = w.var("mean").getVal()
	nuisanceParameters = w.set('modelCfg_NuisParams')
	observables = w.set('modelCfg_Observables')
	POI = w.set('modelCfg_POI')
	#w.Print()

	mass = w.arg("mass")
	data = r.RooDataSet.read("data/masses_data_"+LeptonType+".txt",r.RooArgList(mass)) 
	frame = mass.frame()
	data.plotOn(frame)

	bkgpdf = w.pdf("model")
	bkgpdf.plotOn(frame)
	frame.Draw()

	modelCfg = r.RooStats.ModelConfig('modelCfg',w)
	modelCfg.SetPdf(w.pdf("model"))
	modelCfg.SetNuisanceParameters(nuisanceParameters)
	modelCfg.SetObservables(observables)
	modelCfg.SetParametersOfInterest(POI)
	
	# Get the result

	prior_S = r.RooUniform("prior_S","prior_S",r.RooArgSet(w.var("S")))
	w.var("S").Print()
	modelCfg.SetPriorPdf(prior_S)
	bc = r.RooStats.BayesianCalculator(data,modelCfg)
	bc.SetConfidenceLevel(0.95)
	bc.SetLeftSideTailFraction(0.0)
	Int = bc.GetInterval()

	print w.var("B").getVal()

	# Use RooStats to get the expected limit
	# Note -- this assumes that the lumi and efficiency errors
	# are set up for lognormal. If you switch GetWorkspace back
	# to Gaussian you'll need to change this code here.

	## Whoa is this a hack!
	if (fileName.find("/") > -1):
		eff = w.obj("SigEff").getVal()
		efferr = (w.obj("SigEffE").getVal()-1)*eff
		lum = w.obj("Lumi").getVal()
		lumerr = (w.obj("LumiE").getVal()-1)*lum
	else:
		eff = 1
		efferr = 0
		lum = 1
		lumerr = 0

	if LeptonType == "Muons":
		# massrange = mass.getMax() - mass.getMin()
		# bkg = w.var("B").getVal()/massrange
		# use values from background fit
		bkg = 0.02
		bkgerr = 0.78
	else:
		# mev = r.RooRealVar("mev", "mev", M, 15, 500)
		# obs = r.RooArgSet(mev)
		# bkg = bkgpdf.getVal(obs)
		# use values from background fit now
		if (M > 82.5 and M < 97.5):
			bkg = 0.70
			bkgerr = 0.89
		else:
			bkg = 0.08
			bkgerr = 0.10
		
	# print "Background value for expected estimation is " + str(bkg)

        r.gROOT.LoadMacro("roostats_cl95.C+")
	limit_result = r.roostats_clm(lum, lumerr, eff, efferr, bkg, bkgerr, 1000, 1)

	q = [0]*7
        q[0] = M
        q[1] = Int.UpperLimit()
	q[2] = limit_result.GetExpectedLimit()
	q[3] = limit_result.GetOneSigmaHighRange()
	q[4] = limit_result.GetOneSigmaLowRange()
	q[5] = limit_result.GetTwoSigmaHighRange()
	q[6] = limit_result.GetTwoSigmaLowRange()

	outf = open(fileName+'.txt','write')
	import pickle
	pickle.dump(q,outf)


