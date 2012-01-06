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
	
	'''
	frame = mass.frame()
	data.plotOn(frame)

	bkgpdf = w.pdf("model")
	bkgpdf.plotOn(frame)
	frame.Draw()
	'''

	modelCfg = r.RooStats.ModelConfig('modelCfg',w)
	modelCfg.SetPdf(w.pdf("model"))
	modelCfg.SetNuisanceParameters(nuisanceParameters)
	modelCfg.SetObservables(observables)
	modelCfg.SetParametersOfInterest(POI)
	
	
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

