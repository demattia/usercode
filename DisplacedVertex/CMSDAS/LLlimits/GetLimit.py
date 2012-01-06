import ROOT as r
import sys,os.path

def GetLimit(w,LeptonType):

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
	

	prior_S = r.RooUniform("prior_S","prior_S",r.RooArgSet(w.var("S")))
	w.var("S").Print()
	modelCfg.SetPriorPdf(prior_S)
	bc = r.RooStats.BayesianCalculator(data,modelCfg)
	bc.SetConfidenceLevel(0.95)
	bc.SetLeftSideTailFraction(0.0)
	Int = bc.GetInterval()
	'''
		
	fc = r.RooStats.FeldmanCousins(data,modelCfg)
	fc.SetConfidenceLevel(0.95)
	fc.SetNBins(300)
	fc.UseAdaptiveSampling(1)
#	fc.AdditionalNToysFactor(0.3)
	Int = fc.GetInterval()
	'''	

	# Get the result

	print w.var("B").getVal()
	outf = open(str(M)+'.txt','write')
	#outf.write(str(M) + " " + str(Int.LowerLimit(w.var('S'))) + " " + str(Int.UpperLimit(w.var('S'))) +'\n')
	outf.write(str(M) + " " + str(Int.LowerLimit()) + " " + str(Int.UpperLimit()) +'\n')

