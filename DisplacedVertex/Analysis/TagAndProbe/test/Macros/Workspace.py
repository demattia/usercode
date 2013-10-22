__author__ = 'demattia'

import ROOT
from ROOT import RooRealVar, RooFormulaVar, RooVoigtian, RooChebychev, RooArgList, RooArgSet, \
    RooAddPdf, RooDataSet, RooCategory, RooSimultaneous, RooGenericPdf, RooWorkspace

def buildPdf(ws, p):

    mass = RooRealVar("mass", "mass", p.minMass, p.maxMass)
    # ws.import(mass)
    getattr(ws,'import')(mass)

    # Construct signal pdf
    mean = RooRealVar("mean", "mean", 90, 80, 100)
    width = RooRealVar("width", "width", 2.4952, 1, 3)
    width.setConstant(ROOT.kTRUE)
    sigma = RooRealVar("sigma", "sigma", 1.2, 0.2, 2)
    signalAll = RooVoigtian("signalAll", "signalAll", mass, mean, width, sigma)
    # Construct background pdf
    # a0_all = RooRealVar("a0_all","a0_all",-0.1,-1,1)
    # a1_all = RooRealVar("a1_all","a1_all",0.004,-1,1)
    # backgroundAll = RooChebychev("backgroundAll","backgroundAll",mass,RooArgList(a0_all,a1_all))

    turnOnAll = RooRealVar("turnOnAll","turnOnAll", 80., 50., 150.)
    widthAll_bkg = RooRealVar("widthAll","widthAll", 2., 0., 50.)
    decayAll_bkg = RooRealVar("decayAll","decayAll", 80., 50., 150.)
    backgroundAll = RooGenericPdf("backgroundAll","backgroundAll",
                                  "exp(-@0/@3)*(TMath::Erf((@0-@1)/@2)+1)",
                                  RooArgList(mass, turnOnAll, widthAll_bkg, decayAll_bkg));

    # Construct composite pdf
    sigAll = RooRealVar("sigAll", "sigAll", 2000, 0, 100000)
    bkgAll = RooRealVar("bkgAll", "bkgAll", 100, 0, 1000)
    modelAll = RooAddPdf("modelAll", "modelAll", RooArgList(signalAll, backgroundAll), RooArgList(sigAll, bkgAll))
    # if MC:
    #     modelAll = RooAddPdf("modelAll", "modelAll", RooArgList(signalAll), RooArgList(sigAll))

    # Define pdf for all probes

    # Construct signal pdf.
    # NOTE that sigma is shared with the signal sample model
    signalPass = RooVoigtian("signalPass","signalPass",mass,mean,width,sigma)
    # Construct the background pdf
    # a0_pass = RooRealVar("a0_pass","a0_pass",-0.1,-1,1)
    # a1_pass = RooRealVar("a1_pass","a1_pass",0.5,-0.1,1)
    # backgroundPass = RooChebychev("backgroundPass","backgroundPass",mass,RooArgList(a0_pass,a1_pass))

    # turnOnPass = RooRealVar("turnOnPass","turnOnPass", 80., 50., 150.)
    # widthPass_bkg = RooRealVar("widthPass","widthPass", 2., 0., 50.)
    # decayPass_bkg = RooRealVar("decayPass","decayPass", 80., 50., 150.)
    # backgroundPass = RooGenericPdf("backgroundPass","backgroundPass",
    #                                "exp(-@0/@3)*(TMath::Erf((@0-@1)/@2)+1)",
    #                                RooArgList(mass, turnOnPass, widthPass_bkg, decayPass_bkg));
    backgroundPass = RooGenericPdf("backgroundPass","backgroundPass",
                                   "exp(-@0/@3)*(TMath::Erf((@0-@1)/@2)+1)",
                                   RooArgList(mass, turnOnAll, widthAll_bkg, decayAll_bkg));

    # Construct the composite model
    efficiency = RooRealVar("efficiency","efficiency",0.8,0.,1.)
    sigPass = RooFormulaVar("sigPass", "@0*@1", RooArgList(sigAll, efficiency))
    bkgPass = RooRealVar("bkgPass", "bkgPass", 100, 0, 1000)
    modelPass = RooAddPdf("modelPass", "modelPass", RooArgList(signalPass, backgroundPass), RooArgList(sigPass, bkgPass))
    # if MC:
    #     modelPass = RooAddPdf("modelPass", "modelPass", RooArgList(signalPass), RooArgList(sigPass))

    frac = RooRealVar("frac", "frac", 0.8, 0., 1.)

    # Define combined pdf for simultaneous fit

    # Define category to distinguish physics and control samples events
    sample = RooCategory("sample","sample")
    sample.defineType("all")
    sample.defineType("pass")

    simPdf = RooSimultaneous("simPdf","simultaneous pdf",sample)

    # Associate model with the physics state and model_ctl with the control state
    simPdf.addPdf(modelAll,"all")
    simPdf.addPdf(modelPass,"pass")
    # ws.import(simPdf)
    getattr(ws,'import')(simPdf)

