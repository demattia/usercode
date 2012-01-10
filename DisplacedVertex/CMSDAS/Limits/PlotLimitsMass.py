import ROOT as r
import sys,os
import pickle

r.gROOT.LoadMacro("tdrstyle.C")
r.setTDRStyle()
	
# This is a somewhat hacky script to plot the
# limits vs. mass. It's even hackier for the CLs
# limits -- we read in the Bayesian files to set
# the mass but the actual numbers are hardcoded.


if len(sys.argv) < 2:
	sys.exit('Usage: python PlotLimitsMass.py LeptonType [CLs]')

LeptonType = sys.argv[1]
FileDir = LeptonType + "BayesianMass"

if len(sys.argv) > 2 and sys.argv[2] == "CLs":
	clmode = 1
	FileInfix = "CLs"
else:
	clmode = 0
	FileInfix = "Bayesian"


flist = []
for filename in os.listdir(FileDir):
        if filename.find('.txt')>-1:
                flist.append(float(filename[:-4]))
flist.sort()
print flist
flist = [str(a)+'.txt' for a in flist]

mass = r.RooRealVar("mass","mass",20,500)
data = r.RooDataSet.read("data/masses_data_"+LeptonType+".txt",r.RooArgList(mass))

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
        bkgPdf.fitTo(bkgdata)

if LeptonType=='Muons':
        bkgPdf = r.RooUniform("bkgPdf","bkgPdf",r.RooArgSet(mass))

# graph for limits:
g = r.TGraphErrors(len(flist))
gest = r.TGraphErrors(len(flist))
i = 0
for filename in flist:
        f = open(FileDir+"/"+filename)
        #items = [eval(a) for a in f.readline().split()]
        #m = items[0]
        #low = items[1]
        #up = items[2]

	q = pickle.load(f)
	m = q[0]
	up = q[1]
	est = q[2]
	est1p = q[3]
	est1m = q[4]
	est2p = q[5]
	est2m = q[6]

	# Now fix the numbers for the CLs:
	if (clmode == 1):
		if (LeptonType == "Muons"):
			up = 3.03029
			est1p = 6.43634
			est1m = 3.03029
		if (LeptonType == "Electrons"):
			if (m == 90):
				up = 4.75039
			else:
				up = 3.03029
			if (m > 82.5 and m < 97.5):
				est1p = 6.43634
				est1m = 3.03029
			else:
				est1p = 3.03029
				est1m = 3.03029
	
        g.SetPoint(i,m,up)
        #g.SetPointError(i,0,0.5*(up-low))

	gest.SetPoint(i, m, (est1p+est1m)/2)
	gest.SetPointError(i, 0, (est1p-est1m)/2)

        i+=1

# ROOT bug in plotting TGraphErrors (thank you ROOT)
g.SetPoint(i,m,up)
gest.SetPoint(i, m, (est1p+est1m)/2)
gest.SetPointError(i, 0, (est1p-est1m)/2)

if (LeptonType == "Electrons"):
	gest.SetTitle("Limit on signal events, electron channel")
elif (LeptonType == "Muons"):
	gest.SetTitle("Limit on signal events, muon channel")
else:
	gest.SetTitle("Limit")

g.SetLineColor(r.kBlue)
g.SetFillColor(0)
g.SetMarkerColor(r.kBlue)
g.SetMarkerStyle(r.kFullSquare)
g.SetMarkerSize(1)
g.SetLineWidth(2)

gest.GetXaxis().SetRangeUser(20,500)
gest.GetXaxis().SetTitle("Mass of X boson [GeV/c^{2}]")
gest.GetYaxis().SetTitle("Number of signal events (95% CL)")

#r.gROOT.SetStyle("Plain")
r.gStyle.SetErrorX(0)
# have to restore the title though
r.gStyle.SetOptTitle(1)
r.gStyle.SetTitleBorderSize(0)
#r.gStyle.SetTitleFont(42, "")

gest.SetFillColor(r.kRed)
gest.SetFillStyle(3005)

# plotting
#c = r.TCanvas("c","c",600,1200)
#c.Divide(1,2)
#c.cd(1)
#frame = mass.frame()
#frame.SetTitle("Data")
#data.plotOn(frame)
#bkgPdf.plotOn(frame)
#frame.Draw()
#c.cd(2)
c = r.TCanvas("c","c",600,600)
r.gPad.SetTopMargin(0.10)
r.gPad.SetRightMargin(0.04)

gest.SetMinimum(0)
gest.SetMaximum(6)
gest.Draw("a3")
g.Draw("lp same")

leg = r.TLegend(0.4, 0.17, 0.88, 0.35)
leg.SetFillColor(0)
leg.SetBorderSize(0)
leg.AddEntry(g, "Observed limit")
leg.AddEntry(gest, "Expected limit (#pm 1 #sigma)")
leg.SetTextFont(42)
leg.Draw()

if (LeptonType == "Electrons"):
	gest.SetMinimum(1.5)
	gest.SetTitle("CMS Preliminary #sqrt{s}=7 TeV L=1.1 fb^{-1}")
	if (clmode == 0):
		t1 = r.TLatex(30, 5.5, "e^{+}e^{-}")
	else:
		gest.SetMinimum(1)
		gest.SetMaximum(7)
		t1 = r.TLatex(30, 6.5, "e^{+}e^{-}")
	t1.Draw()

	p = r.TPad("p", "p", 0.40, 0.40, 0.88, 0.88)
	p.Draw()
	p.cd()
	ge1 = gest.Clone()
	ge1.SetMinimum(2)
	ge1.SetMaximum(6.5+0.5*clmode)
	ge1.SetTitle("")
	ge1.GetXaxis().SetTitle("")
	ge1.GetYaxis().SetTitle("")
	ge1.GetXaxis().SetRangeUser(85,95)
	ge1.Draw("a3")
	g1 = g.Clone()
	g1.Draw("lp same")
	
elif (LeptonType == "Muons"):
	gest.SetTitle("CMS Preliminary #sqrt{s}=7 TeV L=1.2 fb^{-1}")
	if (clmode == 0):
		t1 = r.TLatex(30, 5.5, "#mu^{+}#mu^{-}")
	else:
		gest.SetMaximum(8)
		t1 = r.TLatex(30, 7.5, "#mu^{+}#mu^{-}")
	t1.Draw()
t1.SetTextFont(42)

c.Print("limits"+FileInfix+"Events"+LeptonType+"Exp2.png")
c.Print("limits"+FileInfix+"Events"+LeptonType+"Exp2.pdf")

#raw_input("Press ENTER to finish...")
