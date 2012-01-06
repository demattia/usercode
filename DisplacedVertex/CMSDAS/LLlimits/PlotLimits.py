import ROOT as r
import sys,os

flist = []
for filename in os.listdir(sys.argv[1]):
        if filename.find('.txt')>-1:
                flist.append(float(filename[:-4]))
flist.sort()
print flist
flist = [str(a)+'.txt' for a in flist]

LeptonType = sys.argv[2]

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
i = 0
for filename in flist:
        f = open(sys.argv[1]+"/"+filename)
        items = [eval(a) for a in f.readline().split()]
        m = items[0]
        low = items[1]
        up = items[2]
        g.SetPoint(i,m,0.5*(low+up))
        g.SetPointError(i,0,0.5*(up-low))
        i+=1

# ROOT bug in plotting TGraphErrors (thank you ROOT)
g.SetPoint(i,m,0.5*(low+up))
g.SetPointError(i,0,0.5*(up-low))

r.gROOT.SetStyle("Plain")

if (LeptonType == "Electrons"):
	g.SetTitle("Limit on signal events, electron channel")
elif (LeptonType == "Muons"):
	g.SetTitle("Limit on signal events, muon channel")
else:
	g.SetTitle("Limit")
g.SetFillColor(4)
g.SetFillStyle(3005)
g.GetXaxis().SetRangeUser(20,500)
g.GetXaxis().SetTitle("Mass of X boson [GeV/c^{2}]")
g.GetYaxis().SetTitle("Number of signal events (95% CL)")


r.gStyle.SetErrorX(0)
r.gStyle.SetTitleBorderSize(0)
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

g.SetMaximum(5.0)
g.Draw("a4")
if (LeptonType == "Electrons"):
	g.SetTitle("CMS Preliminary #sqrt{s}=7 TeV L=1.1 fb^{-1}")
	t1 = r.TLatex(30, 4.7, "e^{+}e^{-}")
	t1.Draw()
elif (LeptonType == "Muons"):
	g.SetTitle("CMS Preliminary #sqrt{s}=7 TeV L=1.2 fb^{-1}")
	t1 = r.TLatex(30, 4.7, "#mu^{+}#mu^{-}")
	t1.Draw()
c.Print("limitsBayesianEvents"+LeptonType+".png")

raw_input("Press ENTER to finish...")
