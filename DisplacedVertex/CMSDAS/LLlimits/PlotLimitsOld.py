import ROOT as r
import sys,os
from math import *

# Old version of limit plotter -- this takes the
# signal event limits and converts them directly to
# a sig*BR by dividing. This doesn't handle the uncertainties
# on efficiency and L properly, hence its deprecation.

flist = []
for filename in os.listdir(sys.argv[1]):
        if filename.find('.txt')>-1:
                flist.append(float(filename[:-4]))
flist.sort()
print flist
flist = [str(a)+'.txt' for a in flist]

type = int(sys.argv[2])
br = 0.01
eff_syst_rel = 0.2
if ((type % 10) <= 3):
	lumi = 1116
else:
	lumi = 1085
lumi_rel_err = 0.06
maxmass = 99999
if (type % 10) == 2 or (type % 10) == 6:
	maxmass = 175
elif (type % 10) == 3:
	maxmass = 165
elif (type % 10) == 7:
	maxmass = 70

# Expected limit -- for muons this is just 3.0
if ((type % 10) <= 3):
	exp_lim = 2.99
else:
	tr = r.TRandom3()
	tfc = r.TFeldmanCousins()
	tfc.SetCL(0.95)
	
# We fake the error on the expected limit
# for the muons -- since essentially each mass point
# is its own toy, we can just look at their spread
# (and we have to do it for the electrons too)
exp_lim_err = 0.25*2

rf = r.TFile(sys.argv[1]+"/data.root","READ")
data = rf.Get("modelData")
rf.Close()

npts=0
sum=0
sum2=0

# graph for limits:
numpoints = 0
# first loop to see how many points we want to use
for filename in flist:
        f = open(sys.argv[1]+"/"+filename)
        items = [eval(a) for a in f.readline().split()]
        mass = items[0]
	if mass < maxmass:
		numpoints+=1

g = r.TGraphErrors(numpoints)
gex = r.TGraph(numpoints)
gexp = r.TGraph(numpoints)
gexm = r.TGraph(numpoints)
i = 0
for filename in flist:
        f = open(sys.argv[1]+"/"+filename)
        items = [eval(a) for a in f.readline().split()]
        mass = items[0]
	if mass > maxmass:
		continue

        low = items[1]
        up = items[2]

	# Because of binning effects, the "low" and "high" points are
	# the center of their respective bins, not the edge. Since the
	# low edge should always be zero, the "low" point will thus be half
	# the bin width, and so the upper edge should just be low+up.

	lim = low+up

	# expected limit: for muons constant, for electrons calculate it now
	if ((type % 10) > 3):
		exp_n_sig = 0.169421*exp(-0.5*((mass-89.1762)/2.20176)**2)
		limsum = 0
		limsum2 = 0
		# run "pseudo experiments"
		for iex in range(100):
			nevts = tr.Poisson(exp_n_sig)
			lim_i = tfc.CalculateUpperLimit(nevts, exp_n_sig)
			limsum += lim_i
			limsum2 += lim_i*lim_i
		exp_lim = limsum/100
		exp_lim_err = sqrt(abs(limsum2/100 - exp_lim*exp_lim))
		
		# 2 sigma
		exp_lim_err *= 2
		if (exp_lim_err < 0.50):
			exp_lim_err = 0.50

	this_exp_lim = exp_lim
	this_exp_lim_err = exp_lim_err

	if (type != 0 and type != 4):
		# Correct the efficiency

		if type == 1:
			eff1 = 0.11408 + 0.00356652*mass -8.74025e-06*mass*mass
			eff2 = 0.323721 + 0.00156523*mass -4.63571e-06*mass*mass
			eff1e = 0.003009 + 6.43776e-05*mass + 1.69108e-07*mass*mass
			eff2e = 0.00370533 + 7.0838e-05*mass + 1.82419e-07*mass*mass
		elif type == 2:
			eff1 = -0.00260774 + 0.0102262*mass -5.78105e-05*mass*mass
			eff2 = -0.00144389 + 0.0124665*mass -7.20461e-05*mass*mass
			eff1e = 0.000308325 + 8.96677e-05*mass + 6.18919e-07*mass*mass
			eff2e = 0.000410652 + 0.000108473*mass + 7.36942e-07*mass*mass
		elif type == 3:
			eff1 = 0.163 -0.00093*mass
			eff2 = 0.188733 -0.00116667*mass
			eff1e = 0.00392485 + 0.000101379*mass
			eff2e = 0.00428382 + 0.000110805*mass
		elif type == 5:
			eff1 = 0.159192 + 0.000699302*mass -2.43743e-06*mass*mass
			eff2 = 0.173722 + 0.000539425*mass -2.10161e-06*mass*mass
			eff1e = 0.0026725 + 4.59296e-05*mass + 1.17825e-07*mass*mass
			eff2e = 0.00277939 + 4.53527e-05*mass + 1.15278e-07*mass*mass
		elif type == 6:
			eff1 = 0.0888692 + 0.000539154*mass -5.13077e-06*mass*mass
			eff2 = 0.096 + 0.000488334*mass -5.16667e-06*mass*mass
			eff1e = 0.00395039 + 0.000142049*mass + 7.92113e-07*mass*mass
			eff2e = 0.00429141 + 0.000148919*mass + 8.23049e-07*mass*mass
		elif type == 7:
			eff1 = 0.0253333 -0.000316667*mass
			eff2 = 0.0251 -0.00029*mass
			eff1e = 0.00155242 + 3.60555e-05*mass
			eff2e = 0.00155242 + 3.60555e-05*mass

		# and systematics: 10+x = down lifetime, 20+x = up lifetime

		elif type == 11:
			eff1 = 0.169885 + 0.00451204*mass -1.01795e-05*mass*mass
			eff2 = 0.445134 + 0.00200939*mass -5.10414e-06*mass*mass
			eff1e = 0.00302569 + 5.40459e-05*mass + 1.41119e-07*mass*mass
			eff2e = 0.00339361 + 5.70559e-05*mass + 1.47007e-07*mass*mass
		elif type == 12:
			eff1 = -0.00353672 + 0.0135141*mass -7.21097e-05*mass*mass
			eff2 = 0.000469538 + 0.0153261*mass -8.19953e-05*mass*mass
			eff1e = 0.000307387 + 8.11026e-05*mass + 5.69069e-07*mass*mass
			eff2e = 0.000508459 + 8.33765e-05*mass + 5.81304e-07*mass*mass
		elif type == 13:
			eff1 = 0.2083 -0.00038*mass
			eff2 = 0.244567 -0.000683333*mass
			eff1e = 0.00446331 + 0.000115518*mass
			eff2e = 0.00466714 + 0.000122565*mass
		elif type == 15:
			eff1 = 0.21254 + 0.00140118*mass -3.95117e-06*mass*mass
			eff2 = 0.225071 + 0.00150487*mass -4.36692e-06*mass*mass
			eff1e = 0.00296524 + 5.14757e-05*mass + 1.32461e-07*mass*mass
			eff2e = 0.00305974 + 5.32786e-05*mass + 1.36775e-07*mass*mass
		elif type == 16:
			eff1 = 0.121762 + 0.00103803*mass -7.30513e-06*mass*mass
			eff2 = 0.148223 + 0.000437384*mass -4.67692e-06*mass*mass
			eff1e = 0.0046014 + 0.000165092*mass + 9.22706e-07*mass*mass
			eff2e = 0.00481776 + 0.000172684*mass + 9.6299e-07*mass*mass
		elif type == 17:
			eff1 = 0.0424667 -0.000493333*mass
			eff2 = 0.0336667 -0.000233333*mass
			eff1e = 0.00190933 + 4.53382e-05*mass
			eff2e = 0.00192902 + 4.73756e-05*mass
		elif type == 21:
			eff1 = 0.0628431 + 0.00158201*mass -4.06656e-06*mass*mass
			eff2 = 0.164409 + 0.000633272*mass -2.15358e-06*mass*mass
			eff1e = 0.00206921 + 3.9775e-05*mass + 1.0429e-07*mass*mass
			eff2e = 0.00261729 + 4.4687e-05*mass + 1.14185e-07*mass*mass
		elif type == 22:
			eff1 = -0.00145639 + 0.0051654*mass -3.02208e-05*mass*mass
			eff2 = -0.000146017 + 0.00514694*mass -2.98508e-05*mass*mass
			eff1e = 0.0002055 + 6.00639e-05*mass + 4.12482e-07*mass*mass
			eff2e = 0.000305852 + 6.00627e-05*mass + 4.13374e-07*mass*mass
		elif type == 23:
			eff1 = 0.0741333 -0.000406667*mass
			eff2 = 0.0914 -0.00082*mass
			eff1e = 0.00282528 + 7.08676e-05*mass
			eff2e = 0.0029831 + 7.34091e-05*mass
		elif type == 25:
			eff1 = 0.0695794 + 0.000478943*mass -1.62944e-06*mass*mass
			eff2 = 0.0829491 + 0.000163774*mass -7.84839e-07*mass*mass
			eff1e = 0.00185783 + 3.29663e-05*mass + 8.39743e-08*mass*mass
			eff2e = 0.00190788 + 3.12022e-05*mass + 7.88862e-08*mass*mass
		elif type == 26:
			eff1 = 0.0329154 + 0.00044959*mass -3.51795e-06*mass*mass
			eff2 = 0.0338231 + 0.000572385*mass -4.17692e-06*mass*mass
			eff1e = 0.00264957 + 9.61106e-05*mass + 5.35953e-07*mass*mass
			eff2e = 0.00286634 + 0.000103765*mass + 5.78624e-07*mass*mass
		elif type == 27:
			eff1 = -0.00383333 + 0.000536667*mass
			eff2 = 0.0150333 -0.000236667*mass
			eff1e = 0.00102686 + 3.43188e-05*mass
			eff2e = 0.00103495 + 2.4037e-05*mass
			
		eff_factor = 2*(eff1*br*(1-br) + eff2*br*br)
		eff_err = 2*sqrt(eff1e*eff1e*br*br*(1-br)*(1-br) + eff2e*eff2e*br*br*br*br)

		eff_syst_abs = eff_syst_rel*eff_factor
		eff_err_tot = sqrt(eff_err*eff_err + eff_syst_abs*eff_syst_abs)

		print "eff_factor = " + str(eff_factor) + " +/- " + str(eff_err) + " (stat. only) +/- " + str(eff_err_tot)
		
		eff_rel_err = eff_err_tot/eff_factor

		lim = lim/eff_factor
		lim = lim/lumi
		lim = lim*br

		# Add in errors
		# Note that this is not the correct way of doing things, but
		# it will at least give us some final numbers and should
		# be conservative relative to the actual error.
		lim = lim*(1+lumi_rel_err)

		this_exp_lim = exp_lim*br*(1+lumi_rel_err)/(eff_factor*lumi)
		this_exp_lim_err = exp_lim_err*br*(1+lumi_rel_err)/(eff_factor*lumi)
		
		# print str(eff_factor)+" "+str(low)+" "+str(up)+" "+str(lim)

	npts += 1
	sum += lim
	sum2 += lim*lim
	
	# protect against cases where the efficiency goes negative (or very small)
	if lim > 0 and (lim < 0.5 or type == 0 or type == 4):
		g.SetPoint(i,mass,0.5*lim)
		g.SetPointError(i,0,0.5*lim)
		gex.SetPoint(i, mass, this_exp_lim)
		gexp.SetPoint(i, mass, this_exp_lim + this_exp_lim_err)
		gexm.SetPoint(i, mass, this_exp_lim - this_exp_lim_err)
		print str(i) + " " + str(mass) + " " + str(lim)
		i+=1
		lastposmass = mass
		lastposlim = lim
	# end of points loop

# ROOT bug in plotting TGraphErrors (thank you ROOT)
if lim > 0:
	print str(i) + " " + str(lastposmass) + " " + str(lastposlim)
	g.SetPoint(i,lastposmass,0.5*lastposlim)
	g.SetPointError(i,0,0.5*lastposlim)

mean = sum/npts
stdev = sqrt(sum2/npts - mean*mean)
print str(mean) + " +/- " + str(stdev)

r.gROOT.SetStyle("Plain")

if type==0:
	g.SetTitle("Limit on signal events, muon channel")
elif type==1:
	g.SetTitle("Limit on #sigma #times BR, m_{H} = 1000 GeV/c^{2}, muon channel")
elif type==2:
	g.SetTitle("Limit on #sigma #times BR, m_{H} = 400 GeV/c^{2}, muon channel")
elif type==3:
	g.SetTitle("Limit on #sigma #times BR, m_{H} = 200 GeV/c^{2}, muon channel")
elif type==4:
	g.SetTitle("Limit on signal events, electron channel")
elif type==5:
	g.SetTitle("Limit on #sigma #times BR, m_{H} = 1000 GeV/c^{2}, electron channel")
elif type==6:
	g.SetTitle("Limit on #sigma #times BR, m_{H} = 400 GeV/c^{2}, electron channel")
elif type==7:
	g.SetTitle("Limit on #sigma #times BR, m_{H} = 200 GeV/c^{2}, electron channel")

elif type==11:
	g.SetTitle("Limit on #sigma #times BR, m_{H} = 1000 GeV/c^{2}, muon channel, 1/3 c#tau")
elif type==12:
	g.SetTitle("Limit on #sigma #times BR, m_{H} = 400 GeV/c^{2}, muon channel, 1/3 c#tau")
elif type==13:
	g.SetTitle("Limit on #sigma #times BR, m_{H} = 200 GeV/c^{2}, muon channel, 1/3 c#tau")
elif type==15:
	g.SetTitle("Limit on #sigma #times BR, m_{H} = 1000 GeV/c^{2}, electron channel, 1/3 c#tau")
elif type==16:
	g.SetTitle("Limit on #sigma #times BR, m_{H} = 400 GeV/c^{2}, electron channel, 1/3 c#tau")
elif type==17:
	g.SetTitle("Limit on #sigma #times BR, m_{H} = 200 GeV/c^{2}, electron channel, 1/3 c#tau")

elif type==21:
	g.SetTitle("Limit on #sigma #times BR, m_{H} = 1000 GeV/c^{2}, muon channel, 3 c#tau")
elif type==22:
	g.SetTitle("Limit on #sigma #times BR, m_{H} = 400 GeV/c^{2}, muon channel, 3 c#tau")
elif type==23:
	g.SetTitle("Limit on #sigma #times BR, m_{H} = 200 GeV/c^{2}, muon channel, 3 c#tau")
elif type==25:
	g.SetTitle("Limit on #sigma #times BR, m_{H} = 1000 GeV/c^{2}, electron channel, 3 c#tau")
elif type==26:
	g.SetTitle("Limit on #sigma #times BR, m_{H} = 400 GeV/c^{2}, electron channel, 3 c#tau")
elif type==27:
	g.SetTitle("Limit on #sigma #times BR, m_{H} = 200 GeV/c^{2}, electron channel, 3 c#tau")

g.SetFillColor(4)
g.SetFillStyle(3005)
if (maxmass < 350):
	g.GetXaxis().SetRangeUser(20, maxmass)
else:
	g.GetXaxis().SetRangeUser(20, 350)
g.GetXaxis().SetTitle("Mass of X boson [GeV/c^{2}]")
if type==0 or type==4:
	g.GetYaxis().SetTitle("Number of signal events (95% CL)")
else:
	g.GetYaxis().SetTitle("Cross Section*BR [pb^{-1}] (95% CL)")
g.GetYaxis().SetTitleOffset(1.8)

if type==4:
	g.SetMaximum(3.9)

r.gStyle.SetErrorX(0)
# plotting
c = r.TCanvas("c","c",600,600)
#c.Divide(1,2)
#c.cd(1)
#m = r.RooRealVar("mass","mass",20,350)
#frame = m.frame()
#frame.SetTitle("Data")
#data.plotOn(frame)
#frame.Draw()
#c.cd(2)
r.gPad.SetLeftMargin(0.15)
g.Draw("a4")

gex.SetLineColor(r.kRed)
gex.SetLineWidth(2)
gex.Draw("l same")

gexp.SetLineColor(r.kRed)
gexp.SetLineWidth(2)
gexp.SetLineStyle(r.kDashed)
gexp.Draw("l same")

gexm.SetLineColor(r.kRed)
gexm.SetLineWidth(2)
gexm.SetLineStyle(r.kDashed)
gexm.Draw("l same")

if type==0:
	c.Print("LimitMuonsSigExp.png")
elif type==1:
	c.Print("LimitMuons1000Exp.png")
elif type==2:
	c.Print("LimitMuons400Exp.png")
elif type==3:
	c.Print("LimitMuons200Exp.png")
elif type==4:
	c.Print("LimitElectronsSigExp.png")
elif type==5:
	c.Print("LimitElectrons1000Exp.png")
elif type==6:
	c.Print("LimitElectrons400Exp.png")
elif type==7:
	c.Print("LimitElectrons200Exp.png")
elif type==11:
	c.Print("LimitMuons1000CTau0.33Exp.png")
elif type==12:
	c.Print("LimitMuons400CTau0.33Exp.png")
elif type==13:
	c.Print("LimitMuons200CTau0.33Exp.png")
elif type==15:
	c.Print("LimitElectrons1000CTau0.33Exp.png")
elif type==16:
	c.Print("LimitElectrons400CTau0.33Exp.png")
elif type==17:
	c.Print("LimitElectrons200CTau0.33Exp.png")
elif type==21:
	c.Print("LimitMuons1000CTau3.0Exp.png")
elif type==22:
	c.Print("LimitMuons400CTau3.0Exp.png")
elif type==23:
	c.Print("LimitMuons200CTau3.0Exp.png")
elif type==25:
	c.Print("LimitElectrons1000CTau3.0Exp.png")
elif type==26:
	c.Print("LimitElectrons400CTau3.0Exp.png")
elif type==27:
	c.Print("LimitElectrons200CTau3.0Exp.png")
else:
	print "Oops, didn't understand you"

raw_input("Press ENTER to finish...")
