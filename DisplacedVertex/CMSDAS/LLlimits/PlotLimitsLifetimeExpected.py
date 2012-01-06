import ROOT as r
import sys,os
import pickle

r.gROOT.LoadMacro("tdrstyle.C")
r.setTDRStyle()

# Configurable parameters
# Note that the output of mainSigma.py is simply sigma,
# so we have to add the BR back in

br = 0.01
logx = 1
logy = 1
drawExpected = 1

if len(sys.argv) < 3:
	sys.exit('Usage: python PlotLimitsLifetime.py LeptonType HiggsMass')

LeptonType = sys.argv[1]
hmass = int(sys.argv[2])
FileDir = LeptonType

if len(sys.argv) > 3 and sys.argv[3] == "CLs":
	clmode = 1
	FileDir += "CLs"
	FileInfix = "CLs"
	# br = 1 # don't need to add the BR back in in this case
else:
	clmode = 0
	FileInfix = "Bayesian"
	FileDir += "Bayesian"
# clmode == -1 is kept to deal with old behavior. Not supported any more.
	
if (LeptonType == "Electrons"):
	drawExpected = 0

ymin = 0
if (hmass == 1000):
	xmasses = ["350", "150", "050", "020"]
	xmin = 0.04
	xmax = 1100
	ymax = 0.08
	if (logy == 1):
		if (LeptonType == "Electrons"):
			ymin = 0.001
			ymax = 500.0
		else:
			ymin = 0.001
			ymax = 500.0
elif (hmass == 400):
	xmasses = ["150", "050", "020"]
	xmin = 0.1
	xmax = 1500
	ymax = 0.08
	if (logy == 1):
		if (LeptonType == "Electrons"):
			ymin = 0.004
			ymax = 500.0
		else:
			ymin = 0.0005
			ymax = 50.0
elif (hmass == 200):
	xmasses = ["050", "020"]
	xmin = 0.1
	xmax = 1000
	ymax = 1
	if (logy == 1):
		if (LeptonType == "Electrons"):
			ymin = 0.005
			ymax = 800.0
		else:
			ymin = 0.001
			ymax = 80.0
else:
	print "Oops, mass not recognized"

g = {}
gest1 = {}
gest2 = {}
j = 1

#r.gROOT.SetStyle("Plain")
r.gStyle.SetErrorX(0)
# have to restore the title though
r.gStyle.SetOptTitle(1)
r.gStyle.SetTitleBorderSize(0)
#r.gStyle.SetTitleFont(42, "")

c = r.TCanvas("c","c",600,600)
r.gPad.SetLeftMargin(0.15)
r.gPad.SetTopMargin(0.10)
r.gPad.SetRightMargin(0.04)
fr = c.DrawFrame(xmin, ymin, xmax, ymax)
if (LeptonType == "Electrons"):
	fr.SetTitle("CMS Preliminary #sqrt{s}=7 TeV L=1.1 fb^{-1}")
	t1 = r.TLatex(xmin*2, ymax/3, "e^{+}e^{-}")
	t1.Draw()
elif (LeptonType == "Muons"):
	fr.SetTitle("CMS Preliminary #sqrt{s}=7 TeV L=1.2 fb^{-1}")
	t1 = r.TLatex(xmin*2, ymax/3, "#mu^{+}#mu^{-}")
	t1.Draw()
t1.SetTextFont(42)

fr.GetXaxis().SetTitle("c#tau [cm]")
fr.GetYaxis().SetTitle("Cross Section #times BR [pb] (95% CL)")
#fr.GetYaxis().SetTitleOffset(1.8)
fr.GetXaxis().SetTitleOffset(1.05)
if (logx == 1):
	r.gPad.SetLogx()
if (logy == 1):
	r.gPad.SetLogy()

leg = r.TLegend(0.4, 0.6, 0.88, 0.88)
leg.SetFillColor(0)
leg.SetBorderSize(0)
dummy = r.TH1F("dummy", "dummy", 1, 1, 1)
dummy.SetLineColor(0)
dummy.SetMarkerColor(0)
dummy.SetFillColor(0)
leg.AddEntry(dummy, "m_{H} = "+str(hmass)+ " GeV/c^{2}")

for xmass in xmasses:
	if clmode == -1 and xmass[0] == "0":
		xmass = xmass[1:]
	flist = []
	for filename in os.listdir(FileDir):
		if filename.find('.txt')>-1 and filename.find(str(hmass)+"_"+str(xmass))==0 and filename.find("BC")==-1:
			flist.append(filename[:-4])

	flist.sort(key=lambda str: float(str.split("_")[2]))
	print "Files found for " + str(xmass) + ": " + str(flist)
	flist = [str(a)+'.txt' for a in flist]

	# graph for limits:
	g[xmass] = r.TGraph(len(flist))
	gest1[xmass] = r.TGraphErrors(len(flist))
	gest2[xmass] = r.TGraphErrors(len(flist))
	lmax=-1
	lmin=9999
	i = 0
	for filename in flist:
		lt=float(filename[:-4].split("_")[2])
		# print "Lifetime in "+filename+" is " + str(lt)
		if lt > lmax:
			lmax = lt
		if lt < lmin:
			lmin = lt

		f = open(FileDir+"/"+filename)

		if clmode != -1:
			#q = [eval(a) for a in f.readline().split()]
			#m = q[0]
			#up = q[2]
			#est = q[2]
			#est1 = 0
			#est2 = 0

			q = pickle.load(f)
			m = q[0]
			up = q[1]
			est = q[2]
			est1p = q[3]
			est1m = q[4]
			est2p = q[5]
			est2m = q[6]

			# Note: sometimes the Bayesian limit
			# fails to compute for whatever reason.
			# In this case we just use the ratio from
			# the events plot since it should be the same.

			if (drawExpected):
				if (clmode == 0):
					est = up
					est1p = up*1.4636
					est1m = up
				else:
					# On the other hand, the CLs limits are hardcoded here
					est = up
					est1p = up*6.43634/3.03029
					est1m = up

		else:
			# old CLs parsing code
			q = [eval(a) for a in f.readline().split()]
			m = q[0]
			up = q[1]
			est = q[1]   # a little bit of cheating here
			est1p = q[2] - q[1]
			est1m = q[1] - q[3]
			est1 = (est1p + est1m)/2

		# don't know why this is necessary, but it is
		if (est1m > up): est1m = up
		
		g[xmass].SetPoint(i,lt,up*br)

		gest1[xmass].SetPoint(i,lt,br*(est1p+est1m)/2)
		gest1[xmass].SetPointError(i,0,br*(est1p-est1m)/2)
		gest2[xmass].SetPoint(i,lt,br*(est2p+est2m)/2)
		gest2[xmass].SetPointError(i,0,br*(est2p-est2m)/2)
		i+=1

	# Probably don't need this any more
	g[xmass].SetPoint(i,lt,up*br)
	gest1[xmass].SetPoint(i,lt,br*(est1p+est1m)/2)
	gest1[xmass].SetPointError(i,0,br*(est1p-est1m)/2)
	gest2[xmass].SetPoint(i,lt,br*(est2p+est2m)/2)
	gest2[xmass].SetPointError(i,0,br*(est2p-est2m)/2)

	# Draw in order: +/- 2 sigma band, +/- 1 sigma band, actual

	#gest2[xmass].SetFillColor(r.kYellow)
	#gest2[xmass].SetFillStyle(3005)
	#gest2[xmass].Draw("4 same")
	gest1[xmass].SetFillColor(r.kRed)
	gest1[xmass].SetFillStyle(3005)
	gest1[xmass].SetMarkerStyle(0)
	gest1[xmass].SetMarkerColor(0)
	gest1[xmass].SetLineColor(0)
	if (drawExpected):
		gest1[xmass].Draw("4 same")
				
	g[xmass].SetLineColor(j)
	g[xmass].SetFillColor(0)
	g[xmass].SetMarkerColor(j)
	g[xmass].SetMarkerStyle(r.kFullSquare)
	g[xmass].SetMarkerSize(1)
	g[xmass].SetLineWidth(2)
		
	# plotting
	g[xmass].Draw("lp same")

	leg.AddEntry(g[xmass], "m_{X} = " + str(int(xmass)) + " GeV/c^{2}")

	#g[xmass].Fit("pol3", "", "", lmin, lmax)
	#f1 = g[xmass].GetFunction("pol3")
	#f1 = r.TF1("f"+str(xmass), "1/([0]+[1]/x + [2]/(x*x) + [3]/(x*x*x))", lmin, lmax)
	#g[xmass].Fit(f1, "", "", lmin, lmax)
	#f1.SetLineColor(j)

	j+=1

leg.AddEntry(gest1[xmass], "Expected limit (#pm 1 #sigma)")
leg.SetTextFont(42)
leg.Draw()
c.Update()

suffix = ""
if (logx == 1):
	if (logy == 1):
		suffix = "LogXY"
	else:
		suffix = "LogX"
else:
	if (logy == 1):
		suffix = "LogY"

#f = r.TFile("limits"+FileInfix+"Lifetime"+LeptonType+str(hmass)+suffix+"Exp2.png", "RECREATE")
#for xmass in xmasses:
#	g[xmass].Write()
#	gest1[xmass].Write()
#f.Close()

c.Print("limits"+FileInfix+"Lifetime"+LeptonType+str(hmass)+suffix+"Exp2.png")
c.Print("limits"+FileInfix+"Lifetime"+LeptonType+str(hmass)+suffix+"Exp2.pdf")
#raw_input("Press ENTER to finish...")
