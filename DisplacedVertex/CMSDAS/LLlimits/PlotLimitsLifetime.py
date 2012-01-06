import ROOT as r
import sys,os

if len(sys.argv) < 3:
	sys.exit('Usage: python PlotLimitsLifetime.py LeptonType HiggsMass')

LeptonType = sys.argv[1]
hmass = int(sys.argv[2])

# Configurable parameters
# Note that the output of mainSigma.py is simply sigma,
# so we have to add the BR back in

br = 0.01
logx = 1
logy = 1

ymin = 0
if (hmass == 1000):
	xmasses = ["350", "150", "050", "020"]
	xmin = 0.2
	xmax = 120
	ymax = 0.08
	if (logy == 1):
		if (LeptonType == "Electrons"):
			ymin = 0.003
			ymax = 0.3
		else:
			ymin = 0.001
			ymax = 5.0
elif (hmass == 400):
	xmasses = ["150", "050", "020"]
	xmin = 0
	xmax = 130
	ymax = 0.08
	if (logy == 1):
		xmin = 0.8
		if (LeptonType == "Electrons"):
			ymin = 0.008
			ymax = 0.5
		else:
			ymin = 0.002
			ymax = 0.09
elif (hmass == 200):
	xmasses = ["050", "020"]
	xmin = 0
	xmax = 70
	ymax = 1
	if (logy == 1):
		xmin = 2
		if (LeptonType == "Electrons"):
			ymin = 0.03
		else:
			ymin = 0.006
			ymax = 0.1
else:
	print "Oops, mass not recognized"

g = {}
j = 1

r.gROOT.SetStyle("Plain")
r.gStyle.SetErrorX(0)
c = r.TCanvas("c","c",600,600)
r.gPad.SetLeftMargin(0.15)
fr = c.DrawFrame(xmin, ymin, xmax, ymax)
fr.SetTitle("Limit, "+LeptonType[:-1]+" channel, m_{H} = "+str(hmass)+" GeV/c^{2}")
fr.GetXaxis().SetTitle("c#tau [cm]")
fr.GetYaxis().SetTitle("Cross Section*BR [pb] (95% CL)")
fr.GetYaxis().SetTitleOffset(1.8)
if (logx == 1):
	r.gPad.SetLogx()
if (logy == 1):
	r.gPad.SetLogy()

leg = r.TLegend(0.4, 0.6, 0.88, 0.88)
leg.SetFillColor(0)
leg.SetBorderSize(0)

for xmass in xmasses:
	flist = []
	for filename in os.listdir(LeptonType):
		if filename.find('.txt')>-1 and filename.find(str(hmass))>-1 and filename.find(str(xmass))>-1:
			flist.append(filename[:-4])

	flist.sort()
	print "Files found for " + str(xmass) + ": " + str(flist)
	flist = [str(a)+'.txt' for a in flist]

	# graph for limits:
	g[xmass] = r.TGraph(len(flist))
	lmax=-1
	lmin=9999
	i = 0
	for filename in flist:
		lt=float(filename[:-4].split("_")[2])
		f = open(LeptonType+"/"+filename)
		# print "Lifetime in "+filename+" is " + str(lt)
		if lt > lmax:
			lmax = lt
		if lt < lmin:
			lmin = lt
		items = [eval(a) for a in f.readline().split()]
		m = items[0]
		low = items[1]
		up = items[2]
		g[xmass].SetPoint(i,lt,up*br)
		i+=1

	# Probably don't need this any more
	g[xmass].SetPoint(i,lt,up*br)
				
	g[xmass].SetLineColor(j)
	g[xmass].SetFillColor(0)
	g[xmass].SetMarkerColor(j)
	g[xmass].SetMarkerStyle(r.kFullSquare)
	g[xmass].SetMarkerSize(1)
		
	# plotting
	g[xmass].Draw("p same")

	leg.AddEntry(g[xmass], "m_{X} = " + str(int(xmass)) + " GeV/c^{2}")

	g[xmass].Fit("pol2", "", "", lmin, lmax)
	f1 = g[xmass].GetFunction("pol2")
	f1.SetLineColor(j)

	j+=1

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

c.Print("limitsLifetime"+LeptonType+str(hmass)+suffix+".png")
raw_input("Press ENTER to finish...")
