#from setdirs import *
from ROOT import TROOT, TLine, TGraphErrors, TLatex
from ROOT import TFile, TH1F, TH2F, TCanvas, THStack, TLegend
from ROOT import gBenchmark, gStyle, gROOT
#from array import array

plotsDir="plots/"

#arr0 = [ ["a","a1"], ["c"] ]
#arr1 = [ ["a","b1"], ["a"] ]
#
#for a in arr0:
#    for b in arr1:
#        if a[0]==b[0]:
#            print a,"-", b

#arr.append(["d","g"])
#print arr[0][1]
#exit(55)

evts_A = []
evts_B = []

mainBdtCut = { "barrel":0.360, "endcaps":0.368 }

    
file_A = TFile("larger-SgData.root")
tree_A = file_A.Get("t")
hmass_A = TH1F("hm", 'mass', 40, 4.9, 5.9) 
hbdt_A = TH1F("hbdt", 'bdt', 50, 0, 1) 
#tree_A.Show()

file_B = TFile("rootfiles/barrel_mainBDTonMainData.root")
tree_B = file_B.Get("events")
#tree_B.Show()
#exit(55)


hmass_B = TH1F("hmass", 'mass', 40, 4.9, 5.9) 
hbdt_B = TH1F("hbdt", 'bdt', 50, 0, 1) 


for i,event in enumerate(tree_A):
    #if i > 10: continue
    passed = 4.9<event.m and event.m<5.9 and abs(event.m1eta)<1.4 and abs(event.m2eta)<1.4 and event.hlt and event.hltm2 and event.bdt>0.36 and event.muid
    if passed:
        evt = event.evt
        mass = event.m
        bdt = event.bdt
        hmass_A.Fill( mass )
        hbdt_A.Fill( bdt )
        print "A m:",mass," bdt:", bdt
        evts_A.append([evt,mass,bdt])


#print "A:", evts_A[1], evts_A[2] 
#print "B:", evts_B[1], evts_B[2] 
#exit(56)

for event in tree_B:
    passed = True
    if passed:
        evt = event.evt
        mass = event.m
        bdt = event.bdt_v
        hmass_B.Fill( mass )
        hbdt_B.Fill( bdt )
        print "B m:",mass," bdt:", bdt
        evts_B.append([evt,mass,bdt])

#print evts_A, evts_B

notfound_A = []
mismatch_A = []
matched_A  = []

for a in evts_A:
    found = False
    for b in evts_B:
        if a[0]==b[0]:
            found = True
            if abs(a[1]-b[1])>0.00001 or abs(a[2]-b[2])>0.00001:
                print "mismatch:", a,"<->",b
                mismatch_A.append(a)
            else:
                print "matched!:", a,"<->",b
                matched_A.append(a)
                continue
    if not found:
        notfound_A.append(a)


print "notfound_A(",len(notfound_A),"):", sorted(notfound_A, key=lambda up: up[1])
print "mismatch_A(",len(mismatch_A),"):", sorted(mismatch_A, key=lambda up: up[1])
print "matched_A(",len(matched_A),"):", sorted(matched_A, key=lambda up: up[1])

#exit(44)
gStyle.SetOptStat(0)
gStyle.SetLegendBorderSize(0); 
gStyle.SetFillColor(0);
gStyle.SetPadBorderSize(0);

hmass_A.SetLineColor(2)
hmass_B.SetLineColor(3)
hbdt_A.SetLineColor(2)
hbdt_B.SetLineColor(3)
#tl = TLatex()
#tl.SetNDC();
#tl.SetTextSize( 0.033 );
#tl.DrawLatex(0.38,yy,txt)
leg = TLegend(0.35,0.6,0.5,0.75,"","brNDC")
leg.AddEntry(hmass_A, "Main", "f")
leg.AddEntry(hmass_B, "Xcheck", "f")
#leg.SetBorderSize(0)
#leg.SetShadowColor(0)
tl = TLatex()
tl.SetNDC();
tl.SetTextSize( 0.033 );

#mass.SetTitle(region)
            
canvasMass = TCanvas()
hmass_A.Draw()
hmass_B.Draw("same")
hmass_A.Draw("same")
#leg.Draw("same")
tl.SetTextColor(2)
tl.DrawLatex(0.35,0.65,"BDT from file ["+str(hmass_A.GetEntries())+"]")
tl.SetTextColor(3)
tl.DrawLatex(0.35,0.6,"BDT run ["+str(hmass_B.GetEntries())+"]")
canvasMass.SaveAs(plotsDir + "ab_mass.gif")


canvasBdt = TCanvas()
hbdt_B.Draw()
hbdt_A.Draw("same")
hbdt_B.Draw("same")
tl.SetTextColor(2)
tl.DrawLatex(0.6,0.65,"BDT from file ["+str(hbdt_A.GetEntries())+"]")
tl.SetTextColor(3)
tl.DrawLatex(0.6,0.6,"BDT run ["+str(hbdt_B.GetEntries())+"]")
canvasBdt.SaveAs(plotsDir + "ab_bdt.gif")
