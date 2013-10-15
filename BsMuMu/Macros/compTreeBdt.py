from setdirs import *
from ROOT import TROOT, TLine, TGraphErrors, TLatex
from ROOT import TFile, TH1F, TH2F, TCanvas, THStack, TLegend
from ROOT import gBenchmark, gStyle, gROOT
#from array import array


## script knobs
region = "barrel"
#region = "endcaps"

unblind = True

compareMainVsXcheck = True


looseBdtCut = 0 #1:no cut, #2:zero

mainBdtCut = { "barrel":0.36, "endcaps":0.38 }
#mainBdtCut = { "barrel":0., "endcaps":0. }

if len(sys.argv) > 1: 
    if sys.argv[1] == "unblind": unblind = True
    else : unblind = False
if len(sys.argv) > 2:
    if sys.argv[2] == "compareMainVsXcheck": compareMainVsXcheck = True
    else : compareMainVsXcheck = False
if len(sys.argv) > 3:
    region = sys.argv[3]
if len(sys.argv) > 4:
    looseBdtCut = int(sys.argv[4])

figuresDird = "BsMuMuLatex/Figures/mainSel/"
plotsDir="plots/"
logsDir="logs/"

print "applying bdt>", mainBdtCut[region], "for", region
fn = "MainVsOnfile"
if compareMainVsXcheck:
    fn = "XcheckVsOnfile"

evts_A = [] #container for main (in-file) bdt output 
evts_B = [] #container for main (applied) bdt output OR xcheck output with main bdt applied

smallMainTree = False
#unblinded = True
if unblind:
    smallMainTree = True    


file_A = TFile("larger-SgData.root")
tree_A = file_A.Get("t")
if smallMainTree:
    file_A = TFile("small-test.root")
    tree_A = file_A.Get("SgData_bdt")
    #SgData_bdt->Draw("m","bdt>0.36&&m>4.9&&m<5.9&&m1eta<1.4&&m2eta<1.4&&m1bdt&&m2bdt&&hlt&&muid&&me<0.2")
hmass_A = TH1F("hm", 'mass', 40, 4.9, 5.9) 
hbdt_A = TH1F("hbdt", 'bdt', 50, 0, 1) 
#tree_A.Show()


if compareMainVsXcheck:
    strname_ = "_mainBDTonXcheckData"
    if   looseBdtCut == 1: strname_ += "-bdtnone"
    elif looseBdtCut == 2: strname_ += "-bdtzero"
    if unblind: strname_ += "_unblinded"
    strname_ += ".root"
    file_B = TFile("rootfiles/"+region+strname_)
    tree_B = file_B.Get("probe_tree")
    hmass_B = TH1F("hm", 'mass', 40, 4.9, 5.9) 
    hbdt_B = TH1F("hbdt", 'bdt', 50, 0, 1) 
else:
    file_B = TFile("rootfiles/"+region+"_mainBDTonMainData.root")
    tree_B = file_B.Get("events")
    hmass_B = TH1F("hmass", 'mass', 40, 4.9, 5.9) 
    hbdt_B = TH1F("hbdt", 'bdt', 50, 0, 1) 

#tree_B.Show()


for i,event in enumerate(tree_A):
    #continue
    #if i > 10: continue
    passed = abs(event.m1eta)<1.4 and abs(event.m2eta)<1.4
    if "barrel" not in region:
        passed = not passed and abs(event.m1eta < 2.0) and abs(event.m2eta<2.0)
    passed = passed and 4.9<event.m and event.m<5.9
    passed = passed and event.hlt and event.muid
    if not smallMainTree:
        passed = passed and event.hltm2 and event.pvw8>0.7  ##not defined in small
    passed = passed and event.bdt > mainBdtCut[region]
    if not passed: continue
    evt = event.evt
    mass = event.m
    bdt = event.bdt
    hmass_A.Fill( mass )
    hbdt_A.Fill( bdt )
    #print "A m:",mass," bdt:", bdt
    evts_A.append([evt,mass,bdt])


#print "A:", evts_A[1], evts_A[2] 
#print "B:", evts_B[1], evts_B[2] 
#exit(56)
signalRegion_B = []
bsig = open(logsDir+"bdt_signal_"+fn+"_"+region+".txt","w")


for event in tree_B:
    evt = 0
    passed = True
    if compareMainVsXcheck:
        evt = event.event
        passed = passed and abs(event.m1eta)<2.0 and abs(event.m2eta)<2.0
    else:
        evt = event.evt
        passed = passed and abs(event.m1eta)<2.0 and abs(event.m2eta)<2.0
    passed = passed and event.pvw8>0.7
    passed = event.m>4.9 and event.m<5.9
    passed = passed and event.bdt_v > mainBdtCut[region]
    if not passed: continue
    mass = event.m
    bdt = event.bdt_v
    hmass_B.Fill( mass )
    hbdt_B.Fill( bdt )
    #print "B m:",mass," bdt:", bdt
    evts_B.append([evt,mass,bdt])
    if mass>5.2 and mass<5.45:
        signalRegion_B.append([evt,mass,bdt])
        print >> bsig, \
              "\t[x-check unblind] run:%d, evt:%d, bdt:%5.3f, m:%5.3f, pt:%5.3f, m1pt:%5.3f, m2pt:%5.3f, m1eta:%5.3f, m2eta:%5.3f, fl3d:%5.3f, fls3d:%5.3f, chi2dof:%5.3f, pvip:%5.3f, pvips:%5.3f, maxdoca:%5.3f, alpha:%5.3f, docatrk:%5.3f, iso:%5.3f, pvw8:%5.3f \n" \
              % (event.run, event.event, event.bdt_v, event.m, event.pt, event.m1pt, event.m2pt, event.m1eta, event.m2eta, event.fl3d, event.fls3d, event.chi2dof,event.pvip, event.pvips,event.maxdoca,event.alpha, event.docatrk,event.iso, event.pvw8)

#exit(55)


          
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

#hmass.SetTitle(region)

yield_A = "%d" % hmass_A.GetEntries()
yield_B = "%d" % hmass_B.GetEntries()

xl = 0.35
if unblind: xl = 0.55 

canvasMass = TCanvas()
hmass_B.Draw()
hmass_A.Draw("same")
#hmass_B.Draw("same")
#leg.Draw("same")
tl.SetTextColor(1)
tl.DrawLatex(xl,0.72,region.upper())
tl.SetTextColor(2)
tl.DrawLatex(xl,0.65,"BDT onFile ["+yield_A+"]")
tl.SetTextColor(3)
tl.DrawLatex(xl,0.6,"BDT applied ["+yield_B+"]")
canvasMass.SaveAs(plotsDir + "mass_compare_"+fn+"_"+region+".gif")
canvasMass.SaveAs(figuresDird + "mass_compare_"+fn+"_"+region+".pdf")

canvasBdt = TCanvas()
hbdt_B.Draw()
hbdt_A.Draw("same")
#hbdt_B.Draw("same")
tl.SetTextColor(1)
tl.DrawLatex(0.6,0.72,region.upper())
tl.SetTextColor(2)
tl.DrawLatex(0.6,0.65,"BDT onFile ["+yield_A+"]")
tl.SetTextColor(3)
tl.DrawLatex(0.6,0.6,"BDT applied ["+yield_B+"]")
canvasBdt.SaveAs(plotsDir + "bdt_compare_"+fn+"_"+region+".gif")
canvasBdt.SaveAs(figuresDird + "bdt_compare_"+fn+"_"+region+".pdf")

#print evts_A, evts_B

### Do per-event comaprison analysis

notfound_A = []
mismatch_A = []
matched_A  = []

for a in evts_A:
    found = False
    for b in evts_B:
        if a[0]==b[0]:
            found = True
            if abs(a[1]-b[1])>0.0001 or abs(a[2]-b[2])>0.0001:
                print "mismatch:", a,"<->",b
                mismatch_A.append(a)
            else:
                print "A/matched!:", a,"<->",b
                matched_A.append(a)
                continue
    if not found:
        notfound_A.append(a)


ff = open(logsDir+"bdt_compare_"+fn+"_"+region+".txt","w")
print >> ff, "notfound_A(",len(notfound_A),"):", sorted(notfound_A, key=lambda up: up[1]),"\n"
print >> ff, "mismatch_A(",len(mismatch_A),"):", sorted(mismatch_A, key=lambda up: up[1]),"\n"
print >> ff, "matched_A(",len(matched_A),"):", sorted(matched_A, key=lambda up: up[1]),"\n"

notfound_B = []
mismatch_B = []
matched_B  = []

for a in evts_B:
    found = False
    for b in evts_A:
        if a[0]==b[0]:
            found = True
            if abs(a[1]-b[1])>0.00001 or abs(a[2]-b[2])>0.00001:
                print "mismatch:", a,"<->",b
                mismatch_B.append(a)
            else:
                print "B/matched!:", a,"<->",b
                matched_B.append(a)
                continue
    if not found:
        notfound_B.append(a)

print >> ff, "\n"
print >> ff, "notfound_B(",len(notfound_B),"):", sorted(notfound_B, key=lambda up: up[1]),"\n"
print >> ff, "mismatch_B(",len(mismatch_B),"):", sorted(mismatch_B, key=lambda up: up[1]),"\n"
print >> ff, "matched_B(",len(matched_B),"):", sorted(matched_B, key=lambda up: up[1]),"\n"

print >> ff, "\n\n"

print >> ff, "Only in ",file_A.GetName(),"("+str(len(notfound_A))+" candidates):\n"
for n in notfound_A:
    for e in tree_A:
        if e.evt != n[0]: continue
        if not smallMainTree:
            print >> ff, \
                  "\t [A] run:%d, evt:%d, bdt:%5.3f, m:%5.3f, pt:%5.3f, m1pt:%5.3f, m2pt:%5.3f, m1eta:%5.3f, m2eta:%5.3f, fl3d:%5.3f, fls3d:%5.3f, chi2dof:%5.3f, pvip:%5.3f, pvips:%5.3f, maxdoca:%5.3f, alpha:%5.3f, docatrk:%5.3f, iso:%5.3f, pvw8:%5.3f, pvw8Refit:%5.3f \n" \
                  % (e.run, e.evt, e.bdt, e.m, e.pt, e.m1pt, e.m2pt, e.m1eta, e.m2eta, e.fl3d, e.fls3d, e.chi2dof,e.pvip, e.pvips,e.maxdoca,e.alpha, e.docatrk,e.iso, e.pvw8, e.pvw8)
        else:
            print >> ff, \
                  "\t [A] run:%d, evt:%d, bdt:%5.3f, m:%5.3f, me:%5.3f pt:%5.3f, m1bdt:%5.3f, m2bdt:%5.3f, m1eta:%5.3f, m2eta:%5.3f \n" \
                  % (e.run, e.evt, e.bdt, e.m, e.me, e.pt, e.m1bdt, e.m2bdt, e.m1eta, e.m2eta)
### e.closetrk, closetrk:%8.7f,  -- check tree!!
print >> ff, "Only in ",file_B.GetName(),"("+str(len(notfound_B))+" candidates):\n"
for n in notfound_B:
    for e in tree_B:
        evt=0
        if not compareMainVsXcheck:
            if e.evt != n[0]: continue
            print >> ff, \
                  "\t [B] run:%d, evt:%d, bdt:%5.3f, m:%5.3f, pt:%5.3f, m1pt:%5.3f, m2pt:%5.3f, m1eta:%5.3f, m2eta:%5.3f, fl3d:%5.3f, fls3d:%5.3f, chi2dof:%5.3f, pvip:%5.3f, pvips:%5.3f, maxdoca:%5.3f, alpha:%5.3f, docatrk:%5.3f, iso:%5.3f, pvw8:%5.3f, pvw8Refit:%5.3f \n" \
                  % (e.run, e.evt, e.bdt_v, e.m, e.pt, e.m1pt, e.m2pt, e.m1eta, e.m2eta, e.fl3d, e.fls3d, e.chi2dof,e.pvip, e.pvips,e.maxdoca,e.alpha, e.docatrk,e.iso, e.pvw8, e.pvw8)
#            print >> ff, "evt:",e.evt, "bdt:",e.bdt_v, "m:", e.m, "pt:",e.pt, "m1pt:", e.m1pt, "m2pt:",e.m2pt, "m1eta:",e.m1eta, "m2eta:",e.m2eta,"fl3d:", e.fl3d, "fls3d:", e.fls3d, "chi2dof:",e.chi2dof, "pvip:",e.pvip, "pvips:",e.pvips, "maxdoca:",e.maxdoca, "alpha:",e.alpha,"closetrk:",e.closetrk, "docatrk:", e.docatrk, "iso:",e.iso, "\n"
        else:
            if e.event != n[0]: continue
            print >> ff, \
                  "\t [B] run:%d, evt:%d, bdt:%5.3f, m:%5.3f, pt:%5.3f, m1pt:%5.3f, m2pt:%5.3f, m1eta:%5.3f, m2eta:%5.3f, fl3d:%5.3f, fls3d:%5.3f, chi2dof:%5.3f, pvip:%5.3f, pvips:%5.3f, maxdoca:%5.3f, alpha:%5.3f, docatrk:%5.3f, iso:%5.3f, pvw8:%5.3f, pvw8Refit:%5.3f \n" \
                  % (e.run, e.event, e.bdt_v, e.m, e.pt, e.m1pt, e.m2pt, e.m1eta, e.m2eta, e.fl3d, e.fls3d, e.chi2dof,e.pvip, e.pvips,e.maxdoca,e.alpha, e.docatrk,e.iso, e.pvw8, e.pvw8Refit)
            #            print >> ff, "evt:",e.event, "bdt:",e.bdt_v, "m:", e.m, "pt:",e.pt, "m1pt:", e.m1pt, "m2pt:",e.m2pt, "m1eta:",e.m1eta, "m2eta:",e.m2eta,"fl3d:", e.fl3d, "fls3d:", e.fls3d, "chi2dof:",e.chi2dof, "pvip:",e.pvip, "pvips:",e.pvips, "maxdoca:",e.maxdoca, "alpha:",e.alpha,"closetrk:",e.closetrk, "docatrk:", e.docatrk, "iso:",e.iso, "\n"

print >> ff, "Mismatched candidates ("+str(len(mismatch_A+mismatch_B))+"):\n"
for n in mismatch_A:
    evt = n[0]
    print >> ff, "\tEvent: %d" % evt
    for e in tree_A:
        evta=e.evt
        if evt!=evta: continue
        if not smallMainTree:
            print >> ff, \
                  "\t[A] run:%d, evt:%d, bdt:%5.3f, m:%5.3f, pt:%5.3f, m1pt:%5.3f, m2pt:%5.3f, m1eta:%5.3f, m2eta:%5.3f, fl3d:%5.3f, fls3d:%5.3f, chi2dof:%5.3f, pvip:%5.3f, pvips:%5.3f, maxdoca:%5.3f, alpha:%5.3f, docatrk:%5.3f, iso:%5.3f, pvw8:%5.3f, pvw8Refit:%5.3f \n" \
                  % (e.run, e.evt, e.bdt, e.m, e.pt, e.m1pt, e.m2pt, e.m1eta, e.m2eta, e.fl3d, e.fls3d, e.chi2dof,e.pvip, e.pvips,e.maxdoca,e.alpha, e.docatrk,e.iso, e.pvw8, e.pvw8)
        else:
            print >> ff, \
                  "\t [A] run:%d, evt:%d, bdt:%5.3f, m:%5.3f, me:%5.3f pt:%5.3f, m1bdt:%5.3f, m2bdt:%5.3f, m1eta:%5.3f, m2eta:%5.3f \n" \
                  % (e.run, e.evt, e.bdt, e.m, e.me, e.pt, e.m1bdt, e.m2bdt, e.m1eta, e.m2eta)
            
    for e in tree_B:
        evtb=0
        if not compareMainVsXcheck:
            evtb=e.evt
        else:
            evtb=e.event
        if evt!=evtb: continue
        if not compareMainVsXcheck:
            print >> ff, \
                  "\t[B] run:%d, evt:%d, bdt:%5.3f, m:%5.3f, pt:%5.3f, m1pt:%5.3f, m2pt:%5.3f, m1eta:%5.3f, m2eta:%5.3f, fl3d:%5.3f, fls3d:%5.3f, chi2dof:%5.3f, pvip:%5.3f, pvips:%5.3f, maxdoca:%5.3f, alpha:%5.3f, docatrk:%5.3f, iso:%5.3f, pvw8:%5.3f, pvw8Refit:%5.3f \n" \
                  % (e.run, e.evt, e.bdt_v, e.m, e.pt, e.m1pt, e.m2pt, e.m1eta, e.m2eta, e.fl3d, e.fls3d, e.chi2dof,e.pvip, e.pvips,e.maxdoca,e.alpha, e.docatrk,e.iso, e.pvw8, e.pvw8Refit)
        else:
            print >> ff, \
                  "\t[B] run:%d, evt:%d, bdt:%5.3f, m:%5.3f, pt:%5.3f, m1pt:%5.3f, m2pt:%5.3f, m1eta:%5.3f, m2eta:%5.3f, fl3d:%5.3f, fls3d:%5.3f, chi2dof:%5.3f, pvip:%5.3f, pvips:%5.3f, maxdoca:%5.3f, alpha:%5.3f, docatrk:%5.3f, iso:%5.3f, pvw8:%5.3f, pvw8Refit:%5.3f \n" \
                  % (e.run, e.event, e.bdt_v, e.m, e.pt, e.m1pt, e.m2pt, e.m1eta, e.m2eta, e.fl3d, e.fls3d, e.chi2dof,e.pvip, e.pvips,e.maxdoca,e.alpha, e.docatrk,e.iso, e.pvw8, e.pvw8)


  
