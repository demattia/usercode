#############################################
### MAKE REWEIGHTED EFFICIENCY TEXT FILES
#############################################

import os
import ROOT

def dump_reweighted_efficiencies(lepton_name,ptCut,limitfolder,eAnalysis,muAnalysis,workdirs_signal):

    # ctau factors for lifetime reweighting
    # here we use k=3 and k=1./3
    ctfact_list=[1./30.,1./20.,1./10.,1./3.,0.5,2.0,3.,10.,20.,30.]


    # corrected efficiencies, systematics
    leptrack_effi1=[]
    leptrack_effi2=[]
    leptrack_effi1_pileupuncertainty=[]
    leptrack_effi2_pileupuncertainty=[]
    leptrack_effi1_statisticaluncertainty=[]
    leptrack_effi2_statisticaluncertainty=[]
    leptrack_effi1_uncorrected=[]
    leptrack_effi2_uncorrected=[]
    leptrack_effi1_lifetime=[]
    leptrack_effi2_lifetime=[]
    leptrack_effi1_lifetime_statuncert=[]
    leptrack_effi2_lifetime_statuncert=[]
    for entry in ctfact_list:
        leptrack_effi1_lifetime.append([])
        leptrack_effi2_lifetime.append([])
        leptrack_effi1_lifetime_statuncert.append([])
        leptrack_effi2_lifetime_statuncert.append([])
        pass

    if lepton_name.find("mu")>=0:
        treename=muAnalysis+"/dileptons_signal_MC_L2DoubleMu30_NoVertex_v1/bigTree"
    else:
        treename=eAnalysis+"/dileptons_signal_MC_DoublePhoton38_v2/bigTree"
        pass

    for numDecays in [1,2]:
        for workdir in workdirs_signal:
            # extract ctau from DBS name
            ctau=0
            dbsname=os.popen("grep sampleDataSet "+workdir+"/*_cff.py").readlines()
            for line in dbsname:
                for field in line.split("_"):
                    if field.find("CTau-")==0:
                        ctau=int(field.split("-")[1])
                        pass
                    pass
                pass
            ctau/=10. # convert from Pythia units (mm) to cm
            # extract efficiency
            sampleID=workdir.split("/")[-1].replace("_analysis","")
            enumfile=ROOT.TFile.Open(workdir+"/histograms.root")
            enumtree=enumfile.Get(treename)
            # uncorrected, and central value with pile-up correction
            enumtree.Project("valhist","phi_cand",\
                             #"weight*(passesAllCuts"+\
                             "(0.5*(_weight_up+_weight_down))*(passesAllCuts"+\
                             " && _numDecays==%i && leptonPtL>%f)"%(numDecays,ptCut))
            enum_uncorrected=enumfile.Get("valhist").GetEntries()
            enum_central=enumfile.Get("valhist").Integral()
            try:
                enum_central_relerr=1./math.sqrt(enumfile.Get("valhist").GetEffectiveEntries())
            except:
                enum_central_relerr=0
                pass
            # pile-up variation up and down
            enumtree.Project("valhist","phi_cand",
                             "_weight_up*(passesAllCuts"+\
                             " && _numDecays==%i && leptonPtL>%f)"%(numDecays,ptCut))
            enum_up=enumfile.Get("valhist").Integral()
            enumtree.Project("valhist","phi_cand","_weight_down*(passesAllCuts"+\
                             " && _numDecays==%i && leptonPtL>%f)"%(numDecays,ptCut))
            enum_down=enumfile.Get("valhist").Integral()
            # lifetime variation up and down
            # weight (from Ian): (1/newCtau)*exp(-ctau/newCtau)/(1/oldCtau)/exp(-ctau/oldCtau)
            # if we reweight by a specific factor k, this gives
            # (1./k)*exp(-ctau/k/oldCtau)/exp(-ctau/oldCtau)
            # =(1./k)*exp(ctau*(1./oldCtau-1./k/oldCtau))
            # =(1./k)*exp(ctau*(1./k)*(k-1.)/oldCtau)
            enum_list=[]
            enum_relerr_list=[]
            for ctfact in ctfact_list:
                project="(0.5*(_weight_up+_weight_down)"+\
                         "*1./%f*exp(-_ctau1/%f/%f)"%(ctfact,ctfact,ctau)+\
                         "/exp(-_ctau1/%f)"%ctau+\
                         "*1./%f*exp(-_ctau2/%f/%f)"%(ctfact,ctfact,ctau)+\
                         "/exp(-_ctau2/%f))"%ctau+\
                         "*(passesAllCuts"+\
                         " && _numDecays==%i && leptonPtL>%f)"%(numDecays,ptCut)
                enumtree.Project("valhist","phi_cand",project)
                enum_list.append(enumfile.Get("valhist").Integral())
                try:
                    enum_relerr_list.append(1./math.sqrt(enumfile.Get("valhist").GetEffectiveEntries()))
                except:
                    enum_relerr_list.append(0.0)
                    pass
                pass
            enumfile.Close()
            
            # now open prefilter file to get normalization
            # (numSignal in proper decay channel)
            denomfile=ROOT.TFile.Open(workdir+"/prefilter.root")
            if treename.find("eTrack")>=0:
                denomhist=denomfile.Get("isoTrackPrefilter/numSignalE")
            elif treename.find("muTrack")>=0:
                denomhist=denomfile.Get("isoTrackPrefilter/numSignalMu")
                pass
            denom_uncorrected=denomhist.GetBinContent(numDecays+1)*numDecays
            
            # and get overall normalization correction due to lumi weights
            # from prefilter.root
            pileuphist=denomfile.Get("isoTrackPrefilter/pileup_3bx")
                
            pileup_weights_MC=[ 0.104109,
                                0.0703573,
                                0.0698445,
                                0.0698254,
                                0.0697054,
                                0.0697907,
                                0.0696751,
                                0.0694486,
                                0.0680332,
                                0.0651044,
                                0.0598036,
                                0.0527395,
                                0.0439513,
                                0.0352202,
                                0.0266714,
                                0.019411,
                                0.0133974,
                                0.00898536,
                                0.0057516,
                                0.00351493,
                                0.00212087,
                                0.00122891,
                                0.00070592,
                                0.000384744,
                                0.000219377 ]

            # Values used in the PAS
            # pileup_data = [ 0.019091,
            #                 0.0293974,
            #                 0.0667931,
            #                 0.108859,
            #                 0.139533,
            #                 0.149342,
            #                 0.138629,
            #                 0.114582,
            #                 0.0859364,
            #                 0.059324,
            #                 0.0381123,
            #                 0.0229881,
            #                 0.0131129,
            #                 0.00711764,
            #                 0.00369635,
            #                 0.00184543,
            #                 0.000889604,
            #                 0.000415683,
            #                 0.000188921,
            #                 0.000146288,
            #                 0.0,
            #                 0.0,
            #                 0.0,
            #                 0.0,
            #                 0.0 ]

            # Recomputed from data for the A1+A2 period on 07/01/2012
            pileup_data = [ 0.00706972,
                            0.0307666,
                            0.0711507,
                            0.115312,
                            0.146412,
                            0.154726,
                            0.141359,
                            0.114584,
                            0.0839703,
                            0.0564282,
                            0.0351613,
                            0.0204986,
                            0.011264,
                            0.0058706,
                            0.00291759,
                            0.00138913,
                            0.000636229,
                            0.000281323,
                            0.000120478,
                            5.01141e-05,
                            2.02982e-05,
                            8.02373e-06,
                            3.10155e-06,
                            1.17442e-06,
                            4.36284e-07,
                            1.59215e-07,
                            5.71409e-08,
                            2.01857e-08,
                            7.02382e-09,
                            2.40856e-09,
                            8.14197e-10,
                            2.71366e-10,
                            8.91705e-11,
                            2.88915e-11,
                            1.34134e-11,
                            6.96827e-12 ]

            # Computed on 07/01/2012 for all periods: A1+A2+A3+A4+B1
# pileup_data = [ 0.00285942,
#                             0.0125603,
#                             0.0299631,
#                             0.051313,
#                             0.0709713,
#                             0.0847864,
#                             0.0914627,
#                             0.0919255,
#                             0.0879994,
#                             0.0814127,
#                             0.0733995,
#                             0.0647191,
#                             0.0558327,
#                             0.0470663,
#                             0.0386988,
#                             0.0309811,
#                             0.0241175,
#                             0.018241,
#                             0.0133997,
#                             0.00956071,
#                             0.00662814,
#                             0.00446735,
#                             0.00292946,
#                             0.00187057,
#                             0.00116414,
#                             0.000706805,
#                             0.000419059,
#                             0.000242856,
#                             0.0001377,
#                             7.64582e-05,
#                             4.16101e-05,
#                             2.22135e-05,
#                             1.16416e-05,
#                             5.9937e-06,
#                             5.95542e-06,
#                             1.70121e-12 ]

            # basis of calculation for pile-up shifting weights
            p1_minus = [
                -0.677786,
                -0.619614,
                -0.49465,
                -0.357963,
                -0.238359,
                -0.110002,
                0.0348629,
                0.191263,
                0.347648,
                0.516615,
                0.679646,
                0.836673,
                0.97764,
                1.135,
                1.29922,
                1.42467,
                1.55901,
                1.61762,
                1.67275,
                1.96008
                ]
            p2_minus = [
                0.526164,
                0.251816,
                0.11049,
                0.026917,
                -0.0464692,
                -0.087022,
                -0.0931581,
                -0.0714295,
                -0.0331772,
                0.0347473,
                0.108658,
                0.193048,
                0.272314,
                0.376357,
                0.4964,
                0.58854,
                0.684959,
                0.731063,
                0.760044,
                1.02386
                ]
            p1_plus = [
                -0.739059,
                -0.594445,
                -0.477276,
                -0.359707,
                -0.233573,
                -0.103458,
                0.0373401,
                0.176571,
                0.337617,
                0.499074,
                0.675126,
                0.840522,
                1.00917,
                1.15847,
                1.23816,
                1.44271,
                1.52982,
                1.46385,
                1.5802,
                0.988689
                ]
            p2_plus = [
                0.208068,
                0.130033,
                0.0850356,
                0.0448344,
                0.000749832,
                -0.0331347,
                -0.0653281,
                -0.0746009,
                -0.0800667,
                -0.0527636,
                -0.00402649,
                0.103338,
                0.261261,
                0.491084,
                0.857966,
                1.19495,
                1.75071,
                2.65559,
                3.35433,
                5.48835
                ]
            shift=0.6
            pweight_up=[]
            pweight_down=[]
            for i in range(len(p1_plus)):
                pweight_up.append(1+p1_plus[i]*shift+p2_plus[i]*shift*shift)
                pweight_down.append(1-p1_minus[i]*shift+p2_minus[i]*shift*shift)
                pass
            for i in range(5):
                # not quite correct, but these bins are hardly populated anyway
                pweight_up.append(1.0)
                pweight_down.append(1.0)
                pass
            
            # note: denominator for lifetime reweighting is simply whatever
            # we get with central lumi weights. ctau weights should preserve
            # normalization.
            
            # now calculate the weighted number of events
            sum=0.0
            sumw=0.0
            sumw2=0.0
            sumw_up=0.0
            sumw_down=0.0
            for i in range(pileuphist.GetNbinsX()):
                if i<len(pileup_data) and i<len(pileup_weights_MC):
                    w=pileup_data[i]/pileup_weights_MC[i]
                    n=pileuphist.GetBinContent(i+1)
                    sumw+=w*n
                    sumw2+=w*w*n
                    sumw_up+=pweight_up[i]*w*n
                    sumw_down+=pweight_down[i]*w*n
                    sum+=n
                    pass
                pass
            denom_central=denom_uncorrected*sumw/sum
            denom_up=denom_uncorrected*sumw_up/sum
            denom_down=denom_uncorrected*sumw_down/sum
            effi_uncorrected="%-70s %10.4f"%(sampleID,enum_uncorrected/denom_uncorrected)
            effi=enum_central/denom_central
            effi_corrected="%-70s %10.4f"%(sampleID,enum_central/denom_central)
            effi_e1=(enum_up/denom_up)-(enum_central/denom_central)
            effi_e2=(enum_down/denom_down)-(enum_central/denom_central)
            effi_pileuperror="%-50s %10.4f %10.4f"%(sampleID,
                                                    abs(max(effi_e1,effi_e2)),
                                                    -abs(min(effi_e1,effi_e2)))
            effi_statisticalerror="%-50s %10.4f %10.4f"%(sampleID,
                                                         effi*enum_central_relerr,
                                                         -effi*enum_central_relerr)
            
            effi_lifetime=[]
            effi_lifetime_staterror=[]
            for i in range(len(enum_list)):
                enum=enum_list[i]
                enum_relerr=enum_relerr_list[i]
                effi=enum/denom_central
                effi_lifetime.append("%-70s %10.4f"%(sampleID,effi))
                effi_lifetime_staterror.append("%-50s %10.4f %10.4f"\
                                               %(sampleID,
                                                 effi*enum_relerr,
                                                 -effi*enum_relerr))
                pass
            
            if numDecays==1:
                leptrack_effi1.append(effi_corrected)
                leptrack_effi1_uncorrected.append(effi_uncorrected)
                for i in range(len(effi_lifetime)):
                    leptrack_effi1_lifetime[i].append(effi_lifetime[i])
                    leptrack_effi1_lifetime_statuncert[i].append(effi_lifetime_staterror[i])
                    pass
                leptrack_effi1_pileupuncertainty.append(effi_pileuperror)
                leptrack_effi1_statisticaluncertainty.append(effi_statisticalerror)
            elif numDecays==2:
                leptrack_effi2.append(effi_corrected)
                leptrack_effi2_uncorrected.append(effi_uncorrected)
                for i in range(len(effi_lifetime)):
                    leptrack_effi2_lifetime[i].append(effi_lifetime[i])
                    leptrack_effi2_lifetime_statuncert[i].append(effi_lifetime_staterror[i])
                    pass
                leptrack_effi2_pileupuncertainty.append(effi_pileuperror)
                leptrack_effi2_statisticaluncertainty.append(effi_statisticalerror)
                pass
            pass
        pass
    
    lep1file=open(limitfolder+"/efficiencies_di"+lepton_name+"1.txt","w")
    for entry in leptrack_effi1: lep1file.write("%s\n"%entry)
    lep1file.close()
    lep2file=open(limitfolder+"/efficiencies_di"+lepton_name+"2.txt","w")
    for entry in leptrack_effi2: lep2file.write("%s\n"%entry)
    lep2file.close()

    for i in range(len(ctfact_list)):
        filespec="ctaufact_%.3f"%ctfact_list[i]
        lep1file=open(limitfolder+"/efficiencies_di"+lepton_name+"1_"+filespec+".txt","w")
        for entry in leptrack_effi1_lifetime[i]: lep1file.write("%s\n"%entry)
        lep1file.close()
        lep2file=open(limitfolder+"/efficiencies_di"+lepton_name+"2_"+filespec+".txt","w")
        for entry in leptrack_effi2_lifetime[i]: lep2file.write("%s\n"%entry)
        lep2file.close()
        
        lep1file=open(limitfolder+"/efficiencies_di"+lepton_name+"1_"+filespec+"_statistical_uncertainty.txt","w")
        for entry in leptrack_effi1_lifetime_statuncert[i]: lep1file.write("%s\n"%entry)
        lep1file.close()
        lep2file=open(limitfolder+"/efficiencies_di"+lepton_name+"2_"+filespec+"_statistical_uncertainty.txt","w")
        for entry in leptrack_effi2_lifetime_statuncert[i]: lep2file.write("%s\n"%entry)
        lep2file.close()
        pass

    lep1file=open(limitfolder+"/efficiencies_di"+lepton_name+"1_uncorrected.txt","w")
    for entry in leptrack_effi1_uncorrected: lep1file.write("%s\n"%entry)
    lep1file.close()
    lep2file=open(limitfolder+"/efficiencies_di"+lepton_name+"2_uncorrected.txt","w")
    for entry in leptrack_effi2_uncorrected: lep2file.write("%s\n"%entry)
    lep2file.close()

    lep1file=open(limitfolder+"/efficiencies_di"+lepton_name+"1_pileup_uncertainty.txt","w")
    for entry in leptrack_effi1_pileupuncertainty: lep1file.write("%s\n"%entry)
    lep1file.close()
    lep2file=open(limitfolder+"/efficiencies_di"+lepton_name+"2_pileup_uncertainty.txt","w")
    for entry in leptrack_effi2_pileupuncertainty: lep2file.write("%s\n"%entry)
    lep2file.close()

    lep1file=open(limitfolder+"/efficiencies_di"+lepton_name+"1_statistical_uncertainty.txt","w")
    for entry in leptrack_effi1_statisticaluncertainty: lep1file.write("%s\n"%entry)
    lep1file.close()
    lep2file=open(limitfolder+"/efficiencies_di"+lepton_name+"2_statistical_uncertainty.txt","w")
    for entry in leptrack_effi2_statisticaluncertainty: lep2file.write("%s\n"%entry)
    lep2file.close()

    return
