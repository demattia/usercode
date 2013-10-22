import math
import ROOT
from ROOT import RooArgSet
import bisect

muMass = 0.1057


class Properties:
    """Stores the min and max mass and the binning in pt"""
    def __init__(self, minMass, maxMass, ptBins, triggerMatchDeltaR):
        self.minMass = minMass
        self.maxMass = maxMass
        self.ptBins = ptBins
        self.triggerMatchDeltaR = triggerMatchDeltaR


def deltaPhi(phi1, phi2):
    result = phi1 - phi2;
    if result > math.pi: result -= 2*math.pi
    if result <= -math.pi: result += 2*math.pi
    return result

def deltaR(v1, v2):
    return math.sqrt(deltaPhi(v1.phi(), v2.phi())**2 + (v1.eta()-v2.eta())**2)

def deltaR(phi1, eta1, phi2, eta2):
    return math.sqrt(deltaPhi(phi1, phi2)**2 + (eta1-eta2)**2)

def passSelection(track):
    if ( track.pt > 26 and abs(track.eta) < 2 and track.isolation/track.pt < 0.1 and track.trackerLayersWithMeasurement >= 6
         and track.trackQuality and track.dxy < 30. and track.dz < 30. ):
        return True

def passDileptonSelection(track1, track2, cosine):
    if track1.charge == track2.charge: return False
    if deltaR(track1.phi, track1.eta, track2.phi, track2.eta) < 0.2: return False
    if cosine <= -0.79: return False
    return True

def computeCosineAndMass(mu1, mu2):
    muon1 = ROOT.TLorentzVector()
    muon2 = ROOT.TLorentzVector()
    muon1.SetPtEtaPhiM(mu1.pt, mu1.eta, mu1.phi, muMass)
    muon2.SetPtEtaPhiM(mu2.pt, mu2.eta, mu2.phi, muMass)
    cosine = math.cos(muon1.Angle(muon2.Vect()))
    mass = (muon1+muon2).M()
    return (cosine, mass)

def fillTriggerMatchedTrack(track, triggerObjects, matchedTracks, p):
    for triggerMuon in triggerObjects:
        if deltaR(triggerMuon.phi, triggerMuon.eta, track.phi, track.eta) < p.triggerMatchDeltaR and passSelection(track):
            matchedTracks.append(track)

def fillSingleCandidate(mass, p, track1, track2, histoMap, datasetMap):
    cosineAndMass = computeCosineAndMass(track1, track2)
    if not passDileptonSelection(track1, track2, cosineAndMass[0]): return False
    if cosineAndMass[1] < p.minMass or cosineAndMass[1] > p.maxMass: return False
    mass.setVal(cosineAndMass[1])
    bins = find_bins(track1.pt, track2.pt, p.ptBins)
    if bins[0] == -1 or bins[1] == -1: return False
    histoMap[bins].Fill(mass.getVal())
    datasetMap[bins].add(RooArgSet(mass))
    return True

def fillCandidates(mass, properties, matchedTracks, histoMap, datasetMap):
    if len(matchedTracks) == 2:
        fillSingleCandidate(mass, properties, matchedTracks[0], matchedTracks[1], histoMap, datasetMap)
    elif len(matchedTracks) > 2:
        for i in range(0, len(matchedTracks)):
            for j in range(i+1, len(matchedTracks)):
                fillSingleCandidate(mass, properties, matchedTracks[i], matchedTracks[j], histoMap, datasetMap)

def find_bins(pt1, pt2, ptBins):
    return (bisect.bisect_right(ptBins, pt1)-1, bisect.bisect_right(ptBins, pt2)-1)

def find_position(bin1, bin2, ptBins):
    return bin1+len(ptBins)*bin2

def buildName(baseName, ptBin1, ptBin2, ptBins):
    return baseName+str(ptBins[ptBin1])+"_"+str(ptBins[ptBin2])

def buildNamePars(baseName, ptBin1, ptBin12, ptBin2, ptBin22, ptBins):
    return baseName+str(ptBins[ptBin1])+"_"+str(ptBins[ptBin12])+"_"+str(ptBins[ptBin2])+"_"+str(ptBins[ptBin22])

def progressCounter(progress):
    output = ""
    print "\r"
    output += "["
    for i in range (0, 10):
        if i < progress/10: output += "====="
        else: output += "-----"
    output += "] " + str(progress) + "%"
    print output
