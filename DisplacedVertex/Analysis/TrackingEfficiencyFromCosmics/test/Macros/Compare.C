#include "MuonAnalyzer.C"
#include <iostream>

void Compare()
{
  // MuonAnalyzer mu1("globalMuons_reco.root");
  // MuonAnalyzer mu2("globalMuons_Refit.root");
  MuonAnalyzer mu1("globalMuons_SegmentBased.root");
  // MuonAnalyzer mu2("globalMuons_RefitSegmentBased.root");
  MuonAnalyzer mu2("globalMuons_RefitOneIterationSegmentBased.root");


  Long64_t nentries1 = mu1.fChain->GetEntriesFast();
  Long64_t nentries2 = mu2.fChain->GetEntriesFast();

  if( nentries1 != nentries2 ) {
    std::cout << "Different entries for the two trees. Exiting." << std::endl;
    return;
  }

  // Filled for all tracks (there are the cases where multiple tracks enter multiple times giving a broadening...)
  TH1F * hDeltaR = new TH1F("deltaR", "#Delta R", 1000, 0, 10);
  TH1F * hDeltaPt = new TH1F("deltaPt", "#Delta P_{T}", 1000, -100, 100);

  // Filled for matched tracks
  TH1F * hDeltaNhits = new TH1F("deltaNhits", "difference in the number of hits", 20, -10, 10);
  TH1F * hDeltaChi2NDOF = new TH1F("deltaChi2NDOF", "difference in the #chi^{2}/ndof", 100, -10, 10);

  // Filled for unmatched tracks
  TH1F * hUnmatchedMu1Chi2NDOF = new TH1F("unmatchedMu1Chi2NDOF", "#chi^{2}/ndof for unmatched mu1", 100, 0, 10);
  TH1F * hUnmatchedMu2Chi2NDOF = new TH1F("unmatchedMu2Chi2NDOF", "#chi^{2}/ndof for unmatched mu2", 100, 0, 10);

  int nDifferent = 0;
  int nDifferentLess = 0;
  int nDifferentMore = 0;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries1; jentry++) {
    Long64_t ientry = mu1.LoadTree(jentry);
    mu2.LoadTree(jentry);
    if (ientry < 0) break;
    nb = mu1.fChain->GetEntry(jentry);   nbytes += nb;
    mu2.fChain->GetEntry(jentry);

    if( mu1.event != mu2.event ) {
      std::cout << "Event numbers do not match. Exiting" << std::endl;
      exit(1);
    }

    // if (Cut(ientry) < 0) continue;
    if( ientry > 0 ) {
      // std::cout << "number of tracks (tree1) = " << mu1.tracks_ << std::endl;
      // std::cout << "number of tracks (tree2) = " << mu2.tracks_ << std::endl;

      // std::cout << "chi2 (tree1) = " << mu1.tracks_normalizedChi2[0] << std::endl;
      // std::cout << "chi2 (tree2) = " << mu2.tracks_normalizedChi2[0] << std::endl;

      if( mu1.tracks_ != mu2.tracks_ ) {
	// std::cout << "number of tracks (tree1) = " << mu1.tracks_ << ", (tree2) = " << mu2.tracks_ << std::endl;
	// for( int i=0; i<mu1.tracks_; ++i ) {
	//   std::cout << "track["<<i<<"] (tree1) [pt, charge, chi2/ndof] = " << mu2.tracks_pt[i] << ", " << mu2.tracks_charge[i] << ", " << mu2.tracks_normalizedChi2[i] << std::endl;
	// }
	// for( int i=0; i<mu2.tracks_; ++i ) {
	//   std::cout << "track["<<i<<"] (tree2) [pt, charge, chi2/ndof] = " << mu2.tracks_pt[i] << ", " << mu2.tracks_charge[i] << ", " << mu2.tracks_normalizedChi2[i] << std::endl;
	// }
	if( mu1.tracks_ > mu2.tracks_ ) ++nDifferentLess;
	else if( mu1.tracks_ < mu2.tracks_ ) ++nDifferentMore;
	++nDifferent;
      }
    }

    // Find matching tracks
    std::vector<int> matchedMu2;
    for( int i=0; i<mu1.tracks_; ++i ) {
      bool found = false;
      for( int j=0; j<mu2.tracks_; ++j ) {

	if( mu1.tracks_charge[i] == mu2.tracks_charge[j] ) {
	  double deltaRValue = deltaR(mu1.tracks_eta[i], mu1.tracks_phi[i], mu2.tracks_eta[j], mu2.tracks_phi[j]);
	  double deltaPtValue = mu1.tracks_pt[i] - mu2.tracks_pt[j];
	  hDeltaR->Fill(deltaRValue);
	  hDeltaPt->Fill(deltaPtValue);

	  // 	if( mu1.tracks_charge[i] == mu2.tracks_charge[j] ) {
	  // 	  std::cout << "deltaPt = " << mu1.tracks_charge[i] - mu2.tracks_pt[j] << std::endl;
	  // 	  std::cout << "deltaR = " << deltaR(mu1.tracks_eta[i], mu1.tracks_phi[i], mu2.tracks_eta[j], mu2.tracks_phi[j]) << std::endl;
	  // 	}

	  if( (deltaPtValue < 1) &&
	      (deltaRValue < 0.1) &&
	      (std::find(matchedMu2.begin(), matchedMu2.end(), j) == matchedMu2.end()) ) {
	    matchedMu2.push_back(j);
	    found = true;
	    // Fill the matched tracks plots
	    hDeltaNhits->Fill(mu1.tracks_nHits[i] - mu2.tracks_nHits[j]);
	    hDeltaChi2NDOF->Fill(mu1.tracks_normalizedChi2[i] - mu2.tracks_normalizedChi2[j]);
	    break;
	  }
	}
      }
      if( !found ) {
	// Fill the unmatched tracks from the first tree
	hUnmatchedMu1Chi2NDOF->Fill(mu1.tracks_normalizedChi2[i]);
      }
    }
    // Fill the unmatched tracks from the second tree
    for( int j=0; j<mu2.tracks_; ++j) {
      if( std::find(matchedMu2.begin(), matchedMu2.end(), j) == matchedMu2.end() ) {
	hUnmatchedMu2Chi2NDOF->Fill(mu2.tracks_normalizedChi2[j]);
      }
    }
  }

  std::cout << "total entries = " << nentries1 << std::endl;
  std::cout << "entries with different number of tracks = " << nDifferent << std::endl;
  std::cout << "entries with less tracks in the new tree = " << nDifferentLess << std::endl;
  std::cout << "entries with more tracks in the new tree = " << nDifferentMore << std::endl;

  savePlot("deltaR", hDeltaR, "#Delta R", true);
  savePlot("deltaPt", hDeltaPt, "#Delta P_{T}", true);
  savePlot("deltaNhits", hDeltaNhits, "#Delta N(hits)", true);
  savePlot("deltaChi2NDOF", hDeltaChi2NDOF, "#Delta #chi^{2}/ndof", true);
  savePlot("unmatchedMu1Chi2NDOF", hUnmatchedMu1Chi2NDOF, "#chi^{2}/ndof");
  savePlot("unmatchedMu2Chi2NDOF", hUnmatchedMu2Chi2NDOF, "#chi^{2}/ndof");
}
