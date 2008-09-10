/**
 * This small c++ program can be compiled directly with the g++.
 * It merges the tag-matrices produced by the QCDbTagProbability
 * analyzer. It first scales the counts in each matrix bin
 * multiplying it by the cross section of the corresponding qcd
 * bin and dividing it by the number of events in the sample (taken
 * recorded in the txt file). It then sums the different matrix bins
 * for the various qcd bins and normalizes the b-tag and nob-tag
 * so that they can be considered a probability.
 */

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

using namespace std;

void fillProbabilityMatrices(const string & probabilityFileName, unsigned int * binNum, double * binSize, unsigned int ***& trueArray, unsigned int ***& falseArray, int * eventsNumber ) {

  ifstream probabilityFile(probabilityFileName.c_str());
  if ( !probabilityFile.is_open() ) {
    cout << "file: " << probabilityFileName << "not found or unable to open" << endl;
    exit(1);
  }

  string line;
  // Read the first line with the bin numbers and sizes
  getline(probabilityFile, line);
  stringstream probabilityCounts(line);
  // This loop is based on the text file structure
  string skip;
  int wordCount=0;
  int binNumCount=0;
  int binSizeCount=0;
  while ( !probabilityCounts.eof() ) {
    ++wordCount;
    // Take every three words
    if ( wordCount%21 != 0 ) {
      if ( wordCount%3 == 0 ) {
        // Every six words there is the size of the bins, in the other case is the number of bins
        if ( wordCount%6 != 0 ) {
          probabilityCounts >> binNum[binNumCount];
          //cout << "binNum_["<<binNumCount<<"] = " << binNum[binNumCount] << endl;
          ++binNumCount;
        }
        else {
          probabilityCounts >> binSize[binSizeCount];
          //cout << "binSize_["<<binSizeCount<<"] = " << binSize[binSizeCount] << endl;
          ++binSizeCount;
        }
      }
      // Skip the rest of the words
      else probabilityCounts >> skip;
    }
    else {
      probabilityCounts >> *eventsNumber;
      cout << "*eventsNumber = " << *eventsNumber << endl;
    }
  }

  // Create and fill the arrays
  trueArray = new unsigned int**[binNum[0]];
  falseArray = new unsigned int**[binNum[0]];
  for(unsigned int i=0; i < binNum[0]; ++i) {
    trueArray[i] = new unsigned int*[binNum[1]];
    falseArray[i] = new unsigned int*[binNum[1]];
    for(unsigned int j=0; j < binNum[1]; ++j) {
      trueArray[i][j] = new unsigned int[binNum[2]];
      falseArray[i][j] = new unsigned int[binNum[2]];
      for (unsigned int k=0; k < binNum[2]; ++k) {
        if ( probabilityFile.eof() ) { 
          cout << "ERROR: not enough lines in the file for the required bins" << endl;
          exit(1);
        }

        getline(probabilityFile, line);
        stringstream probabilityCounts(line);
        for ( int w=0; w<6; ++w ) {
          // The second word is for the true case probability
          if( w == 2 ) {
            probabilityCounts >> trueArray[i][j][k];
            // cout << "trueArray["<<i<<"]["<<j<<"]["<<k<<"] = " << trueArray[i][j][k];
            // The fourth word is for the false case probability
          }
          if( w == 4 ) {
            probabilityCounts >> falseArray[i][j][k];
            // cout << " falseArray["<<i<<"]["<<j<<"]["<<k<<"] = " << falseArray[i][j][k] << endl;;
          }
          // Skip the other words
          else probabilityCounts >> skip;
        }
      }
    }
  }
}

int main() {

//   cout << "Merging all the files" << endl;

  unsigned int qcdBins = 8;

  double *** mergedTaggedJet;
  double *** mergedNotTaggedJet;
  double *** totNorm;
  unsigned int totQcdBinNum[3];
  double totQcdBinSize[3];
  int eventsNumber[] = { 0,0,0,0,
                         0,0,0,0 };

  const int crossSectionNum = 8;
  double crossSection[crossSectionNum] = { 163000000.,
                            21600000.,
                            3080000.,
                            494000.,
                            101000.,
                            24500.,
                            6240.,
                            2860. };
  string qcdNames[] = { "30-50", "50-80", "80-120",
                        "120-170", "170-230", "230-300",
                        "300-380", "380-incl" };

  double totCrossSection = 0.;
  for( int i=0; i<crossSectionNum; ++i ) totCrossSection += crossSection[i];

  for(unsigned int bin=0; bin < qcdBins; ++bin) {

    unsigned int *** taggedJet;
    unsigned int *** notTaggedJet;
    unsigned int qcdBinNum[3];
    double qcdBinSize[3];

    string fileName = "QCDbTagProbability_QCD_";
    fileName += qcdNames[bin];
    fileName += ".txt";

    fillProbabilityMatrices( fileName, qcdBinNum, qcdBinSize, taggedJet, notTaggedJet, &eventsNumber[bin] );
    // fillProbabilityMatrices( "HiggsPairProbability.txt", qcdBinNum, qcdBinSize, taggedJet, notTaggedJet );

    // Only the first time create also the merged arrays and save the qcdBinNum and qcdBinSize
    if ( bin == 0 ) {

      for( int i=0; i<3; ++i ) {
        totQcdBinNum[i] = qcdBinNum[i];
        totQcdBinSize[i] = qcdBinNum[i];
      }
      mergedTaggedJet = new double**[qcdBinNum[0]];
      mergedNotTaggedJet = new double**[qcdBinNum[0]];
      totNorm = new double**[qcdBinNum[0]];
      for(unsigned int i=0; i < qcdBinNum[0]; ++i) {
        mergedTaggedJet[i] = new double*[qcdBinNum[1]];
        mergedNotTaggedJet[i] = new double*[qcdBinNum[1]];
        totNorm[i] = new double*[qcdBinNum[1]];
        for(unsigned int j=0; j < qcdBinNum[1]; ++j) {
          mergedTaggedJet[i][j] = new double[qcdBinNum[2]];
          mergedNotTaggedJet[i][j] = new double[qcdBinNum[2]];
          totNorm[i][j] = new double[qcdBinNum[2]];
          for (unsigned int k=0; k < qcdBinNum[2]; ++k) {
            mergedTaggedJet[i][j][k] = 0;
            mergedNotTaggedJet[i][j][k] = 0;
            totNorm[i][j][k] = 0;
          }
        }
      }
    }

    double taggedCount = 0.;
    double notTaggedCount = 0.;
//     double taggedProb = 0.;
//     double notTaggedProb = 0.;
    double norm = 0.;
    // Fill the merged arrays with the weighted mean by cross sections of the probability in each bin
    // If the count == 1, do consider it = 0 (it is initialized to 1)
    // Use the probability as it is normalized. Otherwise we could use the counts only if rescaled by
    // the total number of events (and cross section).
    for(unsigned int i=0; i < qcdBinNum[0]; ++i) {
      for(unsigned int j=0; j < qcdBinNum[1]; ++j) {
        for (unsigned int k=0; k < qcdBinNum[2]; ++k) {

          //           if( taggedJet[i][j][k] != 1 ) taggedCount = taggedJet[i][j][k];
          //           else taggedCount = 0.;
          //           if( notTaggedJet[i][j][k] != 1 ) notTaggedCount = notTaggedJet[i][j][k];
          //           else taggedCount = 0.;

//           cout << "taggedJet["<<i<<"]["<<j<<"]["<<k<<"] = " << taggedJet[i][j][k] << endl;
//           cout << "notTaggedJet["<<i<<"]["<<j<<"]["<<k<<"] = " << notTaggedJet[i][j][k] << endl;

          taggedCount = taggedJet[i][j][k]-1;
          notTaggedCount = notTaggedJet[i][j][k]-1;

          // weigh also with the error

          if ( taggedCount != 0. || notTaggedCount != 0 ) norm = taggedCount + notTaggedCount;
          else norm = 1.;
//           taggedProb = taggedCount/norm;
//           notTaggedProb = notTaggedCount/norm;
//           cout << "taggedCount = " << taggedCount << ", taggedProb = " << taggedProb << endl;
//           cout << "notTaggedCount = " << notTaggedCount << ", notTaggedProb = " << notTaggedProb << endl;

//           double normFactor = taggedProb*notTaggedProb;
//           if ( normFactor == 0 ) normFactor = 1.;
//           if ( taggedProb != 0 ) mergedTaggedJet[i][j][k] += taggedProb*crossSection[bin]/(normFactor/norm);
//           if ( notTaggedProb != 0 ) mergedNotTaggedJet[i][j][k] += notTaggedProb*crossSection[bin]/(normFactor/norm);

          if ( eventsNumber[bin] != 0 ) {
            mergedTaggedJet[i][j][k] += taggedCount*crossSection[bin]/eventsNumber[bin];
            mergedNotTaggedJet[i][j][k] += notTaggedCount*crossSection[bin]/eventsNumber[bin];
            // totNorm[i][j][k] += crossSection[bin]/eventsNumber[bin];
          }
          //totNorm[i][j][k] += crossSection[bin]/(normFactor/norm);
        }
      }
    }


    // delete the multidimensional arrays
    for(unsigned int i=0; i != qcdBinNum[0]; ++i) {
      for(unsigned int j=0; j != qcdBinNum[1]; ++j) {
        delete[] taggedJet[i][j];
        delete[] notTaggedJet[i][j];
      }
      delete[] taggedJet[i];
      delete[] notTaggedJet[i];
    }

  } // end loop on bins

  double totTaggedCount = 0.;
  double totNotTaggedCount = 0.;
  double normTot = 0.;

  // Loop one last time on the merged arrays and put to 0.5 all the 0 bins (needed to avoid problems when computing the ratio)
  for(unsigned int i=0; i < totQcdBinNum[0]; ++i) {
    for(unsigned int j=0; j < totQcdBinNum[1]; ++j) {
      for (unsigned int k=0; k < totQcdBinNum[2]; ++k) {
        // If only one of the two
        // if ( !( mergedTaggedJet[i][j][k] == 0 && mergedNotTaggedJet[i][j][k] != 0 ) || !( mergedTaggedJet[i][j][k] != 0 &&  mergedNotTaggedJet[i][j][k] == 0 ) ) cout << "Error: should be both zero or" << endl;

        // Normalize
//         cout << "totNorm = " << totNorm << endl;



//         mergedTaggedJet[i][j][k] = mergedTaggedJet[i][j][k]/totNorm[i][j][k];
//         mergedNotTaggedJet[i][j][k] = mergedNotTaggedJet[i][j][k]/totNorm[i][j][k];

        totTaggedCount = mergedTaggedJet[i][j][k];
        totNotTaggedCount = mergedNotTaggedJet[i][j][k];
        normTot = totTaggedCount+totNotTaggedCount;
        if ( normTot != 0 ) {
          mergedTaggedJet[i][j][k] = totTaggedCount/normTot;
          mergedNotTaggedJet[i][j][k] = totNotTaggedCount/normTot;
        }

        if( mergedTaggedJet[i][j][k] == 0 && mergedNotTaggedJet[i][j][k] == 0 ) {
          mergedTaggedJet[i][j][k] = 0.5;
          mergedNotTaggedJet[i][j][k] = 0.5;
        }

        cout << "mergedTaggedJet["<<i<<"]["<<j<<"]["<<k<<"] = " << mergedTaggedJet[i][j][k];
        cout << " mergedNotTaggedJet["<<i<<"]["<<j<<"]["<<k<<"] = " << mergedNotTaggedJet[i][j][k] << endl;

      }
    }
  }


}
