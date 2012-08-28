#ifndef GETCATEGORY_H
#define GETCATEGORY_H

TString getCategory(const TString & dirName)
{
  // std::cout << "dirName = " << dirName << std::endl;
  TString fileName(dirName);
  if( std::string(fileName).find("/") != std::string::npos ) {
    fileName = std::string(dirName).substr(std::string(dirName).find_last_of("/")+1, std::string::npos);
  }
  if( fileName.BeginsWith("Signal") ) return "Signal";
  else if( fileName.BeginsWith("QCD") ) return "QCD";
  else if( fileName.BeginsWith("ZZ") ) return "ZZ";
  else if( fileName.BeginsWith("WW") ) return "WW";
  else if( fileName.BeginsWith("WZ") ) return "WZ";
  else if( fileName.BeginsWith("TTJets") ) return "TTJets";
  else if( fileName.BeginsWith("Zee") ) return "Zee";
  else if( fileName.BeginsWith("Zmumu") ) return "Zmumu";
  else if( fileName.BeginsWith("Ztautau") ) return "Ztautau";
  else if( fileName.BeginsWith("DYJets") ) return "DYJets";
  else if( fileName.BeginsWith("Data") ) return "Data";
  std::cout << "Unknown cathegory" << std::endl;
  return "";
}

#endif
