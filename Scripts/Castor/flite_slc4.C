// Use this to set the environment to compile FWLite macros for CMSSW versions prior to 1_5_0
{
  gROOT->Reset();

  gSystem->Load("libFWCoreFWLite.so"); 
  AutoLibraryLoader::enable();

  const char* env = gSystem->Getenv("CMSSW_BASE");
  if( 0 != env) {
    string dir(env);
    dir += "/src";
    gROOT->GetInterpreter()->AddIncludePath(dir.c_str());
  }

  env = gSystem->Getenv("CMSSW_RELEASE_BASE");
  if( 0 != env) {
    string dir(env);
    dir += "/src";
    gROOT->GetInterpreter()->AddIncludePath(dir.c_str());
  }

  const char* env = gSystem->Getenv("CLHEP_PARAM_PATH");
  if( 0 != env) {
    string dir(env);
    dir += "/include";
    gROOT->GetInterpreter()->AddIncludePath(dir.c_str());

    //can find boost relative to CLHEP
    string boostDir(env);
    boostDir +="/../../boost/1.33.1-CMS3q/include";
    gROOT->GetInterpreter()->AddIncludePath(boostDir.c_str());
  }
}
