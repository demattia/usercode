{
  gROOT->ProcessLine(".L checkOpenHLT.C+");
  gROOT->ProcessLine("checkOpenHLT a(\"SAMPLENAME/openhlt_merge.root\", PARALLELISMDIFFCUT)");
  gROOT->ProcessLine("a.Loop()");
}
