{
  gROOT->ProcessLine(".L checkOpenHLT.C+");
  gROOT->ProcessLine("checkOpenHLT a(\"part2/openhlt_merge.root\", 3.30000000000000000000)");
  gROOT->ProcessLine("a.Loop()");
}
