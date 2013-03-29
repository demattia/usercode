root -b -q -l compileTMVA.C
root -b -q -l TMVAClassification.C+\(\"barrel\",\"\",\"\",\"barrelWeights\"\) > TMVALog_barrel.txt &
root -b -q -l TMVAClassification.C+\(\"endcaps\",\"\",\"\",\"endcapsWeights\"\) > TMVALog_endcaps.txt &

root -b -q -l TMVAClassification.C+\(\"barrel\",\"\",\"\",\"barrel0Weights\"\) > TMVALog_barrel_0.txt &
root -b -q -l TMVAClassification.C+\(\"barrel\",\"\",\"\",\"barrel1Weights\"\) > TMVALog_barrel_1.txt &
root -b -q -l TMVAClassification.C+\(\"barrel\",\"\",\"\",\"barrel2Weights\"\) > TMVALog_barrel_2.txt &

root -b -q -l TMVAClassification.C+\(\"endcaps\",\"\",\"\",\"endcaps0Weights\"\) > TMVALog_endcaps_0.txt &
root -b -q -l TMVAClassification.C+\(\"endcaps\",\"\",\"\",\"endcaps1Weights\"\) > TMVALog_endcaps_1.txt &
root -b -q -l TMVAClassification.C+\(\"endcaps\",\"\",\"\",\"endcaps2Weights\"\) > TMVALog_endcaps_2.txt &
