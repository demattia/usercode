root -b -q -l compileTMVA.C
root -b -q -l TMVAClassification.C+\(\"barrel\",\"\",\"\"\) > TMVALog_barrel.txt &
root -b -q -l TMVAClassification.C+\(\"endcaps\",\"\",\"\"\) > TMVALog_endcaps.txt &

root -b -q -l TMVAClassification.C+\(\"barrel\",\"\",\"\"\) > TMVALog_barrel_0.txt &
root -b -q -l TMVAClassification.C+\(\"barrel\",\"\",\"\"\) > TMVALog_barrel_1.txt &
root -b -q -l TMVAClassification.C+\(\"barrel\",\"\",\"\"\) > TMVALog_barrel_2.txt &

root -b -q -l TMVAClassification.C+\(\"endcaps\",\"\",\"\"\) > TMVALog_endcaps_0.txt &
root -b -q -l TMVAClassification.C+\(\"endcaps\",\"\",\"\"\) > TMVALog_endcaps_1.txt &
root -b -q -l TMVAClassification.C+\(\"endcaps\",\"\",\"\"\) > TMVALog_endcaps_2.txt &
