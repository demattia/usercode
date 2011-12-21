#!/bin/sh

root -l -b -q acceptance.C+\(\"scaleUp\"\)
root -l -b -q "divide.C(\"_ptup\", \"recreate\")"
root -l -b -q acceptance.C+\(\"scaleDown\"\)
root -l -b -q "divide.C(\"_ptdown\", \"update\")"
mv acceptance0_50_newAccCuts.root acceptance_syst_ptspec.root