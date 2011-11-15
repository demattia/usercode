#!/bin/sh

echo "title = $1"

echo "Saving plots to png"

root -l -b -q TakePlots.C
mv WebPage OldWebPage
mkdir -p WebPage
mv *.png WebPage

echo "Preparing html page"

cat page.html | sed s/TITLE/"${1}"/ > WebPage/page.html

