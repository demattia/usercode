#!/bin/sh

mkdir -p WebPage
mv *.png WebPage

echo "Producing the plots"
root -l -b -q drawCanvases.C

echo "Appending the entries on the web page"
cat template_page.html entriesFile.txt > WebPage/page.html
echo "</body></html>" >> WebPage/page.html

echo "Copy this to the web browser to see the page:"
echo `pwd`"/WebPage/page.html"