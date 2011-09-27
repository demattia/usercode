#!/bin/sh

echo "Producing the plots"
root -l -b -q drawCanvases.C

echo "Appending the entries on the web page"
cat template_page.html entriesFile.txt > page.html
echo "</body></html>" >> page.html

echo "Copy this to the web browser to see the page:"
echo `pwd`"/page.html"