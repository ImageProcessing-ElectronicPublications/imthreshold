#!/bin/sh

# gs -o noIMG.pdf   -sDEVICE=pdfwrite -dFILTERIMAGE                input.pdf
# gs -o noTXT.pdf   -sDEVICE=pdfwrite -dFILTERTEXT                 input.pdf
# gs -o noVCT.pdf   -sDEVICE=pdfwrite -dFILTERVECTOR               input.pdf
#
# gs -o onlyTXT.pdf -sDEVICE=pdfwrite -dFILTERVECTOR -dFILTERIMAGE input.pdf 
# gs -o onlyIMG.pdf -sDEVICE=pdfwrite -dFILTERVECTOR -dFILTERTEXT  input.pdf
# gs -o onlyVCT.pdf -sDEVICE=pdfwrite -dFILTERIMAGE  -dFILTERTEXT  input.pdf

gs -sDEVICE=pdfwrite -dFILTERIMAGE -o "$1-noimages.pdf" "$1"
