#!/bin/sh

# mfbpdf-jbig2jpge.sh
# depends: imthreshold, libtiff, jbig2enc, jpge, bc
# imthreshold: https://github.com/ImageProcessing-ElectronicPublications/imthreshold
# libtiff: http://libtiff.maptools.org/
# jbig2enc: https://github.com/agl/jbig2enc
# jpge: https://github.com/ImageProcessing-ElectronicPublications/jpge

fin=$1
dpi=$2
quality=$3
if [ -z $1 ]
then
    printf "Usage: $0 image [DPI=600] [quality=75]\n"
    exit 1
fi
if [ -z $2 ]
then
    dpi=600
fi
if [ -z $3 ]
then
    quality=75
fi
fout=$fin.pdf
fmask=$fin.m.tif
ffg=$fin.fg.png
fbg=$fin.bg.png
fjb2=$fin.m.jbig2
ffjpg=$fin.fg.jpg
fbjpg=$fin.bg.jpg

tauthor="$(whoami)"
ttitle="$fin"
tdate="$(date +%Y%m%d%H%M%S%z)"

imthreshold-tdjvul -c 3 $fin $fmask $ffg $fbg
tiffset -s 282 $dpi $fmask
tiffset -s 283 $dpi $fmask
jbig2 -p -v $fmask > $fjb2
jpge $ffg $ffjpg $quality
jpge $fbg $fbjpg $quality

sfjb2=$(cat $fjb2 | wc -c)
sffjpg=$(cat $ffjpg | wc -c)
sfbjpg=$(cat $fbjpg | wc -c)

heightm=$(imthreshold-finfo $fin | grep Height | cut -d' ' -f2)
widthm=$(imthreshold-finfo $fin | grep Width | cut -d' ' -f2)
heightf=$(imthreshold-finfo $ffg | grep Height | cut -d' ' -f2)
widthf=$(imthreshold-finfo $ffg | grep Width | cut -d' ' -f2)
heightb=$(imthreshold-finfo $fbg | grep Height | cut -d' ' -f2)
widthb=$(imthreshold-finfo $fbg | grep Width | cut -d' ' -f2)
hpointm=$(echo "$heightm*7200/$dpi" | bc)
wpointm=$(echo "$widthm*7200/$dpi" | bc)
hpoint=$(echo "$hpointm*0.01" | bc)
wpoint=$(echo "$wpointm*0.01" | bc)
spoint=$(printf "%s%s" $wpointm $hpointm | wc -c)
sobj3=$(echo "161+4*$spoint" | bc)

printf "%%PDF-1.4\n" > $fout
printf "%%âãÏÓ\n" >> $fout
printf "\n" >> $fout

xref1=$(cat $fout | wc -c)

printf "1 0 obj\n" >> $fout
printf "<<\n" >> $fout
printf "/Kids [2 0 R]\n" >> $fout
printf "/Type /Pages\n" >> $fout
printf "/Count 1\n" >> $fout
printf ">>\n" >> $fout
printf "endobj\n" >> $fout
printf "\n" >> $fout

xref2=$(cat $fout | wc -c)

printf "2 0 obj\n" >> $fout
printf "<<\n" >> $fout
printf "/pdftk_PageNum 1\n" >> $fout
printf "/Contents 3 0 R\n" >> $fout
printf "/Type /Page\n" >> $fout
printf "/Resources \n" >> $fout
printf "<<\n" >> $fout
printf "/ProcSet [/PDF /ImageC]\n" >> $fout
printf "/ExtGState 4 0 R\n" >> $fout
printf "/XObject 5 0 R\n" >> $fout
printf ">>\n" >> $fout
printf "/Parent 1 0 R\n" >> $fout
printf "/MediaBox [0 0 %s %s]\n" $wpoint $hpoint >> $fout
printf ">>\n" >> $fout
printf "endobj\n" >> $fout
printf "\n" >> $fout

xref3=$(cat $fout | wc -c)

printf "3 0 obj\n" >> $fout
printf "<<\n" >> $fout
printf "/Length %d\n" $sobj3 >> $fout
printf ">>\n" >> $fout
printf "stream\n" >> $fout
printf "q 0.01 0 0 0.01 0 0 cm\n" >> $fout
printf "q\n" >> $fout
printf "1 0 m\n" >> $fout
printf "1 %s l\n" $hpointm >> $fout
printf "%s %s l\n" $wpointm $hpointm >> $fout
printf "%s 0 l\n" $wpointm >> $fout
printf "h\n" >> $fout
printf "W n\n" >> $fout
printf "q %s 0 0 %s 0 0 cm\n" $wpointm $hpointm >> $fout
printf "/R7 Do\n" >> $fout
printf "Q\n" >> $fout
printf "Q\n" >> $fout
printf "q\n" >> $fout
printf "1 0 m\n" >> $fout
printf "1 %s l\n" $hpointm >> $fout
printf "%s %s l\n" $wpointm $hpointm >> $fout
printf "%s 0 l\n" $wpointm >> $fout
printf "h\n" >> $fout
printf "W n\n" >> $fout
printf "q %s 0 0 %s 0 0 cm\n" $wpointm $hpointm >> $fout
printf "/R9 Do\n" >> $fout
printf "Q\n" >> $fout
printf "Q\n" >> $fout
printf "/R10 gs\n" >> $fout
printf "Q\n" >> $fout
printf "endstream\n" >> $fout
printf "endobj\n" >> $fout
printf "\n" >> $fout

xref4=$(cat $fout | wc -c)

printf "4 0 obj\n" >> $fout
printf "<<\n" >> $fout
printf "/R10 6 0 R\n" >> $fout
printf ">>\n" >> $fout
printf "endobj\n" >> $fout
printf "\n" >> $fout

xref5=$(cat $fout | wc -c)

printf "5 0 obj\n" >> $fout
printf "<<\n" >> $fout
printf "/R7 7 0 R\n" >> $fout
printf "/R8 8 0 R\n" >> $fout
printf "/R9 9 0 R\n" >> $fout
printf ">>\n" >> $fout
printf "endobj\n" >> $fout
printf "\n" >> $fout

xref6=$(cat $fout | wc -c)

printf "6 0 obj\n" >> $fout
printf "<<\n" >> $fout
printf "/Type /ExtGState\n" >> $fout
printf "/TR /Identity\n" >> $fout
printf ">>\n" >> $fout
printf "endobj\n" >> $fout
printf "\n" >> $fout

xref7=$(cat $fout | wc -c)

printf "7 0 obj\n" >> $fout
printf "<<\n" >> $fout
printf "/Type /XObject\n" >> $fout
printf "/Subtype /Image\n" >> $fout
printf "/BitsPerComponent 8\n" >> $fout
printf "/ColorSpace /DeviceRGB\n" >> $fout
printf "/Filter /DCTDecode\n" >> $fout
printf "/Height %d\n" $heightf >> $fout
printf "/Width %d\n" $widthf >> $fout
printf "/Length %d\n" $sffjpg >> $fout
printf ">>\n" >> $fout
printf "stream\n" >> $fout

cat $ffjpg  >> $fout

printf "\nendstream\n" >> $fout
printf "endobj\n" >> $fout
printf "\n" >> $fout

xref8=$(cat $fout | wc -c)

printf "8 0 obj\n" >> $fout
printf "<< /Type /XObject\n" >> $fout
printf "/Subtype /Image\n" >> $fout
printf "/ColorSpace /DeviceGray\n" >> $fout
printf "/BitsPerComponent 1\n" >> $fout
printf "/Filter /JBIG2Decode\n" >> $fout
printf "/Height %d\n" $heightm >> $fout
printf "/Width %d\n" $widthm >> $fout
printf "/Length %d\n" $sfjb2 >> $fout
printf ">>\n" >> $fout
printf "stream\n" >> $fout

cat $fjb2 >> $fout

printf "\nendstream\n" >> $fout
printf "endobj\n" >> $fout
printf "\n" >> $fout

xref9=$(cat $fout | wc -c)

printf "9 0 obj\n" >> $fout
printf "<<\n" >> $fout
printf "/Type /XObject\n" >> $fout
printf "/Subtype /Image\n" >> $fout
printf "/BitsPerComponent 8\n" >> $fout
printf "/ColorSpace /DeviceRGB\n" >> $fout
printf "/Filter /DCTDecode\n" >> $fout
printf "/SMask 8 0 R\n" >> $fout
printf "/Height %d\n" $heightb >> $fout
printf "/Width %d\n" $widthb >> $fout
printf "/Length %d\n" $sfbjpg >> $fout
printf ">>\n" >> $fout
printf "stream\n" >> $fout

cat $fbjpg  >> $fout

printf "\nendstream\n" >> $fout
printf "endobj\n" >> $fout
printf "\n" >> $fout

xrefA=$(cat $fout | wc -c)

printf "10 0 obj\n" >> $fout
printf "<<\n" >> $fout
printf "/Type /Catalog\n" >> $fout
printf "/Pages 1 0 R\n" >> $fout
printf ">>\n" >> $fout
printf "endobj\n" >> $fout
printf "\n" >> $fout

xrefB=$(cat $fout | wc -c)

printf "11 0 obj\n" >> $fout
printf "<<\n" >> $fout
printf "/Author ($tauthor)\n" >> $fout
printf "/Title ($ttitle)\n" >> $fout
printf "/Subject (MFB (Mask+Fg+Bg))\n" >> $fout
printf "/Creator (imthreshold: https://github.com/ImageProcessing-ElectronicPublications/imthreshold)\n" >> $fout
printf "/Producer (imthreshold: https://github.com/ImageProcessing-ElectronicPublications/imthreshold, libtiff: http://libtiff.maptools.org/, jbig2enc: https://github.com/agl/jbig2enc, jpge: https://github.com/ImageProcessing-ElectronicPublications/jpge)\n" >> $fout
printf "/CreationDate (D:%s)\n" $tdate >> $fout
printf "/ModDate (D:%s)\n" $tdate >> $fout
printf ">>\n" >> $fout
printf "endobj\n" >> $fout
printf "\n" >> $fout

xrefC=$(cat $fout | wc -c)

printf "xref\n" >> $fout
printf "0 12\n" >> $fout
printf "%010d 65535 f\n" "0" >> $fout
printf "%010d 00000 n\n" "$xref1" >> $fout
printf "%010d 00000 n\n" "$xref2" >> $fout
printf "%010d 00000 n\n" "$xref3" >> $fout
printf "%010d 00000 n\n" "$xref4" >> $fout
printf "%010d 00000 n\n" "$xref5" >> $fout
printf "%010d 00000 n\n" "$xref6" >> $fout
printf "%010d 00000 n\n" "$xref7" >> $fout
printf "%010d 00000 n\n" "$xref8" >> $fout
printf "%010d 00000 n\n" "$xref9" >> $fout
printf "%010d 00000 n\n" "$xrefA" >> $fout
printf "%010d 00000 n\n" "$xrefB" >> $fout
printf "trailer\n" >> $fout
printf "<<\n" >> $fout
printf "/Info 11 0 R\n" >> $fout
printf "/Root 10 0 R\n" >> $fout
printf "/Size 12\n" >> $fout
printf ">>\n" >> $fout
printf "startxref\n" >> $fout
printf "%d\n" "$xrefC" >> $fout
printf "%%%%EOF\n" >> $fout

xrefZ=$(cat $fout | wc -c)
printf "%s: %s\n" $fout $xrefZ
