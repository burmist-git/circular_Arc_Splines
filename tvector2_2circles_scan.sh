rm -rf tvector2_2circles_C.d
rm -rf tvector2_2circles_C_ACLiC_dict_rdict.pcm
rm -rf tvector2_2circles_C.so

outGifName="tvector2_2circles_fit.gif"
rm -rf $outGifName
outNamePreff="tvector2_2circles_fit_"

for xA in 0.75 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0.0 -0.1 -0.2 -0.3 -0.4 -0.5 -0.6 -0.7 -0.8 -0.9 -1.0 -1.1 -1.2 -1.3 -1.4 -1.5 -1.6 -1.7 ; do
    echo "$xA"
    outName=$outNamePreff$xA'.gif'
    root -l -b -q tvector2_2circles.C+($xA,\"$outName\")
    outNameList="$outNameList $outName"
done

convert -delay 40 -loop 0 $outNameList $outGifName

#for i in $(seq 1 1 10); do
#    echo "$i"
#done



