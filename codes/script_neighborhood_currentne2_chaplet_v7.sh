#!/bin/bash
#$ -cwd

rm script_neighborhood_currentne2_chaplet_v7.sh.*

#Set arguments
NIND=$1
BW=$2
SW=$3
SWR=$4
GENS=$5
GENSMULTI=$6
REPS=$7

rm AAA*
rm repfile.dat
rm genfile*
rm dfilename.dat
rm *data*
rm Ne.*
rm OUTFILE*

###################### run neighborhood #########################

for r in $(seq 1 $REPS)
do

num=$(date +%s)$RANDOM
echo "$num" > seedfile

../NEIGHBORHOOD_v7<<@
1
-99
$NIND	NIND
$BW	BW
$SW	SW
$SWR	SWR
$GENS	GENS
$GENSMULTI	GENSMULTI
@

cp genfile.dat genfile.$NIND.$BW.$SW.rep$r.dat

../plink --file data --make-bed --maf 0.000001 --chr-set 20 --recode --out kk
sleep 5
cp kk.map $r.data.map
cp kk.ped $r.data.ped

cp kk.map $r.datafull.$NIND.$BW.$SW.map
cp kk.ped $r.datafull.$NIND.$BW.$SW.ped

../plink --file data2000 --make-bed --maf 0.000001 --chr-set 20 --recode --out kk2000
sleep 5
cp kk2000.map $r.data2000.map
cp kk2000.ped $r.data2000.ped

fichero="$r.data.ped"
lineas=$(wc -l < "$fichero")

if [ "$lineas" -gt 2000 ]; then
        shuf -n 2000 "$fichero" > tmp && mv tmp "$fichero"
fi
sleep 5

../plink --file dataR --make-bed --maf 0.000001 --chr-set 20 --recode --out kkR
sleep 5
cp kkR.map $r.dataR.map
cp kkR.ped $r.dataR.ped

ficheroR="$r.dataR.ped"
lineas=$(wc -l < "$ficheroR")

if [ "$lineas" -gt 2000 ]; then
        shuf -n 2000 "$ficheroR" > tmp && mv tmp "$ficheroR"
fi
sleep 5

rm data.map
rm data.ped
rm dataR.map
rm dataR.ped
rm kk*
rm ekk*

../currentne2 -t 1 -r 1 $r.data.ped
cp $r.data_currentNe2_OUTPUT.txt data.$NIND.$BW.$SW.rep$r.currentNe2_OUTPUT.txt
awk 'NR==26 || NR==30 || NR==48 || NR==53 || NR==59 || NR==61 {printf "%s ", $0} END {print ""}' $r.data_currentNe2_OUTPUT.txt >> OUTFILE_currentne2_$NIND.$BW.$SW

../currentne2 -t 1 -r 1 -x $r.data.ped
cp $r.data_currentNe2_mix_OUTPUT.txt data.$NIND.$BW.$SW.rep$r.currentNe2_mix_OUTPUT.txt
awk 'NR==53 || NR==58 || NR==61 || NR==64 || NR==70 || NR==72 {printf "%s ", $0} END {print ""}' $r.data_currentNe2_mix_OUTPUT.txt >> OUTFILE_currentne2_mix_$NIND.$BW.$SW

../currentne2 -t 1 -r 1 $r.data2000.ped
cp $r.data2000_currentNe2_OUTPUT.txt data2000.$NIND.$BW.$SW.rep$r.currentNe2_OUTPUT.txt
awk 'NR==26 || NR==30 || NR==48 || NR==53 || NR==59 || NR==61 {printf "%s ", $0} END {print ""}' $r.data2000_currentNe2_OUTPUT.txt >> OUTFILE2000_currentne2_$NIND.$BW.$SW

../currentne2 -t 1 -r 1 -x $r.data2000.ped
cp $r.data2000_currentNe2_mix_OUTPUT.txt data2000.$NIND.$BW.$SW.rep$r.currentNe2_mix_OUTPUT.txt
awk 'NR==53 || NR==58 || NR==61 || NR==64 || NR==70 || NR==72 {printf "%s ", $0} END {print ""}' $r.data2000_currentNe2_mix_OUTPUT.txt >> OUTFILE2000_currentne2_mix_$NIND.$BW.$SW

../currentne2 -t 1 -r 1 $r.dataR.ped
cp $r.dataR_currentNe2_OUTPUT.txt dataR.$NIND.$BW.$SW.rep$r.currentNe2_OUTPUT.txt
awk 'NR==26 || NR==30 || NR==48 || NR==53 || NR==59 || NR==61 {printf "%s ", $0} END {print ""}' $r.dataR_currentNe2_OUTPUT.txt >> OUTFILER_currentne2_$NIND.$BW.$SW

../currentne2 -t 1 -r 1 -x $r.dataR.ped
cp $r.dataR_currentNe2_mix_OUTPUT.txt dataR.$NIND.$BW.$SW.rep$r.currentNe2_mix_OUTPUT.txt
awk 'NR==53 || NR==58 || NR==61 || NR==64 || NR==70 || NR==72 {printf "%s ", $0} END {print ""}' $r.dataR_currentNe2_mix_OUTPUT.txt >> OUTFILER_currentne2_mix_$NIND.$BW.$SW

done

#################################################

