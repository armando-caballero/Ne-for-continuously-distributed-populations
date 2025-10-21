#!/bin/bash
#$ -cwd

rm script_neighborhood.sh.*

#Set arguments
NIND=$1
Vk=$2
BW=$3
SW=$4
REPS=$5

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

./NEIGHBORHOOD<<@
0
-99
$NIND	NIND
$Vk	Vk
$BW	BW
$SW	SW
10000	GENS
@

cp genfile.dat genfile.$NIND.$BW.$SW.rep$r.dat

../plink --file data --make-bed --maf 0.0001 --chr-set 20 --recode --out kk
sleep 5
cp kk.map $r.data.map
cp kk.ped $r.data.ped
rm kk*

fichero="$r.data.ped"
lineas=$(wc -l < "$fichero")

if [ "$lineas" -gt 2000 ]; then
        shuf -n 2000 "$fichero" > tmp && mv tmp "$fichero"
fi
sleep 5

./currentne2 $r.data.ped
cp $r.data_currentNe2_OUTPUT.txt data.$NIND.$BW.$SW.rep$r.currentNe2_OUTPUT.txt
awk 'NR==53 || NR==58 || NR==61 || NR==75 || NR==81 || NR==83 {printf "%s ", $0} END {print ""}' $r.data_currentNe2_OUTPUT.txt >> OUTFILE_currentne2_$NIND.$BW.$SW

./currentne2 -x $r.data.ped
cp $r.data_currentNe2_mix_OUTPUT.txt data.$NIND.$BW.$SW.rep$r.currentNe2_mix_OUTPUT.txt
awk 'NR==26 || NR==30 || NR==48 || NR==53 || NR==59 || NR==61 {printf "%s ", $0} END {print ""}' $r.data_currentNe2_mix_OUTPUT.txt >> OUTFILE_currentne2_mix_$NIND.$BW.$SW

done

paste Ne.$NIND.$BW.$SW.rep1 Ne.$NIND.$BW.$SW.rep2 Ne.$NIND.$BW.$SW.rep3 > OUTFILE_gone2_$NIND.$BW.$SW
paste Ne_mix.$NIND.$BW.$SW.rep1 Ne_mix.$NIND.$BW.$SW.rep2 Ne_mix.$NIND.$BW.$SW.rep3 > OUTFILE_gone2_$NIND.$BW.$SW

#################################################


