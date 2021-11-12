#!/usr/bin/bash
rm -f *.exercise.out
goodlines=$(grep -v "#" $inputfile | wc -l | cut -d ' ' -f1)
unset IFS
unset dataline
shortHSP=0; 
hspcounter=0;
echo -e "We have $goodlines data lines for processing...\n"

dupS_acc=()
group1cut=150 
group2cut=250 
group3cut=350 
outfile1="HSPscore.$group1cut.exercise.out" 
outfile2="HSPscore.$group2cut.exercise.out" 
outfile3="HSPscore.$group3cut.exercise.out" 
outfile4="HSPscore. morethan$group3cut.exercise.out”

rm -f $outfile1 $outfile2 $outfile3 $outfile4
while read wholeline
do
# echo "Line is $wholeline"
# Don't want the lines starting with #
# Uses a substring notation: offset from first character, length wanted

if test ${wholeline:0:1} != "#"
   then
   dataline=$((dataline+1))
   # echo "Line $dataline starts with ${wholeline:0:1}"

# Split the line into the fields we want
read Q_acc S_acc pc_identity alignment_length mismatches gap_opens Q_start Q_end S_start S_end evalue bitscore <<< ${wholeline}

# We need to ensure that $bitscore is an integer number not a real/float!
# Bash CAN'T deal properly with non-integers, so make it an integer
# how to format outputs using printf
bitscore=$(printf "%.0f\n" $bitscore)
# list the subject accession for all HSPs
echo -e "$dataline\t$Q_acc\t$S_acc" >> Subject_accessions.exercise.out

# List the alignment length and percent ID for all HSPs
echo -e "$dataline\t$alignment_length\t$pc_identity" >> al_leng_pcID.exercise.out

# Show the HSPs with more than 20 mismatches
if test $mismatches -gt 20 
    then
    echo -e "$dataline\tmore than 20 mismatches:\t$Q_acc $S_acc $mismatches"
fi    

# show the HSPs shorter than 100 amino acids and with more than 20 mismatches
# Are the alignment length amino acids or nucleotides!? aa
if test $alignment_length -lt 100 && test $mismatches -gt 20
     then
     echo -e "$dataline\tHSP shorter than 100aa, more than 20 mismatches:\t$alignment_length\t$mismatches"
fi

# list the first 20 HSPs that have fewer than 20 mismatches
if test $mismatches -lt 20 
    then
    hspcounter=$((hspcounter+1))
    if test $hspcounter -le 20
      then
      hsp_array+=$wholeline
      echo -e "$dataline\t$hspcounter\t$wholeline"   >> Fewer.than20MM.exercise.out
      fi
fi

# how many HSPs are shorter than 100 amino acids?
if test $alignment_length -lt 100 
     then
     shortHSP=$((shortHSP+1))
fi

# list the top ten highest (best) HSPs.
# BLAST output default is to give the best first
if test $dataline -le 10
     then
     echo -e "$dataline\t$wholeline" >> Top10.HSPs.exercise.out
fi

# list the start positions of all matches where the HSP Subject accession includes the letters string AEI.
if [[ $S_acc == *"AEI"* ]]; then
     echo -e "$dataline\t$S_acc contains AEI: Subject starts at $S_start, Query starts at $Q_start" >> AEIinSubjectAcc.starts.exercise.out 
fi     

# how many subject sequences have more than one HSP?
if test $S_acc == $pre_acc  
    then
    dupecount=$((dupecount+1))
    if [[ dupecount == 1 ]]; then
       dupS_acc=$S_acc
    fi

# Some might have more than 2: use wildcard pattern matching to see if we have it already
if [[ $dupS_acc == *$S_acc* ]]; then 
      echo ""
      else
      dupS_acc+=($S_acc)
    fi  
fi    
pre_acc=$S_acc

# what percentage of each HSP is made up of mismatches? This is a rounded
# number as bash doesn't do floating point maths... =-(
MMpercent=$((100*$mismatches/$alignment_length))
echo -e "$dataline\t$alignment_length\t$mismatches\t$MMpercent" >> Mismatchpercent.exercise.out

# allocate HSPs into different groups based on their scores
scorebin=1
if [ $bitscore -gt $group3cut ]; then 
     scorebin=4
fi
if [ $bitscore -le $group3cut ] && [ $bitscore -gt $group2cut ]; then 
     scorebin=3
fi
if [ $bitscore -le $group2cut ] && [ $bitscore -gt $group1cut ]; then 
     scorebin=2
fi
# Use a case statement to send the output to the relevant file
scoregroupdetails=$(echo -e "$dataline\t$Q_acc\t$S_acc\t$bitscore")
case $scorebin in
  4) 
    echo -e "$scoregroupdetails" >> $outfile4
    ;;
  3) 
    echo -e "$scoregroupdetails" >> $outfile3
    ;;
  2) 
    echo -e "$scoregroupdetails" >> $outfile2
    ;;
  1) 
    echo -e "$scoregroupdetails" >> $outfile1
    ;;
esac
# END BLOCK EQUIVALENT
if test $dataline -eq $goodlines
  then
  echo -e "\n\nENDBLOCK\n\nThere were $shortHSP HSPs shorter than 100 amino acids"
  echo -e "There were ${#dupS_acc[@]} Subjects with multiple HSPs"
fi
fi  # was not a commented line in the blast data
done < $inputfile

