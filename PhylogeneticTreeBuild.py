MAIN=`pwd`
BASE=`basename $1 .fasta`
RAW=$MAIN/RAW
GERMLINES=$MAIN/GERMLINES/IGHV4-34.fasta
PAVLOFASTA=$MAIN/MONOCLONALS/EBOV-864.fasta  

if [ ! -d $MAIN/PHYLIP ]; then 
  mkdir $MAIN/PHYLIP
  cd $MAIN/PHYLIP
  mkdir $MAIN/PHYLIP/GENETIC-DISTANCE
  mkdir $MAIN/PHYLIP/NJ-TREE
  mkdir $MAIN/PHYLIP/ML-TREE
  mkdir $MAIN/PHYLIP/DISTANCE-MATRIX
else
  echo "##########################################"
  echo "# CLEAN UP DIRECTORIES BEFORE PROCEEDING #"
  echo "##########################################"
  exit
fi 

echo "##########################################################"
echo "# CREATE FASTA (DEDUPLICATING ON FULL LENGTH SEQUENCES)  #"
echo "##########################################################"
#zcat $RAW/*.csv.gz | grep -v Sequence| gawk -F, '{ table[$NF]++}END{counter=1; for(key in table){printf(">%10.10d\n%s\n",counter,key); counter=counter+1}}'   > TEMP
#cat $GERMLINES $PAVLOFASTA  TEMP > $BASE".fasta"
#rm TEMP

echo "#####################################"
echo "# CREATE PHYLIP FILE USING CLUSTALO #"
echo "#####################################"
if [ -e $BASE".phy" ]; then
  rm $BASE".phy"
fi


##################################
###COPY IN THE FASTAFILE YOU PUT TOGETHER##
##################################
cp $MAIN/864.fasta .
CLUSTAL=/home/elaine/Desktop/Bioinformatics/bin/clustalo-1.2.0-Ubuntu-x86_64
$CLUSTAL  --infile=$BASE".fasta" --infmt=fasta --outfmt=clustal --output-order=input --outfile=$BASE".aln"
$CLUSTAL  --infile=$BASE".fasta" --infmt=fasta --outfmt=phylip --output-order=input --outfile=$BASE".phy"
perl -pi -e 's/\-/\?/g' $BASE".phy"


echo "################################################"
echo "# RUN PHYLIP TO GENERATE THE GENETIC DISTANCES #"
echo "################################################"
cd $MAIN/PHYLIP/GENETIC-DISTANCE
DNADIST="phylip dnadist"
echo "D"  > COMMANDS
echo "D" >> COMMANDS
echo "Y" >> COMMANDS
cp $MAIN/PHYLIP/$BASE".phy" infile
$DNADIST < COMMANDS
if [ -e outfile ]; then 
   mv outfile $BASE".genetic-distance"
else 
   exit
fi
rm infile


echo "####################################"
echo "# RUN PHYLIP TO GENERATE A NJ TREE #"
echo "####################################"
NEIGHBOR="phylip neighbor"
cd $MAIN/PHYLIP/NJ-TREE
cp $MAIN/PHYLIP/GENETIC-DISTANCE/$BASE".genetic-distance" infile
echo "Y" > COMMANDS.NJ
$NEIGHBOR < COMMANDS.NJ
if [ -e outfile ]; then 
   mv outfile $BASE".log-nj"
fi 
if [ -e outtree ]; then
  mv outtree $BASE"-NT.newick"
fi
rm infile


echo "####################################"
echo "# RUN PHYLIP TO GENERATE A ML TREE #"
echo "####################################"
#
#Nucleic acid sequence Maximum Likelihood method, version 3.697
#Settings for this run:
#  U                 Search for best tree?  Yes
#  T        Transition/transversion ratio:  2.0000
#  F       Use empirical base frequencies?  Yes
#  C                One category of sites?  Yes
#  R           Rate variation among sites?  constant rate
#  W                       Sites weighted?  No
#  S        Speedier but rougher analysis?  No, not rough
#  G                Global rearrangements?  No
#  J   Randomize input order of sequences?  No. Use input order
#  O                        Outgroup root?  No, use as outgroup species  1
#  M           Analyze multiple data sets?  No
#  I          Input sequences interleaved?  Yes
#  0   Terminal type (IBM PC, ANSI, none)?  ANSI
#  1    Print out the data at start of run  Yes
#  2  Print indications of progress of run  Yes
#  3                        Print out tree  Yes
#  4       Write out trees onto tree file?  Yes
#  5   Reconstruct hypothetical sequences?  Yes




DNAML="phylip dnaml"
cd $MAIN/PHYLIP/ML-TREE
cp $MAIN/PHYLIP/$BASE".phy" infile
#echo "S"  > COMMANDS.DNAML
echo "1" >> COMMANDS.DNAML
echo "5" >> COMMANDS.DNAML
echo "Y" >> COMMANDS.DNAML
$DNAML < COMMANDS.DNAML
if [ -e outfile ]; then
   mv outfile $BASE".log-dnaml"
fi
if [ -e outtree ]; then
  mv outtree $BASE"-NT-dnaml.newick"
fi
rm infile

echo "#################################"
echo "# GENERATE MATRIX OF IDENTITIES #"
echo "#################################"
cd $MAIN/PHYLIP/DISTANCE-MATRIX
cp $MAIN/PHYLIP/$BASE".fasta" .
$CLUSTAL --percent-id --infile=$BASE".fasta" --infmt=fasta --distmat-out=$BASE".matrix" --full --force
gawk 'NR>1{string=""; for(i=2;i<NR;i++){string=string" "sprintf("%2.2lf ",$i)} print string}' $BASE".matrix" > $BASE".matrix-lower-echelon"
