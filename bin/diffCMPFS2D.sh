#!/bin/bash

########################## diffCMPFS2D ###########################################
#
# author: Frédéric Darboux <Frederic.Darboux@orleans.inra.fr> (2012-2017)
# version: 1.06.01
# date: 2017-03-03
#
# License Cecill-V2 <http://www.cecill.info/licences/Licence_CeCILL_V2-en.html>
#
# (c) CNRS - Universite d'Orleans - INRA (France)
#
# This file is part of FullSWOF_2D software.
# <https://sourcesup.renater.fr/projects/fullswof-2d/>
#
# FullSWOF_2D = Full Shallow-Water equations for Overland Flow,
# in two dimensions of space.
# This software is a computer program whose purpose is to compute
# solutions for 2D Shallow-Water equations.
#
# LICENSE
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software. You can use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# <http://www.cecill.info>.
#
# As a counterpart to the access to the source code and rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty and the software's author, the holder of the
# economic rights, and the successive licensors have only limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading, using, modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean that it is complicated to manipulate, and that also
# therefore means that it is reserved for developers and experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and, more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#
##################################################################################

##################################################################################
# The present script "diffCMPFS2D" helps in benchmarking the results of
# FullSWOF_2D. It computes absolute and relative differences between
# files created by the script BenchEval2D.sh
##################################################################################

VERSION="1.06.01, 2017-03-03"

# Define which command for awk
STDERR=2

#to have . as decimal separator
LANG=C

if [ $# -ne 2 ]
then
  echo "Usage: diffCMPFS2D ReferenceCMP_filename EvalCMP_filename" >&$STDERR
  echo "Version $VERSION" >&$STDERR
  exit 1;
fi ;

# Define which command for awk
AWK_CMD=gawk
# Define which command for paste
PASTE_CMD=paste
# Define which command for mktemp
MKTEMP_CMD=mktemp
# Define which command for rm
RM_CMD=rm

### Check for programs
ALLHERE=0
command -v $AWK_CMD >/dev/null 2>&1 || { echo >&2 "Error: program" $AWK_CMD "required but not installed"; ALLHERE=1; }
command -v $PASTE_CMD >/dev/null 2>&1 || { echo >&2 "Error: program" $PASTE_CMD "required but not installed"; ALLHERE=1; }
command -v $MKTEMP_CMD >/dev/null 2>&1 || { echo >&2 "Error: program" $MKTEMP_CMD "required but not installed"; ALLHERE=1; }
command -v $RM_CMD >/dev/null 2>&1 || { echo >&2 "Error: program" $RM_CMD "required but not installed"; ALLHERE=1; }
if [ $ALLHERE -eq 1 ] #at least one program missing...
then
 exit 1 ;
fi

# Below which absolute value to consider a float as equal to zero
FLOATEQUALZERO=1e-300

# number of header lines in the input files
HEADERLENGTH=6

#Define which input file to use
REF_FILE=$1
EVAL_FILE=$2

#Create temporary files
Ref1Tempfile=$($MKTEMP_CMD Ref1_XXXXXXXX)
Eval1Tempfile=$($MKTEMP_CMD Eval1_XXXXXXXX)
Diff1Tempfile=$($MKTEMP_CMD Diff1_XXXXXXXX)
Diff2Tempfile=$($MKTEMP_CMD Diff2_XXXXXXXX)
Diff3Tempfile=$($MKTEMP_CMD Diff3_XXXXXXXX)
Diff4Tempfile=$($MKTEMP_CMD Diff4_XXXXXXXX)

##function to remove temporary files
function rmtmp(){
rm "$Ref1Tempfile" "$Eval1Tempfile" "$Diff1Tempfile" "$Diff2Tempfile" "$Diff3Tempfile" "$Diff4Tempfile"
return 0
}

##Remove unuseful lines
#Remove header from Reference CMP
$AWK_CMD 'FNR>HEADERLENGTH' HEADERLENGTH=$HEADERLENGTH "$REF_FILE" > "$Ref1Tempfile"

#Remove header from Evaluated CMP
$AWK_CMD 'FNR>HEADERLENGTH' HEADERLENGTH=$HEADERLENGTH "$EVAL_FILE" > "$Eval1Tempfile"

#check if same number of lines
if (test "$(wc -l < "$Ref1Tempfile")" -ne "$(wc -l < "$Eval1Tempfile")"); then
	echo "Error: Numbers of lines are different." >&$STDERR;
	rmtmp;
	exit 1;
fi

#create a single file
$PASTE_CMD "$Ref1Tempfile" "$Eval1Tempfile" > "$Diff1Tempfile"

#check if line titles are identical
$AWK_CMD '{if (($1 != $4)||($2 != $5)){print "Error: Line titles are different at line", FNR, ":", $0; exit 1}}' "$Diff1Tempfile" >&$STDERR
if (test $? -ne 0); then
	echo "Error: Problem while comparing line titles in the files"  >&$STDERR
	rmtmp;
	exit 1;
fi

#write header in output file

echo    "##############################################################################"
echo    "# Generated by diffCMPFS2D version $VERSION"
echo    "# on $(date "+%Y-%m-%d %H:%M:%S") ($(id -un)@$(hostname))"
echo -e "# from reference  CMP file: \t$(pwd)/$REF_FILE"
echo -e "# from evaluation CMP file: \t$(pwd)/$EVAL_FILE"
echo    "##############################################################################"
echo    "#category subcategory absolute_diff relative_diff%"

##compute absolute and relative differences
#absolute differences
$AWK_CMD '{print $6-$3}' "$Diff1Tempfile" > "$Diff2Tempfile"
#relative differences
$AWK_CMD '{if($3!="NaN" && $6!="NaN" && ($3<-'$FLOATEQUALZERO' || $3>'$FLOATEQUALZERO')){print 100*($6/$3-1)}else{print "NaN"}}' "$Diff1Tempfile" > "$Diff3Tempfile"

## concatenate results
# extract line titles
$AWK_CMD '{print $1 " " $2}' "$Diff1Tempfile" > "$Diff4Tempfile"
# paste results
$PASTE_CMD "$Diff4Tempfile" "$Diff2Tempfile" "$Diff3Tempfile"

rmtmp;

exit 0
