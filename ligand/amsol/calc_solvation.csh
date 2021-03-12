#!/bin/csh -f
# doamsol.csh

# modified by Trent E Balius Nov, 2013
# modified by T.B. Adler, Shoichet group
# modified by Teague Sterling, Sept. 2014 (Fix missing file when antechamber is missing)

set mol2file = $1
set amsoltlimit = 1m

# obtain the name of the protonated molecule from the mol2-file name:
set mol2_file_FullName = "$1" #e.g. ZINC00000.0.mol2
set constituents = `echo $mol2_file_FullName:q | sed 's/\./ /g'`
set ProtonatedMoleculeName = "$constituents[1].$constituents[2]"
echo "ProtonatedMoleculeName: $ProtonatedMoleculeName"
set MoleculeName = $ProtonatedMoleculeName


# obtain the name of the protonated molecule from the mol2-file name:
set mol2_file_FullName = "$1" #e.g. ZINC00000.0.mol2
set constituents = `echo $mol2_file_FullName:q | sed 's/\./ /g'`
set ProtonatedMoleculeName = "$constituents[1].$constituents[2]"
echo "ProtonatedMoleculeName: $ProtonatedMoleculeName"
set MoleculeName = $ProtonatedMoleculeName


if -e temp.mol2 then
    echo "warning: temp.mol2. Rewriting link."
endif

echo " Preparing AMSOL7.1 input for $mol2file (first: transformation to ZmatMOPAC by openbabel)"
echo " If there is any trouble, make sure that your DOCKBASE and OBABELBASE"
echo " is correctly set in ~/.cshrc or ~/.bashrc"

ln -svfn $mol2file temp.mol2
if ! $?AMSOLEXE then
    set AMSOLEXE = $DOCKBASE/ligand/amsol/amsol7.1
endif
echo "AMSOLEXE is $AMSOLEXE ."

if ! $?OBABELEXE then
   set OBABELEXE=$OBABELBASE/obabel
   if ! -e $OBABELEXE then
       set OBABELEXE=$OBABELBASE/bin/obabel
   endif
echo " OBABELEXE is ${OBABELEXE} ."
if ! -e ${OBABELEXE} then
    echo "Couldn't fine OBABLE at ${OBABELEXE}"
    exit -1
endif

# obtain the entire current path:
set current_path = `pwd`
echo "current_path ::: $current_path"

# omega does not always assign proper SYBYL atom types (e.g., just S, instead of S.O2 for example): 
# This can lead to obabel warnings during conversion from mol2 to ZmatMOPAC format:
# Therefore, if antechamber is available, we use it to avoid these warnings and we trust antechamber's ability
# to assign correct atom types.

set TEMP_FILE = "${current_path}/temp.mol2"
if( `where  antechamber || echo ''` == "" ) then
      echo "antechamber (part of ambertools (downloadable for free!)) is not available on your computer";
      echo "obabel might write out a warning since atom types cannot be translated/interpreted correctly.";
      echo "The obabel warnings can be confidently disregarded. They don't affect the docking.";
else
      echo "antechamber is used to create reliable SYBYL atom types in temp.mol2 file: --> temp_AtomTypesFixed.mol2";
      echo "antechamber -i ${TEMP_FILE} -fi mol2 -o "${current_path}/temp_AtomTypesFixed.mol2" -fo mol2 -at sybyl;"
      time antechamber -i ${TEMP_FILE} -fi mol2 -o "${current_path}/temp_AtomTypesFixed.mol2" -fo mol2 -at sybyl;
      if ( -e ${current_path}/temp_AtomTypesFixed.mol2 ) then
          echo "Antechamber  Success !";
          set TEMP_FILE = "${current_path}/temp_AtomTypesFixed.mol2";
      else
          echo "Antechamber  Failed !";
      endif
endif


echo "obabel -i mol2 temp_AtomTypesFixed.mol2 -o mopin -O temp.ZmatMOPAC"
${OBABELEXE} -i mol2 ${TEMP_FILE} -o mopin -O ${current_path}/temp.ZmatMOPAC

# prepare the use of python:
if ! $?PYTHONPATH then
    setenv PYTHONPATH $DOCKBASE/ligand/common
else
    setenv PYTHONPATH "$DOCKBASE/ligand/common:$PYTHONPATH"
endif

# create AMSOL7.1 input files (SM5.42R calculations in water and hexadecane solvents) using a Z-matrix in MOPAC style:
#

# python $DOCKBASE/ligand/amsol/make_amsol71_input.py ${current_path}/temp.ZmatMOPAC ${MoleculeName}
perl $DOCKBASE/ligand/amsol/make_amsol71_input.pl ${TEMP_FILE} ${MoleculeName}

# run the AMSOL7.1 calculations:

echo " running AMSOL7.1: SM5.42R (in water solvent) "
timeout $amsoltlimit $SHELL -c "$AMSOLEXE < temp.in-wat > temp.o-wat"
if ( $status ) then
   echo "AMSOL water calculation failed or stalled. Aborting"
   exit -1
endif

echo " running AMSOL7.1: SM5.42R (in water solvent) "
timeout $amsoltlimit $SHELL -c "$AMSOLEXE < temp.in-hex > temp.o-hex"
if ( $status ) then
   echo "AMSOL hexadecane calculation failed or stalled. Aborting"
   exit -1
endif

echo " extract data from AMSOL7.1 water and hexadecane output files:"
echo " starting process_amsol_mol2.py :"
echo $mol2file

python3.6 $DOCKBASE/ligand/amsol/process_amsol_mol2.py ${current_path}/temp.o-wat ${current_path}/temp.o-hex ${TEMP_FILE} ${current_path}/output

echo "process_amsol_mol2.py has finished."


if ( -e ${current_path}/output.solv ) then
    echo "AMSOL  Success !";
else

    python $DOCKBASE/ligand/amsol/make_amsol71_input.py ${current_path}/temp.ZmatMOPAC ${MoleculeName}

    echo " running AMSOL7.1: SM5.42R (in water solvent) "
    timeout $amsoltlimit $SHELL -c "$AMSOLEXE < temp.in-wat > temp.o-wat"
    if ( $status ) then
       echo "AMSOL water calculation failed or stalled. Aborting"
       exit -1
    endif

    echo " running AMSOL7.1: SM5.42R (in water solvent) "
    timeout $amsoltlimit $SHELL -c "$AMSOLEXE < temp.in-hex > temp.o-hex"
    if ( $status ) then
       echo "AMSOL hexadecane calculation failed or stalled. Aborting"
       exit -1
    endif

    echo " extract data from AMSOL7.1 water and hexadecane output files:"
    echo " starting process_amsol_mol2.py :"
    echo $mol2file

    python3.6 $DOCKBASE/ligand/amsol/process_amsol_mol2.py ${current_path}/temp.o-wat ${current_path}/temp.o-hex ${TEMP_FILE} ${current_path}/output
endif



# Ensure (temporary) but expected files are correct
mv -v ${TEMP_FILE} temp-working.mol2
cp -v temp-working.mol2 temp.mol2

#cp temp_AtomTypesFixed.mol2 temp.mol2

# Remove intermediate files
rm fort.* 
rm temp.*
rm *.log
rm temp-working.mol2
rm A*

