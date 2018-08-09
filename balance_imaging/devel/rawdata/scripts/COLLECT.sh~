#!/bin/bash

#This file collects all BF data from each run in model_output/
#  and then collates it in model_output/I211_J211.dat, etc
#  where each run is a column in the file

#!/bin/bash

EXEC_FILE='exec/collect.x'
PARAM_FILE='script_model_parameters.dat'
NAMES_FILE='script_model_names.dat'

SOURCE_FOLDER=''
DEST_FOLDER=''
LINES=0
RUNS=0
FILES=0
while read NAME VALUE; do
   if [ ${NAME} == 'SOURCE_FOLDER' ]
   then
	SOURCE_FOLDER=${VALUE}
   elif [ ${NAME} == 'DEST_FOLDER' ]
   then
       DEST_FOLDER=${VALUE}
   elif [ ${NAME} == 'LINES' ]
   then
       LINES=${VALUE}
   elif [ ${NAME} == 'RUNS' ]
   then 
       RUNS=${VALUE}
   elif [ ${NAME} == 'FILES' ]
   then
       FILES=${VALUE}
   else
       echo "Variable unknown: $NAME $VALUE"
   fi
done < <(egrep -v '^(#|$)' $PARAM_FILE)

if [ $1 ]
then
    SOURCE_FOLDER=$1
    if [ $2 ]
    then
	DEST_FOLDER=$2
    fi
fi

while read NAME; do
   for((j=0;j<$RUNS;j++))
   do
       if [ ! -d $NAME ]
       then
	   mkdir -v $NAME
       fi
       RUN=` printf "${SOURCE_FOLDER}run%04d/${NAME}" $j `
       NEWFILE=` printf "${NAME}/run%04d" $j `   
       egrep -v "^(#|$)" $RUN > $NEWFILE
   done
   DEST="$DEST_FOLDER$NAME"
   ./$EXEC_FILE $RUNS $NAME $DEST
   rm -r $NAME
done < <(egrep -v "^(#|$)" $NAMES_FILE)