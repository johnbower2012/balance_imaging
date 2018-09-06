#!/bin/bash

EXEC_FILE='exec/observables_model.x'
PARAM_FILE='script_model_parameters.dat'
NAMES_FILE='script_model_names.dat'

SOURCE_FOLDER=''
EXP_FOLDER=''
DEST_FOLDER=''
LINES=0
RUNS=0
FILES=0

while read NAME VALUE; do
   if [ ${NAME} == 'SOURCE_FOLDER' ]
   then
	SOURCE_FOLDER=${VALUE}
   elif [ ${NAME} == 'EXP_FOLDER' ]
   then
	EXP_FOLDER=${VALUE}
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
	EXP_FOLDER=$2
	if [ $3 ]
	then
	    DEST_FOLDER=$3
	fi
    fi
fi

declare -a array
i=0
while read name; do
   i=$(($i + 1))
   dest="$SOURCE_FOLDER$name"
   array[$i]="$dest"
done < <(egrep -v '^(#|$)' $NAMES_FILE)

echo "OBSERVABLES_MODEL.sh Running:"
echo $EXEC_FILE $SOURCE_FOLDER $EXP_FOLDER $DEST_FOLDER $LINES $RUNS $FILES ${array[1]} ${array[2]} ${array[3]} ${array[4]} ${array[5]} ${array[6]} ${array[7]}

./$EXEC_FILE $SOURCE_FOLDER $EXP_FOLDER $DEST_FOLDER $LINES $RUNS $FILES ${array[1]} ${array[2]} ${array[3]} ${array[4]} ${array[5]} ${array[6]} ${array[7]}
