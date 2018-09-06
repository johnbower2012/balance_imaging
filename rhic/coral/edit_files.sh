#! /bin/sh

oldstring="ZSYM 1"
newstring="ZSYM true"
checkstring=${oldstring}
echo oldstring is ${oldstring}


    for fn in `find . -iname \*.dat`
    do
	if [ -s $fn ]
	then
	    if [ -e ${fn}.tmp ]
	    then
		rm -f ${fn}.tmp
	    fi
	    if [ -e ${fn}.temp ]
	    then
		rm -f ${fn}.temp
	    fi
	    echo $fn
	    cp $fn $fn.tmp
	    sed -e "s/${oldstring}/${newstring}/g" $fn > $fn.temp
	    mv -f $fn.temp $fn
	fi
    done
