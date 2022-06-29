#!/bin/bash

datapath=`ls data`
total_subject=0

for DATAFOLDER in $datapath
do
	SUBJ_ID=${DATAFOLDER::9}
	#echo $SUBJ_ID

	if [ ! -d "/path/to/dir" ]; then
		#if subject directory not exist, create one
		echo "Dir $SUBJ_ID not exist, creating one"
		mkdir $SUBJ_ID
		((total_subject+=1))
	fi
	
	#once have a folder, move subject's session to the folder
	mv data/$DATAFOLDER data/$SUBJ_ID/$DATAFOLDER
done

echo "Done moving data. Total subjects moved is $total_subject"