# move files from main to structure

datapath=/nfs/masi/BLSA
movedpath=~/data_from_NFS

if [ ! -d "$movedpath" ]; then
	mkdir $movedpath
fi


SESSIONS=`ls -d /nfs/masi/BLSA/BLSA_4*`

total_subject=0

for DATAFOLDER in $SESSIONS
do
	dti=($DATAFOLDER/ASSESSORS/*dti*)
	slant=($DATAFOLDER/ASSESSORS/*slant*)
	#echo "dti: $dti"
	#echo "slant: $slant"

	SESSION=${DATAFOLDER##*/nfs/masi/BLSA/}
	SUBJ_ID=${SESSION::9}
	
	# if both folders are there, then proceed
	if [ -d "${dti[0]}" ] && [ -d "${slant[0]}" ]; then

		echo "Found both folders for $SUBJ_ID"

		#if subject directory not exist, create one
		if [ ! -d "$SUBJ_ID" ]; then
			echo "--- Dir $SUBJ_ID not exist, creating one..."
			mkdir $movedpath/$SUBJ_ID/

			((total_subject+=1))
		fi
		
		mkdir $movedpath/$SUBJ_ID/$SESSION

		#once have a folder, move subject's session to the folder
		echo "-- Copying DWI data..."
		#cp -r data/$DATAFOLDER data/$SUBJ_ID/$DATAFOLDER
		cp $dti/PREPROCESSED/dwmri.nii.gz $movedpath/"$SUBJ_ID"/"$SESSION"/data.nii.gz
		cp $dti/PREPROCESSED/dwmri.bval $movedpath/"$SUBJ_ID"/"$SESSION"/bvals
		cp $dti/PREPROCESSED/dwmri.bvec $movedpath/"$SUBJ_ID"/"$SESSION"/bvecs
		cp $dti/PREPROCESSED/mask.nii.gz $movedpath/"$SUBJ_ID"/"$SESSION"/mask.nii.gz

		echo "-- Copying tensor data..."
		cp $dti/SCALARS/dwmri_tensor_ad.nii.gz $movedpath/"$SUBJ_ID"/"$SESSION"/ad.nii.gz
		cp $dti/SCALARS/dwmri_tensor_md.nii.gz $movedpath/"$SUBJ_ID"/"$SESSION"/md.nii.gz
		cp $dti/SCALARS/dwmri_tensor_rd.nii.gz $movedpath/"$SUBJ_ID"/"$SESSION"/rd.nii.gz
		cp $dti/SCALARS/dwmri_tensor_fa.nii.gz $movedpath/"$SUBJ_ID"/"$SESSION"/fa.nii.gz

		echo "-- Copying T1_seg data..."
		cp $slant/SEG/T1_seg.nii.gz $movedpath/"$SUBJ_ID"/"$SESSION"/T1_seg.nii.gz
		cp $slant/STATS/T1_label_volumes.txt $movedpath/"$SUBJ_ID"/"$SESSION"/T1_label_volumes.txt

		# copy MPRAGE
		mprage=($DATAFOLDER/SCANS/MPRAGE/NIFTI/*MPRAGE*)
		echo "-- Copying MPRAGE file..."
		cp $mprage $movedpath/"$SUBJ_ID"/"$SESSION"/MPRAGE.nii.gz


	fi

	if [ -d "${dti[0]}" ] && [ -d "${slant[0]}" ]; then
		echo "!! Zero or only 1 folder found: $DATAFOLDER"
		echo "$DATAFOLDER" >> session_not_included.txt
	fi

done

echo "Done moving data. Total subjects moved is $total_subject"