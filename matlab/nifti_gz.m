function nifti_gz_out = nifti_gz(fname)

filename = ['/tmp/tmp' num2str(round(rand()*1e8)) '.nii'];
system(sprintf('gzip -dc %s > %s',fname,filename));
nifti_gz_out = nifti(filename);
delete(filename);
