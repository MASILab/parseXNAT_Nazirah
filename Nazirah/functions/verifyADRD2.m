function verifyADRD2(ad,rd,fa,D)

if ~isempty(dir(ad))
    return;
end

disp('Re-computing AD&RD')

if not(isfolder('../AD'))
    mkdir('../RD/')
end

if not(isfolder('../AD'))
    mkdir('../AD/')
end

try
    FA = niftiread(fa);
    FAinfo = niftiinfo(fa);
catch
    FA = niftiread([fa '.gz']);
    FAinfo = niftiinfo([fa '.gz']);
end

fp = fopen(D,'rb','ieee-be');
Ddat = fread(fp,inf,'double');
fclose(fp);
Ddat = reshape(Ddat, [8 size(FA)]);
Ddat = permute(Ddat,[2 3 4 1]);
Ddat = Ddat(:,:,:,3:8);
AD = FA; AD = 0*AD;
RD = FA; RD = 0*RD;
for i=1:size(FA,1)
    for j=1:size(FA,2)
        for k=1:size(FA,3)
            if(FA(i,j,k)>0)
                d = Ddat(i,j,k,:);
                dd = [d(1) d(2) d(3);d(2) d(4) d(5); d(3) d(5) d(6)];
                e = sort(eig(dd));
                AD(i,j,k) = e(3);
                RD(i,j,k) = (e(1)+e(2))/2;
            end
        end
    end
end

% save updated AD,RD
% remove extension (not sure why this cause file to be saved as .nii.nii)
ad(strfind(ad,'.nii'):end) = [];
rd(strfind(rd,'.nii'):end) = [];
niftiwrite(AD,ad,FAinfo);
niftiwrite(RD,rd,FAinfo);
