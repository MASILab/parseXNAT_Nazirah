function verifyADRD(ad,rd,fa,D)
% ONLY HERE FOR TESTING, TO BE DELETED ONCE verifyADRD2 IS OK

if(length(dir(ad))>0)
    return;
end
disp('Re-computing AD&RD')
mkdir('../RD/')
mkdir('../AD/')

try
    FA = load_nii(fa);
catch
    FA = load_nii_gz(fa);
end

fp = fopen(D,'rb','ieee-be');
Ddat = fread(fp,inf,'double');
fclose(fp);
Ddat = reshape(Ddat, [8 size(FA.img)]);
Ddat = permute(Ddat,[2 3 4 1]);
Ddat = Ddat(:,:,:,3:8);
AD = FA; AD.img = 0*AD.img;
RD = FA; RD.img = 0*RD.img;
for i=1:size(FA.img,1)
    for j=1:size(FA.img,2)
        for k=1:size(FA.img,3)
            if(FA.img(i,j,k)>0)
                d = Ddat(i,j,k,:);
                dd = [d(1) d(2) d(3);d(2) d(4) d(5); d(3) d(5) d(6)];
                e = sort(eig(dd));
                AD.img(i,j,k) = e(3);
                RD.img(i,j,k) = (e(1)+e(2))/2;
            end
        end
    end
end

save_nii_gz(AD,ad);
save_nii_gz(RD,rd);
