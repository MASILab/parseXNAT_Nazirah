%compare inv Affine matrix

function invCompare = compareInvAffine(inv1,inv2)

ori_path = fopen(inv1);
new_path = fopen(inv2);

aff_ori = fscanf(ori_path,'%f %f %f %f',[4 Inf]);
aff_new = fscanf(new_path,'%f %f %f %f',[4 Inf]);

aff_ori = aff_ori';
aff_new = aff_new';

if ~isequal(size(aff_new),size(aff_ori))
    invCompare = NaN;
    return;
end

for i = 1:size(aff_ori,1)
    for j = 1:size(aff_ori,2)
        diffpct(i,j) = abs((aff_ori(i,j) - aff_new(i,j))/aff_ori(i,j))*100;
    end
end

invCompare = diffpct;

fclose(ori_path);
fclose(new_path);