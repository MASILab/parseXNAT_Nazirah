load('olddata_5802.mat')

% firstcol = olddata(:,1);
% Index = strfind(firstcol,'BLSA_5802');
% 
% for i = 1:length(Index)
%     if ~isempty(Index{i})
%         Index2 = i;
%     end
% end
% combine = [data(1,:); newstats(Index2,:); data(2,:)];
% variable name: combine
% first row: header name
% old value
%%
combines = [ColHeader;olddata(2,:);ColValues];

check = zeros(length(combines),1);

for i = 1:length(combines)
    if isequal(combines{2,i},combines{3,i})
        check(i) = 1;
    elseif isnan(combines{2,i}) && isnan(combines{3,i})
        check(i) = 1;
    else
        check(i) = 0;
    end
end
combineNotEqual = combines(:,check == 0);
headerNotEqual = ColHeader(check == 0)';