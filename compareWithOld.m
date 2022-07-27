load('dataOld_BLSA5800_cut')

% variable name: combine
% first row: header name
% old value

combine(end+1,:) = ColValues;
%%
check = zeros(length(ColValues),1);

for i = 1:length(ColValues)
    if isequal(combine{2,i},ColValues{i})
        check(i) = 1;
    elseif isnan(combine{2,i}) && isnan(ColValues{i})
        check(i) = 1;
    else
        check(i) = 0;
    end
end
combineNotEqual = combine(:,check == 0);
headerNotEqual = ColHeader(check == 0)';