%oldstat = readcell('/home/local/VANDERBILT/mohdkhn/Documents/BLSA_original_run/statsWithADRDVol_04Dec2019.csv');
%newstat = readcell('/home/local/VANDERBILT/mohdkhn/Documents/Test2-newdata/AllStats-HeaderWithADRDVol-20220810T101048.csv');

oldstat = readcell('/Users/nana/masi-42/Documents/BLSA_original_run/statsWithADRDVol_04Dec2019.csv');
newstat = readcell('/Users/nana/masi-42/Documents/Test2-newdata/AllStats-HeaderWithADRDVol-20220810T164447.csv');

%% check col is similar or not
oldheader = oldstat(1,:);
newheader = newstat(1,:);

[old2newcol,missing_pos] = getcolsfromold(oldheader,newheader);
old2stat = oldstat(:,old2newcol);
newstat(:,missing_pos) = [];

% double check now if the header is correct
diffheader = isequal(old2stat(1,:),newstat(1,:));

if ~diffheader
    error('Check back header')
end

if ~isequal(size(old2stat,2),size(newstat,2))
    error('Check back old and new stats')
end

[rows,cols] = size(old2stat);


%% start comparing
clear check %we need to overwrite. make sure to clear to avoid refill-in.
firstcol = old2stat(:,1);
firstcol1 = newstat(:,1);

%check = zeros(166,length(oldstat(1,:)));
combineAll = (newstat(1,:));

for j=1:length(firstcol1)
    
    stringAtFirstCol=firstcol1{j};
    
    % find location of string of newstat in oldstat
    Index = strcmp(firstcol,stringAtFirstCol);
    Index2 = find(Index);
    
    % if it's not in the oldstat check return
    if isempty(Index2)
        check{j,1} = firstcol1{j};
        check{j,2} = 'Not in old STAT';
    else
        
        combine1 = [old2stat(1,:); old2stat(Index2,:); newstat(j,:)];
        combineAll = [combineAll; old2stat(Index2,:); newstat(j,:)];
        
        for i = 1:cols            
            if j == 1 % first row: header title
                check{j,i} = newheader{i};
            elseif i == 1 % first col: session name
                check{j,i} = combine1{2,i};
            elseif isequal(combine1{2,i},combine1{3,i}) 
                % IF OLD = NEW, set to 1
                check{j,i} = 1;
                %%change the isnan(syms) here
                %         elseif ismissing(combine1{2,i}) | ismissing(combine1{3,i})
                %             check(j,i)=0;
                
            elseif isnan(combine1{2,i}) && isnan(combine1{3,i})
                % IF BOTH OLD & NEW ISNAN --> EQUAL --> set to 1
                check{j,i} = 1;             
            else
                % FIND ABS DIFFERENCE of old and new.
                check{j,i} = abs(combine1{2,i}-combine1{3,i});
            end
        end
    end
    
    
end


%% FUNCTIONS
function [old2newcol,missing_pos] = getcolsfromold(oldheader,newheader)

missing_pos = [];
old2newcol = [];

for i = 1:length(newheader)
    if ismissing(newheader{i})
        missing_pos = [missing_pos; i];
    else
        tmp = find(strcmp(oldheader,newheader{i}));
        
        if length(tmp) > 1
            fprintf('found 2 similar for header %s at col %d\n',newheader{i},i)
        elseif isempty(tmp)
            fprintf('found 0 similar for header %s at col %d\n',newheader{i},i)
        else
            old2newcol = [old2newcol;tmp];
        end
    end
end

end