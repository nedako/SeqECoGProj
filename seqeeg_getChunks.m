function CMB = seqeeg_getChunks(MaxChnkLeng,NumDuo,NumTrio,NumQuadro, GenGroupCLADist)
% MaxChnkLeng = 4;
% NumDuo      = 2;
% NumTrio     = 2;
% NumQuadro   = 2;
% GenGroupCLADist = 0 / 1  only set to 1 if you want to regenrate the group
% chunk arrangements, but they already are in the CMB strcuture
baseDir = ('/Users/nedakordjazi/Documents/SeqECoG/analyze');
load([baseDir , '/CMB.mat']);

%% Generate all possible doubles and triples and quadruples form the 1:5 digits

CMB.fing = [1:5]';

% ---------- comb2 all possible 2 finger combinations
CMB.comb2 = nchoosek(CMB.fing , 2);
CMB.comb2 =[CMB.comb2;fliplr(CMB.comb2);[CMB.fing CMB.fing]];

% ---------- comb3 all possible 3 finger combinations
CMB.comb3 = [];
for f = 1:length(CMB.fing)
    CMB.comb3 = [CMB.comb3 ; [CMB.comb2 CMB.fing(f)*ones(length(CMB.comb2) , 1)]];
end

% ---------- comb4 all possible 4 finger combinations
CMB.comb4 = [];
for f = 1:length(CMB.fing)
    CMB.comb4 = [CMB.comb4 ; [CMB.comb3 CMB.fing(f)*ones(length(CMB.comb3) , 1)]];
end

%% for each double list triples and quadruples that don't include that double
for  i = 1:length(CMB.comb2)
    ind1 = ismember(CMB.comb3(:,1:2) , CMB.comb2(i , :) , 'rows');
    ind2 = ismember(CMB.comb3(:,2:3) , CMB.comb2(i , :) , 'rows');
    ind = ~(ind1 + ind2);
    CMB.comb2NotIn3{i} = CMB.comb3(ind , :);
    
    ind1 = ismember(CMB.comb4(:,1:2) , CMB.comb2(i , :) , 'rows');
    ind2 = ismember(CMB.comb4(:,2:3) , CMB.comb2(i , :) , 'rows');
    ind3 = ismember(CMB.comb4(:,3:4) , CMB.comb2(i , :) , 'rows');
    ind = ~(ind1 + ind2 + ind3);
    CMB.comb2NotIn4{i} = CMB.comb4(ind , :);
end

%% for each triple list all quadruples that don't include that triple
for  i = 1:length(CMB.comb3)
    
    ind1 = ismember(CMB.comb4(:,1:3) , CMB.comb3(i , :) , 'rows');
    ind2 = ismember(CMB.comb4(:,2:4) , CMB.comb3(i , :) , 'rows');
    ind = ~(ind1 + ind2);
    CMB.comb3NotIn4{i} = CMB.comb4(ind , :);
    
end
%% generate all possible Chunk Length Arangements with 2, 3 and 4 chunks such that the total number of presses adds up to 14

CMB.ChnkLengthArngmnt(1,:) = [2 2 2 2 2 2 2]; % the most possible number of chunks

% ----------- 6chunks
c      = {} ;
k      = 6 ;
num = [2 : 4];
[c{1 : k}] = ndgrid(num);
c = cat(k + 1, c{k : -1 : 1});
c = reshape(c, [], k);
for i = 1:length(c)
    if sum(c(i,:)) == 14
        CMB.ChnkLengthArngmnt = [CMB.ChnkLengthArngmnt ; [c(i,:) 0 ]];
    end
end

% ----------- 5chunks
c      = {} ;
k      = 5 ;
num = [2 : 4];
[c{1 : k}] = ndgrid(num);
c = cat(k + 1, c{k : -1 : 1});
c = reshape(c, [], k);
for i = 1:length(c)
    if sum(c(i,:)) == 14
        CMB.ChnkLengthArngmnt = [CMB.ChnkLengthArngmnt ; [c(i,:) 0 0]];
    end
end

% ----------- 4chunks
c      = {} ;
k      = 4 ;
num = [2 : 4];
[c{1 : k}] = ndgrid(num);
c = cat(k + 1, c{k : -1 : 1});
c = reshape(c, [], k);
for i = 1:length(c)
    if sum(c(i,:)) == 14
        CMB.ChnkLengthArngmnt = [CMB.ChnkLengthArngmnt ; [c(i,:) 0 0 0]];
    end
end

%% Create the "chunk placement" vector and the "indeices of first chunk presses" vector for all possible CLAs

% -------------------- first press is always first CB
CMB.ChnkPlcmnt = [ones(length(CMB.ChnkLengthArngmnt) , 1) 2*ones(length(CMB.ChnkLengthArngmnt) , 1) zeros(length(CMB.ChnkLengthArngmnt) , 12)];
CMB.FrstChnkPressIndex = [ones(length(CMB.ChnkLengthArngmnt) , 1) zeros(length(CMB.ChnkLengthArngmnt) , 6)]; % indices for the first chunks presses
for i  = 1:length(CMB.ChnkLengthArngmnt)
    for j = 2 : length(find(CMB.ChnkLengthArngmnt(i,:)))
        CMB.FrstChnkPressIndex(i,j) = CMB.FrstChnkPressIndex(i,j-1) + sum(CMB.ChnkLengthArngmnt(i,j));
        CMB.ChnkPlcmnt(i , CMB.FrstChnkPressIndex(i,j)) = 1;
        CMB.ChnkPlcmnt(i , CMB.FrstChnkPressIndex(i,j) + 1) = 2;
    end
end
for i = 1:length(CMB.ChnkPlcmnt)
    for j = 3:length(CMB.ChnkPlcmnt(i,:))
        if (~CMB.ChnkPlcmnt(i,j) & CMB.ChnkPlcmnt(i,j-1) == 2)
            CMB.ChnkPlcmnt(i , j) = 3;
        elseif (~CMB.ChnkPlcmnt(i,j) & CMB.ChnkPlcmnt(i,j-1) == 3)
            CMB.ChnkPlcmnt(i , j) = 4;
        end
    end
end
%% pool together all the CLAs that dont share any chunk boundries (not used)
for i = 1:length(CMB.FrstChnkPressIndex)
    CMB.UniqChnkLengthArngmnt{i}  = CMB.FrstChnkPressIndex(i,:);
    iUCLA{i} = i;
    for j = 1:length(CMB.FrstChnkPressIndex)
        if j~=i
            if sum(sum(ismember(CMB.FrstChnkPressIndex(j,2:length(find(CMB.FrstChnkPressIndex(j,:))))  , CMB.UniqChnkLengthArngmnt{i})))== 0
                CMB.UniqChnkLengthArngmnt{i}  = [CMB.UniqChnkLengthArngmnt{i} ; CMB.FrstChnkPressIndex(j,:)];
                iUCLA{i} = [iUCLA{i} ; j];
            end
        end
    end
end

%% Generate the 6 chunks
%=========================================================================
%=========================================================================
%=========================================================================

% Randomly pick out a couple of doubles - make sure they are not repetition of the same finger
clear cmb2
rng shuffle % creates a different seed each time
if NumDuo
    CMB.Chunks{2,1} = zeros(NumDuo , MaxChnkLeng);
    temp2 = CMB.comb2;
    for n = 1:NumDuo
        cmb2(n) = randi([1,length(temp2)] , 1);
        % dont allow equal chunks, mirror chunks, or repetition of the same finger chunks
        while CMB.comb2(cmb2(n) , 1) == CMB.comb2(cmb2(n) , 2)
            cmb2(n) = randi([1,length(temp2)], 1);
        end
        ind = ~ismember(temp2 , temp2(cmb2(n) , :) , 'rows');
        CMB.Chunks{2,1}(n,1:2) = temp2(cmb2(n) , :) ;
        temp2 = temp2(ind,:);
    end
    
else
    cmb2 = [];
    CMB.Chunks{2,1} = [];
end

%=========================================================================
%=========================================================================
%=========================================================================
% Randomly pick out a couple of triplets - make sure they dont contain duos

if ~isempty(cmb2)
    temp3 = CMB.comb2NotIn3{cmb2(1)};
    for j  = 2:NumDuo
        ind  = ismember(temp3,CMB.comb2NotIn3{cmb2(j)} , 'rows');
        temp3 = temp3(ind , :);
    end
else
    temp3 = CMB.comb3;
end

clear m
if NumTrio
    CMB.Chunks{3,1} = zeros(NumTrio , MaxChnkLeng);
    for i = 1:NumTrio
        m(i) = randi([1,length(temp3)] , 1);
        % dont allow equal chunks, mirror chunks, or repetition of the same finger chunks
        while length(find(temp3(m(i) , :) == temp3(m(i) , 1))) == 3
            m(i) = randi([1,length(temp3)], 1);
        end
        CMB.Chunks{3,1}(i , 1:3) = temp3(m(i) , :);
        cmb2 = [cmb2  find(ismember(CMB.comb2 , temp3(m(i) , 1:2) , 'rows'))];
        cmb2 = [cmb2  find(ismember(CMB.comb2 , temp3(m(i) , 2:3) , 'rows'))];
        clear temp3
        temp3{1} = CMB.comb2NotIn3{cmb2(1)};
        for j  = 2:length(cmb2)
            temp3{2} = CMB.comb2NotIn3{cmb2(j)};
            ind = ismember(temp3{1} , temp3{2} , 'rows');
            temp3{1} = temp3{1}(ind , :);
        end
        temp3 = temp3{1}; % all the trios that dont contain BOTH duos
    end
else
    CMB.Chunks{3,1}  = [];
end

% cmb3 = find(ismember(CMB.comb3 , CMB.Chunks{3,1}(:,1:3) , 'rows')); % indices for the trios of interest in CMB.comb3
%=========================================================================
%=========================================================================
%=========================================================================
% Randomly pick out a couple of Quadro - make sure they dont contain duos or the trios

% All the unique quadros that don't contain neither of the duos
if ~isempty(cmb2)    
    temp4 = CMB.comb2NotIn4{cmb2(1)};
    for j  = 2:length(cmb2)
        ind = ismember(temp4 , CMB.comb2NotIn4{cmb2(j)} , 'rows');
        temp4 = temp4(ind , :);
    end
else
    temp4 = CMB.comb4;
end
% make sure that temp4 also does not contain Both the Duos, as well as the
% double transitions in the trios
% this isthe sufficient condition that it also would not contain any of the
% trios so no need to check for that

% now choose two quadros randomly
if NumQuadro
    CMB.Chunks{4,1} = zeros(NumQuadro , MaxChnkLeng);
    clear m
    for i = 1:NumQuadro
        m(i) = randi([1,length(temp4)] , 1);
        % dont allow equal chunks, mirror chunks, or repetition of the same finger chunks
        while (length(find(temp4(m(i) , :) == temp4(m(i) , 1))) == 4 | sum((temp4(m(i) , :) == mode(temp4(m(i) , :)))) >2 |  ((temp4(m(i) , 1) == temp4(m(i) , 2) & temp4(m(i) , 3)== temp4(m(i) , 4))) | ((temp4(m(i) , 1) == temp4(m(i) , 3) & temp4(m(i) , 2)== temp4(m(i) , 4))))
            m(i) = randi([1,length(temp4)], 1);
        end
        CMB.Chunks{4,1}(i , 1:4) = temp4(m(i) , :);
        cmb2 = [cmb2  find(ismember(CMB.comb2 , temp4(m(i) , 1:2) , 'rows'))];
        cmb2 = [cmb2  find(ismember(CMB.comb2 , temp4(m(i) , 2:3) , 'rows'))];
        cmb2 = [cmb2  find(ismember(CMB.comb2 , temp4(m(i) , 3:4) , 'rows'))];
        clear temp4
        temp4{1} = CMB.comb2NotIn4{cmb2(1)};
        for j  = 2:length(cmb2)
            temp4{2} = CMB.comb2NotIn4{cmb2(j)};
            ind = ismember(temp4{1} , temp4{2} , 'rows');
            temp4{1} = temp4{1}(ind , :);
        end
        temp4 = temp4{1}; % all the trios that dont contain BOTH duos
    end
else
    CMB.Chunks{4,1} = [];
end


cell2mat(CMB.Chunks)
%% Generate a 6-sequences group of CLAs for each of the two groups
% (you only do this once and then the CLAs remain the same, but the sequences change)
% first pick put asubset of CLAs for each of the two groups (we want 6 for each group)
% first exclude any CLA that has more than 3 duos, trios, or qudros from selection
if GenGroupCLADist
    I2 = sum(CMB.ChnkPlcmnt == 2 , 2);
    I3 = sum(CMB.ChnkPlcmnt == 3 , 2);
    I4 = sum(CMB.ChnkPlcmnt == 4 , 2);
    
    DesInd  = ismember(I4 , [1 2]) & ismember(I2 - I3 , [2 3]);
    CMB.DesiredCnhkargmnt = CMB.ChnkLengthArngmnt(DesInd , :);
    
    
    % option 1 : Pick out the CLAs with the least correlation (not used)
    % allind = nchoosek([1:length(DesInd)] , 6);
    % ChnkPlcmnt = CMB.ChnkPlcmnt(DesInd , :);
    % for i  = 1:length(allind)
    %     a = corr(ChnkPlcmnt(allind(i,:) , :)');
    %     % set diag(a) to 0 and get the upper traingle
    %     CLA_corr(i) = sum(squareform(a - diag(ones(length(a) , 1))));
    % end
    % [a , b] = sort(CLA_corr);
    %
    % ind = b((a == min(a)));
    % CMB.Group1Chnks = ChnkPlcmnt(allind(ind(1),:) , :);
    % i  = 2;
    % while i <= length(ind)
    %     if sum(ismember(allind(ind(i),:), allind(ind(1),:))) == 0
    %         CMB.Group2Chnks = ChnkPlcmnt(allind(ind(i),:) , :);
    %         break
    %     end
    %     i = i +1;
    % end
    % CMB.Group1Chnks
    % CMB.Group2Chnks
    % corr(CMB.Group1Chnks')
    
    % option 2 : create a distance matrix for all the possibel CLA pairs
    % create all the possible 6-group CLAs
    % Pick out the 6-g CLAs with the most accumulaitive distance
    allind = nchoosek(find(DesInd) , 6);
    DistInd = nchoosek([1:6] , 2);
    for i  = 1:length(allind)
        a = CMB.ChnkPlcmnt(allind(i,:),:);
        for k = 1:length(DistInd)
            L(k) = (a(DistInd(k,1),:) - a(DistInd(k,2),:)) * (a(DistInd(k,1),:) - a(DistInd(k,2),:))';
        end
        % set diag(a) to 0 and get the upper traingle
        CMB.CLA_dist(i) = sum(L);
        clear L
        i
    end
    
    
    [a , b] = sort(CMB.CLA_dist,'descend');
    ind = b(a == max(a));
    i = 1;
    % while i <= length(b)
    %     if (length(find(reshape(CMB.DesiredCnhkargmnt(allind(b(i),:) , :) , 42,1))) >=28)
    %         CMB.Group1Chnks = ChnkPlcmnt(allind(b(i),:) , :);
    %         break
    %     end
    %     i = i +1;
    % end
    
    CMB.Group1Chnks = CMB.ChnkPlcmnt(allind(b(1),:) , :);
    CMB.Group1 = allind(b(i),:);
    CMB.DesiredCnhkargmnt(allind(b(i),:) , :)
    
    % for group 2 pich the CLAs that chare any of the CLAs in group 1
    i = 1;
    while i <= length(b)
        if (sum(ismember(allind(b(i),:), allind(b(1),:))) == 0)
            CMB.Group2Chnks = CMB.ChnkPlcmnt(allind(b(i),:) , :);
            break
        end
        i = i +1;
    end
    CMB.Group2 = allind(b(i),:);
    CMB.DesiredCnhkargmnt(allind(b(i),:) , :)
    % CMB.Group1Chnks
    % CMB.Group2Chnks
    CMB.DesiredCnhkargmnt(CMB.Group1 , :)
    CMB.DesiredCnhkargmnt(CMB.Group2 , :)
end

%%


