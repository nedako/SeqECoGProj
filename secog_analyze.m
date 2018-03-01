function Dout=secog_analyze (what , varargin)
getdat = 1;
c = 1;
while(c<=length(varargin))
    switch(varargin{c})
        case {'subjCode'}
            % in case of what = 'sing_subj', the subjCode need to be specified
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'getdat'}
            % 1 if you want to read the dat/MOV files (default)
            % 0 if you want to just modify an already provided Dall
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'Dall'}
            % in case of getdat = 0, Dall needs to be provided
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        otherwise
            error(sprintf('Unknown option: %s',varargin{c}));
    end
end
if getdat == 0 & ~exist('Chan2Plot')
    error('Provide Dall')
end
%%
Dout = [];
prefix = 'ecog1_';
baseDir = '/Users/nkordjazi/Documents/SeqECoG/analyze';
%baseDir = '/Users/nkordjazi/Documents/SeqEye/se1/SeqEye1/se1_data/analyze';
% subj_name = {'XW' , 'ML' , 'DS' , 'BM' , 'HK' , 'BW'};
subj_name = {'P2' , 'P4'};
% the blocks per day are different for every participnat
BlPerDay(1).Day = {[1:12] , [13: 25] , [26:39] , [40:44]};
BlPerDay(2).Day = {[1:12] , [13: 26] , [27:38]};
load([baseDir , '/CMB.mat'])
%load([baseDir , '/se1_all.mat'])
% D   = load('/Users/nedakordjazi/Documents/SeqEye/SequenceHierarchical/Analysis/sh3_avrgPattern.mat');
% MTW = 2 - D.MT(126:end)/max(D.MT(126:end));
% possibleDuo = D.Sequence(126:end , 1:2);
% create an emty structure with the same fields as Dall (empty b/c isError can never be 2)
window = 12;
switch what
    case 'all_subj'
        for i=1:length(subj_name)
            clear ANA
            if getdat
                ANA = secog_subj(subj_name{i} , 0);
            else
                ANA = getrow(Dall , Dall.SN == i);
            end
            
            ANA.Day = zeros(size(ANA.BN));
            for d = 1:length(BlPerDay(i).Day)
                ANA.Day(ismember(ANA.BN , BlPerDay(i).Day{d})) = d;
            end
            ANA.SN(1:length(ANA.BN) , :) = i;
            ChnkArrng = zeros(2,7);
            Chnkplcmnt = zeros(2,7);
            
            ChnkArrang = CMB.DesiredCnhkargmnt;
            for j = 1:size(ChnkArrang , 1)
                temp = [];
                temp1 = [];
                for k = 1:length(find(ChnkArrang(j , :)))
                    temp = [temp k*ones(1,ChnkArrang(j,k))];
                    temp1 = [temp1 1:ChnkArrang(j,k)];
                end
                ChnkArrng(j , :) = temp;
                ChnkPlcmnt(j,:)  = temp1;
                
                ICP(j , :) = diff(ChnkPlcmnt(j,:),1,2);
                ICP(j , (ICP(j , :)<0)) = 0;
                for k = 1:length(ICP(j , :))-1
                    if ICP(j,k) & ICP(j,k+1)
                        ICP(j,k+1) = ICP(j,k) + ICP(j,k+1);
                    end
                end
            end
            IPIarrangement = diff(ChnkArrng , 1 , 2);   % between IPIs will have the index of 1
            IPIarrangement(~IPIarrangement) = 2;         % within IPIs will have the index of 2
            
            ANA.IPI = diff(ANA.AllPressTimes,1,2);
            ANA.badPress = ANA.AllPress ~= ANA.AllResponse;
            
            
            ANA.ChnkArrang     = zeros(size(ANA.AllPressTimes));
            ANA.ChnkPlcmnt     = zeros(size(ANA.AllPressTimes));
            ANA.IPIarrangement = zeros(size(ANA.IPI));
            ANA.IPIChnkPlcmntArr  = zeros(size(ANA.IPI));
            for cl = 1:4
                switch cl
                    case {1 2}
                        ANA.ChnkArrang(ANA.seqNumb == cl , :) = repmat(ChnkArrng(1 , :) , length(find(ANA.seqNumb == cl)) , 1);
                        ANA.ChnkPlcmnt(ANA.seqNumb == cl , :) = repmat(ChnkPlcmnt(1 , :) , length(find(ANA.seqNumb == cl)) , 1);
                        ANA.IPIarrangement(ANA.seqNumb == cl , :) = repmat(IPIarrangement(1 , :) , length(find(ANA.seqNumb == cl)) , 1);
                        ANA.IPIChnkPlcmntArr(ANA.seqNumb == cl , :) = repmat(ICP(1 , :) , length(find(ANA.seqNumb == cl)) , 1);
                    case {3 4}
                        ANA.ChnkArrang(ANA.seqNumb == cl , :) = repmat(ChnkArrng(2 , :) , length(find(ANA.seqNumb == cl)) , 1);
                        ANA.ChnkPlcmnt(ANA.seqNumb == cl , :) = repmat(ChnkPlcmnt(2 , :) , length(find(ANA.seqNumb == cl)) , 1);
                        ANA.IPIarrangement(ANA.seqNumb == cl , :) = repmat(IPIarrangement(2 , :) , length(find(ANA.seqNumb == cl)) , 1);
                        ANA.IPIChnkPlcmntArr(ANA.seqNumb == cl , :) = repmat(ICP(2 , :) , length(find(ANA.seqNumb == cl)) , 1);
                end
            end
            
            uBN = unique(ANA.BN);
            ANA.Rep = repmat([1;2] , 0.5* length(ANA.BN ), 1);
            
            ANA.isgood = ones(length(ANA.TN) , 1);
            for tn = 1:length(ANA.TN)
                
                switch ANA.seqNumb(tn,1)
                    case {0,1,2,3,4,11,22,33,44,55}
                        ANA.seqlength(tn,1) = 7;
                    case {103 , 203}
                        ANA.seqlength(tn,1) = 3;
                    case {104 , 204}
                        ANA.seqlength(tn,1) = 4;
                    otherwise
                        ANA.seqlength(tn,1) = NaN;
                end
                
                % create a smooth press time series
                presses = [0 : length(find(ANA.AllPress(tn , :)))+1];
                chunks  = [1 , ANA.ChnkPlcmnt(tn , :) , 1];
                ANA.PressPressVelocity(tn , :) = NaN * ones(1,14);
                
                
                
                % Create a smooth imposed-chunk timeseries
                ANA.IPIChnkPlcmnt(tn , :) = [NaN NaN NaN NaN];
                if ismember(ANA.seqNumb(tn), [1:4])   % Intermixed or CLAT
                    
                    ANA.IPIwithin(tn,1)  = nanmean(ANA.IPI(tn,ANA.IPIarrangement(tn,:) == 2));
                    ANA.IPIbetween(tn,1) = nanmean(ANA.IPI(tn,ANA.IPIarrangement(tn,:) == 1));
                    
                    ANA.IPIChnkPlcmnt(tn , 1)= nanmean(ANA.IPI(tn,ANA.IPIChnkPlcmntArr(tn,:) == 0)); % basically between
                    ANA.IPIChnkPlcmnt(tn , 2)= nanmean(ANA.IPI(tn,ANA.IPIChnkPlcmntArr(tn,:) == 1));
                    ANA.IPIChnkPlcmnt(tn , 3)= nanmean(ANA.IPI(tn,ANA.IPIChnkPlcmntArr(tn,:) == 2));
                    ANA.IPIChnkPlcmnt(tn , 4)= nanmean(ANA.IPI(tn,ANA.IPIChnkPlcmntArr(tn,:) == 3));
                    
                else
                    ANA.IPIwithin(tn,1)  = NaN;
                    ANA.IPIbetween(tn,1) = NaN;
                end
                
                
                % estimate chunk bounries as the ones that are 20% of a std above the rest of the IPIs
                
                thresh = .3 * std(ANA.IPI(tn , :));
                [dum , estChnkBndry] = findpeaks(ANA.IPI(tn , :));% ,'MinPeakProminence', thresh); % slope of IPIs
                
                if ~isempty(estChnkBndry)
                    goodpeak = ones(1, length(estChnkBndry));
                    for cb = 1:length(estChnkBndry)
                        if ANA.IPI(tn , estChnkBndry(cb)) < nanmean(ANA.IPI(tn  ,:)) + thresh
                            goodpeak(cb) = 0;
                        end
                    end
                    if sum(goodpeak)
                        estChnkBndry = estChnkBndry(logical(goodpeak));
                    else
                        estChnkBndry = [];
                    end
                end
                
                
                ANA.estChnkBndry(tn , :) = zeros(1, 14);  % first presses of chunks will be 1
                ANA.estChnkPlcmnt(tn , :) = zeros(1, 14);  % chunk placements of digits
                ANA.estChnkBndry(tn , estChnkBndry+1) = 1;
                ANA.estChnkBndry(tn , 1) = 1;
                dum = find(~ANA.estChnkBndry(tn,:));
                
                if sum(ANA.estChnkBndry(tn,:))>0
                    ANA.estChnkPlcmnt(tn,logical(ANA.estChnkBndry(tn,:))) = 1;
                    ANA.estChnkPlcmnt(tn,1) = 1;
                    dum = find(~ANA.estChnkPlcmnt(tn,:));
                    for h = 1:length(dum)
                        ANA.estChnkPlcmnt(tn,dum(h)) = ANA.estChnkPlcmnt(tn,dum(h)-1) + 1;
                    end
                end
            end
            
            Dout = addstruct(Dout , ANA);
            
        end
        
    case 'sing_subj'
        i  = find(strcmp(subj_name , subjCode));
        clear ANA
        if getdat
            ANA = secog_subj(subj_name{i} , 0);
        else
            ANA = getrow(Dall , Dall.SN == i);
        end
        
        ANA.Day = zeros(size(ANA.BN));
        for d = 1:length(BlPerDay(i).Day)
            ANA.Day(ismember(ANA.BN , BlPerDay(i).Day{d})) = d;
        end
        ANA.SN(1:length(ANA.BN) , :) = i;
        ChnkArrng = zeros(2,7);
        Chnkplcmnt = zeros(2,7);
        
        ChnkArrang = CMB.DesiredCnhkargmnt;
        for j = 1:size(ChnkArrang , 1)
            temp = [];
            temp1 = [];
            for k = 1:length(find(ChnkArrang(j , :)))
                temp = [temp k*ones(1,ChnkArrang(j,k))];
                temp1 = [temp1 1:ChnkArrang(j,k)];
            end
            ChnkArrng(j , :) = temp;
            ChnkPlcmnt(j,:)  = temp1;
            
            ICP(j , :) = diff(ChnkPlcmnt(j,:),1,2);
            ICP(j , (ICP(j , :)<0)) = 0;
            for k = 1:length(ICP(j , :))-1
                if ICP(j,k) & ICP(j,k+1)
                    ICP(j,k+1) = ICP(j,k) + ICP(j,k+1);
                end
            end
        end
        IPIarrangement = diff(ChnkArrng , 1 , 2);   % between IPIs will have the index of 1
        IPIarrangement(~IPIarrangement) = 2;         % within IPIs will have the index of 2
        
        ANA.IPI = diff(ANA.AllPressTimes,1,2);
        ANA.badPress = ANA.AllPress ~= ANA.AllResponse;
        
        
        ANA.ChnkArrang     = zeros(size(ANA.AllPressTimes));
        ANA.ChnkPlcmnt     = zeros(size(ANA.AllPressTimes));
        ANA.IPIarrangement = zeros(size(ANA.IPI));
        ANA.IPIChnkPlcmntArr  = zeros(size(ANA.IPI));
        for cl = 1:4
            switch cl
                case {1 2}
                    ANA.ChnkArrang(ANA.seqNumb == cl , :) = repmat(ChnkArrng(1 , :) , length(find(ANA.seqNumb == cl)) , 1);
                    ANA.ChnkPlcmnt(ANA.seqNumb == cl , :) = repmat(ChnkPlcmnt(1 , :) , length(find(ANA.seqNumb == cl)) , 1);
                    ANA.IPIarrangement(ANA.seqNumb == cl , :) = repmat(IPIarrangement(1 , :) , length(find(ANA.seqNumb == cl)) , 1);
                    ANA.IPIChnkPlcmntArr(ANA.seqNumb == cl , :) = repmat(ICP(1 , :) , length(find(ANA.seqNumb == cl)) , 1);
                case {3 4}
                    ANA.ChnkArrang(ANA.seqNumb == cl , :) = repmat(ChnkArrng(2 , :) , length(find(ANA.seqNumb == cl)) , 1);
                    ANA.ChnkPlcmnt(ANA.seqNumb == cl , :) = repmat(ChnkPlcmnt(2 , :) , length(find(ANA.seqNumb == cl)) , 1);
                    ANA.IPIarrangement(ANA.seqNumb == cl , :) = repmat(IPIarrangement(2 , :) , length(find(ANA.seqNumb == cl)) , 1);
                    ANA.IPIChnkPlcmntArr(ANA.seqNumb == cl , :) = repmat(ICP(2 , :) , length(find(ANA.seqNumb == cl)) , 1);
            end
        end
        
        uBN = unique(ANA.BN);
        ANA.Rep = repmat([1;2] , 0.5* length(ANA.BN ), 1);
        
        ANA.isgood = ones(length(ANA.TN) , 1);
        for tn = 1:length(ANA.TN)
            
            switch ANA.seqNumb(tn,1)
                case {0,1,2,3,4,11,22,33,44,55}
                    ANA.seqlength(tn,1) = 7;
                case {103 , 203}
                    ANA.seqlength(tn,1) = 3;
                case {104 , 204}
                    ANA.seqlength(tn,1) = 4;
                otherwise
                    ANA.seqlength(tn,1) = NaN;
            end
            
            % create a smooth press time series
            presses = [0 : length(find(ANA.AllPress(tn , :)))+1];
            chunks  = [1 , ANA.ChnkPlcmnt(tn , :) , 1];
            ANA.PressPressVelocity(tn , :) = NaN * ones(1,14);
            
            
            
            % Create a smooth imposed-chunk timeseries
            ANA.IPIChnkPlcmnt(tn , :) = [NaN NaN NaN NaN];
            if ismember(ANA.seqNumb(tn), [1:4])   % Intermixed or CLAT
                
                ANA.IPIwithin(tn,1)  = nanmean(ANA.IPI(tn,ANA.IPIarrangement(tn,:) == 2));
                ANA.IPIbetween(tn,1) = nanmean(ANA.IPI(tn,ANA.IPIarrangement(tn,:) == 1));
                
                ANA.IPIChnkPlcmnt(tn , 1)= nanmean(ANA.IPI(tn,ANA.IPIChnkPlcmntArr(tn,:) == 0)); % basically between
                ANA.IPIChnkPlcmnt(tn , 2)= nanmean(ANA.IPI(tn,ANA.IPIChnkPlcmntArr(tn,:) == 1));
                ANA.IPIChnkPlcmnt(tn , 3)= nanmean(ANA.IPI(tn,ANA.IPIChnkPlcmntArr(tn,:) == 2));
                ANA.IPIChnkPlcmnt(tn , 4)= nanmean(ANA.IPI(tn,ANA.IPIChnkPlcmntArr(tn,:) == 3));
                
            else
                ANA.IPIwithin(tn,1)  = NaN;
                ANA.IPIbetween(tn,1) = NaN;
            end
            
            
            % estimate chunk bounries as the ones that are 20% of a std above the rest of the IPIs
            
            thresh = .3 * std(ANA.IPI(tn , :));
            [dum , estChnkBndry] = findpeaks(ANA.IPI(tn , :));% ,'MinPeakProminence', thresh); % slope of IPIs
            
            if ~isempty(estChnkBndry)
                goodpeak = ones(1, length(estChnkBndry));
                for cb = 1:length(estChnkBndry)
                    if ANA.IPI(tn , estChnkBndry(cb)) < nanmean(ANA.IPI(tn  ,:)) + thresh
                        goodpeak(cb) = 0;
                    end
                end
                if sum(goodpeak)
                    estChnkBndry = estChnkBndry(logical(goodpeak));
                else
                    estChnkBndry = [];
                end
            end
            
            
            ANA.estChnkBndry(tn , :) = zeros(1, 14);  % first presses of chunks will be 1
            ANA.estChnkPlcmnt(tn , :) = zeros(1, 14);  % chunk placements of digits
            ANA.estChnkBndry(tn , estChnkBndry+1) = 1;
            ANA.estChnkBndry(tn , 1) = 1;
            dum = find(~ANA.estChnkBndry(tn,:));
            
            if sum(ANA.estChnkBndry(tn,:))>0
                ANA.estChnkPlcmnt(tn,logical(ANA.estChnkBndry(tn,:))) = 1;
                ANA.estChnkPlcmnt(tn,1) = 1;
                dum = find(~ANA.estChnkPlcmnt(tn,:));
                for h = 1:length(dum)
                    ANA.estChnkPlcmnt(tn,dum(h)) = ANA.estChnkPlcmnt(tn,dum(h)-1) + 1;
                end
            end
        end
        
        Dout = addstruct(Dout , ANA);
        
end
