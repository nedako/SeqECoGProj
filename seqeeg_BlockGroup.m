function seqeeg_BlockGroup(Pall , subjNum, what,varargin)


%% setup the defaults and deal with the varargin
DownsampleRate = 10;
NormType = 'stim';
NumWarpSampFast = 150;
NumWarpSampSlow = 300;
TimeDelay = 0.5; % sec
FreqRange = [4 184];
numFreqBins = 45;
subjname = {'P2', 'P4' , 'P5'};
c = 1;
saveDir = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' subjname{subjNum} , '/'] ;
while(c<=length(varargin))
    switch(varargin{c})
        case {'NumWarpSampChunk'}
            % Number of warping sample for chunks, fast single finger  Default = 200
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'NumWarpSampSeq'}
            % Number of warping sample for sequences, slow single fingers  Default = 500
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'BlockGroup'}
            % Block Group to plot
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'DownsampleRate'}
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'numFreqBins'}
            % number of frequency bins to consider inside the FreqRange
            % default  = 90
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'FreqRange'}
            % min frequency default = [2 150]
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'Channels'}
            % channels of interest Default : everythig
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        otherwise
            error(sprintf('Unknown option: %s',varargin{c}));
    end
end

Fs = 1024;
Fs_ds = floor(Fs/DownsampleRate);

%% HouseKeeping  : load required data

mainDir = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' subjname{subjNum} ,'/'];
load([mainDir , 'ChanLabels.mat'])
ChanLabels = ChanLabels(Channels);

load([mainDir , 'AllData_Behav.mat'])
load([mainDir , 'AllData_AvgMarker.mat'])
Events  = seqeeg_addEventMarker(Dall,subjNum, Fs_ds, 'addEvent' , 'NumWarpSampFast' , NumWarpSampFast, 'NumWarpSampSlow'  ,NumWarpSampSlow);
E  = seqeeg_addEventMarker(Events, subjNum, Fs_ds , 'CalcAveragePattern_seqType' , 'NumWarpSampFast' , NumWarpSampFast, 'NumWarpSampSlow'  ,NumWarpSampSlow)';
save([mainDir , 'AllData_AvgMarker_SeqType.mat'] , 'E');

E  = seqeeg_addEventMarker(Events, subjNum, Fs_ds , 'CalcAveragePattern' , 'NumWarpSampFast' , NumWarpSampFast, 'NumWarpSampSlow'  ,NumWarpSampSlow)';
save([mainDir , 'AllData_AvgMarker.mat'] , 'E');


[~, ~, BLockGroups] = xlsread([saveDir , 'BLockGroups.xlsx'],'Sheet1');
BLockGroups = BLockGroups(1:end,:);
BLockGroups(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),BLockGroups)) = {''};
idx = cellfun(@ischar, BLockGroups);
BLockGroups(idx) = cellfun(@(x) string(x), BLockGroups(idx), 'UniformOutput', false);
clearvars idx;

for bg = 1:length(BLockGroups)
    if ~isnumeric(BLockGroups{bg,1})
        blockGroups{bg} = str2num(char(BLockGroups{bg,1}));
    else
        blockGroups{bg} = BLockGroups{bg,1};
    end
end
blockGroupNames = BLockGroups(:,2);
fastBlock = horzcat(blockGroups{1} ,blockGroups{6} , blockGroups{7}, blockGroups{8}, blockGroups{9},...
    blockGroups{12}, blockGroups{16}, blockGroups{20}, blockGroups{23});
%% Define freq bands
min_freq =  FreqRange(1);
max_freq = FreqRange(2);
frex = linspace(min_freq, max_freq,numFreqBins);
BandInfo.bandsLab = {'Delta <4Hz' , 'Theta 5-8Hz' , 'Alpha 9-16Hz' , 'Beta 17-36Hz' , 'L-Gamma 37-80Hz' , 'H-Gamma 80-100HZ' , 'HIGH 100-180HZ'};
BandInfo.bands = {[0 4], [5 8] [9 16] [17 36] [37 80] [80 100] [100 180]};
for b = 1:length(BandInfo.bands)
    BandInfo.bandid{b} = [find(frex>BandInfo.bands{b}(1) ,1, 'first') , find(frex<BandInfo.bands{b}(2) ,1, 'last')];
end

%%
switch what
    case 'binned_BlockGroup'
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped.mat

        Pall.Fast = zeros(size(Pall.TN));
        E.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
        E.NumWarpSamp([1 6 7 8 9 12 16 20 23]) = NumWarpSampFast;
        Pall.Fast(ismember(Pall.BN , fastBlock)) = 1;
        clear Dall
        for BG = 1:length(blockGroups)
            clear tempPow binnedPow tcount F Rep
            saveName = [mainDir,'AverageBinnedPSD',num2str(BG) , '.mat'];
            E1 = getrow(E ,  BG);
            if ~isempty(E1.SN{1})
                for sn = 1:length(E1.SN{1})
                    id = ismember(Pall.BN , E1.blockGroups{1}) & ismember(Pall.seqNumb , E1.SN{1}(sn));
                    F = getrow(Pall , id);
                    % sum the warped PSDs insife the F structure
                    tcount = 1;
                    for tn = 1:length(F.PSD)
                        if isequal(size(F.PSD{tn}) ,[length(ChanLabels) , length(BandInfo.bandsLab) , E1.NumWarpSamp])
                            tempPow(tcount , :,:,:) = F.PSD{tn};
                            R(tcount , 1) = F.Rep(tn);
                            tcount = tcount +1;
                        end
                    end
                    P.binnedPow{sn,1} = tempPow;
                    P.Rep{sn,1}       = R;
                    clear tempPow R
                end
            else
                P.binnedPow = [];
                P.Rep = [];
            end
            save(saveName , '-struct','P');
        end
    case 'raw_BlockGroup'
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped.mat
        % patterned PSDs so the output  of Pall  = seqeeg_parseEEG_PSD('TimeWarpPSD' , Dall, subjNum);
        
        Pall.Fast = zeros(size(Pall.TN));
        E.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
        E.NumWarpSamp([1 6 7 8 9 12 16 20 23]) = NumWarpSampFast;
        Pall.Fast(ismember(Pall.BN , fastBlock)) = 1;
        P = [];
        for BG = 1:length(blockGroups)
            E1 = getrow(E ,  BG);
            saveName = [mainDir,'AverageRawPSD',num2str(BG) , '.mat'];
            if length(E1.SN{1})>0
                for bn = 1:length(E.blockGroups{BG})
                    BN = E.blockGroups{BG}(bn);
                    filename = [mainDir ,  'warped_PSD_B' , num2str(BN),'.mat'];
                    Pall = getrow(Dall , Dall.BN == BN);
                    for tn = 1:length(Pall.TN)
                        trialName = ['PSD',num2str(tn)];
                        A = load(filename , trialName);
                        Pall.PSD_stim{tn,1} = eval(['A.',trialName,';']);
                    end
                    P = addstruct(P,Pall);
                end
                clear tempPow rawPow tcount F Rep
                for sn = 1:length(E1.SN{1})
                    id = ismember(P.BN , E1.blockGroups{1}) & ismember(P.seqNumb , E1.SN{1}(sn));
                    F = getrow(P , id);
                    % sum the warped PSDs insife the F structure
                    tcount = 1;
                    for tn = 1:length(F.PSD_stim)
                        if isequal(size(F.PSD_stim{tn}) ,[length(ChanLabels) , numFreqBins , E1.NumWarpSamp])
                            tempPow(tcount , :,:,:) = F.PSD_stim{tn};
                            R(tcount , 1) = F.Rep(tn);
                            tcount = tcount +1;
                        end
                    end
                    rawPow{sn,1} = tempPow;
                    Rep{sn , 1}  = R;
                    clear tempPow R
                end
            else
                rawPow  = [];
                Rep = [];
            end
            E1.rawPow = rawPow;
            E1.Rep    = Rep;
            Pall = E1;
            save(saveName , 'Pall');
        end
    case 'raw_SeqType'
        load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped_SeqType.mat
        % patterned PSDs so the output  of Pall  = seqeeg_parseEEG_PSD('TimeWarpPSD_Raw_Binned_seqType' , Dall, subjNum);
        
        % Define sequence numbers and their transformations:
        SeqTrans = [5 11 22 33 44 55 0 1 2 3 4 103 104 203 204;...
            100 10 10 10 10 10 20 30 30 30 30 40 50 40 50];
        for sn = 1:length(SeqTrans)
            id  = Dall.seqNumb == SeqTrans(1 , sn);
            Dall.seqNumb(id) = SeqTrans(2 , sn);
        end
        
        Dall.Fast = zeros(size(Dall.TN));
        E.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
        E.NumWarpSamp([1 6 7 8 9 12 16 20 23]) = NumWarpSampFast;
        Dall.Fast(ismember(Dall.BN , fastBlock)) = 1;
        P = [];
        for BG = 1:length(blockGroups)
            saveName = [mainDir,'AverageRawPSD_SeqType',num2str(BG) , '.mat'];
            E1 = getrow(E ,  BG);
            if length(E1.SN{1})>0
                for bn = 1:length(E.blockGroups{BG})
                    BN = E.blockGroups{BG}(bn);
                    filename = [mainDir ,  'warped_PSD_B_SeqType' , num2str(BN),'.mat'];
                    Pall = getrow(Dall , Dall.BN == BN);
                    for tn = 1:length(Pall.TN)
                        trialName = ['PSD',num2str(tn)];
                        A = load(filename , trialName);
                        Pall.PSD_stim{tn,1} = eval(['A.',trialName,';']);
                    end
                    P = addstruct(P,Pall);
                end
                
                clear tempPow AvgPow tcount F Rep
                
                
                for sn = 1:length(E1.SN{1})
                    %             NEM = E.NEM{1}(sn , :);
                    id = ismember(P.BN , E1.blockGroups{1}) & ismember(P.seqNumb , E1.SN{1}(sn));
                    F = getrow(P , id);
                    % sum the warped PSDs insife the F structure
                    tcount = 1;
                    for tn = 1:length(F.PSD_stim)
                        if isequal(size(F.PSD_stim{tn}) ,[length(ChanLabels) , numFreqBins , E1.NumWarpSamp])
                            tempPow(tcount , :,:,:) = F.PSD_stim{tn};
                            R(tcount , 1) = F.Rep(tn);
                            tcount = tcount +1;
                        end
                    end
                    AvgPow{sn,1} = tempPow;
                    Rep{sn,1}    = R;
                    clear tempPow R
                end
            else
                AvgPow = [];
                Rep = [];
            end
            E1.AvgPow = AvgPow;
            E1.Rep    = Rep;
            Pall = E1;
            save(saveName , 'Pall');
        end
    case 'binned_SeqType'
        % Pall is the AllData_PSD_Warped_SeqType.mat
        % patterned PSDs so the output of Pall  = seqeeg_parseEEG_PSD('TimeWarpPSD_Raw_Binned_seqType' , Pall, subjNum);
        load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
        SeqTrans = [5 11 22 33 44 55 0 1 2 3 4 103 104 203 204;...
            100 10 10 10 10 10 20 30 30 30 30 40 50 40 50];
        for sn = 1:length(SeqTrans)
            id  = Pall.seqNumb == SeqTrans(1 , sn);
            Pall.seqNumb(id) = SeqTrans(2 , sn);
        end
        
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped_SeqType.mat
        % patterned PSDs so the output of Pall  = seqeeg_parseEEG_PSD('TimeWarpPSD_Raw_Binned_seqType' , Pall, subjNum);
        Pall.Fast = zeros(size(Pall.TN));
        E.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
        E.NumWarpSamp([1 6 7 8 9 12 16 20 23]) = NumWarpSampFast;
        Pall.Fast(ismember(Pall.BN , fastBlock)) = 1;
        clear Dall
        for BG = 23:length(blockGroups)
            clear tempPow binnedPow tcount F Rep
            saveName = [mainDir,'AverageBinnedPSD_SeqType',num2str(BG) , '.mat'];
            E1 = getrow(E ,  BG);
            if length(E1.SN{1})>0
                for sn = 1:length(E1.SN{1})
                    [sn BG]
                    %             NEM = E.NEM{1}(sn , :);
                    id = ismember(Pall.BN , E1.blockGroups{1}) & ismember(Pall.seqNumb , E1.SN{1}(sn));
                    F = getrow(Pall , id);
                    % sum the warped PSDs insife the F structure
                    tcount = 1;
                    for tn = 1:length(F.PSD)
                        if isequal(size(F.PSD{tn}) ,[length(ChanLabels) , length(BandInfo.bandsLab) , E1.NumWarpSamp])
                            tempPow(tcount , :,:,:) = F.PSD{tn};
                            R(tcount , 1) = F.Rep(tn);
                            tcount = tcount +1;
                        end
                    end
                    P.binnedPow{sn} = tempPow;
                    P.Rep{sn} = R;
                    clear tempPow R
                end
            else
                P.binnedPow = [];
                P.Rep = [];
            end
            save(saveName , '-struct','P');
        end
    case 'raw_AvgPower_BlockGroup'
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped.mat
        % patterned PSDs so the output  of Pall  = seqeeg_parseEEG_PSD('TimeWarpPSD_Raw_Binned' , Dall, subjNum);
        
        Pall = Dall;
        Pall.Fast = zeros(size(Pall.TN));
        E.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
        E.NumWarpSamp([1 6 7 8 9 12 16 20 23]) = NumWarpSampFast;
        Pall.Fast(ismember(Pall.BN , fastBlock)) = 1;
        P = [];
        repTypes = {[1] , [2] , [1 2]};
        for BG = 1:length(blockGroups)
            saveName = [mainDir,'AverageSpectBG',num2str(BG) , '.mat'];
            E1 = getrow(E ,  BG);
            if length(E1.SN{1})>0
                for bn = 1:length(E.blockGroups{BG})
                    BN = E.blockGroups{BG}(bn);
                    filename = [mainDir ,  'Raw_PSD_B' , num2str(BN),'.mat'];
                    Pall = getrow(Dall , Dall.BN == BN );
                    for tn = 1:length(Pall.TN)
                        T = load(filename , ['PSD',num2str(tn)] , ['BaseLine',num2str(tn)]);
                        psd = eval(['T.PSD', num2str(tn) , '(:,:,floor(Fs_ds*TimeDelay):end - floor(Fs_ds*2*TimeDelay))']);
                        bl = repmat(eval(['T.BaseLine', num2str(tn)]) ,1,1, size(psd , 3));
                        Pall.AvgPowTR{tn,1} = nanmean(100*(psd - bl)./bl  , 3); % percent change from baseline
                    end
                    P = addstruct(P,Pall);
                end
                clear tempPow AvgPow tcount F
                for sn = 1:length(E1.SN{1})
                    id = ismember(P.BN , E1.blockGroups{1}) & ismember(P.seqNumb , E1.SN{1}(sn));
                    F = getrow(P , id);
                    % sum the warped PSDs insife the F structure
                    tcount = 1;
                    for tn = 1:length(F.AvgPowTR)
                        if isequal(length(F.AvgPowTR{tn}) , numFreqBins)
                            tempPowTR(tcount,:,:) = F.AvgPowTR{tn};
                            Rep(tcount , 1)  = F.Rep(tn);
                            tcount = tcount +1;
                        end
                    end
                    E1.PowTR{sn} = tempPowTR;
                    E1.Rep{sn} = Rep;
                    clear tempPowTR tempPowBL tempPowAF
                end
                
            else
                E1.PowTR = {};
                E1.Rep= {};
            end
            Pall = E1;
            save(saveName , 'Pall');
        end        
    case 'raw_AvgPower_SeqType'
        % Pall is the AllData_PSD_Warped_SeqType.mat
        % patterned PSDs so the output of Pall  = seqeeg_parseEEG_PSD('TimeWarpPSD_Raw_Binned_seqType' , Pall, subjNum);
        load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
        Dall;
        SeqTrans = [5 11 22 33 44 55 0 1 2 3 4 103 104 203 204;...
            100 10 10 10 10 10 20 30 30 30 30 40 50 40 50];
        for sn = 1:length(SeqTrans)
            id  = Dall.seqNumb == SeqTrans(1 , sn);
            Dall.seqNumb(id) = SeqTrans(2 , sn);
        end
        
        
        
        P = [];
        for BG = 1:length(blockGroups)
            saveName = [mainDir,'PowerSpectrum_SeqType',num2str(BG) , '.mat'];
            E1 = getrow(E ,  BG);
            if length(E1.SN{1})>0
                for bn = 1:length(E.blockGroups{BG})
                    BN = E.blockGroups{BG}(bn);
                    filename = [mainDir ,  'Raw_PSD_B' , num2str(BN),'.mat'];
                    Pall = getrow(Dall , Dall.BN == BN );
                    for tn = 1:length(Pall.TN)
                        T = load(filename , ['PSD',num2str(tn)] , ['BaseLine',num2str(tn)]);
                        psd = eval(['T.PSD', num2str(tn) , '(:,:,floor(Fs_ds*TimeDelay):end - floor(Fs_ds*2*TimeDelay))']);
                        bl = repmat(eval(['T.BaseLine', num2str(tn)]) ,1,1, size(psd , 3));
                        Pall.AvgPowTR{tn,1} = nanmean(100*(psd - bl)./bl  , 3); % percent change from baseline
                    end
                    P = addstruct(P,Pall);
                end
                
                clear tempPow AvgPow tcount F
                
                
                for sn = 1:length(E1.SN{1})
                    id = ismember(P.BN , E1.blockGroups{1}) & ismember(P.seqNumb , E1.SN{1}(sn));
                    F = getrow(P , id);
                    % sum the warped PSDs insife the F structure
                    tcount = 1;
                    for tn = 1:length(F.AvgPowTR)
                        if isequal(length(F.AvgPowTR{tn}) , numFreqBins)
                            tempPowTR(tcount,:,:) = F.AvgPowTR{tn};
                            Rep(tcount , 1)  = F.Rep(tn);
                            tcount = tcount +1;
                        end
                    end
                    E1.PowTR{sn} = tempPowTR;
                    E1.Rep{sn} = Rep;
                    clear tempPowTR tempPowBL tempPowAF
                end
            else
                E1.PowTR = {};
                E1.Rep = {};
            end
            Pall = E1;
            save(saveName , 'Pall');
        end        
    case 'raw_Power_SeqType'
        
        % Pall is the AllData_PSD_Warped_SeqType.mat
        % patterned PSDs so the output of Pall  = seqeeg_parseEEG_PSD('TimeWarpPSD_Raw_Binned_seqType' , Pall, subjNum);
        load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
        SeqTrans = [5 11 22 33 44 55 0 1 2 3 4 103 104 203 204;...
            100 10 10 10 10 10 20 30 30 30 30 40 50 40 50];
        for sn = 1:length(SeqTrans)
            id  = Dall.seqNumb == SeqTrans(1 , sn);
            Dall.seqNumb(id) = SeqTrans(2 , sn);
        end
        
        
        
        P = [];
        for BG = 1:length(blockGroups)
            saveName = [mainDir,'PowerSpectrum_TimeKept_SeqType',num2str(BG) , '.mat'];
                E1 = getrow(E ,  BG);
            if length(E1.SN{1})>0
                for bn = 1:length(E.blockGroups{BG})
                    BN = E.blockGroups{BG}(bn);
                    filename = [mainDir ,  'warped_PSD_B_SeqType' , num2str(BN),'.mat'];
                    Pall = getrow(Dall , Dall.BN == BN);
                    for tn = 1:length(Pall.TN)
                        trialName = ['PSD',num2str(tn)];
                        A = load(filename , trialName);
                        Pall.PSD{tn,1} = eval(['A.',trialName,';']);
                    end
                    P = addstruct(P,Pall);
                end

                for sn = 1:length(E1.SN{1})
                    %             NEM = E.NEM{1}(sn , :);
                    id = ismember(P.BN , E1.blockGroups{1}) & ismember(P.seqNumb , E1.SN{1}(sn));
                    F = getrow(P , id);
                    % sum the warped PSDs insife the F structure
                    tcount = 1;
                    for tn = 1:length(F.PSD)
                        if isequal(size(F.PSD{tn},2) , numFreqBins)
                            tempPowTR(tcount,:,:,:) = F.PSD{tn};
                            Rep(tcount ,1)          = F.Rep(tn);
                            tcount = tcount +1;
                        end
                    end
                    E1.Pow{sn} = tempPowTR;
                    E1.Rep{sn} = Rep;
                    clear tempPowTR 
                end
            else
                E1.Pow = {};
                E1.Rep = {};
            end
            Pall = E1;
            save(saveName , 'Pall');
        end
    case 'binned_Power_SeqType'
        % Pall is the AllData_PSD_Warped_SeqType.mat
        % patterned PSDs so the output of Pall  = seqeeg_parseEEG_PSD('TimeWarpPSD_Raw_Binned_seqType' , Pall, subjNum);
        load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
        SeqTrans = [5 11 22 33 44 55 0 1 2 3 4 103 104 203 204;...
            100 10 10 10 10 10 20 30 30 30 30 40 50 40 50];
        for sn = 1:length(SeqTrans)
            id  = Dall.seqNumb == SeqTrans(1 , sn);
            Dall.seqNumb(id) = SeqTrans(2 , sn);
        end
        for bg = 1:length(blockGroups)
            E1 = getrow(E , bg);
            filename = [mainDir ,  'PowerSpectrum_TimeKept_SeqType' , num2str(bg),'.mat'];
            saveName = [mainDir,'PowerSpectrumBinned_TimeKept_SeqType',num2str(bg) , '.mat'];
            if length(E1.SN{1})>0
                E1.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
                E1.NumWarpSamp([1 6 7 8 9 12 16 20 23]) = NumWarpSampFast;
                load(filename)
                %  apply banding to make it more managable
                for sn = 1:length(Pall.Pow)
                    for b = 1:length(BandInfo.bands)
                        ind = find(ismember(frex , [BandInfo.bands{b}(1) : BandInfo.bands{b}(2)]));
                        Tr = squeeze(nanmean(Pall.Pow{sn}(:,:,ind,:) , 3));
                        eval(['P.PSD_band' , num2str(b) , '{sn} = Tr;'])
                    end
                    P.Rep{sn,1} = Pall.Rep{sn};
                end
            else
                for b = 1:length(BandInfo.bands)
                    eval(['P.PSD_band' , num2str(b) , '={};'])
                end
                P.Rep = {};
            end
            save(saveName , '-struct','P', '-v7.3');
        end
    case 'Aligned_SeqType'

        load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
        SeqTrans = [5 11 22 33 44 55 0 1 2 3 4 103 104 203 204;...
            100 10 10 10 10 10 20 30 30 30 30 40 50 40 50];
        for sn = 1:length(SeqTrans)
            id  = Dall.seqNumb == SeqTrans(1 , sn);
            Dall.seqNumb(id) = SeqTrans(2 , sn);
        end
        
        for BG = 1:length(blockGroups)
            P = [];
            saveName = [mainDir,'Group_Aligned_PSD',num2str(BG) , '.mat'];
            E1 = getrow(E , BG);
            for bn = 1:length(blockGroups{BG})
                BN = blockGroups{BG}(bn);
                filename = [mainDir ,  'EventAligned_PSD_B' , num2str(BN),'.mat'];
                Pall = getrow(Dall , Dall.BN == BN);
                
                if length(E1.SN{1})>0 
                    for tn = 1:length(Pall.TN)
                        trialName1 = ['PSD_Aligned',num2str(tn)];
                        trialName2 = ['PSD_betweenEvnt',num2str(tn)];
                        A1 = load(filename , trialName1);
                        tempPSD1 = eval(['A1.',trialName1,';']);
                        A2 = load(filename , trialName2);
                        tempPSD2 = eval(['A2.',trialName2,';']);
                        %% bin the PSDs
                        for bin = 1:length(BandInfo.bands)
%                             ind = find(ismember(frex , [BandInfo.bands{bin}(1) : BandInfo.bands{bin}(2)]));
                            ind = find(frex>= BandInfo.bands{bin}(1) &  frex<=BandInfo.bands{bin}(2));
                            for event = 1:length(tempPSD1)
                                if ~isempty(tempPSD1{event})
                                    Pall.PSD1{tn , event}(:,bin, :)  = squeeze(nanmean(tempPSD1{event}(:,ind, :) , 2));
                                else
                                    Pall.PSD1{tn , event} = [];
                                end
                            end
                            for event = 1:length(tempPSD2)
                                if ~isempty(tempPSD2{event})
                                    Pall.PSD2{tn , event}(:,bin, :)  = squeeze(nanmean(tempPSD2{event}(:,ind, :) , 2));
                                else
                                    Pall.PSD2{tn , event} = [];
                                end
                            end
                        end
                    end
                else
                    Pall.PSD1 = {};
                    Pall.PSD2 = {};
                end
                if size(Pall.PSD1 , 2)>1
                    P = addstruct(P,Pall);
                end
            end
            
            %% bin the PSDs and pad them with nans
            Pall = P;
            save(saveName , 'Pall');
            clear P Pall
        end        
    case 'AlignedWarped_SeqType'
        % Pall is the AllData_PSD_Warped_SeqType.mat
        % patterned PSDs so the output of Pall  = seqeeg_parseEEG_PSD('TimeWarpPSD_Raw_Binned_seqType' , Pall, subjNum);
        load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
        SeqTrans = [5 11 22 33 44 55 0 1 2 3 4 103 104 203 204;...
            100 10 10 10 10 10 20 30 30 30 30 40 50 40 50];
        for sn = 1:length(SeqTrans)
            id  = Dall.seqNumb == SeqTrans(1 , sn);
            Dall.seqNumb(id) = SeqTrans(2 , sn);
        end

         for BG = 1:length(blockGroups)
            P = [];
            saveName = [mainDir,'Group_WarpedAligned_PSD',num2str(BG) , '.mat'];
            E1 = getrow(E , BG);
            for bn = 1:length(blockGroups{BG})
                BN = blockGroups{BG}(bn);
                filename = [mainDir ,  'EventAligned_WarpedPSD_B' , num2str(BN),'.mat'];
                Pall = getrow(Dall , Dall.BN == BN);
                
                if length(E1.SN{1})>0 
                    for tn = 1:length(Pall.TN)
                        trialName = ['PSD_Aligned',num2str(tn)];
                        A = load(filename , trialName);
                        tempPSD = eval(['A.',trialName,';']);
                        %% bin the PSDs
                        for bin = 1:length(BandInfo.bands)
%                             ind = find(ismember(frex , [BandInfo.bands{bin}(1) : BandInfo.bands{bin}(2)]));
                            ind = find(frex>= BandInfo.bands{bin}(1) &  frex<=BandInfo.bands{bin}(2));
                            for event = 1:length(tempPSD)
                                if ~isempty(tempPSD{event})
                                    Pall.PSD{tn , event}(:,bin, :)  = squeeze(nanmean(tempPSD{event}(:,ind, :) , 2));
                                else
                                    Pall.PSD{tn , event} = [];
                                end
                            end
                        end
                    end
                else
                    Pall.PSD = {};
                    
                end
                if size(Pall.PSD , 2)>1
                    P = addstruct(P,Pall);
                end
            end
            
            %% bin the PSDs and pad them with nans
            
            Pall = P;
            save(saveName , 'Pall');
            clear P Pall
        end        
        
end



