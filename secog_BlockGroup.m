function secog_BlockGroup(Pall , subjNum, what,varargin)


%% setup the defaults and deal with the varargin
DownsampleRate = 10;
NormType = 'stim';
NumWarpSampFast = 200;
NumWarpSampSlow = 500;
TimeDelay = 0.5; % sec
FreqRange = [2 180];
numFreqBins = 90;
Channels = [1:129];
subjname = {'P2'};
c = 1;
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
min_freq =  FreqRange(1);
max_freq = FreqRange(2);
frex = linspace(min_freq, max_freq,numFreqBins);
Fs = 1024;
Fs_ds = floor(Fs/DownsampleRate);

%% HouseKeeping  : load required data
mainDir = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' subjname{subjNum} ,'/'];
load([mainDir , 'ChanLabels.mat'])
ChanLabels = ChanLabels(Channels);

load([mainDir , 'AllData_Behav.mat'])
load([mainDir , 'AllData_AvgMarker.mat'])

% block groupings for subject 1
BG(1).blockGroups =  {[1 2] , [3], [13], [26], [40] , [4], [14], [27] [41] , [5:7] , [9:11] , [8 12] , [15:17] , [19:21] , [23:25],...
    [18 22] , [28:30] , [32:34] , [36:38], [31 35 39],[42:44]}';
% block groupings for subject 2
BG(2).blockGroups = {[ ] , [2 8], [14 20 26], [29 38], [], [1 7],[13 19 25], [28 37], [] , [3:5] , [9:11] , [6 12] , [15:17] , [21:23] , [],...
    [18 24] , [30:32] , [34:36] , [], [27 33],[]}';

% define block types
blockGroups = BG(subjNum).blockGroups;
blockGroupNames = {'SingleFingNat' , 'SingleFingSlow1' , 'SingleFingSlow2'  , 'SingleFingSlow3' ,'SingleFingSlow4',...
    'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3', 'SingleFingFast4' , 'Intermixed1' , 'Intermixed2' , ...
    'ChunkDay1' , 'Intermixed3' , 'Intermixed4' , 'Intermixed5', 'ChunkDay2' , 'Intermixed6' , ...
    'Intermixed7' , 'Intermixed8', 'ChunkDay3', 'Intermixed9'}';
fastBlock = horzcat(blockGroups{1} ,blockGroups{6} , blockGroups{7}, blockGroups{8}, blockGroups{9},...
    blockGroups{12}, blockGroups{16}, blockGroups{20});


BandInfo.bandsLab = {'Delta <4Hz' , 'Theta 5-8Hz' , 'Alpha 9-16Hz' , 'Beta 17-36Hz' , 'L-Gamma 37-70Hz' , 'H-Gamma 70-130' , 'NoBandLand 130-180'};
BandInfo.bands = {[0 4], [5 8] [9 16] [17 36] [37 70] [70 110] [110 180]};


switch what
    case 'binned_BlockGroup'
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped.mat

        Pall.Fast = zeros(size(Pall.TN));
        E.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
        E.NumWarpSamp([1 6 7 8 9 12 16 20]) = NumWarpSampFast;
        Pall.Fast(ismember(Pall.BN , fastBlock)) = 1;
        clear Dall
        for BG = 1:length(blockGroups)
            clear tempPow AvgPow tcount F
            saveName = [mainDir,'AverageBinnedPSD',num2str(BG) , '.mat'];
            
            E1 = getrow(E ,  BG);
            
            for sn = 1:length(E1.SN{1})
                [sn BG]
                %             NEM = E.NEM{1}(sn , :);
                id = ismember(Pall.BN , E1.blockGroups{1}) & ismember(Pall.seqNumb , E1.SN{1}(sn));
                F = getrow(Pall , id);
                % sum the warped PSDs insife the F structure
                tcount = 1;
                for tn = 1:length(F.Pow_Norm_stim)
                    if isequal(size(F.Pow_Norm_stim{tn}) ,[length(ChanLabels) , length(BandInfo.bandsLab) , E1.NumWarpSamp])
                        tempPow(tcount , :,:,:) = F.Pow_Norm_stim{tn};
                        tcount = tcount +1;
                    end
                end
                AvgPow{sn} = squeeze(nanmean(tempPow , 1));
                clear tempPow
            end
            save(saveName , 'AvgPow');
        end
    case 'raw_BlockGroup'
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped.mat
        % patterned PSDs so the output  of Pall  = secog_parseEEG_PSD('TimeWarpPSD' , Dall, subjNum);

        Pall.Fast = zeros(size(Pall.TN));
        E.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
        E.NumWarpSamp([1 6 7 8 9 12 16 20]) = NumWarpSampFast;
        Pall.Fast(ismember(Pall.BN , fastBlock)) = 1;
        P = [];
        for BG = 1:length(blockGroups)
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
            
            clear tempPow AvgPow tcount F
            saveName = [mainDir,'AverageRawPSD',num2str(BG) , '.mat'];
            
            E1 = getrow(E ,  BG);
            
            for sn = 1:length(E1.SN{1})
                %             NEM = E.NEM{1}(sn , :);
                id = ismember(P.BN , E1.blockGroups{1}) & ismember(P.seqNumb , E1.SN{1}(sn));
                F = getrow(P , id);
                % sum the warped PSDs insife the F structure
                tcount = 1;
                for tn = 1:length(F.PSD_stim)
                    if isequal(size(F.PSD_stim{tn}) ,[length(ChanLabels) , numFreqBins , E1.NumWarpSamp])
                        tempPow(tcount , :,:,:) = F.PSD_stim{tn};
                        tcount = tcount +1;
                    end
                end
                AvgPow{sn,1} = squeeze(nanmean(tempPow , 1));
                clear tempPow
            end
            
            E1.AvgPow = AvgPow;
            Pall = E1;
            save(saveName , 'Pall');
        end
    case 'raw_SeqType'
        load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped_SeqType.mat
        % patterned PSDs so the output  of Pall  = secog_parseEEG_PSD('TimeWarpPSD_Raw_Binned_seqType' , Dall, subjNum);
        
        % Define sequence numbers and their transformations:
        SeqTrans = [5 11 22 33 44 55 0 1 2 3 4 103 104 203 204;...
            100 10 10 10 10 10 20 30 30 30 30 40 50 40 50];
        for sn = 1:length(SeqTrans)
            id  = Dall.seqNumb == SeqTrans(1 , sn);
            Dall.seqNumb(id) = SeqTrans(2 , sn);
        end
        
        Dall.Fast = zeros(size(Dall.TN));
        E.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
        E.NumWarpSamp([1 6 7 8 9 12 16 20]) = NumWarpSampFast;
        Dall.Fast(ismember(Dall.BN , fastBlock)) = 1;
        P = [];
        for BG = 1:length(blockGroups)
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
            
            clear tempPow AvgPow tcount F
            saveName = [mainDir,'AverageRawPSD_SeqType',num2str(BG) , '.mat'];
            
            E1 = getrow(E ,  BG);
            
            for sn = 1:length(E1.SN{1})
                %             NEM = E.NEM{1}(sn , :);
                id = ismember(P.BN , E1.blockGroups{1}) & ismember(P.seqNumb , E1.SN{1}(sn));
                F = getrow(P , id);
                % sum the warped PSDs insife the F structure
                tcount = 1;
                for tn = 1:length(F.PSD_stim)
                    if isequal(size(F.PSD_stim{tn}) ,[length(ChanLabels) , numFreqBins , E1.NumWarpSamp])
                        tempPow(tcount , :,:,:) = F.PSD_stim{tn};
                        tcount = tcount +1;
                    end
                end
                AvgPow{sn,1} = squeeze(nanmean(tempPow , 1));
                clear tempPow
            end
            
            E1.AvgPow = AvgPow;
            Pall = E1;
            save(saveName , 'Pall');
        end
    case 'binned_SeqType'
        % Pall is the AllData_PSD_Warped_SeqType.mat
        % patterned PSDs so the output of Pall  = secog_parseEEG_PSD('TimeWarpPSD_Raw_Binned_seqType' , Pall, subjNum);
        load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
        SeqTrans = [5 11 22 33 44 55 0 1 2 3 4 103 104 203 204;...
            100 10 10 10 10 10 20 30 30 30 30 40 50 40 50];
        for sn = 1:length(SeqTrans)
            id  = Pall.seqNumb == SeqTrans(1 , sn);
            Pall.seqNumb(id) = SeqTrans(2 , sn);
        end
        
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped_SeqType.mat
        % patterned PSDs so the output of Pall  = secog_parseEEG_PSD('TimeWarpPSD_Raw_Binned_seqType' , Pall, subjNum);
        Pall.Fast = zeros(size(Pall.TN));
        E.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
        E.NumWarpSamp([1 6 7 8 9 12 16 20]) = NumWarpSampFast;
        Pall.Fast(ismember(Pall.BN , fastBlock)) = 1;
        clear Dall
        for BG = 1:length(blockGroups)
            clear tempPow AvgPow tcount F
            saveName = [mainDir,'AverageBinnedPSD_SeqType',num2str(BG) , '.mat'];
            E1 = getrow(E ,  BG);
            
            for sn = 1:length(E1.SN{1})
                [sn BG]
                %             NEM = E.NEM{1}(sn , :);
                id = ismember(Pall.BN , E1.blockGroups{1}) & ismember(Pall.seqNumb , E1.SN{1}(sn));
                F = getrow(Pall , id);
                % sum the warped PSDs insife the F structure
                tcount = 1;
                for tn = 1:length(F.Pow_Norm_stim)
                    if isequal(size(F.Pow_Norm_stim{tn}) ,[length(ChanLabels) , length(BandInfo.bandsLab) , E1.NumWarpSamp])
                        tempPow(tcount , :,:,:) = F.Pow_Norm_stim{tn};
                        tcount = tcount +1;
                    end
                end
                AvgPow{sn} = squeeze(nanmean(tempPow , 1));
                clear tempPow
            end
            save(saveName , 'AvgPow');
        end
    case 'raw_AvgPower_BlockGroup'
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped.mat
        % patterned PSDs so the output  of Pall  = secog_parseEEG_PSD('TimeWarpPSD_Raw_Binned' , Dall, subjNum);
     
        Pall = Dall;
        Pall.Fast = zeros(size(Pall.TN));
        E.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
        E.NumWarpSamp([1 6 7 8 9 12 16 20]) = NumWarpSampFast;
        Pall.Fast(ismember(Pall.BN , fastBlock)) = 1;
        P = [];
        for BG = 1:length(blockGroups)
            for bn = 1:length(E.blockGroups{BG})
                BN = E.blockGroups{BG}(bn);
                filename = [mainDir ,  'Raw_Decomp_B' , num2str(BN),'.mat'];
                Pall = getrow(Dall , Dall.BN == BN);
                for tn = 1:length(Pall.TN)
                    trialName = ['DEC',num2str(tn)];
                    A = load(filename , trialName);
                    temp = eval(['10*log10(abs(A.',trialName,'.decompTR));;']);
                    Pall.AvgPowTR{tn,1} = squeeze(mean(temp  , 3));
                    temp = eval(['10*log10(abs(A.',trialName,'.decompBefTR));;']);
                    Pall.AvgPowBL{tn,1} = squeeze(mean(temp  , 3));
                    temp = eval(['10*log10(abs(A.',trialName,'.decompAftTR));;']);
                    Pall.AvgPowAF{tn,1} = squeeze(mean(temp  , 3));
                end
                P = addstruct(P,Pall);
            end
            
            clear tempPow AvgPow tcount F
            saveName = [mainDir,'AverageSpectBG',num2str(BG) , '.mat'];
            
            E1 = getrow(E ,  BG);
            
            for sn = 1:length(E1.SN{1})
                [BG sn]
                %             NEM = E.NEM{1}(sn , :);
                id = ismember(P.BN , E1.blockGroups{1}) & ismember(P.seqNumb , E1.SN{1}(sn));
                F = getrow(P , id);
                % sum the warped PSDs insife the F structure
                tcount = 1;
                for tn = 1:length(F.AvgPowTR)
                    if isequal(length(F.AvgPowTR{tn}) , numFreqBins)
                        tempPowTR(tcount,:,:) = F.AvgPowTR{tn};
                        tempPowBL(tcount,:,:) = F.AvgPowBL{tn};
                        tempPowAF(tcount,:,:) = F.AvgPowAF{tn};
                        tcount = tcount +1;
                    end
                end
                AvgPowTR{sn,1} = squeeze(nanmean(tempPowTR , 1));
                SePowTR{sn , 1} = squeeze(std(tempPowTR , 1))/sqrt(size(tempPowTR , 1));
                AvgPowBL{sn,1} = squeeze(nanmean(tempPowBL , 1));
                SePowBL{sn , 1} = squeeze(std(tempPowBL , 1))/sqrt(size(tempPowBL , 1));
                AvgPowAF{sn,1} = squeeze(nanmean(tempPowAF , 1));
                SePowAF{sn , 1} = squeeze(std(tempPowAF , 1))/sqrt(size(tempPowAF , 1));
                clear tempPowTR tempPowBL tempPowAF
            end
            
            E1.AvgPowTR = AvgPowTR;
            E1.AvgPowBL = AvgPowBL;
            E1.AvgPowAF = AvgPowAF;
            E1.SePowTR = SePowTR;
            E1.SePowBL = SePowBL;
            E1.SePowAF = SePowAF;
            Pall = E1;
            save(saveName , 'Pall');
        end        
    case 'raw_AvgPower_SeqType'
        % Pall is the AllData_PSD_Warped_SeqType.mat
        % patterned PSDs so the output of Pall  = secog_parseEEG_PSD('TimeWarpPSD_Raw_Binned_seqType' , Pall, subjNum);
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
            for bn = 1:length(E.blockGroups{BG})
                BN = E.blockGroups{BG}(bn);
                filename = [mainDir ,  'Raw_Decomp_B' , num2str(BN),'.mat'];
                Pall = getrow(Dall , Dall.BN == BN);
                for tn = 1:length(Pall.TN)
                    trialName = ['DEC',num2str(tn)];
                    A = load(filename , trialName);
                    temp = eval(['10*log10(abs(A.',trialName,'.decompTR).^2);']);
                    Pall.AvgPowTR{tn,1} = squeeze(mean(temp  , 3));
                    temp = eval(['10*log10(abs(A.',trialName,'.decompBefTR).^2);']);
                    Pall.AvgPowBL{tn,1} = squeeze(mean(temp  , 3));
                    temp = eval(['10*log10(abs(A.',trialName,'.decompAftTR).^2);']);
                    Pall.AvgPowAF{tn,1} = squeeze(mean(temp  , 3));
                end
                P = addstruct(P,Pall);
            end
            
            clear tempPow AvgPow tcount F
            saveName = [mainDir,'PowerSpectrum_SeqType',num2str(BG) , '.mat'];
            
            E1 = getrow(E ,  BG);
            
            for sn = 1:length(E1.SN{1})
                %             NEM = E.NEM{1}(sn , :);
                id = ismember(P.BN , E1.blockGroups{1}) & ismember(P.seqNumb , E1.SN{1}(sn));
                F = getrow(P , id);
                % sum the warped PSDs insife the F structure
                tcount = 1;
                for tn = 1:length(F.AvgPowTR)
                    if isequal(length(F.AvgPowTR{tn}) , numFreqBins)
                        tempPowTR(tcount,:,:) = F.AvgPowTR{tn};
                        tempPowBL(tcount,:,:) = F.AvgPowBL{tn};
                        tempPowAF(tcount,:,:) = F.AvgPowAF{tn};
                        tcount = tcount +1;
                    end
                end
                E1.PowTR{sn} = tempPowTR;
                E1.PowBL{sn} = tempPowBL;
                E1.PowAF{sn} = tempPowAF;
                clear tempPowTR tempPowBL tempPowAF
            end
            Pall = E1;
            save(saveName , 'Pall');
        end        
    case 'raw_Power_SeqType'
        % Pall is the AllData_PSD_Warped_SeqType.mat
        % patterned PSDs so the output of Pall  = secog_parseEEG_PSD('TimeWarpPSD_Raw_Binned_seqType' , Pall, subjNum);
        load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
        SeqTrans = [5 11 22 33 44 55 0 1 2 3 4 103 104 203 204;...
            100 10 10 10 10 10 20 30 30 30 30 40 50 40 50];
        for sn = 1:length(SeqTrans)
            id  = Dall.seqNumb == SeqTrans(1 , sn);
            Dall.seqNumb(id) = SeqTrans(2 , sn);
        end
        
        
        
        P = [];
        for BG = 1:length(blockGroups)
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
            
            saveName = [mainDir,'PowerSpectrum_TimeKept_SeqType',num2str(BG) , '.mat'];
            
            E1 = getrow(E ,  BG);
            
            for sn = 1:length(E1.SN{1})
                %             NEM = E.NEM{1}(sn , :);
                id = ismember(P.BN , E1.blockGroups{1}) & ismember(P.seqNumb , E1.SN{1}(sn));
                F = getrow(P , id);
                % sum the warped PSDs insife the F structure
                tcount = 1;
                for tn = 1:length(F.PSD)
                    if isequal(size(F.PSD{tn},2) , numFreqBins)
                        tempPowTR(tcount,:,:,:) = F.PSD{tn};
                        tcount = tcount +1;
                    end
                end
                E1.Pow{sn} = tempPowTR;
                clear tempPowTR tempPowBL tempPowAF
            end
            
            
            
            Pall = E1;
            save(saveName , 'Pall');
        end        
    case 'binned_Power_SeqType'
        % Pall is the AllData_PSD_Warped_SeqType.mat
        % patterned PSDs so the output of Pall  = secog_parseEEG_PSD('TimeWarpPSD_Raw_Binned_seqType' , Pall, subjNum);
        load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
        Dall;
        SeqTrans = [5 11 22 33 44 55 0 1 2 3 4 103 104 203 204;...
            100 10 10 10 10 10 20 30 30 30 30 40 50 40 50];
        for sn = 1:length(SeqTrans)
            id  = Dall.seqNumb == SeqTrans(1 , sn);
            Dall.seqNumb(id) = SeqTrans(2 , sn);
        end
        for bg = 1:length(blockGroups)
            saveName = [mainDir,'PowerSpectrumBinned_TimeKept_SeqType',num2str(bg) , '.mat'];
            load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
            E.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
            E.NumWarpSamp([1 6 7 8 9 12 16 20]) = NumWarpSampFast;
            
            BG(bg) = find(strcmp(E.blockGroupNames , blockGroupNames{bg}));
            E = getrow(E , BG(bg));
            
            filename = [mainDir ,  'PowerSpectrum_TimeKept_SeqType' , num2str(BG(bg)),'.mat'];
            load(filename)
            %  apply banding to make it more managable
            for sn = 1:length(Pall.Pow)
                for b = 1:length(BandInfo.bands)
                    ind = find(ismember(frex , [BandInfo.bands{b}(1) : BandInfo.bands{b}(2)]));
                    Tr = squeeze(nanmean(Pall.Pow{sn}(:,:,ind,:) , 3));
                    eval(['P.PSD_band' , num2str(b) , '{sn} = Tr;'])
                end
            end
            save(saveName , '-struct','P', '-v7.3');
        end
    case 'Aligned_SeqType'
        % Pall is the AllData_PSD_Warped_SeqType.mat
        % patterned PSDs so the output of Pall  = secog_parseEEG_PSD('TimeWarpPSD_Raw_Binned_seqType' , Pall, subjNum);
        load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
        SeqTrans = [5 11 22 33 44 55 0 1 2 3 4 103 104 203 204;...
            100 10 10 10 10 10 20 30 30 30 30 40 50 40 50];
        for sn = 1:length(SeqTrans)
            id  = Dall.seqNumb == SeqTrans(1 , sn);
            Dall.seqNumb(id) = SeqTrans(2 , sn);
        end

        for BG = 1:length(blockGroups)
            P = [];
            for bn = 1:length(blockGroups{BG})
                BN = blockGroups{BG}(bn);
                filename = [mainDir ,  'EventAligned_PSD_B' , num2str(BN),'.mat'];
                Pall = getrow(Dall , Dall.BN == BN);
                for tn = 1:length(Pall.TN)
                    tn
                    trialName = ['PSD_Aligned',num2str(tn)];
                    A = load(filename , trialName);
                    tempPSD = eval(['A.',trialName,';']);
                    %% bin the PSDs
                    for bin = 1:length(BandInfo.bands)
                        ind = find(ismember(frex , [BandInfo.bands{bin}(1) : BandInfo.bands{bin}(2)]));
                        for event = 1:length(tempPSD)
                            if ~isempty(tempPSD{event})
                                Pall.PSD{tn , event}(:,bin, :)  = squeeze(nanmean(tempPSD{event}(:,ind, :) , 2));
                            else
                                Pall.PSD{tn , event} = [];
                            end
                        end
                    end
                end
                P = addstruct(P,Pall);
            end
            %% bin the PSDs and pad them with nans
            saveName = [mainDir,'Group_Aligned_PSD',num2str(BG) , '.mat'];
            Pall = P;
            save(saveName , 'Pall');
            clear P Pall
        end
        
    case 'AlignedWarped_SeqType'
        % Pall is the AllData_PSD_Warped_SeqType.mat
        % patterned PSDs so the output of Pall  = secog_parseEEG_PSD('TimeWarpPSD_Raw_Binned_seqType' , Pall, subjNum);
        load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
        SeqTrans = [5 11 22 33 44 55 0 1 2 3 4 103 104 203 204;...
            100 10 10 10 10 10 20 30 30 30 30 40 50 40 50];
        for sn = 1:length(SeqTrans)
            id  = Dall.seqNumb == SeqTrans(1 , sn);
            Dall.seqNumb(id) = SeqTrans(2 , sn);
        end

        for BG = 1:length(blockGroups)
            P = [];
            for bn = 1:length(blockGroups{BG})
                BN = blockGroups{BG}(bn);
                filename = [mainDir ,  'EventAligned_WarpedPSD_B' , num2str(BN),'.mat'];
                Pall = getrow(Dall , Dall.BN == BN);
                for tn = 1:length(Pall.TN)
                    tn
                    trialName = ['PSD_Aligned',num2str(tn)];
                    A = load(filename , trialName);
                    tempPSD = eval(['A.',trialName,';']);
                    %% bin the PSDs
                    for bin = 1:length(BandInfo.bands)
                        ind = find(ismember(frex , [BandInfo.bands{bin}(1) : BandInfo.bands{bin}(2)]));
                        for event = 1:length(tempPSD)
                            if ~isempty(tempPSD{event})
                                Pall.PSD{tn , event}(:,bin, :)  = squeeze(nanmean(tempPSD{event}(:,ind, :) , 2));
                            else
                                Pall.PSD{tn , event} = [];
                            end
                        end
                    end
                end
                P = addstruct(P,Pall);
            end
            %% bin the PSDs and pad them with nans
            saveName = [mainDir,'Group_WarpedAligned_PSD',num2str(BG) , '.mat'];
            Pall = P;
            save(saveName , 'Pall');
            clear P Pall
        end
        
end



