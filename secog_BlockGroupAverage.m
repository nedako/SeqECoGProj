function secog_BlockGroupAverage(Pall , subjNum, what,varargin)


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

% define block types
blockGroups = {[1 2] , [3 13 26 40] , [4 14 27 41] , [5:7] , [9:11] , [8 12] , [15:17] , [19:21] , [23:25],...
    [18 22] , [28:30] , [32:34] , [36:38], [31 35 39],[42:44]}';
blockGroupNames = {'SingleFingNat' , 'SingleFingSlow' , 'SingleFingFast' , 'Intermixed1' , 'Intermixed2' , 'ChunkDay1' , 'Intermixed3' , 'Intermixed4' , 'Intermixed5',...
    'ChunkDay2' , 'Intermixed6' , 'Intermixed7' , 'Intermixed8', 'ChunkDay3', 'Intermixed9'}';


bandsLab = {'Delta <4Hz' , 'Theta 5-8Hz' , 'Alpha 9-16Hz' , 'Beta 17-366Hz' , 'L-Gamma 37-70Hz' , 'H-Gamma 70-130' , 'NoBandLand 130-180'};



switch what
    case 'binned_BlockGroup'
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped.mat
        % patterned PSDs so the output of Pall  = secog_parseEEG_PSD('TimeWarpPSD_Raw_Binned' , Pall, subjNum);
        Pall.Fast = zeros(size(Pall.TN));
        E.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
        E.NumWarpSamp([1 3 6 10 14]) = NumWarpSampFast;
        fastBlock = horzcat(blockGroups{1} , blockGroups{2} , blockGroups{3} , blockGroups{6}, blockGroups{10},blockGroups{14});
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
                    if isequal(size(F.Pow_Norm_stim{tn}) ,[length(ChanLabels) , length(bandsLab) , E1.NumWarpSamp])
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
        
        Pall = Dall;
        Pall.Fast = zeros(size(Pall.TN));
        E.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
        E.NumWarpSamp([1 3 6 10 14]) = NumWarpSampFast;
        fastBlock = horzcat(blockGroups{1} , blockGroups{3} , blockGroups{6}, blockGroups{10},blockGroups{14});
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
        E.NumWarpSamp([1 3 6 10 14]) = NumWarpSampFast;
        fastBlock = horzcat(blockGroups{1} , blockGroups{3} , blockGroups{6}, blockGroups{10},blockGroups{14});
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
        E.NumWarpSamp([1 3 6 10 14]) = NumWarpSampFast;
        fastBlock = horzcat(blockGroups{1} , blockGroups{3} , blockGroups{6}, blockGroups{10},blockGroups{14});
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
                    if isequal(size(F.Pow_Norm_stim{tn}) ,[length(ChanLabels) , length(bandsLab) , E1.NumWarpSamp])
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
        E.NumWarpSamp([1 3 6 10 14]) = NumWarpSampFast;
        fastBlock = horzcat(blockGroups{1} , blockGroups{3} , blockGroups{6}, blockGroups{10},blockGroups{14});
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
            saveName = [mainDir,'AverageSpect_SeqType',num2str(BG) , '.mat'];
            
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
end


