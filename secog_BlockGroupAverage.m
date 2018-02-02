function secog_BlockGroupAverage(Pall , subjNum, what,varargin)

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
        otherwise
            error(sprintf('Unknown option: %s',varargin{c}));
    end
end
if ~exist('NumWarpSampFast')
    NumWarpSampFast = 200;
end
if ~exist('NumWarpSampSlow')
    NumWarpSampSlow = 500;
end

if ~exist('DownsampleRate')
    DownsampleRate = 10;
end
if ~exist('numFreqBins')
    numFreqBins = 75;
end
Fs = 1024;
Fs_ds = floor(Fs/DownsampleRate);

%% load required data
mainDir = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' subjname{subjNum} ,'/'];
load([mainDir , 'ChanLabels.mat'])
load([mainDir , 'AllData_Behav.mat'])
load([mainDir , 'AllData_AvgMarker.mat'])

%% HouseKeeping
blockGroups = {[1 2] , [3 13 26 40] , [4 14 27 41] , [5:7] , [9:11] , [8 12] , [15:17] , [19:21] , [23:25],...
    [18 22] , [28:30] , [32:34] , [36:38], [31 35 39],[42:44]}';
blockGroupNames = {'SingleFingNat' , 'SingleFingSlow' , 'SingleFingFast' , 'Intermixed1' , 'Intermixed2' , 'ChunkDay1' , 'Intermixed3' , 'Intermixed4' , 'Intermixed5',...
    'ChunkDay2' , 'Intermixed6' , 'Intermixed7' , 'Intermixed8', 'ChunkDay3', 'Intermixed9'}';


bandsLab = {'Delta <4Hz' , 'Theta 4-8Hz' , 'Alpha 8-13Hz' , 'L-Beta 13-24Hz' , 'H-Beta 24-36Hz' , 'L-Gamma 36-48Hz' , 'H-Gamma >48Hz'};


switch what
    case 'binned_BlockGroup'
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped.mat
        % patterned PSDs so the output of Pall  = secog_parseEEG_PSD('TimeWarpPSD' , Pall, subjNum);
        Dall.PSD_stim = Pall.PSD_stim;
        Pall = Dall;
        Pall.Fast = zeros(size(Pall.TN));
        E.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
        E.NumWarpSamp([1 2 3 6 10 14]) = NumWarpSampFast;
        fastBlock = horzcat(blockGroups{1} , blockGroups{2} , blockGroups{3} , blockGroups{6}, blockGroups{10},blockGroups{14});
        Pall.Fast(ismember(Pall.BN , fastBlock)) = 1;
        clear Dall
        for BG = 1:length(blockGroups)
            clear tempPow AvgPow tcount F
            saveName = [mainDir,'AverageBinnedPSD',num2str(BG) , '.mat'];
            
            E1 = getrow(E ,  BG);
            
            for sn = 1:length(E1.SN{1})
                %             NEM = E.NEM{1}(sn , :);
                id = ismember(Pall.BN , E1.blockGroups{1}) & ismember(Pall.seqNumb , E1.SN{1}(sn));
                F = getrow(Pall , id);
                % sum the warped PSDs insife the F structure
                tcount = 1;
                for tn = 1:length(F.PSD_stim)
                    if isequal(size(F.PSD_stim{tn}) ,[length(ChanLabels) , length(bandsLab) , E1.NumWarpSamp])
                        tempPow(tcount , :,:,:) = F.PSD_stim{tn};
                        tcount = tcount +1;
                    end
                end
                AvgPow{sn} = abs(squeeze(nanmean(tempPow , 1)));
            end
            save(saveName , 'AvgPow');
        end
    case 'raw_BlockGroup'
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped.mat
        % patterned PSDs so the output of Pall  = secog_parseEEG_PSD('TimeWarpPSD' , Dall, subjNum);
        
        Pall = Dall;
        Pall.Fast = zeros(size(Pall.TN));
        E.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
        E.NumWarpSamp([1 3 6 10 14]) = NumWarpSampFast;
        fastBlock = horzcat(blockGroups{1} , blockGroups{3} , blockGroups{6}, blockGroups{10},blockGroups{14});
        Pall.Fast(ismember(Pall.BN , fastBlock)) = 1;
        P = [];
        for BG = 2%3:length(blockGroups)
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
                AvgPow{sn,1} = abs(squeeze(nanmean(tempPow , 1)));
            end
            
            E1.AvgPow = AvgPow; 
            Pall = E1;
            save(saveName , 'Pall');
        end
end
