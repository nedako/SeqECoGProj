function C = seqeeg_ForcePower(what ,Pall, subjNum, varargin)

c = 1;
%% setup the defaults and deal with the varargin
subjname = {'P2' , 'P4' , 'P5'};
DownsampleRate = 10;
NumWarpSampFast = 150;
NumWarpSampSlow = 300;
TimeDelay = 0.5; % sec
FreqRange = [2 180];
numFreqBins = 90;

mainDir = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' subjname{subjNum} ,'/'];
load([mainDir , 'ChanLabels.mat']);
load([mainDir , 'AllData_Behav_Force.mat'])
load([mainDir , 'AllData_AvgMarker.mat'])
Channels = [1:length(ChanLabels)];   
%%  control for too short IPIs that the keys get accidentally pressed
if subjNum==1
    for i = 1:length(Dall.TN)
        if sum(Dall.IPI(i,:)<120)
            Dall.isError(i) =1;
        end
    end
end
%%
while(c<=length(varargin))
    switch(varargin{c})
        case {'DownsampleRate'}
            % folds by which you wnat the PSD to be downsampled
            % default 10
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'NumWarpSampChunk'}
            % Number of warping sample for chunks, fast single finger Default = 200
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'NumWarpSampSeq'}
            % Number of warping sample for sequences, slow single fingers  Default = 500
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'TimeDelay'}
            % The time delay before the stimulus comes on to consider for baseline normalization
            % Default  = 0.5sec
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'FreqRange'}
            % min frequency default = [2 150]
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
        case {'BlockGroup'}
            % Blockgroup within which you want to calculate  force-power
            % relationship - should be cell
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        otherwise
            error(sprintf('Unknown option: %s',varargin{c}));
    end
end
  
% block groupings for subjects
BG(1).blockGroups =  {[1 2] , [3], [13], [26], [40] , [4], [14], [27] [41] , [5:7] , [9:11] , [8 12] , [15:17] , [19:21] , [23:25],...
    [18 22] , [28:30] , [32:34] , [36:38], [31 35 39],[42:44],[],[]}';
% block groupings for subject 2
BG(2).blockGroups = {[ ] , [2 8], [14 20 26], [29 38], [], [1 7],[13 19 25], [28 37], [] , [3:5] , [9:11] , [6 12] , [15:17] , [21:23] , [],...
    [18 24] , [30:32] , [34:36] , [], [27 33],[],[],[]}';
% block groupings for subject 3
BG(3).blockGroups = {[ ] , [1 7], [13 19], [25 31], [37 43], [2 8],[14 20], [26 32], [38 44] , [3:5] , [9:11] , [6 12] , [15:17] , [21:23] , [],...
    [18 24] , [27:29] , [33:35] , [], [30 36],[39:41] , [45:47] [42 48]}';

% define block types
blockGroups = BG(subjNum).blockGroups;
blockGroupNames = {'SingleFingNat' , 'SingleFingSlow1' , 'SingleFingSlow2'  , 'SingleFingSlow3' ,'SingleFingSlow4',...
    'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3', 'SingleFingFast4' , 'Intermixed1' , 'Intermixed2' , ...
    'ChunkDay1' , 'Intermixed3' , 'Intermixed4' , 'Intermixed5', 'ChunkDay2' , 'Intermixed6' , ...
    'Intermixed7' , 'Intermixed8', 'ChunkDay3', 'Intermixed9','Intermixed10','ChunkDay4'}';
fastBlock = horzcat(blockGroups{1} ,blockGroups{6} , blockGroups{7}, blockGroups{8}, blockGroups{9},...
    blockGroups{12}, blockGroups{16}, blockGroups{20}, blockGroups{23});



load([mainDir , 'ChanLabels.mat'])
ChanLabels = ChanLabels(Channels);

BandInfo.bandsLab = {'Delta <4Hz' , 'Theta 5-8Hz' , 'Alpha 9-16Hz' , 'Beta 17-36Hz' , 'L-Gamma 37-80Hz' , 'H-Gamma 80-100HZ' , 'HIGH 100-180HZ'};
BandInfo.bands = {[0 4], [5 8] [9 16] [17 36] [37 80] [80 110] [100 180]};
switch what
    case 'corr_ForePow'
        if isnumeric(BlockGroup)
            blocks = BlockGroup;
        else
            blocks = [];
            for bg = 1:length(BlockGroup)
                E.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
                BG = find(strcmp(E.blockGroupNames , BlockGroup{bg}));
                blocks = [blocks E.blockGroups{BG}];
            end
        end
        E.NumWarpSamp([1 6 7 8 9 12 16 20 23]) = NumWarpSampFast;
        P = getrow(Pall , ismember(Pall.BN , blocks) & ~ Dall.isError & Dall.seqNumb~=5);
        D = getrow(Dall , ismember(Dall.BN , blocks) & ~ Dall.isError & Dall.seqNumb~=5);
        % downsample the forces to 10Hz
        for tn = 1:length(D.TN)
            clear ftemp
            for f = 1:size(D.F{tn} , 2)
                ftemp(:,f) = downsample(D.F{tn}(:,f) , 5);
            end
            D.F{tn} = ftemp;
            for ch = 1:size(P.Pow_Norm_stim{tn} , 1)
                for band = 1:size(P.Pow_Norm_stim{tn} , 2)
                    A = sum(D.F{tn},2);
                    B = squeeze(P.Pow_Norm_stim{tn}(ch , band , :));
                    L = min(length(A) , length(B));
                    %                     temp = corrcoef(A(1:L) , B(1:L));
                    [acor,lag] = xcorr(A(1:L) , B(1:L),'coeff' , 10); % limit the maximum lag to 200 ms
                    C.acor{tn,ch , band} = acor;
                    C.lag{tn,ch , band}  = lag;
%                     [C.pks{tn,ch , band},C.locs{tn,ch , band},w,p]= findpeaks(acor,lag,'MinPeakHeight',0.1);
                    [C.maxCorr(tn, ch , band) , i] = max(abs(acor));
                    C.maxCorr(tn, ch , band) = acor(i);
                    C.maxLag(tn, ch , band) = lag(i);
                end
            end
        end
        C.maxCorr = squeeze(nanmean(C.maxCorr , 1));
        C.maxLag = squeeze(nanmedian(C.maxLag ,1));
        for ch = 1:size(C.maxCorr , 1)
            [~ , C.bestBandPerChan(ch,:)] = sort(abs(C.maxCorr(ch,:)) , 'descend');
        end
        for b = 1:size(C.maxCorr , 2)
            [~ , C.bestChanPerBand(b,:)] = sort(abs(C.maxCorr(1:length(Channels)-1,b)) , 'descend');
        end
        [C.SortedBands(2,:),C.SortedBands(1,:)] = sort(nansum(abs(C.maxCorr) , 1) , 'descend');
        [C.SortedChans{1}(2,:),C.SortedChans{1}(1,:)] = sort(nansum(abs(C.maxCorr) , 2) , 'descend');
        C.SortedChans{2} = ChanLabels(C.SortedChans{1}(1,:));
    case 'AlignPlot'
        if isnumeric(BlockGroup)
            blocks = BlockGroup;
        else
            blocks = [];
            for bg = 1:length(BlockGroup)
                E.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
                BG = find(strcmp(E.blockGroupNames , BlockGroup{bg}));
                blocks = [blocks E.blockGroups{BG}];
            end
        end
        E.NumWarpSamp([1 6 7 8 9 12 16 20 23]) = NumWarpSampFast;
        C = seqeeg_ForcePower('corr_ForePow' ,Pall, subjNum, 'Channels' , Channels , 'BlockGroup' , BlockGroup);
        Corr.cor = reshape(C.maxCorr , numel(C.maxCorr) , 1);
        Corr.band = repmat([1:length(BandInfo.bandsLab)] , length(Channels) , 1);
        Corr.band = reshape(Corr.band ,numel(Corr.band) , 1);
        
        Corr.ch = repmat([1:length(Channels)] , 1 ,  length(BandInfo.bandsLab));
        Corr.ch = reshape(Corr.ch ,numel(Corr.ch) , 1);
        % get rid of the TTL channel
        Corr = getrow(Corr,Corr.ch ~= length(Channels));
        [Corr.Sortedcor(:,2),Corr.Sortedcor(:,1)] = sort(Corr.cor  , 'descend');
        Corr.Sortedband = Corr.band(Corr.Sortedcor(:,1));
        Corr.Sortedch = Corr.ch(Corr.Sortedcor(:,1));
        Corr.SortedchLabels = ChanLabels(Corr.Sortedch);
%         seqeeg_ForcePower('binned_AlignTrials' ,Pall,subjNum, 'Channels' , Channels,'BlockGroup',BlockGroup)
%         seqeeg_ForcePower('binned_AlignTrials' ,Pall,subjNum, 'Chan2Plot' , Corr.Sortedch(2),'Channels' , Channels,'BlockGroup',BlockGroup ,...
%             'BandPowerForceCorr' , C.maxCorr(Corr.Sortedch(1) , :))
%         seqeeg_ForcePower('binned_AlignTrials' ,Pall,subjNum, 'Chan2Plot' , Corr.Sortedch(3),'Channels' , Channels,'BlockGroup',BlockGroup ,...
%             'BandPowerForceCorr' , C.maxCorr(Corr.Sortedch(1) , :))
    case 'binned_AlignTrials'
        C = seqeeg_ForcePower('AlignPlot' ,Pall, subjNum, 'Channels' , Channels , 'BlockGroup' , BlockGroup);
        % pall --> load([mainDir , 'AllData_PSD_StimNorm.mat'])
        N = input('How many channels to average over?');
        if isnumeric(BlockGroup)
            blocks = BlockGroup;
        else
            blocks = [];
            for bg = 1:length(BlockGroup)
                E.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
                BG = find(strcmp(E.blockGroupNames , BlockGroup{bg}));
                blocks = [blocks E.blockGroups{BG}];
            end
        end
        E.NumWarpSamp([1 6 7 8 9 12 16 20 23]) = NumWarpSampFast;
        P = getrow(Pall , ismember(Pall.BN , blocks) & ~Pall.isError & Pall.seqNumb~=5);
        D = getrow(Dall , ismember(Dall.BN , blocks) & ~Dall.isError & Dall.seqNumb~=5);  
        % obtain the average for the first -500ms to 1000msec after stim onset
        % obtain the average for the last 1000ms to 1000msec after stim offset
        % plot the average of the first 3 channels in every band
        PSD.psd_s_time = repmat([0:10:1499] , length(P.TN) , 1);
        PSD.psd_s_time = reshape(PSD.psd_s_time , numel(PSD.psd_s_time) , 1);
        PSD.psd_e_time = repmat([0:10:1999] , length(P.TN) , 1);
        PSD.psd_e_time = reshape(PSD.psd_e_time , numel(PSD.psd_e_time) , 1);
        h1 = figure;
        for bnd = 1:length(BandInfo.bandsLab)
            Chan2Plot = C.bestChanPerBand(bnd , 1:N);
            for tn = 1:length(P.TN)
                PSD.psd_s{bnd}(tn,:) = squeeze(nanmean(P.Pow_Norm_stim{tn}(Chan2Plot , bnd , 1:150),1)); % align the start
                PSD.psd_e{bnd}(tn,:) = squeeze(nanmean(P.Pow_Norm_stim{tn}(Chan2Plot , bnd , end-199:end),1)); % align the start
            end
            [xBand_s{bnd} , pBand_s{bnd} , eBand_s{bnd}] = lineplot(PSD.psd_s_time, reshape(PSD.psd_s{bnd},numel(PSD.psd_s{bnd}),1) , 'plotfcn' , 'nanmean');
            [xBand_e{bnd} , pBand_e{bnd} , eBand_e{bnd}] = lineplot(PSD.psd_e_time, reshape(PSD.psd_e{bnd},numel(PSD.psd_e{bnd}),1) , 'plotfcn' , 'nanmean');
        end
        close(h1)
        figure('color' , 'white')
        figCount = 1;
        for b =1:length(BandInfo.bandsLab)
            subplot(length(BandInfo.bandsLab),2, figCount)
            Chan2Plot = C.bestChanPerBand(b , 1:N);
            hold on
            plotshade(xBand_s{b}' , pBand_s{b} , eBand_s{b},'patchcolor',[.8 0 .3] , 'linecolor' , [.8 0 .3])
            corTag = num2str(C.maxCorr(Chan2Plot(1),b));
            namTag = ChanLabels{1};
            for n = 2:N
                corTag = [corTag,', ' ,num2str(C.maxCorr(Chan2Plot(n),b))];
                namTag = [namTag,', ' , ChanLabels{Chan2Plot(n)}];
            end
            title (['500ms pre stim onset - ch(s) ',namTag , ' Cor_F = ', corTag])
            line([500 500] , [min(pBand_s{b}) max(pBand_s{b})] , 'color' , 'b' , 'LineStyle' , ':' , 'LineWidth' , 3)
            line([1 1500]  , [0 0], 'color' , 'k' , 'LineWidth' , 1)
            ylabel(['%C ',BandInfo.bandsLab{b}])
            xlabel('Time (ms)')
            set(gca ,'FontSize' , 14,'Box' , 'off');
            figCount = figCount + 1;
            
            
            
            subplot(length(BandInfo.bandsLab),2, figCount)
            hold on
            plotshade(xBand_e{b}' , pBand_e{b} , eBand_e{b},'patchcolor',[.8 0 .3] , 'linecolor' , [.8 0 .3])
            corTag = num2str(C.maxCorr(Chan2Plot(1),b));
            namTag = ChanLabels{1};
            for n = 2:N
                corTag = [corTag,', ' ,num2str(C.maxCorr(Chan2Plot(n),b))];
                namTag = [namTag,', ' , ChanLabels{Chan2Plot(n)}];
            end
            title (['1000ms post stim offset - ch(s) ',namTag , ' Cor_F = ', corTag])
            line([1000 1000] , [min(pBand_e{b}) max(pBand_e{b})] , 'color' , 'b' , 'LineStyle' , ':' , 'LineWidth' , 3)
            line([1 2000]  , [0 0], 'color' , 'k' , 'LineWidth' , 1)
            ylabel(['%C ',BandInfo.bandsLab{b}])
            xlabel('Time (ms)')
            set(gca ,'FontSize' , 14,'Box' , 'off');
            figCount = figCount + 1;
        end
        
end

