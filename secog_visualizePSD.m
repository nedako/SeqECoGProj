function secog_visualizePSD(Pall , subjNum, what,varargin)
% for patient 2 Peter:
%       SMA channels:     [123:126 ]
%       preSMA channels : [128 , 129]
%       PMd channels:     [12 14 , 104:107 , 109]
%       PMv channel :     [6]
%% plotting case summary
%      'binned_SingleTrial'
%      'binned_SingleTrial_AvgChann'
%      'binned_BlockGroup_AvgChann'
%      'binned_BlockGroup'
%      'raw_SingleTrial'
%      'raw_BlockGroup'
%      'raw_BlockGroup_AvgChann'

%% set defaults and deal with varargin
subjname = {'P2' , 'P4'};
mainDir = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' subjname{subjNum} ,'/'];
NumWarpSampFast = 200;
NumWarpSampSlow = 500;
DownsampleRate = 10;
FreqRange = [2 180];
numFreqBins = 90;
load([mainDir , 'ChanLabels.mat']);
Channels = length(ChanLabels);

c = 1;
while(c<=length(varargin))
    switch(varargin{c})
        case {'Chan2Plot'}
            % Channel(s) to plot
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'Bands2Plot'}
            % Band(s) to plot
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
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
        case {'TBNum'} % number of row in the struture
            % 1*2 vector including block and trail number to plot [trail  block]
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'DownsampleRate'}
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
        case {'bandofInterest'}
            % the band you want to plot
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        otherwise
            error(sprintf('Unknown option: %s',varargin{c}));
    end
end
if ~exist('Chan2Plot')
    error('Define Channels --> Chan2Plot')
end
if sum(strcmp(what , {'binned_BlockGroup_AvgChann','binned_BlockGroup','raw_BlockGroup' , 'raw_BlockGroup_AvgChann'})) & ~exist('BlockGroup')
    error('Define Block Group to plot')
end
if  sum(strcmp(what , {'binned_SingleTrial','raw_SingleTrial'}))  & ~exist('TBNum')
    error('Define Block and trial number to plot')
end


Fs = 1024;
Fs_ds = floor(Fs/DownsampleRate);
min_freq =  FreqRange(1);
max_freq = FreqRange(2);
frex = linspace(min_freq, max_freq,numFreqBins);
%% load required data
mainDir = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' subjname{subjNum} ,'/'];
load([mainDir , 'ChanLabels.mat'])
ChanLabels = ChanLabels(Channels);
load([mainDir , 'AllData_Behav.mat'])
load([mainDir , 'AllData_Events.mat'])
load([mainDir , 'AllData_AvgMarker.mat'])

%% HouseKeeping

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

blockGroupTags = {'Natural Speed' , 'Slow Speed Day 1' , 'Slow Speed Day 2' ,'Slow Speed Day 3' ,'Slow Speed Day 4' ,...
    'Fast Speed Day 1' ,'Fast Speed Day 2' ,'Fast Speed Day 3' ,'Fast Speed Day 4',...
    ' Day 1' , ' Day 1' , ' Day 1' , ' Day 2' , ...
    ' Day 2' , ' Day 2' , ' Day 2' , ' Day 3' , ...
    ' Day 3' , ' Day 3' , ' Day 3' , ' Day 4'};
SeqTrans = [5 11 22 33 44 55 0 1 2 3 4 103 104 203 204;...
            100 10 10 10 10 10 20 30 30 30 30 40 50 40 50];
        
sequenceType.nums = [100 10 20 30 40 50];
sequenceType.tags = {'Null Trials' , 'Single Finger Sequences Trials' , 'Random Sequences' , 'Trained Sequences' , 'Triplet Segments' , 'Quadruple Segments'}';

BandInfo.bandsLab = {'Delta <4Hz' , 'Theta 5-8Hz' , 'Alpha 9-16Hz' , 'Beta 17-36Hz' , 'L-Gamma 37-70Hz' , 'H-Gamma 70-130' , 'NoBandLand 130-180'};
BandInfo.bands = {[0 4], [5 8] [9 16] [17 36] [37 70] [70 110] [110 180]};

switch what
    %% PLOT average normalized binned power
    case 'binned_SingleTrial'
        % input to this has to be AllData_PSD_StimNorm.mat
        Pall  = secog_addEventMarker(Pall, subjNum, Fs_ds , 'addEvent' , 'NumWarpSampFast' , NumWarpSampFast, 'NumWarpSampSlow'  ,NumWarpSampSlow)';
        D = getrow(Pall , Pall.BN == TBNum(2) & Pall.TN == TBNum(1));
        time = [1:size(D.Pow_Norm_stim{1} , 3)]/Fs_ds;
        figure('color' , 'white')
        Pow = real(squeeze(D.Pow_Norm_stim{1}(Chan2Plot,:,:)));
        EM = find(D.EventMarker{1})/Fs_ds;
        for ch = 1:length(Chan2Plot)
            subplot(length(Chan2Plot),1, ch)
            for b =1:length(BandInfo.bandsLab)
                B = real(squeeze(Pow(Chan2Plot(ch) , b , :))+8*b);
                plot(time , B , 'LineWidth' , 3)
                hold on
                T(b) = nanmean(B);
            end
            title (['Raw PSD for ',num2str(D.AllPress(1,:)) ', in Channel ' , ChanLabels{Chan2Plot(ch)}])
            for lin = 1:length(EM)
                line([EM(lin) EM(lin)] , [0 max(B)] , 'color' , 'red' , 'LineStyle' , ':' , 'LineWidth' , 3)
            end
            xlabel('Norm Time')
            
            set(gca , 'YTickLabels' , BandInfo.bandsLab, 'YTick' , T );
            
            line([median(EM) median(EM)] , [max(B)-5 max(B)],'color' , 'black' , 'LineWidth' , 5)
            text(median(EM),max(B),'5','FontSize' , 16 )
            set(gca ,'FontSize' , 16,'Box' , 'off')
        end
    case 'raw_SingleTrial'
        % input to this has to be AllData_PSD_StimNorm.mat
        % find the block number for that trial number to load up the
        % raw_PSD file
        Pall  = secog_addEventMarker(Pall, subjNum, Fs_ds , 'addEvent' , 'NumWarpSampFast' , NumWarpSampFast, 'NumWarpSampSlow'  ,NumWarpSampSlow)';
        D = getrow(Pall , Dall.BN == TBNum(2) & Dall.TN == TBNum(1));
        filename = [mainDir ,  'Raw_PSD_B' , num2str(TBNum(2)),'.mat'];
        trialName = ['PSD',num2str(TBNum(1))];
        load(filename , trialName)
        eval(['D.Pow_Norm_stim = ',num2str(trialName)]);
        time = [1:size(D.Pow_Norm_stim , 3)]/Fs_ds;
        figure('color' , 'white')
        Pow = real(squeeze(D.Pow_Norm_stim(Chan2Plot,:,:)));
        EM = find(D.EventMarker{1})/Fs_ds;
        for ch = 1:length(Chan2Plot)
            subplot(length(Chan2Plot),1, ch)
            B = real(squeeze(D.Pow_Norm_stim(Chan2Plot(ch),:,:)));
            contourf(time,frex,B,60,'linecolor','none')
            colorbar
%             imagesc(time,frex,B)
            % caxis([-15 7])
            title (['Raw PSD for ',num2str(D.AllPress(1,:)) ', in Channel ' , ChanLabels{Chan2Plot(ch)}])
            for lin = 1:length(EM)
                line([EM(lin) EM(lin)] , [0 max(frex)] , 'color' , 'red' , 'LineStyle' , ':' , 'LineWidth' , 3)
            end
            xlabel('SEC')
            ylabel('Frequency (Hz)')
            
            set(gca ,'FontSize' , 16,'Box' , 'off', 'YTick' , [5:5:150]);
        end
        %% PLOT average normalized binned power
    case 'raw_BlockGroup_AvgChann'
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped.mat
        % patterned PSDs so the output of Pall  = secog_parseEEG_PSD('TimeWarpPSD' , Dall, subjNum);
        E.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
        BG = find(strcmp(E.blockGroupNames , BlockGroup));
        E.NumWarpSamp([1 3 6 10 14]) = NumWarpSampFast;
        E = getrow(E , BG);
        if isempty(E.EM)
            error('Nothing to plot!')
        end
        
        filename = [mainDir ,  'AverageRawPSD' , num2str(BG),'.mat'];
        load(filename)
        figure('color' , 'white')
        figCount = 1;
        for sn = 1:length(E.SN{1})
            NEM = E.NEM{1}(sn , :);
            subplot(length(E.SN{1}) ,1, figCount)
            B = squeeze(real(nanmean(Pall.AvgPow{sn}(Chan2Plot,:,:) , 1)));
            % contourf([1:E.NumWarpSamp],frex,B,60,'linecolor','none')
            imagesc([1:E.NumWarpSamp],frex,B);
            % caxis([0 4])
            colorbar
            title (['Raw PSD for SeqNumb =  ', num2str(E.SN{1}(sn)),' in Blocks of ' , BlockGroup])
            hold on
            for lin = 1:length(NEM)
                line([NEM(lin) NEM(lin)] , [0  max(frex)] , 'color' , 'red' , 'LineStyle' , ':' , 'LineWidth' , 3)
            end
            xlabel('Norm Time')
            
            ylabel('Frequency (Hz)')
            
            set(gca ,'FontSize' , 16,'Box' , 'off')
            
            figCount = figCount + 1;
        end
    case 'binned_BlockGroup_AvgChann'
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped.mat
        % patterned PSDs so the output of Pall  = secog_parseEEG_PSD('TimeWarpPSD' , Dall, subjNum);
        
        Dall.PSD_stim = Pall.Pow_Norm_stim;
        Pall = Dall;
        Pall.Fast = zeros(size(Pall.TN));
        E.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
        E.NumWarpSamp([1 3 6 10 14]) = NumWarpSampFast;
        E = getrow(E , strcmp(E.blockGroupNames , BlockGroup));
        if isempty(E.EM)
            error('Nothing to plot!')
        end
        fastBlock = horzcat(blockGroups{1} ,blockGroups{3} , blockGroups{6}, blockGroups{10},blockGroups{14});
        Pall.Fast(ismember(Pall.BN , fastBlock)) = 1;
        clear Dall
        figure('color' , 'white')
        figCount = 1;
        for sn = 1:length(E.SN{1})
            NEM = E.NEM{1}(sn , :);
            id = ismember(Pall.BN , E.blockGroups{1}) & ismember(Pall.seqNumb , E.SN{1}(sn));
            F = getrow(Pall , id);
            % sum the warped PSDs insife the F structure
            tcount = 1;
            for tn = 1:length(F.PSD_stim)
                if isequal(size(F.PSD_stim{tn}) ,[length(ChanLabels) , length(BandInfo.bandsLab) , E.NumWarpSamp])
                    tempPow(tcount , :,:,:) = F.PSD_stim{tn};
                    tcount = tcount +1;
                end
            end
            AvgPow{sn} = real(squeeze(nanmean(tempPow , 1)));
            
            subplot(length(E.SN{1}) ,1, figCount)
            for b =1:length(BandInfo.bandsLab)
                B  = squeeze(nanmean(AvgPow{sn}(Chan2Plot , b , :) , 1))+7*b;
                plot([1:E.NumWarpSamp] , B , 'LineWidth' , 3)
                T(b) = nanmedian(B);
                hold on
            end
            hold off
            title (['Average Time-Warped PSD for ' ,E.blockGroupNames{1} , ' SeqNumb =  ', num2str(E.SN{1}(sn))])
            for lin = 1:length(NEM)
                line([NEM(lin) NEM(lin)] , [0 max(B)] , 'color' , 'red' , 'LineStyle' , ':' , 'LineWidth' , 3)
            end
            T = linspace(min(T) , max(T) , length(T));
            xlabel('Norm Time')
            set(gca ,'YTickLabels' , BandInfo.bandsLab, 'YTick' , T );
            
            line([median(NEM) median(NEM)] , [max(B)-5 max(B)],'color' , 'black' , 'LineWidth' , 5)
            text(median(NEM),max(B),'5','FontSize' , 16 )
            set(gca , 'XLim' , [1 E.NumWarpSamp] ,'FontSize' , 16,'Box' , 'off')
            figCount = figCount + 1;
        end
    case 'raw_BlockGroup'
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped.mat
        % patterned PSDs so the output of Pall  = secog_parseEEG_PSD('TimeWarpPSD' , Dall, subjNum);
        E.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
        BG = find(strcmp(E.blockGroupNames , BlockGroup));
        E.NumWarpSamp([1 6 7 8 9 12 16 20]) = NumWarpSampFast;
        E = getrow(E , BG);
        if isempty(E.EM)
            error('Nothing to plot!')
        end
        % return the sequence numbers to the original
        Pall.seqNumb = Dall.seqNumb;
        filename = [mainDir ,  'AverageRawPSD' , num2str(BG),'.mat'];
        load(filename)
        figure('color' , 'white')
        figCount = 1;
        for sn = 1:length(E.SN{1})
            if E.SN{1}(sn) ~= 5
                NEM = E.NEM{1}(sn ,  ~isnan(E.EM{1}(sn,:)));
                EM  = E.EM{1}(sn ,  ~isnan(E.EM{1}(sn,:)))/100;
                for ch = 1:length(Chan2Plot)
                    subplot(length(E.SN{1})-1 ,length(Chan2Plot), figCount)
                    B = real(squeeze(Pall.AvgPow{sn}(ch,:,:)));
                    contourf([1:E.NumWarpSamp],frex,100*B,60,'linecolor','none')
                    colorbar
                    caxis([-8 10])
                    title (['Raw PSD for SeqNumb =  ', num2str(E.SN{1}(sn)),', in Channel ' , ChanLabels{Chan2Plot(ch)}])
                    hold on
                    xtikx = {};
                    for lin = 1:length(NEM)
                        xtikx = [xtikx  , num2str(EM(lin))];
                        line([NEM(lin) NEM(lin)] , [0  max(frex)] , 'color' , 'red' , 'LineStyle' , ':' , 'LineWidth' , 3)
                    end
                    
                    xlabel('Norm Time')
                    
                    ylabel('Frequency (Hz)')
                    
                    set(gca ,'FontSize' , 10,'Box' , 'off' , 'YTick' , [10:10:max(frex)],'XTickLabels' , xtikx, 'XTick' , NEM )
                    
                    figCount = figCount + 1;
                end
            end
        end
    case 'binned_BlockGroup'
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped.mat
        % patterned PSDs so the output of Pall  = secog_parseEEG_PSD('TimeWarpPSD' , Dall, subjNum);
        
        Dall.PSD_stim = Pall.Pow_Norm_stim;
        Pall = Dall;
        Pall.Fast = zeros(size(Pall.TN));
        E.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
        E.NumWarpSamp([1 6 7 8 9 12 16 20]) = NumWarpSampFast;
        E = getrow(E , strcmp(E.blockGroupNames , BlockGroup));
        if isempty(E.EM)
            error('Nothing to plot!')
        end
        fastBlock = horzcat(blockGroups{1} , blockGroups{3} , blockGroups{6}, blockGroups{10},blockGroups{14});
        Pall.Fast(ismember(Pall.BN , fastBlock)) = 1;
        clear Dall
        figure('color' , 'white')
        figCount = 1;
        for sn = 1:length(E.SN{1})
            if E.SN{1}(sn) ~= 5
%                 stim = unique(Events.AllPress(Events.seqNumb == E.SN{1}(sn) , :) , 'rows');
                NEM = E.NEM{1}(sn , ~isnan(E.NEM{1}(sn,:)));
                % actual time in seconds to use for labeling the time axis
                EM  = E.EM{1}(sn ,  ~isnan(E.EM{1}(sn,:)))/100;
                
                id = ismember(Pall.BN , E.blockGroups{1}) & ismember(Pall.seqNumb , E.SN{1}(sn));
                F = getrow(Pall , id);
                % sum the warped PSDs insife the F structure
                tcount = 1;
                for tn = 1:length(F.PSD_stim)
                    if isequal(size(F.PSD_stim{tn}) ,[length(ChanLabels) , length(BandInfo.bandsLab) , E.NumWarpSamp])
                        tempPow(tcount , :,:,:) = F.PSD_stim{tn};
                        tcount = tcount +1;
                    end
                end
                AvgPow{sn} = real(squeeze(nanmean(tempPow , 1)));
                
                for ch = 1:length(Chan2Plot)
                    subplot(length(E.SN{1})-1 ,length(Chan2Plot), figCount)
                    for b =1:length(BandInfo.bandsLab)
                        B  = 100*squeeze(AvgPow{sn}(ch , b , :))+7*b;
                        plot([1:E.NumWarpSamp] , B , 'LineWidth' , 3)
                        T(b) = nanmedian(B);
                        hold on
                    end
                    hold off
                    title (['Average PSD for ' ,E.blockGroupNames{1}, ' Channel ' , ChanLabels{Chan2Plot(ch)} , 'Seqyence Type ' , num2str(E.SN{1}(sn))])
                    xtikx = {};
                    for lin = 1:length(NEM)
                        xtikx = [xtikx  , num2str(EM(lin))];
                        line([NEM(lin) NEM(lin)] , [0 max(B)] , 'color' , 'red' , 'LineStyle' , ':' , 'LineWidth' , 3)
                    end
                    T = linspace(min(T) , max(T) , length(T));
                    xlabel('Normalized Time (ms)')
                    set(gca ,'YTick' , T ,'YTickLabels' , BandInfo.bandsLab,'XTickLabels' , xtikx, 'XTick' , NEM ,...
                        'XLim' , [1 E.NumWarpSamp] ,'FontSize' , 10,'Box' , 'off');
                    line([max(NEM)-std(NEM) max(NEM)-std(NEM)] , [max(B)-1- max(B)],'color' , 'black' , 'LineWidth' , 5)
                    text(max(NEM),max(B),'1%','FontSize' , 16 )
                    figCount = figCount + 1;
                end
            end
        end
    case 'raw_BlockGroup_SeqType'
        load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped.mat
        % patterned PSDs so the output of Pall  = secog_parseEEG_PSD('TimeWarpPSD' , Dall, subjNum);
        E.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
        BG = find(strcmp(E.blockGroupNames , BlockGroup));
        E.NumWarpSamp([1 6 7 8 9 12 16 20]) = NumWarpSampFast;
        E = getrow(E , BG);
        if isempty(E.EM)
            error('Nothing to plot!')
        end
        filename = [mainDir ,  'AverageRawPSD_SeqType' , num2str(BG),'.mat'];
        load(filename)
        figure('color' , 'white')
        figCount = 1;
        for sn = 1:length(E.SN{1})
            if E.SN{1}(sn) ~= 100
                NEM = E.NEM{1}(sn ,  ~isnan(E.EM{1}(sn,:)));
                EM  = E.EM{1}(sn ,  ~isnan(E.EM{1}(sn,:)))/100;
                for ch = 1:length(Chan2Plot)
                    subplot(length(E.SN{1})-1 ,length(Chan2Plot), figCount)
                    B = squeeze(Pall.AvgPow{sn}(Chan2Plot(ch),:,:));
                    contourf([1:E.NumWarpSamp],frex,B,60,'linecolor','none')
                    colorbar
                    caxis([-0.08 0.07])
                    title (['Raw PSD for SeqNumb =  ', num2str(E.SN{1}(sn)),', in Channel ' , ChanLabels{Chan2Plot(ch)}])
                    hold on
                    xtikx = {};
                    for lin = 1:length(NEM)
                        xtikx = [xtikx  , num2str(EM(lin))];
                        line([NEM(lin) NEM(lin)] , [0  max(frex)] , 'color' , 'red' , 'LineStyle' , ':' , 'LineWidth' , 3)
                    end
                    xlabel('Norm Time')
                    ylabel('Frequency (Hz)')
                    set(gca ,'FontSize' , 10,'Box' , 'off' , 'YTick' , [10:10:max(frex)],'XTickLabels' , xtikx, 'XTick' , NEM )
                    figCount = figCount + 1;
                end
            end
        end
    case 'binned_BlockGroup_SeqType'
        load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
        % Pall is the AllData_PSD_Warped_SeqType.mat
        % patterned PSDs so the output of Pall  = secog_parseEEG_PSD('TimeWarpPSD_Raw_Binned_seqType' , Pall, subjNum);
        
        SeqTrans = [5 11 22 33 44 55 0 1 2 3 4 103 104 203 204;...
            100 10 10 10 10 10 20 30 30 30 30 40 50 40 50];
        for sn = 1:length(SeqTrans)
            id  = Pall.seqNumb == SeqTrans(1 , sn);
            Pall.seqNumb(id) = SeqTrans(2 , sn);
        end
        Pall.Fast = zeros(size(Pall.TN));
        E.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
        E.NumWarpSamp([1 6 7 8 9 12 16 20]) = NumWarpSampFast;
        E = getrow(E , strcmp(E.blockGroupNames , BlockGroup));
        if isempty(E.EM)
            error('Nothing to plot!')
        end
        fastBlock = horzcat(blockGroups{1} , blockGroups{3} , blockGroups{6}, blockGroups{10},blockGroups{14});
        Pall.Fast(ismember(Pall.BN , fastBlock)) = 1;
        figure('color' , 'white')
        figCount = 1;
        for sn = 1:length(E.SN{1})
            if E.SN{1}(sn) ~= 100
                stim = unique(Events.AllPress(Events.seqNumb == E.SN{1}(sn) , :) , 'rows');
                NEM = E.NEM{1}(sn , ~isnan(E.NEM{1}(sn,:)));
                % actual time in seconds to use for labeling the time axis
                EM  = E.EM{1}(sn ,  ~isnan(E.EM{1}(sn,:)))/100;
                
                id = ismember(Pall.BN , E.blockGroups{1}) & ismember(Pall.seqNumb , E.SN{1}(sn));
                F = getrow(Pall , id);
                % sum the warped PSDs insife the F structure
                tcount = 1;
                for tn = 1:length(F.Pow_Norm_stim)
                    if isequal(size(F.Pow_Norm_stim{tn}) ,[length(ChanLabels) , length(BandInfo.bandsLab) , E.NumWarpSamp])
                        tempPow(tcount , :,:,:) = F.Pow_Norm_stim{tn};
                        tcount = tcount +1;
                    end
                end
                AvgPow{sn} = real(squeeze(nanmean(tempPow , 1)));
                
                for ch = 1:length(Chan2Plot)
                    subplot(length(E.SN{1})-1 ,length(Chan2Plot), figCount)
                    for b =1:length(BandInfo.bandsLab)
                        B  = 10*squeeze(AvgPow{sn}(Chan2Plot(ch) , b , :))+.7*b;
                        plot([1:E.NumWarpSamp] , B , 'LineWidth' , 3)
                        T(b) = nanmedian(B);
                        hold on
                    end
                    hold off
                    title (['Average PSD for ' ,E.blockGroupNames{1}, ' Channel ' , ChanLabels{Chan2Plot(ch)} ,', ' , num2str(stim)])
                    xtikx = {};
                    for lin = 1:length(NEM)
                        xtikx = [xtikx  , num2str(EM(lin))];
                        line([NEM(lin) NEM(lin)] , [0 max(B)] , 'color' , 'red' , 'LineStyle' , ':' , 'LineWidth' , 3)
                    end
                    T = linspace(min(T) , max(T) , length(T));
                    xlabel('Normalized Time (ms)')
                    set(gca ,'YTick' , T ,'YTickLabels' , BandInfo.bandsLab,'XTickLabels' , xtikx, 'XTick' , NEM );
                    
                    line([max(NEM)-std(NEM) max(NEM)-std(NEM)] , [max(B)-1- max(B)],'color' , 'black' , 'LineWidth' , 5)
                    text(max(NEM),max(B),'1%','FontSize' , 16 )
                    set(gca , 'XLim' , [1 E.NumWarpSamp] ,'FontSize' , 10,'Box' , 'off')
                    figCount = figCount + 1;
                end
            end
        end
    case 'AvgPower_SeqType'
        load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped.mat
        % patterned PSDs so the output of Pall  = secog_parseEEG_PSD('TimeWarpPSD' , Dall, subjNum);
        % for this case because we are doigna comparisn, the BlockGroup is
        % a 1xN cell containnig the groups to be compared
        E.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
        BG = find(strcmp(E.blockGroupNames , BlockGroup));
        E.NumWarpSamp([1 6 7 8 9 12 16 20]) = NumWarpSampFast;
        E = getrow(E , BG);
        if isempty(E.EM)
            error('Nothing to plot!')
        end
        filename = [mainDir ,  'AverageSpect_SeqType' , num2str(BG),'.mat'];
        load(filename)
        figure('color' , 'white')
        chanGroup = {Chan2Plot};
        fCont = 1;
        for sn = 1:length(Pall.AvgPowTR)
            Tr = Pall.AvgPowTR{sn};
            Tre = Pall.SePowTR{sn};
            
            Bl = Pall.AvgPowBL{sn};
            Ble = Pall.SePowBL{sn};
            
            for cg = 1:length(chanGroup)
                APtr{cg}(sn,:) = nanmean(Tr(chanGroup{cg}, :) , 1);
                AEtr{cg}(sn,:) = nanmean(Tre(chanGroup{cg}, :) , 1);
                APbl{cg}(sn,:) = nanmean(Bl(chanGroup{cg}, :) , 1);
                AEbl{cg}(sn,:) = nanmean(Ble(chanGroup{cg}, :) , 1);
                subplot(length(chanGroup) , length(Pall.AvgPowTR) ,fCont)
                title(['sn = ' , num2str(Pall.SN{1}(sn))])
                h1 = plotshade(frex , APtr{cg}(sn,:) ,AEtr{cg}(sn,:),'transp' , .2 , 'patchcolor' , 'b' , 'linecolor' , 'b' , 'linewidth' , 3 , 'linestyle' , ':')
                hold on
                h2 = plotshade(frex , APbl{cg}(sn,:) ,AEbl{cg}(sn,:),'transp' , .2 , 'patchcolor' , 'r' , 'linecolor' , 'r' , 'linewidth' , 3 , 'linestyle' , ':')
                title(['Average Powerfor sequence type ' , num2str(E.SN{1}(sn)) ])
                fCont = fCont +1;
                legend([h1,h2] , {'Trial' , 'Baseline'})
                set(gca , 'YLim' , [30 80]);
            end
        end
    case 'AvgPower_SeqType_comp'
        colorz = {[0 0 1],[1 0 0],[0 1 0],[1 0 1],[0 1 1],[0.7 0.7 0.7],[1 1 0],[.3 .3 .3]};
        load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped.mat
        % patterned PSDs so the output of Pall  = secog_parseEEG_PSD('TimeWarpPSD' , Dall, subjNum);
        % for this case because we are doigna comparisn, the BlockGroup is
        % a 1xN cell containnig the groups to be compared
        
        for n = 1:length(BlockGroup)
            bgid = find(strcmp(blockGroupNames , BlockGroup{n}));
 
            BlG(n) = find(strcmp(blockGroupNames , BlockGroup{n}));
        
            filename = [mainDir ,  'PowerSpectrum_SeqType' , num2str(BlG(n)),'.mat'];
            load(filename)
            chanGroup = {Chan2Plot};
            fig = figure;
            for sn = 1:length(Pall.PowTR)
                Tr = Pall.PowTR{sn};
                Bl = Pall.PowBL{sn};
                for cg = 1:length(chanGroup)
                    Tcg = Tr(:,chanGroup{cg} , :);
                    Bcg = Bl(:,chanGroup{cg} , :);
                    diff  = 100*squeeze(mean((Tcg - Bcg)./Bcg,2)); % percent change from baseline
                    bandLabel = repmat(frex' , size(diff,1) , 1);
                    data{cg , n}{sn} = reshape(diff' , numel(diff) , 1);
                    [x{cg,n}{sn} , p{cg,n}{sn} , e{cg,n}{sn}] = lineplot(bandLabel , data{cg , n}{sn} , 'plotfcn' , 'nanmean');
%                     diff{sn}(n, :) = 100*(APtr{cg,n}(sn,:) - APbl{cg,n}(sn,:))./ APbl{cg,n}(sn,:);
                end
            end
            close(fig)
        end
        figure('color' , 'white')
        fCont = 1;
        
        for sn = 1:length(Pall.PowTR)
            subplot(1,length(Pall.PowTR),sn)
            leg = [];
            legen = {};
            rectangle('position' , [17 -4 19 10] , 'EdgeColor' , 'none' , 'FaceColor' , [.9 .9 .9])
            text(19 ,5 , 'Beta' , 'FontSize' , 30)
            hold on
            rectangle('position' , [70 -4 60 10] , 'EdgeColor' , 'none' , 'FaceColor' , [.9 .9 .9])
            text(80 ,5 , 'High Gamma' , 'FontSize' , 30)
            for n = 1:length(BlockGroup)
                [n sn]
                snid = find(sequenceType.nums == E.SN{BlG(n)}(sn));
                bgid = find(strcmp(blockGroupNames , BlockGroup{n}));
                legen = [legen , blockGroupTags{bgid}];
                h(n).fig = plotshade(x{cg,n}{sn}' , p{cg,n}{sn} , e{cg,n}{sn},'transp' , .2 , 'patchcolor' , colorz{n} , 'linecolor' , colorz{n} , 'linewidth' , 3 , 'linestyle' , ':');
                leg = [leg , h(n).fig];
                hold on
                %                 grid on
            end
            legend(leg , legen , 'Box' , 'off')
            title(['Percent Change in Power from Baseline in ' , sequenceType.tags{snid} ])
            set(gca , 'YLim',[-8 8] , 'FontSize' , 16 , 'Box' , 'off' , 'YTick' , [-4 0 4]);
            ylabel('(%)')
            xlabel('Frequency (Hz)')
            fCont = fCont +1;
        end
    case 'Compare_SeqType_TemporalPattern'
        colorz = {[0 0  1],[1 0 0],[0 1 0],[1 0 1],[0 1 1],[0.7 0.7 0.7],[1 1 0],[.3 .3 .3]};
        load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped.mat
        % patterned PSDs so the output of Pall  = secog_parseEEG_PSD('TimeWarpPSD' , Dall, subjNum);
        % for this case because we are doigna comparisn, the BlockGroup is
        % a 1xN cell containnig the groups to be compared
        if ~exist('bandofInterest')
            error('Define band to  plot -- > bandofInterest')
        end

        D = [];
        D.data = [];
        D.blockgroup = [];
        D.band  = [];
        D.seqNumb    = [];
        D.changroup  = [];
        
        
        for bg = 1:length(BlockGroup)
            load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
            E.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
            E.NumWarpSamp([1 6 7 8 9 12 16 20]) = NumWarpSampFast;
            BG(bg) = find(strcmp(E.blockGroupNames , BlockGroup{bg}));
            E = getrow(E , BG(bg));
            if isempty(E.EM)
                error('Nothing to plot!')
            end
            filename = [mainDir ,  'PowerSpectrumBinned_TimeKept_SeqType' , num2str(BG(bg)),'.mat'];
            load(filename , ['PSD_band' , num2str(bandofInterest)]);
            eval(['psd_b = PSD_band' , num2str(bandofInterest)]);
            chanGroup = {Chan2Plot};
            for sn = 1:length(psd_b)
                [sn bg]
                for cg = 1:length(chanGroup)
                    % average over channels in the chagroup
                    temp = squeeze(nanmean(psd_b{sn}(: ,chanGroup{cg},:) , 2));
                    D.data = [D.data ; temp];
                    D.changroup  = [D.changroup; cg*ones(size(temp,1) ,1)];
                    D.band  = [D.band ; bandofInterest*ones(size(temp , 1) ,1)];
                    D.seqNumb    = [D.seqNumb; E.SN{1}(sn)*ones(size(temp , 1),1)];
                    D.blockgroup = [D.blockgroup; bg*ones(size(temp , 1),1)];
                end
            end
        end
        
        %% prepare data for testing Block Groups against eachother
        fig = figure;
        for bg = 1:length(BlockGroup)
            load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
            BG(bg) = find(strcmp(E.blockGroupNames , BlockGroup{bg}));
            E = getrow(E , BG(bg));
            for sn = 1:length(psd_b)
                T = getrow(D , D.seqNumb == E.SN{1}(sn) & D.blockgroup == bg);
                data = reshape(T.data'  , numel(T.data) , 1);
                timeLabel = repmat([1:500]' , size(T.data , 1) , 1);
                [x{sn,bg} , p{sn,bg} , e{sn,bg}] = lineplot(timeLabel , 100*data , 'plotfcn' , 'nanmean' , 'style_shade'); % percent change
                hold on
            end
        end
        close(fig)
        
        figure('color' , 'white')
        figCount = 1;
        for bg = 1:length(BlockGroup)
            load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
            BG(bg) = find(strcmp(E.blockGroupNames , BlockGroup{bg}));
            BlNames{bg} = blockGroupTags{BG(bg)};
            E = getrow(E , BG(bg));
            leg = [];
            legen = {};
            for sn = 1:length(E.SN{1})
                if E.SN{1}(sn)~=100
                    NEM = E.NEM{1}(sn , ~isnan(E.NEM{1}(sn,:)));
                    % actual time in seconds to use for labeling the time axis
                    EM  = E.EM{1}(sn ,  ~isnan(E.EM{1}(sn,:)))/100;
                    subplot(1,length(BlockGroup), figCount)
                    h(sn).fig = plotshade(x{sn,bg}' , p{sn,bg} , e{sn,bg},'transp' , .2 , 'patchcolor' , colorz{sn} , 'linecolor' , colorz{sn} , 'linewidth' , 3 , 'linestyle' , ':');
                    leg = [leg , h(sn).fig];
                    hold on
                    grid on
                    snid = find(sequenceType.nums == E.SN{1}(sn));
                    bgid = find(strcmp(blockGroupNames , BlockGroup{bg}));
                    legen = [legen , sequenceType.tags{snid}];
                    xlabel('Normalized Time (ms)')
                    set(gca ,'YLim' , [-8 +3] , 'FontSize' , 16);
                end
            end
            legend(leg , legen)
            title (['Average PSD on' ,blockGroupTags{bgid}])
            figCount = figCount + 1;
        end
        
        figure('color' , 'white')
        figCount = 1;
        for sn = 1:length(E.SN{1})
            load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
            BG(bg) = find(strcmp(E.blockGroupNames , BlockGroup{bg}));
            E = getrow(E , BG(bg));
            leg = [];
            legen = {};
            for bg = 1:length(BlockGroup)
                if E.SN{1}(sn)~=100
                    NEM = E.NEM{1}(sn , ~isnan(E.NEM{1}(sn,:)));
                    % actual time in seconds to use for labeling the time axis
                    EM  = E.EM{1}(sn ,  ~isnan(E.EM{1}(sn,:)))/100;
                    snid = find(sequenceType.nums == E.SN{1}(sn));
                    bgid = find(strcmp(blockGroupNames , BlockGroup{bg}));
                    subplot(1,length(E.SN{1})-1, figCount)
                    h(sn).fig = plotshade(x{sn,bg}' , p{sn,bg} , e{sn,bg},'transp' , .2 , 'patchcolor' , colorz{bg} , 'linecolor' , colorz{bg} , 'linewidth' , 3 , 'linestyle' , ':');
                    leg = [leg , h(sn).fig];
                    hold on
                    grid on
                    legen = [legen , blockGroupTags{bgid}];

                    xlabel('Normalized Time (ms)')
                    set(gca ,'YLim' , [-8 +3] , 'FontSize' , 16);
                end
            end
            legend(leg , legen)
            figCount = figCount + 1;
            title (['Average PSD for ' , sequenceType.tags{snid}])
        end
    case 'Aligned_SeqType'
        colorz = {[0 0  1],[1 0 0],[0 1 0],[1 0 1],[0 1 1],[0.7 0.7 0.7],[1 1 0],[.3 .3 .3]};
        load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped.mat
        % patterned PSDs so the output of Pall  = secog_parseEEG_PSD('TimeWarpPSD' , Dall, subjNum);
        % for this case because we are doigna comparisn, the BlockGroup is
        % a 1xN cell containnig the groups to be compared
        blockGroups = {[1 2] , [3], [13], [26], [40] , [4], [14], [27] [41] , [5:7] , [9:11] , [8 12] , [15:17] , [19:21] , [23:25],...
            [18 22] , [28:30] , [32:34] , [36:38], [31 35 39],[42:44]}';
        blockGroupNames = {'SingleFingNat' , 'SingleFingSlow1' , 'SingleFingSlow2'  , 'SingleFingSlow3' ,'SingleFingSlow4',...
            'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3', 'SingleFingFast4' , 'Intermixed1' , 'Intermixed2' , ...
            'ChunkDay1' , 'Intermixed3' , 'Intermixed4' , 'Intermixed5', 'ChunkDay2' , 'Intermixed6' , ...
            'Intermixed7' , 'Intermixed8', 'ChunkDay3', 'Intermixed9'}';

        if ~exist('bandofInterest')
            error('Define band to  plot -- > bandofInterest')
        end

        C = [];
        T = [];
        
        if sum(ismember (BlockGroup , blockGroupNames([12 16 20])))
            seq = input('Triplet or Quadruple? (t/q)', 's');
            switch seq
                case {'q'}
                    selectSN = 50;
                case {'t'}
                    selectSN = 40;
            end
        end
        for bg = 1:length(BlockGroup)
            load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
           
            BlG(bg) = find(strcmp(blockGroupNames , BlockGroup{bg}));
            E = getrow(E , BlG(bg));
            if isempty(E.EM)
                error('Nothing to plot!')
            end
            filename = [mainDir ,  'Group_Aligned_PSD' , num2str(BlG(bg)),'.mat'];
            load(filename);
            Pall = getrow(Pall , Pall.seqNumb ~= 100 & ~Pall.isError);
            SN = unique(Pall.seqNumb);
            if exist('selectSN')
              SN =   selectSN;
            end
            for sn = 1:length(SN)    
                T = getrow(Pall , Pall.seqNumb == SN(sn));
                % reduce the band and channel dimensions
                for evnt = 1:T.seqlength(1)+1
                    for tn = 1:length(T.TN)
                        if~isempty(T.PSD{tn , evnt})
                            temp = squeeze(T.PSD{tn , evnt}(:,bandofInterest,:));
                            T.reducPSD{tn , evnt} = nanmean(temp(Chan2Plot , :) , 1);
                        else
                            T.reducPSD{tn , evnt} = nan(1, 5);
                        end
                    end
                end
                T.blockgroup = bg*ones(size(T.BN));
                C = addstruct(C , T);
            end
        end
        for sn = 1:length(SN)
            D = getrow(C , C.seqNumb == SN(sn));
            for evnt = 1:T.seqlength(1)+1
                medEventLen = floor(median(cellfun(@length, D.reducPSD(:,evnt))));
                maxLen = max(cellfun(@length, D.reducPSD(:,evnt)));
                for tn = 1:length(D.TN)
                    D.reducPSD{tn , evnt} = padarray(D.reducPSD{tn , evnt} , [0 maxLen-length(D.reducPSD{tn , evnt})] , NaN , 'post');
                end
                eval(['D.stackPSD' , num2str(evnt),' = cell2mat(D.reducPSD(: , evnt));']);
                eval(['D.stackPSD' , num2str(evnt),' = D.stackPSD' , num2str(evnt),'(:,1:medEventLen)']);
            end
            All(sn)  = D;
        end
        %% prepare data for testing Block Groups against eachother
        
        for sn = 1:length(SN)
            fig = figure;
            D = All(sn);
            for bg = 1:length(BlockGroup)
                T = getrow(D , D.blockgroup == bg);
                for evnt = 1:T.seqlength(1)+1
                    eval(['temp = T.stackPSD' , num2str(evnt),';']);
                    data = reshape(temp' , numel(temp) , 1);
                    timeStamp = repmat([1:size(temp , 2)]' , size(temp , 1) , 1);
                    [x{sn,bg,evnt} , p{sn,bg,evnt} , e{sn,bg,evnt}] = lineplot(timeStamp , data , 'plotfcn' , 'nanmean' , 'style_shade'); % percent change
                end
            end
            
            V = evnt;
            close(fig)
        end
%         P  = cellfun(@nanmean , p);
%         E = cellfun(@nanmean , e);
%         errorbar([1:4] , squeeze(mean(P(2,:,:) , 3)), squeeze(mean(E(2,:,:) , 3)))
        figure('color' , 'white')
        figCount = 1;
        for sn = 1:length(SN)
            snid = find(sequenceType.nums == E.SN{1}(sn));
            for evnt = 1:T.seqlength(1)+1
                subplot(length(SN),V , figCount);
                leg = [];
                legen = {};
                for bg = 1:length(BlockGroup)
                    bgid = find(strcmp(blockGroupNames , BlockGroup{bg}));
                    h(sn).fig = plotshade(x{sn,bg,evnt}' , p{sn,bg,evnt} , e{sn,bg,evnt},'transp' , .2 , 'patchcolor' , colorz{bg} , 'linecolor' , colorz{bg} , 'linewidth' , 3 , 'linestyle' , ':');
                    leg = [leg , h(sn).fig];
                    hold on
                    grid on
                    legen = [legen , blockGroupTags{bgid}];
                end
                figCount = figCount+1;
                
                title (['Press ' ,num2str(evnt-1),' PSD for ' , sequenceType.tags{snid}])
                set(gca , 'Box' , 'off' , 'YLim' , [-8 , 8]);
                line([10,10] , [min(p{sn,bg,evnt}) max(p{sn,bg,evnt})] , 'LineWidth' , 3 , 'color' , 'k')
            end
        end
        legend(leg , legen)
        
        figure('color' , 'white')
        figCount = 1;
        for bg = 1:length(BlockGroup)
            bgid = find(strcmp(blockGroupNames , BlockGroup{bg}));
            for evnt = 1:T.seqlength(1)+1
                subplot(length(BlockGroup),V , figCount);
                leg = [];
                legen = {};
                for sn = 1:length(SN)
                    snid = find(sequenceType.nums == E.SN{1}(sn));
                    h(sn).fig = plotshade(x{sn,bg,evnt}' , p{sn,bg,evnt} , e{sn,bg,evnt},'transp' , .2 , 'patchcolor' , colorz{sn} , 'linecolor' , colorz{sn} , 'linewidth' , 3 , 'linestyle' , ':');
                    leg = [leg , h(sn).fig];
                    hold on
                    grid on
                    legen = [legen , sequenceType.tags{snid}];
                end
                figCount = figCount+1;
                
                title (['Press ' ,num2str(evnt-1),' PSD for ' , blockGroupTags{bgid}])
                set(gca , 'Box' , 'off' , 'YLim' , [-8 , 8])
                line([10,10] , [min(p{sn,bg,evnt}) max(p{sn,bg,evnt})] , 'LineWidth' , 3 , 'color' , 'k')
            end
        end 
        legend(leg , legen)
        
        
        
    case 'AlignedWarped_SeqType'
        colorz = {[0 0  1],[1 0 0],[0 1 0],[1 0 1],[0 1 1],[0.7 0.7 0.7],[1 1 0],[.3 .3 .3]};
        load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped.mat
        % patterned PSDs so the output of Pall  = secog_parseEEG_PSD('TimeWarpPSD' , Dall, subjNum);
        % for this case because we are doigna comparisn, the BlockGroup is
        % a 1xN cell containnig the groups to be compared
        blockGroups = {[1 2] , [3], [13], [26], [40] , [4], [14], [27] [41] , [5:7] , [9:11] , [8 12] , [15:17] , [19:21] , [23:25],...
            [18 22] , [28:30] , [32:34] , [36:38], [31 35 39],[42:44]}';
        blockGroupNames = {'SingleFingNat' , 'SingleFingSlow1' , 'SingleFingSlow2'  , 'SingleFingSlow3' ,'SingleFingSlow4',...
            'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3', 'SingleFingFast4' , 'Intermixed1' , 'Intermixed2' , ...
            'ChunkDay1' , 'Intermixed3' , 'Intermixed4' , 'Intermixed5', 'ChunkDay2' , 'Intermixed6' , ...
            'Intermixed7' , 'Intermixed8', 'ChunkDay3', 'Intermixed9'}';
        
        if ~exist('bandofInterest')
            error('Define band to  plot -- > bandofInterest')
        end

        C = [];
        T = [];
        
        if sum(ismember (BlockGroup , blockGroupNames([12 16 20])))
            seq = input('Triplet or Quadruple? (t/q)', 's');
            switch seq
                case {'q'}
                    selectSN = 50;
                case {'t'}
                    selectSN = 40;
            end
        end
        for bg = 1:length(BlockGroup)
            load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
          
            BlG(bg) = find(strcmp(blockGroupNames , BlockGroup{bg}));
            E = getrow(E , BlG(bg));
            if isempty(E.EM)
                error('Nothing to plot!')
            end
            filename = [mainDir ,  'Group_WarpedAligned_PSD' , num2str(BlG(bg)),'.mat'];
            load(filename);
            for i = 1:length(Pall.TN)
                if sum(Pall.IPI(i,:)<120)
                    Pall.isError(i) =1;
                end
            end
            Pall = getrow(Pall , Pall.seqNumb ~= 100 & ~Pall.isError);
            SN = unique(Pall.seqNumb);
            if exist('selectSN')
                SN =   selectSN;
            end
            for sn = 1:length(SN)    
                T = getrow(Pall , Pall.seqNumb == SN(sn));
                % reduce the band and channel dimensions
                for evnt = 1:T.seqlength(1)+1
                    for tn = 1:length(T.TN)
                        [bg sn evnt tn]
                        if~isempty(T.PSD{tn , evnt})
                            temp = squeeze(T.PSD{tn , evnt}(:,bandofInterest,:));
                            T.reducPSD{tn , evnt} = nanmean(temp(Chan2Plot , :) , 1);
                        else
                            T.reducPSD{tn , evnt} = nan(1, 5);
                        end
                    end
                end
                T.blockgroup = bg*ones(size(T.BN));
                C = addstruct(C , T);
            end
        end
        for sn = 1:length(SN)
            D = getrow(C , C.seqNumb == SN(sn));
            for evnt = 1:T.seqlength(1)+1
                medEventLen = floor(median(cellfun(@length, D.reducPSD(:,evnt))));
                maxLen = max(cellfun(@length, D.reducPSD(:,evnt)));
                for tn = 1:length(D.TN)
                    D.reducPSD{tn , evnt} = padarray(D.reducPSD{tn , evnt} , [0 maxLen-length(D.reducPSD{tn , evnt})] , NaN , 'post');
                end
                eval(['D.stackPSD' , num2str(evnt),' = cell2mat(D.reducPSD(: , evnt));']);
                eval(['D.stackPSD' , num2str(evnt),' = D.stackPSD' , num2str(evnt),'(:,1:medEventLen)']);
            end
            All(sn)  = D;
        end
        %% prepare data for testing Block Groups against eachother
        
        for sn = 1:length(SN)
            fig = figure;
            D = All(sn);
            for bg = 1:length(BlockGroup)
                T = getrow(D , D.blockgroup == bg);
                for evnt = 1:T.seqlength(1)+1
                    eval(['temp = T.stackPSD' , num2str(evnt),';']);
                    data = reshape(temp' , numel(temp) , 1);
                    timeStamp = repmat([1:size(temp , 2)]' , size(temp , 1) , 1);
                    [x{sn,bg,evnt} , p{sn,bg,evnt} , e{sn,bg,evnt}] = lineplot(timeStamp , data , 'plotfcn' , 'nanmean' , 'style_shade'); % percent change
                end
            end
            
            V = evnt;
            close(fig)
        end
        %         P  = cellfun(@nanmean , p);
        %         E = cellfun(@nanmean , e);
        %         errorbar([1:4] , squeeze(mean(P(2,:,:) , 3)), squeeze(mean(E(2,:,:) , 3)))
        figure('color' , 'white')
        figCount = 1;
        for sn = 1:length(SN)
            snid = find(sequenceType.nums == E.SN{1}(sn));
            for evnt = 1:T.seqlength(1)+1
                subplot(length(SN),V , figCount);
                leg = [];
                legen = {};
                for bg = 1:length(BlockGroup)
                    bgid = find(strcmp(blockGroupNames , BlockGroup{bg}));
                    h(sn).fig = plotshade(x{sn,bg,evnt}' , p{sn,bg,evnt} , e{sn,bg,evnt},'transp' , .2 , 'patchcolor' , colorz{bg} , 'linecolor' , colorz{bg} , 'linewidth' , 3 , 'linestyle' , ':');
                    leg = [leg , h(sn).fig];
                    hold on
                    grid on
                    legen = [legen , blockGroupTags{bgid}];
                end
                figCount = figCount+1;
                
                title (['Press ' ,num2str(evnt-1),' PSD for ' , sequenceType.tags{snid}])
                %                 set(gca , 'Box' , 'off' , 'YLim' , [-8 , 8]);
                line([5,5] , [min(p{sn,bg,evnt}) max(p{sn,bg,evnt})] , 'LineWidth' , 3 , 'color' , 'k')
            end
        end
        legend(leg , legen)
        
        figure('color' , 'white')
        figCount = 1;
        for bg = 1:length(BlockGroup)
            bgid = find(strcmp(blockGroupNames , BlockGroup{bg}));
            for evnt = 1:T.seqlength(1)+1
                subplot(length(BlockGroup),V , figCount);
                leg = [];
                legen = {};
                for sn = 1:length(SN)
                    snid = find(sequenceType.nums == E.SN{1}(sn));
                    h(sn).fig = plotshade(x{sn,bg,evnt}' , p{sn,bg,evnt} , e{sn,bg,evnt},'transp' , .2 , 'patchcolor' , colorz{sn} , 'linecolor' , colorz{sn} , 'linewidth' , 3 , 'linestyle' , ':');
                    leg = [leg , h(sn).fig];
                    hold on
                    grid on
                    legen = [legen , sequenceType.tags{snid}];
                end
                figCount = figCount+1;
                
                title (['Press ' ,num2str(evnt-1),' PSD for ' , blockGroupTags{bgid}])
                %                 set(gca , 'Box' , 'off' , 'YLim' , [-8 , 8])
                line([5,5] , [min(p{sn,bg,evnt}) max(p{sn,bg,evnt})] , 'LineWidth' , 3 , 'color' , 'k')
            end
        end
        legend(leg , legen)
        
    case 'ChunkAligned_SeqType'
        colorz = {[0 0  1],[1 0 0],[0 1 0],[1 0 1],[0 1 1],[0.7 0.7 0.7],[1 1 0],[.3 .3 .3] , [.4 .2 .5]};
        load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped.mat
        % patterned PSDs so the output of Pall  = secog_parseEEG_PSD('TimeWarpPSD' , Dall, subjNum);
        % for this case because we are doigna comparisn, the BlockGroup is
        % a 1xN cell containnig the groups to be compared
        blockGroups = {[1 2] , [3], [13], [26], [40] , [4], [14], [27] [41] , [5:7] , [9:11] , [8 12] , [15:17] , [19:21] , [23:25],...
            [18 22] , [28:30] , [32:34] , [36:38], [31 35 39],[42:44]}';
        
        
        if ~exist('bandofInterest')
            error('Define band to  plot -- > bandofInterest')
        end
        
        SNofInterest = sequenceType.nums(find(strcmp(sequenceType.tags , 'Trained Sequences' )));
        
        C = [];
        T = [];
        
        if sum(ismember (BlockGroup , blockGroupNames([12 16 20])))
            seq = input('Triplet or Quadruple? (t/q)', 's');
            switch seq
                case {'q'}
                    selectSN = 50;
                case {'t'}
                    selectSN = 40;
            end
        end
        for bg = 1:length(BlockGroup)
            load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
            
            BlG(bg) = find(strcmp(blockGroupNames , BlockGroup{bg}));
            E = getrow(E , BlG(bg));
            if isempty(E.EM)
                error('Nothing to plot!')
            end
            filename = [mainDir ,  'Group_Aligned_PSD' , num2str(BlG(bg)),'.mat'];
            load(filename);
            Pall = getrow(Pall , Pall.seqNumb ~= 100 & ~Pall.isError);
            SN = unique(Pall.seqNumb);
            if exist('selectSN')
                SN =   selectSN;
            end
            T = getrow(Pall , Pall.seqNumb == SNofInterest);
            % reduce the band and channel dimensions
            for evnt = 1:T.seqlength(1)+1
                for tn = 1:length(T.TN)
                    temp = squeeze(T.PSD{tn , evnt}(:,bandofInterest,:));
                    T.reducPSD{tn , evnt} = nanmean(temp(Chan2Plot , :) , 1);
                end
            end
            T.blockgroup = bg*ones(size(T.BN));
            C = addstruct(C , T);
        end
        D = C;
        for tn = 1:length(D.TN)
            chnkBnd = find(D.ChnkPlcmnt(tn , 2:end) == 1);
            chnkWth  = [chnkBnd - 1 ,  chnkBnd+1];
            chnkInd  = [chnkBnd , chnkWth];
            cntr = 1;
            for evnt = chnkInd
                % the frist element is always the chunk boundry and the
                % next two are the one before and one after
                D.chunkArrPSD{tn , cntr} = D.reducPSD{tn , evnt};
                cntr =  cntr+1;
            end
        end
        titleLab = {'Chunk Boundry Press'};
        titleLab = [titleLab  , repmat({'Before Chunk Boundry' , 'After Chunk Boundry'} , 1 , .5*(length(chnkInd)-1))];
        for evnt = 1:size(D.chunkArrPSD , 2)
            medEventLen = floor(median(cellfun(@length, D.chunkArrPSD(:,evnt))));
            maxLen = max(cellfun(@length, D.chunkArrPSD(:,evnt)));
            for tn = 1:length(D.TN)
                D.chunkArrPSD{tn , evnt} = padarray(D.chunkArrPSD{tn , evnt} , [0 maxLen-length(D.chunkArrPSD{tn , evnt})] , NaN , 'post');
                
            end
            eval(['D.Chunk_stackPSD' , num2str(evnt),' = cell2mat(D.chunkArrPSD(: , evnt));']);
            eval(['D.Chunk_stackPSD' , num2str(evnt),' = D.Chunk_stackPSD' , num2str(evnt),'(:,1:medEventLen)']);
        end

        %% prepare data for testing Block Groups against eachother
        
        fig = figure;
        for bg = 1:length(BlockGroup)
            T = getrow(D , D.blockgroup == bg);
            for evnt = 1:size(D.chunkArrPSD  , 2)
                eval(['temp = T.Chunk_stackPSD' , num2str(evnt),';']);
                data = reshape(temp' , numel(temp) , 1);
                timeStamp = repmat([1:size(temp , 2)]' , size(temp , 1) , 1);
                [x{bg,evnt} , p{bg,evnt} , e{bg,evnt}] = lineplot(timeStamp , data , 'plotfcn' , 'nanmean' , 'style_shade'); % percent change
            end
        end
        
        V = evnt;
        close(fig)
        %         P  = cellfun(@nanmean , p);
        %         E = cellfun(@nanmean , e);
        %         errorbar([1:4] , squeeze(mean(P(2,:,:) , 3)), squeeze(mean(E(2,:,:) , 3)))
        figure('color' , 'white')
        figCount = 1;
            for evnt = 1:size(D.chunkArrPSD  , 2)
                subplot(1,V , figCount);
                leg = [];
                legen = {};
                for bg = 1:length(BlockGroup)
                    bgid = find(strcmp(blockGroupNames , BlockGroup{bg}));
                    h.fig = plotshade(x{bg,evnt}' , p{bg,evnt} , e{bg,evnt},'transp' , .2 , 'patchcolor' , colorz{bg} , 'linecolor' , colorz{bg} , 'linewidth' , 3 , 'linestyle' , ':');
                    leg = [leg , h.fig];
                    hold on
                    grid on
                    legen = [legen , blockGroupTags{bgid}];
                end
                figCount = figCount+1;
                
                title (titleLab{evnt})
                set(gca , 'Box' , 'off' , 'YLim' , [-8 , 8]);
                line([10,10] , [min(p{bg,evnt}) max(p{bg,evnt})] , 'LineWidth' , 3 , 'color' , 'k')
            end
        legend(leg , legen)
        
        figure('color' , 'white')
        figCount = 1;
        for bg = 1:length(BlockGroup)
            bgid = find(strcmp(blockGroupNames , BlockGroup{bg}));
            subplot(length(BlockGroup),1 , figCount);
            hold on
             leg = [];
            for evnt = 1:size(D.chunkArrPSD  , 2)
               
                legen = {};
                h(evnt).fig = plotshade(x{bg,evnt}' , p{bg,evnt} , e{bg,evnt},'transp' , .2 , 'patchcolor' , colorz{evnt} , 'linecolor' , colorz{evnt} , 'linewidth' , 3 , 'linestyle' , ':');
                leg = [leg , h(evnt).fig];
                hold on
                grid on
                
                title (['PSD for ' , blockGroupTags{bgid}])
                set(gca , 'Box' , 'off' , 'YLim' , [-8 , 8])
                line([10,10] , [min(p{bg,evnt}) max(p{bg,evnt})] , 'LineWidth' , 3 , 'color' , 'k')
            end
            figCount = figCount+1;
        end 
        legend(leg , titleLab)
        
        
end

