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
NumWarpSampFast = 200;
NumWarpSampSlow = 500;
DownsampleRate = 10;
FreqRange = [2 180];
numFreqBins = 90;
Channels = 129;
subjname = {'P2'};

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
blockGroups = {[1 2] , [3 13 26 40] , [4 14 27 41] , [5:7] , [9:11] , [8 12] , [15:17] , [19:21] , [23:25],...
    [18 22] , [28:30] , [32:34] , [36:38], [31 35 39],[42:44]}';
blockGroupNames = {'SingleFingNat' , 'SingleFingSlow' , 'SingleFingFast' , 'Intermixed1' , 'Intermixed2' , 'ChunkDay1' , 'Intermixed3' , 'Intermixed4' , 'Intermixed5',...
    'ChunkDay2' , 'Intermixed6' , 'Intermixed7' , 'Intermixed8', 'ChunkDay3', 'Intermixed9'}';


bandsLab = {'Delta <4Hz' , 'Theta 5-8Hz' , 'Alpha 9-16Hz' , 'Beta 17-36Hz' , 'L-Gamma 37-70Hz' , 'H-Gamma 70-130' , 'NoBandLand 130-180'};


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
            for b =1:length(bandsLab)
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
            
            set(gca , 'YTickLabels' , bandsLab, 'YTick' , T );
            
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
            %             contourf(time,frex,B,60,'linecolor','none')
            imagesc(time,frex,B)
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
                if isequal(size(F.PSD_stim{tn}) ,[length(ChanLabels) , length(bandsLab) , E.NumWarpSamp])
                    tempPow(tcount , :,:,:) = F.PSD_stim{tn};
                    tcount = tcount +1;
                end
            end
            AvgPow{sn} = real(squeeze(nanmean(tempPow , 1)));
            
            subplot(length(E.SN{1}) ,1, figCount)
            for b =1:length(bandsLab)
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
            set(gca ,'YTickLabels' , bandsLab, 'YTick' , T );
            
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
            if E.SN{1}(sn) ~= 5
                NEM = E.NEM{1}(sn ,  ~isnan(E.EM{1}(sn,:)));
                EM  = E.EM{1}(sn ,  ~isnan(E.EM{1}(sn,:)))/100;
                for ch = 1:length(Chan2Plot)
                    subplot(length(E.SN{1})-1 ,length(Chan2Plot), figCount)
                    B = real(squeeze(Pall.AvgPow{sn}(ch,:,:)));
                    contourf([1:E.NumWarpSamp],frex,B,60,'linecolor','none')
                    colorbar
                    caxis([-0.08 0.1])
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
        E.NumWarpSamp([1 3 6 10 14]) = NumWarpSampFast;
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
                stim = unique(Events.AllPress(Events.seqNumb == E.SN{1}(sn) , :) , 'rows');
                NEM = E.NEM{1}(sn , ~isnan(E.NEM{1}(sn,:)));
                % actual time in seconds to use for labeling the time axis
                EM  = E.EM{1}(sn ,  ~isnan(E.EM{1}(sn,:)))/100;
                
                id = ismember(Pall.BN , E.blockGroups{1}) & ismember(Pall.seqNumb , E.SN{1}(sn));
                F = getrow(Pall , id);
                % sum the warped PSDs insife the F structure
                tcount = 1;
                for tn = 1:length(F.PSD_stim)
                    if isequal(size(F.PSD_stim{tn}) ,[length(ChanLabels) , length(bandsLab) , E.NumWarpSamp])
                        tempPow(tcount , :,:,:) = F.PSD_stim{tn};
                        tcount = tcount +1;
                    end
                end
                AvgPow{sn} = real(squeeze(nanmean(tempPow , 1)));
                
                for ch = 1:length(Chan2Plot)
                    subplot(length(E.SN{1})-1 ,length(Chan2Plot), figCount)
                    for b =1:length(bandsLab)
                        B  = 100*squeeze(AvgPow{sn}(ch , b , :))+.7*b;
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
                    set(gca ,'YTick' , T ,'YTickLabels' , bandsLab,'XTickLabels' , xtikx, 'XTick' , NEM );
                    
                    line([max(NEM)-std(NEM) max(NEM)-std(NEM)] , [max(B)-1- max(B)],'color' , 'black' , 'LineWidth' , 5)
                    text(max(NEM),max(B),'1%','FontSize' , 16 )
                    set(gca , 'XLim' , [1 E.NumWarpSamp] ,'FontSize' , 10,'Box' , 'off')
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
        E.NumWarpSamp([1 3 6 10 14]) = NumWarpSampFast;
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
                    B = real(squeeze(Pall.AvgPow{sn}(Chan2Plot(ch),:,:)));
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
        E.NumWarpSamp([1 3 6 10 14]) = NumWarpSampFast;
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
                    if isequal(size(F.Pow_Norm_stim{tn}) ,[length(ChanLabels) , length(bandsLab) , E.NumWarpSamp])
                        tempPow(tcount , :,:,:) = F.Pow_Norm_stim{tn};
                        tcount = tcount +1;
                    end
                end
                AvgPow{sn} = real(squeeze(nanmean(tempPow , 1)));
                
                for ch = 1:length(Chan2Plot)
                    subplot(length(E.SN{1})-1 ,length(Chan2Plot), figCount)
                    for b =1:length(bandsLab)
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
                    set(gca ,'YTick' , T ,'YTickLabels' , bandsLab,'XTickLabels' , xtikx, 'XTick' , NEM );
                    
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
        E.NumWarpSamp([1 3 6 10 14]) = NumWarpSampFast;
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
                h1 = plotshade([2:2:180] , APtr{cg}(sn,:) ,AEtr{cg}(sn,:),'transp' , .2 , 'patchcolor' , 'b' , 'linecolor' , 'b' , 'linewidth' , 3 , 'linestyle' , ':')
                hold on
                h2 = plotshade([2:2:180] , APbl{cg}(sn,:) ,AEbl{cg}(sn,:),'transp' , .2 , 'patchcolor' , 'r' , 'linecolor' , 'r' , 'linewidth' , 3 , 'linestyle' , ':')
                title(['Average Powerfor sequence type ' , num2str(E.SN{1}(sn)) ])
                fCont = fCont +1;
                legend([h1,h2] , {'Trial' , 'Baseline'})
                set(gca , 'YLim' , [30 80]);
            end
        end
        
    case 'AvgPower_SeqType_comp'
        colorz = {[0 0  1],[1 0 0],[0 1 0],[1 0 1],[0 1 1],[0.7 0.7 0.7],[1 1 0],[.3 .3 .3]};
        load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped.mat
        % patterned PSDs so the output of Pall  = secog_parseEEG_PSD('TimeWarpPSD' , Dall, subjNum);
        % for this case because we are doigna comparisn, the BlockGroup is
        % a 1xN cell containnig the groups to be compared
        
        for n = 1:length(BlockGroup)
            n
            load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
            E.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
            E.NumWarpSamp([1 3 6 10 14]) = NumWarpSampFast;
            
            BG(n) = find(strcmp(E.blockGroupNames , BlockGroup{n}));
            E = getrow(E , BG(n));
            if isempty(E.EM)
                error('Nothing to plot!')
            end
            filename = [mainDir ,  'AverageSpect_SeqType' , num2str(BG(n)),'.mat'];
            load(filename)
            
            chanGroup = {Chan2Plot};
            fig = figure;
            for sn = 1:length(Pall.AvgPowTR)
                Tr = Pall.AvgPowTR{sn};
                Tre = Pall.SePowTR{sn};
                
                Bl = Pall.AvgPowBL{sn};
                Ble = Pall.SePowBL{sn};
                
                for cg = 1:length(chanGroup)
                    Tcg = Tr(chanGroup{cg} , :);
                    Bcg = Bl(chanGroup{cg} , :);
                    diff  = Tcg - Bcg;
                    diffe = Tre(chanGroup{cg} , :) - Ble(chanGroup{cg} , :);
                    bandLabel = repmat([2:2:180]' , length(chanGroup{cg}) , 1);
                    data{cg , n}(sn , :) = reshape(diff' , numel(diff) , 1);
                    datae{cg , n}(sn , :) = reshape(diffe' , numel(diffe) , 1);
                    [x{cg,n}(sn,:) , p{cg,n}(sn,:) , e{cg,n}(sn,:)] = lineplot(bandLabel , data{cg , n}(sn , :)' , 'plotfcn' , 'nanmean');
%                     diff{sn}(n, :) = 100*(APtr{cg,n}(sn,:) - APbl{cg,n}(sn,:))./ APbl{cg,n}(sn,:);
                end
            end
            close(fig)
        end
        figure('color' , 'white')
        fCont = 1;
        if length(chanGroup{1}) == 1
            for sn = 1:size(x{cg,n} , 1)
                leg = [];
                subplot(length(chanGroup) , size(x{cg,n} , 1) ,fCont)
                for n = 1:length(BlockGroup)
                    h(n).fig = plotshade(x{cg,n}(sn,:) ,data{cg , n}(sn , :) , datae{cg , n}(sn , :),'transp' , .2 , 'patchcolor' , colorz{n} , 'linecolor' , colorz{n} , 'linewidth' , 3 , 'linestyle' , ':')
                    leg = [leg , h(n).fig];
                    hold on
                    grid on
                end
                legend(leg , BlockGroup)
                title(['Power Percent change from baseline for sequence type ' , num2str(E.SN{1}(sn)) ])
                set(gca , 'YLim',[-4 6] , 'FontSize' , 16);
                fCont = fCont +1;
            end
        else
            for sn = 1:size(x{cg,n} , 1)
                leg = [];
                subplot(length(chanGroup) , size(x{cg,n} , 1) ,fCont)
                for n = 1:length(BlockGroup)
                    h(n).fig = plotshade(x{cg,n}(sn,:) , p{cg,n}(sn,:) , e{cg,n}(sn,:),'transp' , .2 , 'patchcolor' , colorz{n} , 'linecolor' , colorz{n} , 'linewidth' , 3 , 'linestyle' , ':')
                    leg = [leg , h(n).fig];
                    hold on
                    grid on
                end
                legend(leg , BlockGroup)
                title(['Power Percent change from baseline for sequence type ' , num2str(E.SN{1}(sn)) ])
                set(gca , 'YLim',[-4 6] , 'FontSize' , 16);
                fCont = fCont +1;
            end
        end
        
    case 'binned_SingleTrial_AvgChann'
        % input to this has to be AllData_PSD_StimNorm.mat
        Pall  = secog_addEventMarker(Pall, subjNum, Fs_ds , 'addEvent' , 'NumWarpSampFast' , NumWarpSampFast, 'NumWarpSampSlow'  ,NumWarpSampSlow)';
        D = getrow(Pall , Pall.BN == TBNum(2) & Pall.TN == TBNum(1));
        time = [1:size(D.Pow_Norm_stim{1} , 3)]/Fs_ds;
        figure('color' , 'white')
        Pow = real(squeeze(D.Pow_Norm_stim{1}));
        EM = find(D.EventMarker{1})/Fs_ds;
        for b =1:length(bandsLab)
            B = real(squeeze(mean(Pow(Chan2Plot , b , :),1))+8*b);
            plot(time , B , 'LineWidth' , 3)
            hold on
            T(b) = nanmean(B);
        end
        title (['Raw PSD for ',num2str(D.AllPress(1,:))])
        for lin = 1:length(EM)
            line([EM(lin) EM(lin)] , [0 max(B)] , 'color' , 'red' , 'LineStyle' , ':' , 'LineWidth' , 3)
        end
        xlabel('Norm Time')
        
        set(gca , 'YTickLabels' , bandsLab, 'YTick' , T );
        
        line([median(EM) median(EM)] , [max(B)-5 max(B)],'color' , 'black' , 'LineWidth' , 5)
        text(median(EM),max(B),'5','FontSize' , 16 )
        set(gca ,'FontSize' , 16,'Box' , 'off')
        %% PLOT average normalized binned power
        
        
end

