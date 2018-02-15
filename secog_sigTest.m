function stats  = secog_sigTest(what ,subjNum,varargin)
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
        case {'Groups'}
            % cell containing the divisions of the blockgroups
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'BlockGroup'}
            % Block Group to plot
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'Chan'}
            % channels of interest Default : everythig
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        otherwise
            error(sprintf('Unknown option: %s',varargin{c}));
    end
end
if ~exist('Chan')
    error('Define Channels --> Chan')
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


BandInfo.bandsLab = {'Delta <4Hz' , 'Theta 5-8Hz' , 'Alpha 9-16Hz' , 'Beta 17-36Hz' , 'L-Gamma 37-70Hz' , 'H-Gamma 70-130' , 'NoBandLand 130-180'};
BandInfo.bands = {[0 4], [5 8] [9 16] [17 36] [37 70] [70 130] [130 180]};

colorz = {[0 0  1],[1 0 0],[0 1 0],[1 0 1],[0 1 1],[0.7 0.7 0.7],[1 1 0],[.3 .3 .3]};
load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
switch what
    case 'seqs_across_days'
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped.mat
        % patterned PSDs so the output of Pall  = secog_parseEEG_PSD('TimeWarpPSD' , Dall, subjNum);
        % for this case because we are doigna comparisn, the BlockGroup is
        % a 1xN cell containnig the groups to be compared
        D = [];
        D.data = [];
        D.blockgroup = [];
        D.frequency = [];
        D.seqNumb = [];
        D.changroup = [];
        D.chanNumb = [];
        for bg = 1:length(BlockGroup)
            bg
            load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
            E.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
            E.NumWarpSamp([1 3 6 10 14]) = NumWarpSampFast;
            
            BG(bg) = find(strcmp(E.blockGroupNames , BlockGroup{bg}));
            E = getrow(E , BG(bg));
            if isempty(E.EM)
                error('Nothing to plot!')
            end
            filename = [mainDir ,  'AverageSpect_SeqType' , num2str(BG(bg)),'.mat'];
            load(filename)
            
            chanGroup = {Chan};
            %             fig = figure;
            for sn = 1:length(Pall.AvgPowTR)
                Tr = Pall.AvgPowTR{sn};
                Tre = Pall.SePowTR{sn};
                
                Bl = Pall.AvgPowBL{sn};
                Ble = Pall.SePowBL{sn};
                
                for cg = 1:length(chanGroup)
                    Tcg = Tr(chanGroup{cg} , :);
                    Bcg = Bl(chanGroup{cg} , :);
                    diff  = Tcg - Bcg;
                    data  = reshape(diff' , numel(diff),1);
                    D.data = [D.data ; data];
                    D.changroup  = [D.changroup; cg*ones(size(data))];
                    for ch = 1:length(chanGroup{1})
                        D.chanNumb = [D.chanNumb ;ch*ones(size(diff , 2) ,1)];
                        D.frequency  = [D.frequency; frex'];
                    end
                    D.seqNumb    = [D.seqNumb; Pall.SN{1}(sn)*ones(size(data))];
                    D.blockgroup = [D.blockgroup; bg*ones(size(data))];
                    bandLabel = repmat(frex' , length(chanGroup{cg}) , 1);
                    %                     [x{cg,sn}(bg,:) , p{cg,sn}(bg,:) , e{cg,sn}(bg,:)] = lineplot(bandLabel , data , 'plotfcn' , 'nanmean');
                end
            end
            %             close(fig)
            
            
        end
        %% prepare data for testing Block Groups against eachother
        
        
        
        SN = length(Pall.AvgPowTR);
        for cg = 1:length(chanGroup)
            for sn = 1 : SN
                for bg = 1:length(BlockGroup)
                    groupNum     = Groups(bg);
                    D.data       = [D.data ;p{cg,sn}(bg,:)'];
                    
                    D.changroup  = [D.changroup; cg*ones(size(p{cg,sn}(bg,:)'))];
                end
            end
        end
        
        % abbplyBanding
        T = []
        for b = 1:length(BandInfo.bands)
            A = getrow(D , ismember(D.frequency , [BandInfo.bands{b}(1) : BandInfo.bands{b}(2)]));
            temp =  tapply(A , {'changroup','blockgroup','seqNumb'},{'data'} )
            temp.band = b*ones(size(temp.seqNumb));
            T = addstruct(T , temp);
        end
        SN = length(Pall.AvgPowTR);
        figure('color' , 'white')
        subplot(311)
        lineplot([D.seqNumb D.blockgroup] , D.data , 'subset' , ismember(D.frequency , [BandInfo.bands{6}(1) : BandInfo.bands{6}(2)]) , 'style_thickline')
        set(gca , 'XTickLabel' , repmat({'Day1' , 'Day2' , 'Day3' , 'Day4'} , 1 , SN))
        title('High Gamma Band Power Across days for Random, Structured and Null Trails respectively left to right')
        grid on
        
        subplot(312)
        lineplot([D.seqNumb D.blockgroup ] , D.data , 'subset' , ismember(D.frequency , [BandInfo.bands{5}(1) : BandInfo.bands{5}(2)]) , 'style_thickline')
        set(gca , 'XTickLabel' , repmat({'Day1' , 'Day2' , 'Day3' , 'Day4'} , 1 , SN))
        title('Low Gamma Band Power Across days for Random, Structured and Null Trails respectively left to right')
        grid on
        
        subplot(313)
        lineplot([D.seqNumb D.blockgroup ] , D.data , 'subset' , ismember(D.frequency , [BandInfo.bands{4}(1) : BandInfo.bands{4}(2)]) , 'style_thickline')
        set(gca , 'XTickLabel' , repmat({'Day1' , 'Day2' , 'Day3' , 'Day4'} , 1 , SN))
        title('Beta Band Power Across days for Random, Structured and Null Trails respectively left to right')
        grid on
        
        %% significance tests
        %% testing the effect of sequence type and days on high gamma
        % P = anovaMixed(D.data , [1:length(D.data)]','between',[D.blockgroup D.seqNumb],...
        %     {'BlockGroup' , 'Sequence Type'},'intercept',1 , 'subset',ismember(D.frequency , [BandInfo.bands{6}(1) : BandInfo.bands{6}(2)]));
        A = getrow(D , ismember(D.frequency , [BandInfo.bands{6}(1) : BandInfo.bands{6}(2)]));
        [p,tbl,stats] =  anovan(A.data , [A.blockgroup A.seqNumb],'model','interaction','varnames',{'Days' , 'Sequence Type'});
        P = anovaMixed(A.data , [1:length(A.data)]','between',[A.blockgroup A.seqNumb],{'BlockGroup' , 'Sequence Type'},'intercept',1);
        %% testing the effect of sequence type and days on low gamma
        A = getrow(D , ismember(D.frequency , [BandInfo.bands{5}(1) : BandInfo.bands{5}(2)]));
        [p,tbl,stats] =  anovan(A.data , [A.blockgroup A.seqNumb],'model','interaction','varnames',{'Days' , 'Sequence Type'});
        P = anovaMixed(A.data , [1:length(A.data)]','between',[A.blockgroup A.seqNumb],{'BlockGroup' , 'Sequence Type'},'intercept',1)
        %% testing the effect of sequence type and days on beta
        A = getrow(D , ismember(D.frequency , [BandInfo.bands{4}(1) : BandInfo.bands{4}(2)]));
        [p,tbl,stats] =  anovan(A.data , [A.blockgroup A.seqNumb],'model','interaction','varnames',{'Days' , 'Sequence Type'});
        P = anovaMixed(D.data , [1:length(D.data)]','between',[D.blockgroup D.seqNumb],{'BlockGroup' , 'Sequence Type'},'intercept',1)
        %% testing the effect of days on random sequences with every day as a single degree of freedom in high gamma
        A = getrow(D , ismember(D.frequency , [BandInfo.bands{6}(1) : BandInfo.bands{6}(2)]) & ismember(D.seqNumb , 20));
        [p,tbl,stats] =  anovan(A.data , [A.blockgroup],'model','interaction','varnames',{'Days'});
        P = anovaMixed(A.data , [1:length(A.data)]','between',[A.blockgroup],{'BlockGroup'},'intercept',1);
        %% testing the effect first vs. all the other days in high gamma in Random sequences
        % P = anovaMixed(D.data , [1:length(D.data)]','between',[D.blockgroup],{'BlockGroup'},'intercept',1 ,...
        %     ismember(D.frequency , [BandInfo.bands{6}(1) : BandInfo.bands{6}(2)]) & ismember(D.seqNumb , 20));
        A = getrow(D , ismember(D.frequency , [BandInfo.bands{6}(1) : BandInfo.bands{6}(2)]) & ismember(D.seqNumb , 20));
        A.blockgroup(A.blockgroup>1) = 0;
        [p,tbl,stats] =  anovan(A.data , [A.blockgroup],'model','interaction','varnames',{'Days'});
        P = anovaMixed(A.data , [1:length(A.data)]','between',[A.blockgroup],{'BlockGroup'},'intercept',1);
        %% testing the effect of days on Structured sequences with every day as a single degree of freedom in high gamma
        % P = anovaMixed(D.data , [1:length(D.data)]','between',[D.blockgroup],{'BlockGroup'},'intercept',1 ,...
        %     ismember(D.frequency , [BandInfo.bands{6}(1) : BandInfo.bands{6}(2)]) & ismember(D.seqNumb , 20));
        A = getrow(D , ismember(D.frequency , [BandInfo.bands{6}(1) : BandInfo.bands{6}(2)]) & ismember(D.seqNumb , 30));
        [p,tbl,stats] =  anovan(A.data , [A.blockgroup],'model','interaction','varnames',{'Days'});
        P = anovaMixed(A.data , [1:length(A.data)]','between',[A.blockgroup],{'BlockGroup'},'intercept',1);
        %% testing the effect first vs. all the other days in high gamma in structured sequences
        % P = anovaMixed(D.data , [1:length(D.data)]','between',[D.blockgroup],{'BlockGroup'},'intercept',1 ,...
        %     ismember(D.frequency , [BandInfo.bands{6}(1) : BandInfo.bands{6}(2)]) & ismember(D.seqNumb , 20));
        A = getrow(D , ismember(D.frequency , [BandInfo.bands{6}(1) : BandInfo.bands{6}(2)]) & ismember(D.seqNumb , 30));
        A.blockgroup(A.blockgroup>1) = 0;
        [p,tbl,stats] =  anovan(A.data , [A.blockgroup],'model','interaction','varnames',{'Days'});
        P = anovaMixed(A.data , [1:length(A.data)]','between',[A.blockgroup],{'BlockGroup'},'intercept',1);
        
        
        
        %%
        % testing the effect of days on structured sequences with every day as a single degree of freedom in high gamma
        P = anovaMixed(D.data , [1:length(D.data)]','between',[D.blockgroup],{'BlockGroup'},'intercept',1 ,...
            ismember(D.frequency , [BandInfo.bands{6}(1) : BandInfo.bands{6}(2)]) & ismember(D.seqNumb , 100));
        
        close all
        
end


