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


subjname = {'P2' , 'P4'};
mainDir = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' subjname{subjNum} , '/'] ;
DownsampleRate = 10;
NumWarpSampFast = 200;
NumWarpSampSlow = 500;
FreqRange = [2 180];
numFreqBins = 90;
load([mainDir , 'ChanLabels.mat']);
Channels = [1:length(ChanLabels)];



c = 1;
while(c<=length(varargin))
    switch(varargin{c})
        case {'Groups'}
            % cell containing the divisions of the blockgroups
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
load([mainDir , 'ChanLabels.mat'])
load([mainDir , 'AllData_Behav.mat'])
load([mainDir , 'AllData_Events.mat'])
load([mainDir , 'AllData_AvgMarker.mat'])

%% HouseKeeping

blockGroups = {[1 2] , [3], [13], [26], [40] , [4], [14], [27] [41] , [5:7] , [9:11] , [8 12] , [15:17] , [19:21] , [23:25],...
    [18 22] , [28:30] , [32:34] , [36:38], [31 35 39],[42:44]}';
blockGroupNames = {'SingleFingNat' , 'SingleFingSlow1' , 'SingleFingSlow2'  , 'SingleFingSlow3' ,'SingleFingSlow4',...
    'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3', 'SingleFingFast4' , 'Intermixed1' , 'Intermixed2' , ...
    'ChunkDay1' , 'Intermixed3' , 'Intermixed4' , 'Intermixed5', 'ChunkDay2' , 'Intermixed6' , ...
    'Intermixed7' , 'Intermixed8', 'ChunkDay3', 'Intermixed9'}';

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
        allSubj_BlockGroup(1).bg = {'Intermixed1','Intermixed5','Intermixed8','Intermixed9'};
        allSubj_BlockGroup(2).bg = {'Intermixed1','Intermixed4','Intermixed7'};
        BlockGroup = allSubj_BlockGroup(subjNum).bg;
        D = [];
        D.data = [];
        D.blockgroup = [];
        D.frequency = [];
        D.seqNumb = [];
        D.changroup = [];
        D.chanNumb = [];
        D. TN = [];
        
        for bg = 1:length(BlockGroup)
            load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
            E.NumWarpSamp = NumWarpSampSlow*ones(size(blockGroups));
            E.NumWarpSamp([1 3 6 10 14]) = NumWarpSampFast;
            
            BG(bg) = find(strcmp(E.blockGroupNames , BlockGroup{bg}));
            E = getrow(E , BG(bg));
            if isempty(E.EM)
                error('Nothing to plot!')
            end
            filename = [mainDir ,  'PowerSpectrum_SeqType' , num2str(BG(bg)),'.mat'];
            load(filename)
            
            chanGroup = {Chan};
            %             fig = figure;
            for sn = 1:length(Pall.PowTR)
                Tr = Pall.PowTR{sn};
                
                Bl = Pall.PowBL{sn};
                
                for cg = 1:length(chanGroup)
                    Tcg = Tr(:,chanGroup{cg} , :);
                    Bcg = Bl(:,chanGroup{cg} , :);
                    diff  = Tcg - Bcg;
                    
                    for tr = 1:size(diff , 1)
                        temp = squeeze(diff(tr , :,:));
                        data  = reshape(temp' , numel(temp),1);
                        D.data = [D.data ; data];
                        D.changroup  = [D.changroup; cg*ones(size(data))];
                        for ch = 1:length(chanGroup{1})
                            D.chanNumb = [D.chanNumb ;ch*ones(length(frex) ,1)];
                            D.frequency  = [D.frequency; frex'];
                        end
                        D.seqNumb    = [D.seqNumb; Pall.SN{1}(sn)*ones(size(data))];
                        D.blockgroup = [D.blockgroup; bg*ones(size(data))];
                        D.TN         = [D.TN ; tr*ones(length(data),1)];
                    end
                end
            end
            %             close(fig)
        end
        D.band = zeros(size(D.frequency));
        for b = 1:length(BandInfo.bands)
            ind = ismember(D.frequency , [BandInfo.bands{b}(1) : BandInfo.bands{b}(2)]);
            D.band(ind) = b;
        end
        %% prepare data for testing Block Groups against eachother
        
        

        % applyBanding
        T = tapply(D , {'changroup','blockgroup','seqNumb','TN','band' },{'data'} );
           
        SN = length(Pall.PowTR);
        figure('color' , 'white')
        subplot(311)
        lineplot([T.seqNumb T.blockgroup] , T.data , 'subset' , ismember(T.band, 6) , 'style_thickline')
        set(gca , 'XTickLabel' , repmat({'Day1' , 'Day2' , 'Day3' , 'Day4'} , 1 , SN))
        title('High Gamma Band Power Across days for Random, Structured and Null Trails respectively left to right')
        grid on
        
        subplot(312)
        lineplot([T.seqNumb T.blockgroup] , T.data , 'subset' , ismember(T.band, 5) , 'style_thickline')
        set(gca , 'XTickLabel' , repmat({'Day1' , 'Day2' , 'Day3' , 'Day4'} , 1 , SN))
        title('Low Gamma Band Power Across days for Random, Structured and Null Trails respectively left to right')
        grid on
        
        subplot(313)
        lineplot([T.seqNumb T.blockgroup] , T.data , 'subset' , ismember(T.band, 4) , 'style_thickline')
        set(gca , 'XTickLabel' , repmat({'Day1' , 'Day2' , 'Day3' , 'Day4'} , 1 , SN))
        title('Beta Band Power Across days for Random, Structured and Null Trails respectively left to right')
        grid on
        
       
        disp('***********************************************************************************')
        disp('****************************** CHOSE THE TEST TO RUN ******************************')
        disp('***********************************************************************************')
        keyboard
        
        %% significance tests
        %% testing the effect of sequence type and days on high gamma
        bandnumb = 7;
        seqNumb = [30 20];
        A = getrow(D , ismember(D.band, bandnumb)  & ismember(D.seqNumb , seqNumb ));
        %         D.blockGroup(D.blockGroup>1) = 0;     % test day one vs all other days
        if length(seqNumb)>1
            [p,tbl,stats] =  anovan(A.data , [A.blockgroup A.seqNumb],'model','interaction','varnames',{'Days' , 'Sequence Type'});
        else
            [p,tbl,stats] =  anovan(A.data , [A.blockgroup],'model','interaction','varnames',{'Days'});
        end
        close all
    case 'chunks_across_days'
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped.mat
        % patterned PSDs so the output of Pall  = secog_parseEEG_PSD('TimeWarpPSD' , Dall, subjNum);
        % for this case because we are doigna comparisn, the BlockGroup is
        % a 1xN cell containnig the groups to be compared
        BlockGroup = {'ChunkDay1','ChunkDay2','ChunkDay3'};
       
        C = [];
        T = [];
        
        if sum(ismember (BlockGroup , blockGroupNames([12 16 20])))
            seq = input('Triplet or Quadruple? (t/q)', 's');
            switch seq
                case {'q'}
                    selectSN = 50;
                    chnkLabel = 'Quadruples';
                case {'t'}
                    selectSN = 40;
                    chnkLabel = 'Triplets';
            end
        end
        for bg = 1:length(BlockGroup)
            load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
           
            BG(bg) = find(strcmp(blockGroupNames , BlockGroup{bg}));
            E = getrow(E , BG(bg));
            if isempty(E.EM)
                error('Nothing to plot!')
            end
            filename = [mainDir ,  'Group_Aligned_PSD' , num2str(BG(bg)),'.mat'];
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
                        for bands = 1:size(T.PSD{tn , evnt} , 2)
                            temp = squeeze(T.PSD{tn , evnt}(:,bands,:));
                            eval(['T.reducPSD' , num2str(bands) , '(tn , evnt) = nanmean(nanmean(temp(Chan , :)));'])
                        end
                    end
                end
                T.blockgroup = bg*ones(size(T.BN));
                C = addstruct(C , T);
            end
        end

        %% prepare data for testing Block Groups against eachother
        D = [];
        for tn = 1:length(C.TN)
            temp = getrow(C , tn);
            T = [];
            T.avgPSD = [];
            T.seqNumb = [];
            T.band = [];
            T.blockGroup = [];
            T.TN = [];
            T.eventNum = [];
            for bands = 1:size(temp.PSD{1} , 2)
               eval(['A = temp.reducPSD' , num2str(bands) , ';']);
               T.avgPSD = [T.avgPSD ; A'];
               T.eventNum = [T.eventNum ; [1:length(A)]'];
               T.seqNumb = [T.seqNumb ; temp.seqNumb*ones(size(A'))];
               T.band = [T.band; bands*ones(size(A'))];
               T.blockGroup = [T.blockGroup ; temp.blockgroup*ones(size(A'))];
               T.TN = [T.TN ; temp.TN*ones(size(A'))];
            end
            D = addstruct(D , T);
        end
        %% prepare data for testing Block Groups against eachother
        
        

        % applyBanding
        T = D;
        
        figure('color' , 'white')
        subplot(311)
        T = getrow(D , D.band == 6);  % high gamma
        lineplot([T.eventNum T.blockGroup] , T.avgPSD , 'style_thickline' , 'linecolor' , 'b' , 'errorcolor' , 'b' , 'leg' , 'Random')
        set(gca , 'XTickLabel' , repmat({'Day1' , 'Day2' , 'Day3' , 'Day4'} , 1 , max(T.eventNum)))
        title(['High Gamma Band Power Across days for ' , chnkLabel])
        grid on
        
        subplot(312)
        T = getrow(D , D.band == 5);  % low gamma
        lineplot([T.eventNum T.blockGroup] , T.avgPSD , 'style_thickline' , 'linecolor' , 'b' , 'errorcolor' , 'b' , 'leg' , 'Random')
        set(gca , 'XTickLabel' , repmat({'Day1' , 'Day2' , 'Day3' , 'Day4'} , 1 , max(T.eventNum)))
        title(['Low Gamma Band Power Across days for ' , chnkLabel])
        grid on
        
        subplot(313)
        T = getrow(D , D.band == 4);  % betta
        lineplot([T.eventNum T.blockGroup] , T.avgPSD , 'style_thickline' , 'linecolor' , 'b' , 'errorcolor' , 'b' , 'leg' , 'Random')
        set(gca , 'XTickLabel' , repmat({'Day1' , 'Day2' , 'Day3' , 'Day4'} , 1 , max(T.eventNum)))
        title(['Beta Band Power Across days for ' , chnkLabel])
        grid on
        disp('***********************************************************************************')
        disp('****************************** CHOSE THE TEST TO RUN ******************************')
        disp('***********************************************************************************')
        keyboard
        %% significance tests
        %% testing the effect of sequence type and days on high gamma
        % P = anovaMixed(D.data , [1:length(D.data)]','between',[D.blockgroup D.seqNumb],...
        %     {'BlockGroup' , 'Sequence Type'},'intercept',1 , 'subset',ismember(D.frequency , [BandInfo.bands{6}(1) : BandInfo.bands{6}(2)]));
        bandnumb = 6; 
        EventsofInterest =1;
        seqNumb = selectSN;
        A = getrow(D , ismember(D.band, bandnumb) & ismember(D.eventNum , EventsofInterest )  & ismember(D.seqNumb , seqNumb ));
%         D.blockGroup(D.blockGroup>1) = 0;     % test day one vs all other days
        if length(EventsofInterest)>1
            if length(seqNumb)>1
                [p,tbl,stats] =  anovan(A.avgPSD , [A.blockGroup A.seqNumb A.eventNum],'model','interaction','varnames',{'Days' , 'Sequence Type' , 'Event'});
            else
                [p,tbl,stats] =  anovan(A.avgPSD , [A.blockGroup A.eventNum],'model','interaction','varnames',{'Days' , 'Event'});
            end
        else
            if length(seqNumb)>1
                [p,tbl,stats] =  anovan(A.avgPSD , [A.blockGroup A.seqNumb],'model','interaction','varnames',{'Days' , 'Sequence Type'});
            else
                [p,tbl,stats] =  anovan(A.avgPSD , [A.blockGroup],'model','interaction','varnames',{'Days'});
            end
        end
        P = anovaMixed(A.data , [1:length(A.data)]','between',[A.blockgroup A.seqNumb],{'BlockGroup' , 'Sequence Type'},'intercept',1);        
    case 'Alined_Sequences'


        allSubj_BlockGroup(1).bg = {'Intermixed1','Intermixed5','Intermixed8','Intermixed9'};
        allSubj_BlockGroup(2).bg = {'Intermixed1','Intermixed4','Intermixed7'};
        BlockGroup = allSubj_BlockGroup(subjNum).bg;


        C = [];
        T = [];

        for bg = 1:length(BlockGroup)
            load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
           
            BG(bg) = find(strcmp(blockGroupNames , BlockGroup{bg}));
            E = getrow(E , BG(bg));
            if isempty(E.EM)
                error('Nothing to plot!')
            end
            filename = [mainDir ,  'Group_Aligned_PSD' , num2str(BG(bg)),'.mat'];
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
                        for bands = 1:size(T.PSD{tn , evnt} , 2)
                            temp = squeeze(T.PSD{tn , evnt}(:,bands,:));
                            eval(['T.reducPSD' , num2str(bands) , '(tn , evnt) = nanmean(nanmean(temp(Chan , :)));'])
                        end
                    end
                end
                T.blockgroup = bg*ones(size(T.BN));
                C = addstruct(C , T);
            end
        end

        %% prepare data for testing Block Groups against eachother
        D = [];
        for tn = 1:length(C.TN)
            temp = getrow(C , tn);
            T = [];
            T.avgPSD = [];
            T.seqNumb = [];
            T.band = [];
            T.blockGroup = [];
            T.TN = [];
            T.eventNum = [];
            for bands = 1:size(temp.PSD{1} , 2)
               eval(['A = temp.reducPSD' , num2str(bands) , ';']);
               T.avgPSD = [T.avgPSD ; A'];
               T.eventNum = [T.eventNum ; [1:length(A)]'];
               T.seqNumb = [T.seqNumb ; temp.seqNumb*ones(size(A'))];
               T.band = [T.band; bands*ones(size(A'))];
               T.blockGroup = [T.blockGroup ; temp.blockgroup*ones(size(A'))];
               T.TN = [T.TN ; temp.TN*ones(size(A'))];
            end
            D = addstruct(D , T);
        end
        %% prepare data for testing Block Groups against eachother
        % applyBanding
        figure('color' , 'white')
        subplot(311)
        T = getrow(D , D.band == 6 & D.seqNumb == 20);  % high gamma
        lineplot([T.eventNum T.blockGroup] , T.avgPSD , 'style_thickline' , 'linecolor' , 'b' , 'errorcolor' , 'b' , 'leg' , 'Random')
        hold on 
        T = getrow(D , D.band == 6 & D.seqNumb == 30);  % high gamma
        lineplot([T.eventNum T.blockGroup] , T.avgPSD , 'style_thickline' , 'linecolor' , 'r' , 'errorcolor' , 'r' , 'leg' , 'Trained')
        set(gca , 'XTickLabel' , repmat({'Day1' , 'Day2' , 'Day3' , 'Day4'} , 1 , max(T.eventNum)))
        title('High Gamma Band Power Across days for Random (blue) and Structured (red)')
        grid on
        
        subplot(312)
        T = getrow(D , D.band == 5 & D.seqNumb == 20);  % low gamma
        lineplot([T.eventNum T.blockGroup] , T.avgPSD , 'style_thickline' , 'linecolor' , 'b' , 'errorcolor' , 'b' , 'leg' , 'Random')
        hold on 
        T = getrow(D , D.band == 5 & D.seqNumb == 30);  % low gamma
        lineplot([T.eventNum T.blockGroup] , T.avgPSD , 'style_thickline' , 'linecolor' , 'r' , 'errorcolor' , 'r' , 'leg' , 'Trained')
        set(gca , 'XTickLabel' , repmat({'Day1' , 'Day2' , 'Day3' , 'Day4'} , 1 , max(T.eventNum)))
        title('Low Gamma Band Power Across days for Random (blue) and Structured (red)')
        grid on
        
        subplot(313)
        T = getrow(D , D.band == 4 & D.seqNumb == 20);  % betta
        lineplot([T.eventNum T.blockGroup] , T.avgPSD , 'style_thickline' , 'linecolor' , 'b' , 'errorcolor' , 'b' , 'leg' , 'Random')
        hold on 
        T = getrow(D , D.band == 4 & D.seqNumb == 30);  % betta
        lineplot([T.eventNum T.blockGroup] , T.avgPSD , 'style_thickline' , 'linecolor' , 'r' , 'errorcolor' , 'r' , 'leg' , 'Trained')
        set(gca , 'XTickLabel' , repmat({'Day1' , 'Day2' , 'Day3' , 'Day4'} , 1 , max(T.eventNum)))
        title('Beta Band Power Across days for Random (blue) and Structured (red)')
        grid on
        disp('***********************************************************************************')
        disp('****************************** CHOSE THE TEST TO RUN ******************************')
        disp('***********************************************************************************')
        keyboard
        %% significance tests
        %% testing the effect of sequence type and days on high gamma (6) lowgamma (5) beta (4)
        % P = anovaMixed(D.data , [1:length(D.data)]','between',[D.blockgroup D.seqNumb],...
        %     {'BlockGroup' , 'Sequence Type'},'intercept',1 , 'subset',ismember(D.frequency , [BandInfo.bands{6}(1) : BandInfo.bands{6}(2)]));
        bandnumb = 4; 
        EventsofInterest =1;
        seqNumb = [20];
        A = getrow(D , ismember(D.band, bandnumb) & ismember(D.eventNum , EventsofInterest )  & ismember(D.seqNumb , seqNumb ));
%         D.blockGroup(D.blockGroup>1) = 0;     % test day one vs all other days
        if length(EventsofInterest)>1
            if length(seqNumb)>1
                [p,tbl,stats] =  anovan(A.avgPSD , [A.blockGroup A.seqNumb A.eventNum],'model','interaction','varnames',{'Days' , 'Sequence Type' , 'Event'});
            else
                [p,tbl,stats] =  anovan(A.avgPSD , [A.blockGroup A.eventNum],'model','interaction','varnames',{'Days' , 'Event'});
            end
        else
            if length(seqNumb)>1
                [p,tbl,stats] =  anovan(A.avgPSD , [A.blockGroup A.seqNumb],'model','interaction','varnames',{'Days' , 'Sequence Type'});
            else
                [p,tbl,stats] =  anovan(A.avgPSD , [A.blockGroup],'model','interaction','varnames',{'Days'});
            end
        end
        P = anovaMixed(A.data , [1:length(A.data)]','between',[A.blockgroup A.seqNumb],{'BlockGroup' , 'Sequence Type'},'intercept',1);        
    case 'Alined_SingleFinger'
        
        blockGroups = {[1 2] , [3], [13], [26], [40] , [4], [14], [27] [41] , [5:7] , [9:11] , [8 12] , [15:17] , [19:21] , [23:25],...
            [18 22] , [28:30] , [32:34] , [36:38], [31 35 39],[42:44]}';
        blockGroupNames = {'SingleFingNat' , 'SingleFingSlow1' , 'SingleFingSlow2'  , 'SingleFingSlow3' ,'SingleFingSlow4',...
            'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3', 'SingleFingFast4' , 'Intermixed1' , 'Intermixed2' , ...
            'ChunkDay1' , 'Intermixed3' , 'Intermixed4' , 'Intermixed5', 'ChunkDay2' , 'Intermixed6' , ...
            'Intermixed7' , 'Intermixed8', 'ChunkDay3', 'Intermixed9'}';
        
        seq = input('Fast or Slow? (f/s)', 's');
            switch seq
                case {'f'}
                    BlockGroup = {'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3', 'SingleFingFast4'};
                case {'s'}
                    BlockGroup = {'SingleFingSlow1' , 'SingleFingSlow2'  , 'SingleFingSlow3' ,'SingleFingSlow4'}; 
            end
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
           
            BG(bg) = find(strcmp(blockGroupNames , BlockGroup{bg}));
            E = getrow(E , BG(bg));
            if isempty(E.EM)
                error('Nothing to plot!')
            end
            filename = [mainDir ,  'Group_Aligned_PSD' , num2str(BG(bg)),'.mat'];
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
                        for bands = 1:size(T.PSD{tn , evnt} , 2)
                            temp = squeeze(T.PSD{tn , evnt}(:,bands,:));
                            eval(['T.reducPSD' , num2str(bands) , '(tn , evnt) = nanmean(nanmean(temp(Chan , :)));'])
                        end
                    end
                end
                T.blockgroup = bg*ones(size(T.BN));
                C = addstruct(C , T);
            end
        end

        %% prepare data for testing Block Groups against eachother
        D = [];
        for tn = 1:length(C.TN)
            temp = getrow(C , tn);
            T = [];
            T.avgPSD = [];
            T.seqNumb = [];
            T.band = [];
            T.blockGroup = [];
            T.TN = [];
            T.eventNum = [];
            for bands = 1:size(temp.PSD{1} , 2)
               eval(['A = temp.reducPSD' , num2str(bands) , ';']);
               T.avgPSD = [T.avgPSD ; A'];
               T.eventNum = [T.eventNum ; [1:length(A)]'];
               T.seqNumb = [T.seqNumb ; temp.seqNumb*ones(size(A'))];
               T.band = [T.band; bands*ones(size(A'))];
               T.blockGroup = [T.blockGroup ; temp.blockgroup*ones(size(A'))];
               T.TN = [T.TN ; temp.TN*ones(size(A'))];
            end
            D = addstruct(D , T);
        end
        %% prepare data for testing Block Groups against eachother
        % applyBanding
        figure('color' , 'white')
        subplot(311)
        T = getrow(D , D.band == 6);  % high gamma
        lineplot([T.eventNum T.blockGroup] , T.avgPSD , 'style_thickline' , 'linecolor' , 'b' , 'errorcolor' , 'b' , 'leg' , 'Random')
        set(gca , 'XTickLabel' , repmat({'Day1' , 'Day2' , 'Day3' , 'Day4'} , 1 , max(T.eventNum)))
        title('High Gamma Band Power Across days for Single finger repetitions')
        grid on
        
        subplot(312)
        T = getrow(D , D.band == 5);  % low gamma
        lineplot([T.eventNum T.blockGroup] , T.avgPSD , 'style_thickline' , 'linecolor' , 'b' , 'errorcolor' , 'b' , 'leg' , 'Random')
        set(gca , 'XTickLabel' , repmat({'Day1' , 'Day2' , 'Day3' , 'Day4'} , 1 , max(T.eventNum)))
        title('Low Gamma Band Power Across days for Single finger repetitions')
        grid on
        
        subplot(313)
        T = getrow(D , D.band == 4 );  % betta
        lineplot([T.eventNum T.blockGroup] , T.avgPSD , 'style_thickline' , 'linecolor' , 'b' , 'errorcolor' , 'b' , 'leg' , 'Random')
        set(gca , 'XTickLabel' , repmat({'Day1' , 'Day2' , 'Day3' , 'Day4'} , 1 , max(T.eventNum)))
        title('Beta Band Power Across days for Single finger repetitions')
        grid on
        disp('***********************************************************************************')
        disp('****************************** CHOSE THE TEST TO RUN ******************************')
        disp('***********************************************************************************')
        keyboard
        %% significance tests
        %% testing the effect of sequence type and days on high gamma (6) lowgamma (5) beta (4)
        % P = anovaMixed(D.data , [1:length(D.data)]','between',[D.blockgroup D.seqNumb],...
        %     {'BlockGroup' , 'Sequence Type'},'intercept',1 , 'subset',ismember(D.frequency , [BandInfo.bands{6}(1) : BandInfo.bands{6}(2)]));
        bandnumb = 6; 
        EventsofInterest =6;
        seqNumb = [10];
        A = getrow(D , ismember(D.band, bandnumb) & ismember(D.eventNum , EventsofInterest )  & ismember(D.seqNumb , seqNumb ));
%         D.blockGroup(D.blockGroup>1) = 0;     % test day one vs all other days
        if length(EventsofInterest)>1
            if length(seqNumb)>1
                [p,tbl,stats] =  anovan(A.avgPSD , [A.blockGroup A.seqNumb A.eventNum],'model','interaction','varnames',{'Days' , 'Sequence Type' , 'Event'});
            else
                [p,tbl,stats] =  anovan(A.avgPSD , [A.blockGroup A.eventNum],'model','interaction','varnames',{'Days' , 'Event'});
            end
        else
            if length(seqNumb)>1
                [p,tbl,stats] =  anovan(A.avgPSD , [A.blockGroup A.seqNumb],'model','interaction','varnames',{'Days' , 'Sequence Type'});
            else
                [p,tbl,stats] =  anovan(A.avgPSD , [A.blockGroup],'model','interaction','varnames',{'Days'});
            end
        end

end


