function stats  = seqeeg_sigTest(what ,subjNum,varargin)
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


subjname = {'P2' , 'P4' , 'P5'};
mainDir = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' subjname{subjNum} , '/'] ;
DownsampleRate = 10;
FreqRange = [2 180];
numFreqBins = 90;
load([mainDir , 'ChanLabels.mat']);
Channels = [1:length(ChanLabels)];

daylab(1).dl = {'Day1' , 'Day2' , 'Day3' , 'Day4'};
daylab(2).dl = {'Day1' , 'Day2' , 'Day3' };
daylab(3).dl = {'Day1' , 'Day2' , 'Day3' , 'Day4'};


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
% block groupings for subject 1

blockGroupNames = {'SingleFingNat' , 'SingleFingSlow1' , 'SingleFingSlow2'  , 'SingleFingSlow3' ,'SingleFingSlow4',...
    'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3', 'SingleFingFast4' , 'Intermixed1' , 'Intermixed2' , ...
    'ChunkDay1' , 'Intermixed3' , 'Intermixed4' , 'Intermixed5', 'ChunkDay2' , 'Intermixed6' , ...
    'Intermixed7' , 'Intermixed8', 'ChunkDay3', 'Intermixed9','Intermixed10','ChunkDay4'}';
BandInfo.bandsLab = {'Delta <4Hz' , 'Theta 5-8Hz' , 'Alpha 9-16Hz' , 'Beta 17-36Hz' , 'L-Gamma 37-80Hz' , 'H-Gamma 80-100HZ' , 'HIGH 100-180HZ'};
BandInfo.bands = {[0 4], [5 8] [9 16] [17 36] [37 80] [80 110] [100 180]};
sequenceType.TypeNums = [100 10 20 30 40 50];
sequenceType.Typetags = {'Null Trials' , 'Single Finger Sequences Trials' , 'Random Sequences' , 'Trained Sequences' , 'Triplet Segments' , 'Quadruple Segments'}';

sequenceType.SeqNums  = [5 11 22 33 44 55 0 1 2 3 4 103 203 104 204];
sequenceType.Numtags = {'Null Trials' , 'Single Finger 1' , 'Single Finger 2' , 'Single Finger 3' , 'Single Finger 4' , 'Single Finger 5',...
    'Random Sequences' , 'Structure 1231234' , 'Structure 1231234' , 'Structure 1234123' , 'Structure 1234123' ,...
    'Triplet Segments','Triplet Segments' , 'Quadruple Segments','Quadruple Segments'}';


colorz = {[0 0  1],[1 0 0],[0 1 0],[1 0 1],[0 1 1],[0.7 0.7 0.7],[1 1 0],[.3 .3 .3]};
load([mainDir , 'AllData_AvgMarker_SeqType.mat'])

switch what
    case 'seqs_across_days'
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped.mat
        % patterned PSDs so the output of Pall  = seqeeg_parseEEG_PSD('TimeWarpPSD' , Dall, subjNum);
        % for this case because we are doigna comparisn, the BlockGroup is
        % a 1xN cell containnig the groups to be compared
        allSubj_BlockGroup(1).bg = {'Intermixed1','Intermixed5','Intermixed8','Intermixed9'};
        allSubj_BlockGroup(2).bg = {'Intermixed1','Intermixed4','Intermixed7'};
        allSubj_BlockGroup(1).bg = {'Intermixed1','Intermixed4','Intermixed7','Intermixed10'};
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
        T = D;
        %% prepare data for testing Block Groups against eachother
        
        

        % applyBanding
%         T = tapply(D , {'changroup','blockgroup','seqNumb','TN','band' },{'data'} );
        
        D = getrow(T , ~ismember(T.seqNumb,100));
        
        SN = length(Pall.PowTR);
        figure('color' , 'white')
        
        subplot(311)
        D1 = getrow(D, D.seqNumb == 20 & ismember(D.band, 6));
        lineplot([D1.seqNumb D1.blockgroup] , D1.data , 'style_thickline' , 'linecolor' , colorz{1} , 'errorcolor' , colorz{1})
        hold on
        D1 = getrow(D, D.seqNumb == 30 & ismember(D.band, 6));
        lineplot([D1.seqNumb D1.blockgroup] , D1.data , 'style_thickline', 'linecolor' , colorz{2} , 'errorcolor' , colorz{2})
        set(gca , 'XTickLabel' , repmat(daylab(subjNum).dl , 1 , SN))
        title('High Gamma Band Power Across days for for Random(blue) and Structured(red)')
        grid on
        
        subplot(312)
        D1 = getrow(D, D.seqNumb == 20 & ismember(D.band, 5));
        lineplot([D1.seqNumb D1.blockgroup] , D1.data , 'style_thickline','subset' , D1.seqNumb == 20)% , 'linecolor' , colorz{1} , 'errorcolor' , colorz{1})
        hold on
        D1 = getrow(D, D.seqNumb == 30 & ismember(D.band, 5));
        lineplot([D1.seqNumb D1.blockgroup] , D1.data , 'style_thickline','subset' , D1.seqNumb == 30 , 'linecolor' , colorz{2} , 'errorcolor' , colorz{2})
        set(gca , 'XTickLabel' , repmat(daylab(subjNum).dl , 1 , SN))
        title('Low Gamma Band Power Across days for Random(blue) and Structured(red)')
        grid on
        
        subplot(313)
        D1 = getrow(D, D.seqNumb == 20 & ismember(D.band, 4));
        lineplot([D1.seqNumb D1.blockgroup] , D1.data , 'style_thickline','subset' ,D1.seqNumb == 20 , 'linecolor' , colorz{1} , 'errorcolor' , colorz{1})
        hold on
        D1 = getrow(D, D.seqNumb == 30 & ismember(D.band, 4));
        lineplot([D1.seqNumb D1.blockgroup] , D1.data , 'style_thickline','subset' ,D1.seqNumb == 30 , 'linecolor' , colorz{2} , 'errorcolor' , colorz{2})
        set(gca , 'XTickLabel' , repmat(daylab(subjNum).dl , 1 , SN))
        title('Beta Band Power Across days for for Random(blue) and Structured(red)')
        grid on
        
       
        disp('***********************************************************************************')
        disp('****************************** CHOSE THE TEST TO RUN ******************************')
        disp('***********************************************************************************')
        keyboard
        
        %% significance tests
        %% testing the effect of sequence type and days on high gamma
        bandnumb = 4;
        seqNumb = [30 20];
        A = getrow(D , ismember(D.band, bandnumb)  & ismember(D.seqNumb , seqNumb ));
        %         D.blockGroup(D.blockGroup>1) = 0;     % test day one vs all other days
        if length(seqNumb)>1
            [p,tbl,stats] =  anovan(A.data , [A.blockgroup A.seqNumb],'model','interaction','varnames',{'Days' , 'Sequence Type'});
            figure;lineplot([A.blockgroup A.seqNumb] , A.data)
        else
            [p,tbl,stats] =  anovan(A.data , [A.blockgroup],'model','interaction','varnames',{'Days'});
            figure;lineplot([A.blockgroup] , A.data)
        end
        close all
    case 'singleFing_across_days'
        seq = input('Slow or Fast? (s/f)', 's');
        switch seq
            case {'s'}
                allSubj_BlockGroup(1).bg = {'SingleFingSlow1','SingleFingSlow2','SingleFingSlow3','SingleFingSlow4'};
                allSubj_BlockGroup(2).bg = {'SingleFingSlow1','SingleFingSlow2','SingleFingSlow3'};
                allSubj_BlockGroup(1).bg = {'SingleFingSlow1','SingleFingSlow2','SingleFingSlow3','SingleFingSlow4'};
                
            case {'f'}
                allSubj_BlockGroup(1).bg = {'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3', 'SingleFingFast4' };
                allSubj_BlockGroup(2).bg = {'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3'};
                allSubj_BlockGroup(1).bg = {'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3', 'SingleFingFast4' };
        end
        
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped.mat
        % patterned PSDs so the output of Pall  = seqeeg_parseEEG_PSD('TimeWarpPSD' , Dall, subjNum);
        % for this case because we are doigna comparisn, the BlockGroup is
        % a 1xN cell containnig the groups to be compared

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
        T = getrow(D , D.seqNumb~=100);   
        SN = length(Pall.PowTR);
        figure('color' , 'white')
        subplot(311)
        lineplot([T.seqNumb T.blockgroup] , T.data , 'subset' , ismember(T.band, 6) , 'style_thickline')
        set(gca , 'XTickLabel' , repmat(daylab(subjNum).dl , 1 , SN))
        title('High Gamma Band Power Across days for finger repetition')
        grid on
        
        subplot(312)
        lineplot([T.seqNumb T.blockgroup] , T.data , 'subset' , ismember(T.band, 5) , 'style_thickline')
        set(gca , 'XTickLabel' , repmat(daylab(subjNum).dl , 1 , SN))
        title('Low Gamma Band Power Across days for finger repetition')
        grid on
        
        subplot(313)
        lineplot([T.seqNumb T.blockgroup] , T.data , 'subset' , ismember(T.band, 4) , 'style_thickline')
        set(gca , 'XTickLabel' , repmat(daylab(subjNum).dl , 1 , SN))
        title('Beta Band Power Across days for finger repetition')
        grid on
        
       
        disp('***********************************************************************************')
        disp('****************************** CHOSE THE TEST TO RUN ******************************')
        disp('***********************************************************************************')
        keyboard
        
        %% significance tests
        %% testing the effect of sequence type and days on high gamma
        bandnumb = 4;
        seqNumb = [10];
        A = getrow(D , ismember(D.band, bandnumb)  & ismember(D.seqNumb , seqNumb ));
        if length(seqNumb)>1
            [p,tbl,stats] =  anovan(A.data , [A.blockgroup A.seqNumb],'model','interaction','varnames',{'Days' , 'Sequence Type'});
            figure
            lineplot([A.blockgroup A.seqNumb],A.data);
        else
            [p,tbl,stats] =  anovan(A.data , [A.blockgroup],'model','interaction','varnames',{'Days'});
            figure
            lineplot([A.blockgroup],A.data);
        end
        close all
    case 'chunks_across_days'
        % Pall needs to be the structure containing time normalized, average
        % Pall is the AllData_PSD_Warped.mat
        % patterned PSDs so the output of Pall  = seqeeg_parseEEG_PSD('TimeWarpPSD' , Dall, subjNum);
        % for this case because we are doigna comparisn, the BlockGroup is
        % a 1xN cell containnig the groups to be compared
        allSubj_BlockGroup(1).bg = {'ChunkDay1','ChunkDay2','ChunkDay3'};
        allSubj_BlockGroup(2).bg = {'ChunkDay1','ChunkDay2','ChunkDay3'};
        BlockGroup = allSubj_BlockGroup(subjNum).bg;
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
        set(gca , 'XTickLabel' , repmat(daylab(subjNum).dl, 1 , max(T.eventNum)))
        title(['High Gamma Band Power Across days for ' , chnkLabel])
        grid on
        
        subplot(312)
        T = getrow(D , D.band == 5);  % low gamma
        lineplot([T.eventNum T.blockGroup] , T.avgPSD , 'style_thickline' , 'linecolor' , 'b' , 'errorcolor' , 'b' , 'leg' , 'Random')
        set(gca , 'XTickLabel' , repmat( daylab(subjNum).dl, 1 , max(T.eventNum)))
        title(['Low Gamma Band Power Across days for ' , chnkLabel])
        
        grid on
        
        subplot(313)
        T = getrow(D , D.band == 4);  % betta
        lineplot([T.eventNum T.blockGroup] , T.avgPSD , 'style_thickline' , 'linecolor' , 'b' , 'errorcolor' , 'b' , 'leg' , 'Random')
        set(gca , 'XTickLabel' , repmat( daylab(subjNum).dl, 1 , max(T.eventNum)))
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
%         A.blockGroup(A.blockGroup>1) = 0;     % test day one vs all other days
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
        chans{1} = [4:12 ,14, 122:129 ,36] ; % channels for patient P2
        chans{2} = [15:21 93:99 101:110 112 28]; % channels for patient P4
        seqeeg_visualizePSD([],subjNum,'Aligned_seqType_average' , 'BlockGroup' , [],'Chan2Plot' , Chan,'Channels' , chans{subjNum});

        allSubj_BlockGroup(1).bg = {'Intermixed1','Intermixed5','Intermixed8','Intermixed9'};
        allSubj_BlockGroup(2).bg = {'Intermixed1','Intermixed4','Intermixed7'};
        allSubj_BlockGroup(3).bg = {'Intermixed1','Intermixed4','Intermixed7','Intermixed10'};
        BlockGroup = allSubj_BlockGroup(subjNum).bg;

        NT = input('Plot according to SeqNumb or SeqType? (n/t)', 's');
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
            if strcmp(NT , 'n')
                D = getrow(Dall , ismember(Dall.BN , unique(Pall.BN)));
                D.seqNumb(D.seqNumb == 2) = 1; % seqNumb 1  and 2 have the same structure
                
                D.seqNumb(D.seqNumb == 4) = 3; % seqNumb 4  and 3 have the same structure
                Pall.seqNumb = D.seqNumb;
                Pall = getrow(Pall , Pall.seqNumb ~= 5 & ~Pall.isError);
            else
                Pall = getrow(Pall , Pall.seqNumb ~= 100 & ~Pall.isError);
                seqInfo.nums = sequenceType.TypeNums;
                seqInfo.tags = sequenceType.Typetags;
            end
            SN = unique(Pall.seqNumb);
            for sn = 1:length(SN)    
                T = getrow(Pall , Pall.seqNumb == SN(sn));
                % reduce the band and channel dimensions
                for evnt = 1:T.seqlength(1)+1
                    for tn = 1:length(T.TN)
                        for bands = 1:size(T.PSD1{tn , evnt} , 2)
                            temp = squeeze(T.PSD1{tn , evnt}(:,bands,:));
                            eval(['T.reducPSD' , num2str(bands) , '(tn , evnt) = nanmean(nanmean(temp(Chan , :)));'])
                        end
                    end
                end
                T.blockgroup = bg*ones(size(T.BN));
                C = addstruct(C , T);
            end
        end
        % account for the first event which is the stim onset
        C.ChnkPlcmnt = [zeros(length(C.ChnkPlcmnt) , 1) , C.ChnkPlcmnt ];
        C.ChnkPlcmnt(:,2) = 0; % make sure the first press does not enter analysis
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
            T.ChnkPlcmnt = [];
            for bands = 1:size(temp.PSD1{1} , 2)
               eval(['A = temp.reducPSD' , num2str(bands) , ';']);
               T.avgPSD = [T.avgPSD ; A'];
               T.eventNum = [T.eventNum ; [1:length(A)]'];
               T.ChnkPlcmnt = [T.ChnkPlcmnt ; temp.ChnkPlcmnt'];
               T.seqNumb = [T.seqNumb ; temp.seqNumb*ones(size(A'))];
               T.band = [T.band; bands*ones(size(A'))];
               T.blockGroup = [T.blockGroup ; temp.blockgroup*ones(size(A'))];
               T.TN = [T.TN ; temp.TN*ones(size(A'))];
            end
            D = addstruct(D , T);
        end
        
        disp('***********************************************************************************')
        disp('****************************** CHOSE THE TEST TO RUN ******************************')
        disp('***********************************************************************************')
        keyboard
        %% significance tests
        %% testing the effect of sequence type and days on high gamma (6) lowgamma (5) beta (4)
        % P = anovaMixed(D.data , [1:length(D.data)]','between',[D.blockgroup D.seqNumb],...
        %     {'BlockGroup' , 'Sequence Type'},'intercept',1 , 'subset',ismember(D.frequency , [BandInfo.bands{6}(1) : BandInfo.bands{6}(2)]));
        bandnumb = 6; 
        EventsofInterest =[2:8];
        seqNumb = [20 30];
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

        seq = input('Fast or Slow? (f/s)', 's');
        switch seq
            case {'f'}
                allSubj_BlockGroup(1).bg = {'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3', 'SingleFingFast4'};
                allSubj_BlockGroup(2).bg = {'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3'};
                allSubj_BlockGroup(3).bg = {'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3', 'SingleFingFast4'};
            case {'s'}
                allSubj_BlockGroup(1).bg = {'SingleFingSlow1' , 'SingleFingSlow2'  , 'SingleFingSlow3' ,'SingleFingSlow4'};
                allSubj_BlockGroup(2).bg = {'SingleFingSlow1' , 'SingleFingSlow2'  , 'SingleFingSlow3'};
                allSubj_BlockGroup(3).bg = {'SingleFingSlow1' , 'SingleFingSlow2'  , 'SingleFingSlow3' ,'SingleFingSlow4'};
        end
        zBlockGroup = allSubj_BlockGroup(subjNum).bg;
        C = [];
        T = [];
        whatToCompare = input('Compare around press-time or between presses?(a/b)', 's');

        for bg = 1:length(BlockGroup)
            load([mainDir , 'AllData_AvgMarker_SeqType.mat'])
            
            BG(bg) = find(strcmp(blockGroupNames , BlockGroup{bg}));
            E = getrow(E , BG(bg));
            if isempty(E.EM)
                error('Nothing to plot!')
            end
            filename = [mainDir ,  'Group_Aligned_PSD' , num2str(BG(bg)),'.mat'];
            load(filename);
            switch whatToCompare
                case {'a'}
                    eval('Pall.PSD = Pall.PSD1;') % around presses
                case {'b'}
                    eval('Pall.PSD = Pall.PSD2;') % between Presses
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
        set(gca , 'XTickLabel' , repmat(daylab(subjNum).dl , 1 , max(T.eventNum)))
        title('High Gamma Band Power Across days for Single finger repetitions')
        grid on
        
        subplot(312)
        T = getrow(D , D.band == 5);  % low gamma
        lineplot([T.eventNum T.blockGroup] , T.avgPSD , 'style_thickline' , 'linecolor' , 'b' , 'errorcolor' , 'b' , 'leg' , 'Random')
        set(gca , 'XTickLabel' , repmat(daylab(subjNum).dl , 1 , max(T.eventNum)))
        title('Low Gamma Band Power Across days for Single finger repetitions')
        grid on
        
        subplot(313)
        T = getrow(D , D.band == 4 );  % betta
        lineplot([T.eventNum T.blockGroup] , T.avgPSD , 'style_thickline' , 'linecolor' , 'b' , 'errorcolor' , 'b' , 'leg' , 'Random')
        set(gca , 'XTickLabel' , repmat(daylab(subjNum).dl , 1 , max(T.eventNum)))
        title('Beta Band Power Across days for Single finger repetitions')
        grid on
        disp('***********************************************************************************')
        disp('****************************** CHOSE THE TEST TO RUN ******************************')
        disp('***********************************************************************************')
        keyboard
        %% significance tests
        %% testing the effect of sequence type and days on high gamma (6) lowgamma (5) beta (4)
        bandnumb = 6; 
        EventsofInterest =[2:8];
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
    case 'ChunkPlacement'
        
        allSubj_BlockGroup(1).bg = {'Intermixed1','Intermixed5','Intermixed8','Intermixed9'};
        allSubj_BlockGroup(2).bg = {'Intermixed1','Intermixed4','Intermixed7'};
        allSubj_BlockGroup(3).bg = {'Intermixed1','Intermixed4','Intermixed7','Intermixed10'};
        BlockGroup = allSubj_BlockGroup(subjNum).bg;
        whatToCompare = input('Compare around press-time or between presses?(a/b)', 's');
        
        NT = 'n';
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
            switch whatToCompare
                case {'a'}
                    eval('Pall.PSD = Pall.PSD1;') % around presses
                case {'b'}
                    eval('Pall.PSD = Pall.PSD2;') % between Presses
            end
            if strcmp(NT , 'n')
                D = getrow(Dall , ismember(Dall.BN , unique(Pall.BN)));
                D.seqNumb(D.seqNumb == 2) = 1; % seqNumb 1  and 2 have the same structure
                D.seqNumb(D.seqNumb == 4) = 3; % seqNumb 4  and 3 have the same structure
                Pall.seqNumb = D.seqNumb;
                Pall = getrow(Pall , Pall.seqNumb ~= 5 & ~Pall.isError);
            else
                Pall = getrow(Pall , Pall.seqNumb ~= 100 & ~Pall.isError);
                seqInfo.nums = sequenceType.TypeNums;
                seqInfo.tags = sequenceType.Typetags;
            end
            SN = unique(Pall.seqNumb);
            for sn = 1:length(SN)    
                T = getrow(Pall , Pall.seqNumb == SN(sn));
                % reduce the band and channel dimensions
                for evnt = 1:T.seqlength(1)+1
                    for tn = 1:length(T.TN)
                        for bands = 1:size(T.PSD{tn , evnt} , 2)
                            temp = squeeze(T.PSD{tn , evnt}(:,bands,:));
                            eval(['T.reducPSD' , num2str(bands) , '(tn , evnt) = nanmean(nanmean(temp(Chan , :)));'])
                            switch T.seqNumb(tn)
                                case {1 2}
                                    T.ChnkPlcmnt(tn , :) = [-1 2 3 1 2 2 -2]; % [null within last first within within null]
                                case {3 4}
                                    T.ChnkPlcmnt(tn , :) = [-1 2 2 3 1 2 -2]; % [null within within last first within null]
                                case {0}
                                    T.ChnkPlcmnt(tn , :) = [-1 0 0 0 0 0 -2]; % [null within within last first within null]
                            end
                        end
                    end
                end
                T.blockgroup = bg*ones(size(T.BN));
                C = addstruct(C , T);
            end
        end
        % account for the first event which is the stim onset
        C.ChnkPlcmnt = [-ones(length(C.ChnkPlcmnt) , 1) , C.ChnkPlcmnt ];
        %% prepare data for testing Block Groups against each other
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
            T.ChnkPlcmnt = [];
            for bands = 1:size(temp.PSD{1} , 2)
               eval(['A = temp.reducPSD' , num2str(bands) , ';']);
               T.avgPSD = [T.avgPSD ; A'];
               T.eventNum = [T.eventNum ; [1:length(A)]'];
               T.ChnkPlcmnt = [T.ChnkPlcmnt ; temp.ChnkPlcmnt'];
               T.seqNumb = [T.seqNumb ; temp.seqNumb*ones(size(A'))];
               T.band = [T.band; bands*ones(size(A'))];
               T.blockGroup = [T.blockGroup ; temp.blockgroup*ones(size(A'))];
               T.TN = [T.TN ; temp.TN*ones(size(A'))];
            end
            D = addstruct(D , T);
        end
        %% prepare data for testing Block Groups against eachother
        includeFirstLast = 1;
        figure('color' , 'white')
        SN = unique(D.seqNumb);
        CP = unique(D.ChnkPlcmnt);
        C = D;
        if includeFirstLast
            id = C.ChnkPlcmnt==-1 & ismember(C.seqNumb , [1 2 3 4]);
            C.ChnkPlcmnt(id) = 1;
            id = C.ChnkPlcmnt==-2 & ismember(C.seqNumb , [1 2 3 4]);
            C.ChnkPlcmnt(id) = 3;
            id = ismember(C.ChnkPlcmnt ,[-1 -2]) & ismember(C.seqNumb , [0]);
            C.ChnkPlcmnt(id) = 0;
        end
        CP = CP(~ismember(CP , [-1 -2]));
        subplot(311)
        for cp = 1:length(CP)
            T = getrow(D , D.band == 6 & D.ChnkPlcmnt == CP(cp));  % high gamma
            lineplot([T.ChnkPlcmnt T.blockGroup] , T.avgPSD , 'style_thickline' , 'linecolor' , colorz{cp} , 'errorcolor' , colorz{cp})
            hold on
        end
        set(gca , 'XTickLabel' , repmat(daylab(subjNum).dl , 1 , max(T.eventNum)))
        title('High Gamma Band Power Across days for        Random (blue),      first (red),      middle (green),       last (cyan)')
        grid on
        
        subplot(312)
        for cp = 1:length(CP)
            T = getrow(D , D.band == 5 & D.ChnkPlcmnt == CP(cp));  % high gamma
            lineplot([T.ChnkPlcmnt T.blockGroup] , T.avgPSD , 'style_thickline' , 'linecolor' , colorz{cp} , 'errorcolor' , colorz{cp})
            hold on
        end
        set(gca , 'XTickLabel' , repmat(daylab(subjNum).dl , 1 , max(T.eventNum)))
         title('Low Gamma Band Power Across days for        Random (blue),      first (red),      middle (green),       last (cyan)')
        grid on
        
        subplot(313)
        for cp = 1:length(CP)
            T = getrow(D , D.band == 4 & D.ChnkPlcmnt == CP(cp));  % high gamma
            lineplot([T.ChnkPlcmnt T.blockGroup] , T.avgPSD , 'style_thickline' , 'linecolor' , colorz{cp} , 'errorcolor' , colorz{cp})
            hold on
        end
        set(gca , 'XTickLabel' , repmat(daylab(subjNum).dl , 1 , max(T.eventNum)))
         title('Beta Band Power Across days for        Random (blue),      first (red),      middle (green),       last (cyan)')
        grid on
        disp('***********************************************************************************')
        disp('****************************** CHOSE THE TEST TO RUN ******************************')
        disp('***********************************************************************************')
        keyboard
        %% significance tests
        %% testing the effect of sequence type and days on high gamma (6) lowgamma (5) beta (4)
        % P = anovaMixed(D.data , [1:length(D.data)]','between',[D.blockgroup D.seqNumb],...
        %     {'BlockGroup' , 'Sequence Type'},'intercept',1 , 'subset',ismember(D.frequency , [BandInfo.bands{6}(1) : BandInfo.bands{6}(2)]));
        bandnumb = 5; 
        seqNumb = [0:4];
        CombineSeqs = 1;     % combine all the seqNums to look at placement effect
        %______________________________________________________________ DEFINE TYPE OF TEST
        %______________________________________________________________
        firstVsmiddleChunk = 0;     % take out the random presses and just look at the difference between within vs first
        lastVsmiddleChunk = 0;     % take out the random presses and just look at the difference between within vs last
        GeneralChunking  = 1;     % take out the random presses and look at the difference between within, middle and last
        RanVsSeq = 0;        % look at the difference between random vs seq
        begVsMiddle = 0;  % look at the difference beginning the seq to the middle
        EndVsMiddle = 0;  % look at the difference end the seq to the middle
        LastVsFirstDay = 0; % just  look at the oveeal effect of taraining
        %______________________________________________________________ 
        %______________________________________________________________ 
        if LastVsFirstDay
            ubg = unique(C.blockGroup);
            bg_f = min(ubg);
            bg_l = max(ubg);
            bgofInterest = [bg_f , bg_l];
        else
            bgofInterest = unique(C.blockGroup);
        end
        if CombineSeqs 
            C.seqNumb = ones(length(C.seqNumb),1); 
            seqNumb = 1;
        end
        if firstVsmiddleChunk
            CofInterest =[1 2];
            A = getrow(C , ismember(C.band, bandnumb) & ismember(C.ChnkPlcmnt , CofInterest )  & ismember(C.seqNumb , seqNumb ) & ismember(C.blockGroup , bgofInterest ));
            A.ChnkPlcmnt(ismember(A.ChnkPlcmnt , [2 3])) = 0;
            CofInterest = [0 1];
        end
        
        if GeneralChunking
            CofInterest =[1  2  3];
            A = getrow(C , ismember(C.band, bandnumb) & ismember(C.ChnkPlcmnt , CofInterest )  & ismember(C.seqNumb , seqNumb ) & ismember(C.blockGroup , bgofInterest ));

        end
        if firstVsmiddleChunk
            CofInterest =[1 2];
            A = getrow(C , ismember(C.band, bandnumb) & ismember(C.ChnkPlcmnt , CofInterest )  & ismember(C.seqNumb , seqNumb ) & ismember(C.blockGroup , bgofInterest ));
            A.ChnkPlcmnt(ismember(A.ChnkPlcmnt , [2 3])) = 0;
            CofInterest = [0 1];
        end
        if begVsMiddle
            C= D;
            CofInterest =[-1 0 1 2 3];
            A = getrow(C , ismember(C.band, bandnumb) & ismember(C.ChnkPlcmnt , CofInterest )  & ismember(C.seqNumb , seqNumb ) & ismember(C.blockGroup , bgofInterest ));
            A.ChnkPlcmnt(ismember(A.ChnkPlcmnt , [1 2 3])) = 0;
            CofInterest = [0 -1];
        end
        if EndVsMiddle
            C= D;
            CofInterest =[-2 0 1 2 3];
            A = getrow(C , ismember(C.band, bandnumb) & ismember(C.ChnkPlcmnt , CofInterest )  & ismember(C.seqNumb , seqNumb ) & ismember(C.blockGroup , bgofInterest ));
            A.ChnkPlcmnt(ismember(A.ChnkPlcmnt , [1 2 3])) = 0;
            CofInterest = [0 -2];
        end
        if RanVsSeq
            CofInterest =[0 1 2 3];
            A = getrow(C , ismember(C.band, bandnumb) & ismember(C.ChnkPlcmnt , CofInterest )  & ismember(C.seqNumb , seqNumb ) & ismember(C.blockGroup , bgofInterest ));
            A.ChnkPlcmnt(ismember(A.ChnkPlcmnt , [1 2 3])) = 1;
            CofInterest = [0 1];
        end
        
        %______________________________________________________________
        %______________________________________________________________
        % test day one vs all other days
        if length(CofInterest)>1
            if length(seqNumb)>1
                [p,tbl,stats] =  anovan(A.avgPSD , [A.blockGroup A.seqNumb A.ChnkPlcmnt],'model','interaction','varnames',{'Days' , 'Sequence Type' , 'ChnkPlcmnt'});
                figure;
                lineplot([A.blockGroup A.seqNumb A.ChnkPlcmnt] , A.avgPSD)
            else
                [p,tbl,stats] =  anovan(A.avgPSD , [A.blockGroup A.ChnkPlcmnt],'model','interaction','varnames',{'Days' , 'ChnkPlcmnt'});
                figure;
                lineplot([A.blockGroup A.ChnkPlcmnt] , A.avgPSD)
            end
        else
            if length(seqNumb)>1
                [p,tbl,stats] =  anovan(A.avgPSD , [A.blockGroup A.seqNumb],'model','interaction','varnames',{'Days' , 'Sequence Type'});
                figure;
                lineplot([A.blockGroup A.seqNumb] , A.avgPSD)
            else
                [p,tbl,stats] =  anovan(A.avgPSD , [A.blockGroup],'model','interaction','varnames',{'Days'});
                figure;
                lineplot([A.blockGroup] , A.avgPSD)
            end
        end
    case 'SeqFingRepetition'
        allSubj_BlockGroup(1).bg = {'SingleFingSlow1' , 'SingleFingSlow2'  , 'SingleFingSlow3' ,'SingleFingSlow4'};
        allSubj_BlockGroup(2).bg = {'SingleFingSlow1' , 'SingleFingSlow2'  , 'SingleFingSlow3'};
        allSubj_BlockGroup(3).bg = {'SingleFingSlow1' , 'SingleFingSlow2'  , 'SingleFingSlow3' ,'SingleFingSlow4'};
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
%         T = tapply(D , {'changroup','blockgroup','seqNumb','TN','band' },{'data'} );
        T = D;
        SN = length(Pall.PowTR);
        figure('color' , 'white')
        subplot(311)
        
        lineplot([D.seqNumb D.blockgroup] , D.data , 'subset' , ismember(T.band, 6) , 'style_thickline' , 'subset' , D.seqNumb == 20 , 'linecolor' , colorz{1} , 'errorcolor' , colorz{1})
        hold on
        lineplot([D.seqNumb D.blockgroup] , D.data , 'subset' , ismember(T.band, 6) , 'style_thickline' , 'subset' , D.seqNumb == 30 , 'linecolor' , colorz{2} , 'errorcolor' , colorz{2})
        set(gca , 'XTickLabel' , repmat(daylab(subjNum).dl , 1 , SN))
        title('High Gamma Band Power Across days for Random, Structured and Null Trails respectively left to right')
        grid on
        
        subplot(312)
        lineplot([D.seqNumb T.blockgroup] , D.data , 'subset' , ismember(D.band, 5) , 'style_thickline','subset' , D.seqNumb == 20)% , 'linecolor' , colorz{1} , 'errorcolor' , colorz{1})
        hold on
        lineplot([D.seqNumb D.blockgroup] , D.data , 'subset' , ismember(T.band, 5) , 'style_thickline','subset' , D.seqNumb == 30 , 'linecolor' , colorz{2} , 'errorcolor' , colorz{2})
        set(gca , 'XTickLabel' , repmat(daylab(subjNum).dl , 1 , SN))
        title('Low Gamma Band Power Across days for Random, Structured and Null Trails respectively left to right')
        grid on
        
        subplot(313)
        lineplot([T.seqNumb T.blockgroup] , T.data , 'subset' , ismember(T.band, 4) , 'style_thickline','subset' ,D.seqNumb == 20 , 'linecolor' , colorz{1} , 'errorcolor' , colorz{1})
        hold on
        lineplot([T.seqNumb T.blockgroup] , T.data , 'subset' , ismember(T.band, 4) , 'style_thickline','subset' ,D.seqNumb == 30 , 'linecolor' , colorz{2} , 'errorcolor' , colorz{2})
        set(gca , 'XTickLabel' , repmat(daylab(subjNum).dl , 1 , SN))
        title('Beta Band Power Across days for Random, Structured and Null Trails respectively left to right')
        grid on
        
       
        disp('***********************************************************************************')
        disp('****************************** CHOSE THE TEST TO RUN ******************************')
        disp('***********************************************************************************')
        keyboard
        
        %% significance tests
        %% testing the effect of sequence type and days on high gamma
        bandnumb = 6;
        seqNumb = [30 20];
        A = getrow(D , ismember(D.band, bandnumb)  & ismember(D.seqNumb , seqNumb ));
        %         D.blockGroup(D.blockGroup>1) = 0;     % test day one vs all other days
        if length(seqNumb)>1
            [p,tbl,stats] =  anovan(A.data , [A.blockgroup A.seqNumb],'model','interaction','varnames',{'Days' , 'Sequence Type'});
            lineplot()
        else
            [p,tbl,stats] =  anovan(A.data , [A.blockgroup],'model','interaction','varnames',{'Days'});
        end
        close all
end


