function Dout  = secog_addEventMarker(Dall, subjNum, Fs , what,varargin)
% adds event markers for the EEG / PSD data pased on press times and the sampling frequency
% make sure to account for downsampling in Fs

subjname = {'P2' , 'P4'};
mainDir = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' subjname{subjNum} ,'/'];
load([mainDir , 'ChanLabels.mat'])
c= 1;
NumWarpSampFast = 200;
NumWarpSampSlow = 500;
TimeDelay = 0.5; % sec
while(c<=length(varargin))
    switch(varargin{c})
        case {'NumWarpSampFast'}
            % Number of warping sample for chunks  Default = 200
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'NumWarpSampSlow'}
            % Number of warping sample for sequences  Default = 500
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'TimeDelay'}
            % The time delay before the stimulus comes on to consider for baseline normalization
            % Default  = 0.5sec
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        otherwise
            error(sprintf('Unknown option: %s',varargin{c}));
    end
end
% block groupings for subject 1
BG(1).blockGroups =  {[1 2] , [3], [13], [26], [40] , [4], [14], [27] [41] , [5:7] , [9:11] , [8 12] , [15:17] , [19:21] , [23:25],...
    [18 22] , [28:30] , [32:34] , [36:38], [31 35 39],[42:44]}';
% block groupings for subject 2
BG(2).blockGroups = {[ ] , [2 8], [14 20 26], [29 38], [], [1 7],[13 19 25], [28 37], [] , [3:5] , [9:11] , [6 12] , [15:17] , [21:23] , [],...
    [18 24] , [30:32] , [34:36] , [], [27 33],[]}';

Dout.blockGroups = BG(subjNum).blockGroups;
Dout.blockGroupNames = {'SingleFingNat' , 'SingleFingSlow1' , 'SingleFingSlow2'  , 'SingleFingSlow3' ,'SingleFingSlow4',...
    'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3', 'SingleFingFast4' , 'Intermixed1' , 'Intermixed2' , ...
    'ChunkDay1' , 'Intermixed3' , 'Intermixed4' , 'Intermixed5', 'ChunkDay2' , 'Intermixed6' , ...
    'Intermixed7' , 'Intermixed8', 'ChunkDay3', 'Intermixed9'}';
Dall.Fast = zeros(size(Dall.TN));
fastBlock = horzcat(Dout.blockGroups{1} ,Dout.blockGroups{6} , Dout.blockGroups{7},Dout.blockGroups{8},Dout.blockGroups{9},...
    Dout.blockGroups{12},Dout.blockGroups{16},Dout.blockGroups{20});
Dall.Fast(ismember(Dall.BN , fastBlock)) = 1;
%%  control for too short IPIs that the keys get accidentally pressed
if subjNum==1
    for i = 1:length(Dall.TN)
        if sum(Dall.IPI(i,:)<120)
            Dall.isError(i) =1;
        end
    end
end
% the length of the Null trials was diferent for P2 and P4
NullTrailTime = [15 7]; % sec
Dall.seqlength(Dall.seqNumb == 5) = 0;
switch what
    case 'addEvent'
        %% find the non-normal and normal event marker and averages
        % can be performed on individual blocks or trails
        
        for i = 1:length(Dall.TN)
            if ~Dall.isError(i) && sum(~isnan(Dall.AllPressTimes(i,:))) == Dall.seqlength(i)
                if ismember(Dall.TN(i) , [1,2])
                    eegTime = floor((NullTrailTime(subjNum)+3*TimeDelay)*Fs);
                else
                    eegTime = floor((.001*(Dall.AllPressTimes(i , Dall.seqlength(i)))+3*TimeDelay)*Fs);
                end
                if Dall.Fast(i)
                    NumWarpSamp = NumWarpSampFast;
                else
                    NumWarpSamp = NumWarpSampSlow;
                end
                
                % add TimeDelay ms to the event times because we chopp EEG TimeDelay ms before the stim comes on'
                % EEGEventTime will contain [0 ,  time the stimulus comes on(TimeDelay ms mark) , event times(if any) ,
                % time the stimulus goes off (simultaneous to the last press) , endtime (TimeDelayms after the stim goes off)]
                E.EventMarker{i , 1} = zeros(1,eegTime);
                E.NormEventMarker{i , 1} = zeros(1,NumWarpSamp);
                
                if Dall.TN(i)>2
                    mark = ((Dall.AllPressTimes(i,~isnan(Dall.AllPressTimes(i,:))))/1000 + TimeDelay)*Fs;
                    EEGEventidx = floor([0 TimeDelay*Fs , mark , eegTime]);
                else
                    mark = (Dall.AllPressTimes(i,~isnan(Dall.AllPressTimes(i,:))))/1000 + TimeDelay;
                    EEGEventidx = floor([0 TimeDelay*Fs , eegTime - 2*TimeDelay*Fs  eegTime]);
                end
                % warp the events to the interval of interest.
                NEEGEventidx = floor([EEGEventidx/max(EEGEventidx)]*NumWarpSamp);
                
                for pr = 3:length(EEGEventidx)-1
                    E.EventMarker{i , 1}(EEGEventidx(pr)) = pr-2;  % First press
                    E.NormEventMarker{i , 1}(NEEGEventidx(pr)) = pr-2;
                end
                E.EventMarker{i , 1}(EEGEventidx(2)) = -1;  % stimCome on
                E.NormEventMarker{i , 1}(NEEGEventidx(2)) = -1;
            else
                E.EventMarker{i , 1} = NaN;  % stimCome on
                E.NormEventMarker{i , 1} = NaN;
            end
        end
        Dall.EventMarker = E.EventMarker;
        Dall.NormEventMarker = E.NormEventMarker;
        Dout = Dall;
    case 'CalcAveragePattern'
        
        %% Calculating the average pattern
        % the input to this case in the completed Bahavioral data with the
        % real and normalized time stamps.
        % this would be the saved structure throught
        % Dout = secog_parseEEG_PSD('ParseEEG-calcPSD' , Dall, 1); undeer AllData_Events.mat
        % we want the average pattern to be calculated across the groups of 3 block of seq training in between the chunk training blocks
        % this should be defined for every participant
        
        for BG = 1:length(Dout.blockGroups)
            D1 = getrow(Dall , ~Dall.isError & ismember(Dall.BN , Dout.blockGroups{BG}));
            Dout.SN{BG,1} = unique(D1.seqNumb);
            Dout.SN{BG,1} = Dout.SN{BG,1};
            if length(Dout.SN{BG,1})>0
                for sn= 1:length(Dout.SN{BG,1})
                    clear EventMarker NormEventMarker EM NEM diffNEM
                    EventMarker = [];
                    NormEventMarker = [];
                    D = getrow(D1 , ismember(D1.seqNumb , Dout.SN{BG,1}(sn)) & ~D1.isError & ismember(D1.BN , Dout.blockGroups{BG,1}));
                    D.seqlength(isnan(D.seqlength)) = 1;
                    SL = unique(D.seqlength);
                    Dout.EM{BG,1}(sn , :)= nan(1,max(Dall.seqlength)+1);
                    Dout.NEM{BG,1}(sn , :)= nan(1,max(Dall.seqlength)+1);
                    if ~isempty(D.TN)
                        for k = 1:length(D.TN)
                            temp = nan(1,max(Dall.seqlength)+1);
                            Ntemp = nan(1,max(Dall.seqlength)+1);
                            events = [-1 , 1:D.seqlength(k)];
                            clear temp Ntemp
                            for jj = 1:length(events)
                                if  ~isempty(find(D.EventMarker{k} == events(jj)))
                                    temp(1, jj) =  find(D.EventMarker{k} == events(jj));
                                else
                                    temp(1, jj) = NaN;
                                end
                                if ~isempty(find(D.NormEventMarker{k} == events(jj)))
                                    Ntemp(1,jj) =  find(D.NormEventMarker{k} == events(jj));
                                else
                                    Ntemp(1,jj) = NaN;
                                end
                            end
                            EventMarker = [EventMarker; temp];
                            NormEventMarker = [NormEventMarker ;Ntemp];
                        end
                    end
                    Dout.EM{BG,1}(sn , 1:SL+1) = floor(nanmean(EventMarker));
                    Dout.NEM{BG,1}(sn ,1:SL+1) = floor(nanmin(NormEventMarker));
                    
                end
            else
                Dout.EM{BG,1} = [];
                Dout.NEM{BG,1} = [];
                Dout.SN{BG,1} = [];
            end
            
        end
        
    case 'CalcAveragePattern_seqType'
        
        %% Calculating the average pattern for every sequence type as a general
        %% the goal here is to get average time pattern for general sequence types regardless of finger
        % so all the slow single fingers, fast singel finger, random, structured, triplets, quadruples
        % so we will change the SeqNumb and average the average times of the SeqNumbs
        % within the same type
        
        
        % the input to this case in the completed Bahavioral data with the
        % real and normalized time stamps.
        % this would be the saved structure throught
        % Dout = secog_parseEEG_PSD('ParseEEG-calcPSD' , Dall, 1); undeer AllData_Events.mat
        % we want the average pattern to be calculated across the groups of 3 block of seq training in between the chunk training blocks
        % this should be defined for every participant
        
        % Define sequence numbers and their transformations:
        SeqTrans = [5 11 22 33 44 55 0 1 2 3 4 103 104 203 204;...
            100 10 10 10 10 10 20 30 30 30 30 40 50 40 50];
        for sn = 1:length(SeqTrans)
            id  = Dall.seqNumb == SeqTrans(1 , sn);
            Dall.seqNumb(id) = SeqTrans(2 , sn);
        end
        for BG = 1:length(Dout.blockGroups)
            
            D1 = getrow(Dall , ~Dall.isError & ismember(Dall.BN , Dout.blockGroups{BG}));
            Dout.SN{BG,1} = unique(D1.seqNumb);
            Dout.SN{BG,1} = Dout.SN{BG,1};
            if length(Dout.SN{BG,1})>1
                for sn= 1:length(Dout.SN{BG,1})
                    clear EventMarker NormEventMarker EM NEM diffNEM
                    EventMarker = [];
                    NormEventMarker = [];
                    D = getrow(D1 , ismember(D1.seqNumb , Dout.SN{BG,1}(sn)) & ~D1.isError & ismember(D1.BN , Dout.blockGroups{BG,1}));
                    D.seqlength(isnan(D.seqlength)) = 1;
                    SL = unique(D.seqlength);
                    Dout.EM{BG,1}(sn , :)= nan(1,max(Dall.seqlength)+1);
                    Dout.NEM{BG,1}(sn , :)= nan(1,max(Dall.seqlength)+1);
                    if ~isempty(D.TN)
                        for k = 1:length(D.TN)
                            temp = nan(1,max(Dall.seqlength)+1);
                            Ntemp = nan(1,max(Dall.seqlength)+1);
                            events = [-1 , 1:D.seqlength(k)];
                            clear temp Ntemp
                            for jj = 1:length(events)
                                if  ~isempty(find(D.EventMarker{k} == events(jj)))
                                    temp(1, jj) =  find(D.EventMarker{k} == events(jj));
                                else
                                    temp(1, jj) = NaN;
                                end
                                if ~isempty(find(D.NormEventMarker{k} == events(jj)))
                                    Ntemp(1,jj) =  find(D.NormEventMarker{k} == events(jj));
                                else
                                    Ntemp(1,jj) = NaN;
                                end
                            end
                            EventMarker = [EventMarker; temp];
                            NormEventMarker = [NormEventMarker ;Ntemp];
                        end
                    end
                    Dout.EM{BG,1}(sn , 1:SL+1) = floor(nanmean(EventMarker));
                    Dout.NEM{BG,1}(sn ,1:SL+1) = floor(nanmin(NormEventMarker));
                    
                end
            else
                Dout.EM{BG,1} = [];
                Dout.NEM{BG,1} = [];
                Dout.SN{BG,1} = [];
            end
        end
end
