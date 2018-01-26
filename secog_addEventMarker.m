function Dout  = secog_addEventMarker(Dall, subjNum, Fs , what,varargin)
% adds event markers for the EEG / PSD data pased on press times and the sampling frequency
% make sure to account for downsampling in Fs

subjname = {'P2'};
mainDir = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' subjname{subjNum} ,'/'];
load([mainDir , 'ChanLabels.mat'])
c= 1;
while(c<=length(varargin))
    switch(varargin{c})
        case {'NumWarpSamp'}
            % number of sine cycles to use in the Morlet wavelet
            % default 10
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        otherwise
            error(sprintf('Unknown option: %s',varargin{c}));
    end
end

if ~exist('NumWarpSamp')
    NumWarpSamp = 300;
end
Dout.blockGroups = {[1 2] , [3 13 26 40] , [4 14 27 41] , [5:7] , [9:11] , [8 12] , [15:17] , [19:21] , [23:25],...
            [18 22] , [28:30] , [32:34] , [36:38], [31 35 39],[42:44]}';
        

Dout.blockGroupNames = {'SingleFingNat' , 'SingleFingSlow' , 'SingleFingFast' , 'Intermixed1' , 'Intermixed2' , 'ChunkDay1' , 'Intermixed3' , 'Intermixed4' , 'Intermixed5',...
            'ChunkDay2' , 'Intermixed6' , 'Intermixed7' , 'Intermixed8', 'ChunkDay3', 'Intermixed9'}';
switch what
    case 'addEvent'
        %% find the non-normal and normal event marker and averages
        % can be performed on individual blocks or trails 
        for i = 1:length(Dall.TN)
            if isfield(Dall , 'EEG')
                eegTime = length(Dall.EEG{i});
            elseif isfield(Dall , 'PSD')
                eegTime = size(Dall.PSD{i} , 3);
            end
            
            if Dall.TN(i)>2 && ~Dall.isError(i)
                Dall.AllPressIdx(i ,:) = (Dall.AllPressTimes(i,:)/1000)*Fs;
                Dall.EEGEventTime(i , :) = (Dall.AllPressTimes(i,:) +500)/1000; % because we chopp EEG 500 ms before the stim comes on'
                EEGEventidx = floor(Dall.EEGEventTime(i , :)* Fs);
                EEGEventidx =  [Fs/2 , EEGEventidx ];
                NEEGEventidx = floor([EEGEventidx/(max(EEGEventidx)+(Fs/2)) , 1]*NumWarpSamp);
                for pr = 1:length(EEGEventidx)
                    E.EventMarker{i , 1}(EEGEventidx(pr)) = pr-1;  % First press
                    E.NormEventMarker{i , 1}(NEEGEventidx) = pr-1;
                end
            else
                trialTime  = 0:1/Fs:(eegTime/Fs) - 1/Fs;
                
                E.EventMarker{i , 1} = zeros(size( trialTime));
                E.NormEventMarker{i , 1} = zeros(NumWarpSamp,1);
                idx = floor(Fs/2);
                E.EventMarker{i , 1}(idx) = -1;       % equivalent to stimulus apears (because we take .5sec before the stim appears)
                Nidx = floor(idx*NumWarpSamp/length(trialTime));
                E.NormEventMarker{i , 1}(Nidx) = -1;
                
                idx = floor(Fs/2);
                Nidx = floor(idx*NumWarpSamp/length(trialTime));
                E.EventMarker{i , 1}(end - idx) = 7;  % equivalent to last press (because we take .5sec after the trail's over)
                E.NormEventMarker{i , 1}(end - Nidx) = 7;
            end
        end
        Dall.EventMarker = E.EventMarker;
        Dall.NormEventMarker = E.NormEventMarker;
        Dout = Dall;
    case 'CalcAveragePattern'
        
        %% Calculating the average pattern
        Dall  = secog_addEventMarker(Dall, subjNum, Fs , 'addEvent')';
        % single finger blocks that participant is instructed to go slow 
        slowSingleFing = [3 13 26 40];
        BNid_slow  = logical(sum(Dall.BN == slowSingleFing , 2));
        Dall.seqNumb(BNid_slow) = Dall.seqNumb(BNid_slow)*10;
        
        % single finger blocks that participant is instructed to go fast as she can
        fastSingleFing = [4 14 27 41];
        BNid_fast  = logical(sum(Dall.BN == fastSingleFing , 2));
        Dall.seqNumb(BNid_fast) = Dall.seqNumb(BNid_fast)*100;
        
        
        % we want the average pattern to be calculated across the groups of 3 block of seq training in between the chunk training blocks
        % this should be defined for every participant
        
        
        for BG = 1:length(Dout.blockGroups)
            D1 = getrow(Dall , ~Dall.isError & ismember(Dall.BN , Dout.blockGroups{BG}));
            Dout.SN{BG,1} = unique(D1.seqNumb);
            Dout.SN{BG,1} = Dout.SN{BG,1}(~ismember(Dout.SN{BG,1} , [5 50 500])); % exclude the stars
            for sn= 1:length(Dout.SN{BG,1})
                clear EventMarker NormEventMarker EM NEM diffNEM
                EventMarker = [];
                NormEventMarker = [];
                Dout.EM{BG,1}(sn , :)= nan(1,max(Dall.seqlength)+1);
                Dout.NEM{BG,1}(sn , :)= nan(1,max(Dall.seqlength)+1);
                D = getrow(D1 , ismember(D1.seqNumb , Dout.SN{BG,1}(sn)) & ~D1.isError & ismember(D1.BN , Dout.blockGroups{BG,1}));
                SL = unique(D.seqlength);
                if ~isempty(D.TN)
                    for k = 1:length(D.TN)
                        temp = nan(1,max(Dall.seqlength)+1);
                        Ntemp = nan(1,max(Dall.seqlength)+1);
                        events = [-1 , 1:D.seqlength(k)];
                        if D.TN(k) >3
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
                    
                end
                Dout.EM{BG,1}(sn , 1:SL+1) = floor(nanmean(EventMarker));
                Dout.NEM{BG,1}(sn ,1:SL+1) = floor(nanmean(NormEventMarker));
            end
        end
end
