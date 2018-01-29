function Dout  = secog_addEventMarker(Dall, subjNum, Fs , what,varargin)
% adds event markers for the EEG / PSD data pased on press times and the sampling frequency
% make sure to account for downsampling in Fs

subjname = {'P2'};
mainDir = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' subjname{subjNum} ,'/'];
load([mainDir , 'ChanLabels.mat'])
c= 1;
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
Dout.blockGroups = {[1 2] , [3 13 26 40] , [4 14 27 41] , [5:7] , [9:11] , [8 12] , [15:17] , [19:21] , [23:25],...
            [18 22] , [28:30] , [32:34] , [36:38], [31 35 39],[42:44]}';
        

Dout.blockGroupNames = {'SingleFingNat' , 'SingleFingSlow' , 'SingleFingFast' , 'Intermixed1' , 'Intermixed2' , 'ChunkDay1' , 'Intermixed3' , 'Intermixed4' , 'Intermixed5',...
            'ChunkDay2' , 'Intermixed6' , 'Intermixed7' , 'Intermixed8', 'ChunkDay3', 'Intermixed9'}';
Dall.Fast = zeros(size(Dall.TN));
fastBlock = horzcat(Dout.blockGroups{1} , Dout.blockGroups{2} , Dout.blockGroups{3} , Dout.blockGroups{6}, Dout.blockGroups{10},Dout.blockGroups{14});
Dall.Fast(ismember(Dall.BN , fastBlock)) = 1;
switch what
    case 'addEvent'
        %% find the non-normal and normal event marker and averages
        % can be performed on individual blocks or trails

        for i = 1:length(Dall.TN)
            if Dall.Fast(i)
                NumWarpSamp = NumWarpSampFast;
            else
                NumWarpSamp = NumWarpSampSlow;
            end
            if isfield(Dall , 'EEG')
                eegTime = length(Dall.EEG{i});
            elseif isfield(Dall , 'PSD')
                eegTime = size(Dall.PSD{i} , 3);
            elseif isfield(Dall , 'Pow_Norm_stim')
                eegTime = size(Dall.Pow_Norm_stim{i} , 3);
            elseif isfield(Dall , 'PSD_stim')
                eegTime = size(Dall.PSD_stim{i} , 3);
            end
            % add 200 ms to the event times because we chopp EEG 200 ms before the stim comes on'
            % EEGEventTime will contain [0 ,  time the stimulus comes on(200 ms mark) , event times(if any) ,
            % time the stimulus goes off (simultaneous to the last press) , endtime (200ms after the stim goes off)]
            if Dall.TN(i)>2
                EEGEventTime = [0 .2 , (Dall.AllPressTimes(i,~isnan(Dall.AllPressTimes(i,:))) +200)/1000 , eegTime/Fs];
            else
                EEGEventTime = [0 .2 , (Dall.AllPressTimes(i,~isnan(Dall.AllPressTimes(i,:))) +200)/1000 , (eegTime/Fs)-.2 eegTime/Fs];
            end
            EEGEventidx = floor(EEGEventTime* Fs);
            % warp the events to the interval of interest.

            NEEGEventidx = floor([EEGEventidx/max(EEGEventidx)]*NumWarpSamp);

            for pr = 3:length(EEGEventidx)-1
                E.EventMarker{i , 1}(EEGEventidx(pr)) = pr-2;  % First press 
                E.NormEventMarker{i , 1}(NEEGEventidx(pr)) = pr-2;
            end
            E.EventMarker{i , 1}(EEGEventidx(2)) = -1;  % First press
            E.NormEventMarker{i , 1}(NEEGEventidx(2)) = -1;
        end
        Dall.EventMarker = E.EventMarker;
        Dall.NormEventMarker = E.NormEventMarker;
        Dout = Dall;
    case 'CalcAveragePattern'
        
        %% Calculating the average pattern
        Dall  = secog_addEventMarker(Dall, subjNum, Fs , 'addEvent')';

        % we want the average pattern to be calculated across the groups of 3 block of seq training in between the chunk training blocks
        % this should be defined for every participant

        for BG = 1:length(Dout.blockGroups)
            D1 = getrow(Dall , ~Dall.isError & ismember(Dall.BN , Dout.blockGroups{BG}));
            Dout.SN{BG,1} = unique(D1.seqNumb);
            Dout.SN{BG,1} = Dout.SN{BG,1};
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
                Dout.NEM{BG,1}(sn ,1:SL+1) = floor(nanmean(NormEventMarker));
            end
        end
end
