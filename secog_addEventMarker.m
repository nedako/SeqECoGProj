function Dout  = secog_addEventMarker(Dall, subjNum, Fs , what)
% adds event markers for the EEG / PSD data pased on press times and the sampling frequency
% make sure to account for downsampling in Fs

subjname = {'P2'};
mainDir = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' subjname{subjNum} ,'/'];
load([mainDir , 'ChanLabels.mat'])
while(c<=length(varargin))
    switch(varargin{c})
        case {'NumNormSamp'}
            % number of sine cycles to use in the Morlet wavelet
            % default 10
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        otherwise
            error(sprintf('Unknown option: %s',varargin{c}));
    end
end

if ~exist('NumNormSamp')
    NumNormSamp = 500;
end

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
            i
            if Dall.TN(i)>2 && ~Dall.isError(i)
                trialTime  = 0:1/Fs:(eegTime/Fs) - 1/Fs;
                E.EventMarker{i , 1} = zeros(size( trialTime));
                E.NormEventMarker{i , 1} = zeros(NumNormSamp,1);
                for pr = 1:Dall.seqlength(i)
                    numTaskSamp = floor(( Dall.AllPressIdx(i , Dall.seqlength(i)) - Dall.AllPressIdx(i , pr)) * (Fs/NumNormSamp));
                    idx = numTaskSamp + Fs/2;
                    Nidx = floor(idx*NumNormSamp/length(trialTime));
                    E.EventMarker{i , 1}(end - idx) = pr;  % First press
                    E.NormEventMarker{i , 1}(end - Nidx) = pr;
                end
                idx = floor(Fs/2);
                E.EventMarker{i , 1}(idx) = -1;       % stimulus apears
                Nidx = floor(idx*NumNormSamp/length(trialTime));
                E.NormEventMarker{i , 1}(Nidx) = -1;
            else
                trialTime  = 0:1/Fs:(eegTime/Fs) - 1/Fs;
                
                E.EventMarker{i , 1} = zeros(size( trialTime));
                E.NormEventMarker{i , 1} = zeros(NumNormSamp,1);
                idx = floor(Fs/2);
                E.EventMarker{i , 1}(idx) = -1;       % equivalent to stimulus apears
                Nidx = floor(idx*NumNormSamp/length(trialTime));
                E.NormEventMarker{i , 1}(Nidx) = -1;
                
                idx = Fs/2;
                Nidx = floor(idx*NumNormSamp/length(trialTime));
                E.EventMarker{i , 1}(end - idx) = 7;  % equivalent to last press
                E.NormEventMarker{i , 1}(end - Nidx) = 7;
            end
        end
        Dall.EventMarker = E.EventMarker;
        Dall.NormEventMarker = E.NormEventMarker;
        Dout = Dall;
    case 'CalcAveragePattern'
        
        %% Calculating the average pattern
        % should be performed on completed data - no need for EEG, but at
        % least the length of data sample in each trial should be included
        
        % single finger blocks that participant is instructed to go slow 
        slowSingleFing = [3 13 26 40];
        BNid_slow  = logical(sum(Dall.BN == slowSingleFing , 2));
        Dall.seqNumb(BNid_slow) = Dall.seqNumb(BNid_slow)*10;
        
        % single finger blocks that participant is instructed to go fast as she can
        fastSingleFing = [4 14 27 41];
        BNid_fast  = logical(sum(Dall.BN == fastSingleFing , 2));
        Dall.seqNumb(BNid_fast) = Dall.seqNumb(BNid_fast)*100;
        
        Dout.SN = unique(Dall.seqNumb);
        Dout.SN = Dout.SN(Dout.SN~=5); % exclude the stars
        
        clear EventMarker NormEventMarker EM NEM diffNEM
        for sn= 1:length(Dout.SN)
            D = getrow(Dall , ismember(Dall.seqNumb , Dout.SN(sn)) & ~Dall.isError);
            EventMarker{sn} = [];
            NormEventMarker{sn} = [];
            for k = 1:length(D.TN)
                events = [-1 , 1:D.seqlength(k)];
                if D.TN(k) >3
                    clear temp Ntemp
                    for jj = 1:length(events)
                        temp(1, jj) =  find(D.EventMarker{k} == events(jj));
                        Ntemp(1,jj) =  find(D.NormEventMarker{k} == events(jj));
                    end
                    EventMarker{sn} = [EventMarker{sn}; temp];
                    NormEventMarker{sn} = [NormEventMarker{sn} ;Ntemp];
                end
            end
            Dout.EM{sn, :} = floor(mean(EventMarker{sn}));
            Dout.NEM{sn, :} = floor(mean(NormEventMarker{sn}));
            Dout.diffNEM{sn, :} = diff([0 , Dout.NEM{sn} , NumNormSamp]);
        end
        
end
