function Dout  = seqeeg_parseEEG_PSD(what , Dall , subjNum, varargin)
%% reads the all channels-packed EEG data, uses the BlockInfo file to parse the EEG into single trials and ammend the bihavioral data structure
%% It also calculates the PSD on the whole block and then parses up the PSD inot trials. This is mainly to avoid any window effect on single trials
c = 1;
%% setup the defaults and deal with the varargin
subjname = {'P2' , 'P4' , 'P5'};
mainDir = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' subjname{subjNum} , '/Packed/'] ;
saveDir = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' subjname{subjNum} , '/'] ;
DownsampleRate = 10;
NumWarpSampFast = 150;
NumWarpSampSlow = 300;
TimeDelay = 0.5; % sec
FreqRange = [4 184];
numFreqBins = 45;
load([saveDir , 'ChanLabels.mat']);
Channels = [1:length(ChanLabels)];
%%  control for too short IPIs that the keys get accidentally pressed - > happens for Subj 1
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
        case {'HilbertOrWavelet'}
            % Hilbert or Wavelet time-freq decomposition , 'H' , 'W'
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        otherwise
            error(sprintf('Unknown option: %s',varargin{c}));
    end
end


%%



min_freq =  FreqRange(1);
max_freq = FreqRange(2);
frex = linspace(min_freq, max_freq,numFreqBins);
BandInfo.bandsLab = {'Delta <4Hz' , 'Theta 5-8Hz' , 'Alpha 9-16Hz' , 'Beta 17-36Hz' , 'L-Gamma 37-80Hz' , 'H-Gamma 80-100HZ' , 'HIGH 100-180HZ'};
BandInfo.bands = {[0 4], [5 8] [9 16] [17 36] [37 80] [80 100] [100 180]};
for b = 1:length(BandInfo.bands)
    BandInfo.bandid{b} = [find(frex>BandInfo.bands{b}(1) ,1, 'first') , find(frex<BandInfo.bands{b}(2) ,1, 'last')];
end
%%

cd(mainDir)
% Import the BLock Info
[~, ~, BlockInfo] = xlsread([mainDir , 'BlockInfo.xlsx'],'Sheet1');
BlockInfo = BlockInfo(2:end,:);
BlockInfo(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),BlockInfo)) = {''};
idx = cellfun(@ischar, BlockInfo);
BlockInfo(idx) = cellfun(@(x) string(x), BlockInfo(idx), 'UniformOutput', false);
Fs = 1024;
Fs_ds = Fs/DownsampleRate;
clearvars idx;
Dout = [];


[~, ~, BLockGroups] = xlsread([saveDir , 'BLockGroups.xlsx'],'Sheet1');
BLockGroups = BLockGroups(1:end,:);
BLockGroups(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),BLockGroups)) = {''};
idx = cellfun(@ischar, BLockGroups);
BLockGroups(idx) = cellfun(@(x) string(x), BLockGroups(idx), 'UniformOutput', false);
clearvars idx;

for bg = 1:length(BLockGroups)
    if ~isnumeric(BLockGroups{bg,1})
        blockGroups{bg} = str2num(char(BLockGroups{bg,1}));
    else
        blockGroups{bg} = BLockGroups{bg,1};
    end
end
blockGroupNames = BLockGroups(:,2);
fastBlock = horzcat(blockGroups{1} ,blockGroups{6} , blockGroups{7}, blockGroups{8}, blockGroups{9},...
    blockGroups{12}, blockGroups{16}, blockGroups{20}, blockGroups{23});
Dall.Fast = zeros(size(Dall.BN));
Dall.Fast(ismember(Dall.BN , fastBlock)) = 1;

%% chop up the EEG data into trials, and filter out the power line noise
switch what
    case 'ParseEEG-freqDecomp'
        
        %% preprocess EEG and filtering
        
        ChanLabels = ChanLabels(Channels);
        % Q : quality factor is the center frequency divided by the bandwidth.
        % Q = 35;
        % BW = Fo/(Fs/2);
        % [b,a] = iircomb(10,BW,'notch');
        Fo = 60;
        Fs = 1024;
        wo = Fo/(Fs/2);  bw = wo/35; % notch to eliminate 60 Hz
        [b,a] = iirnotch(wo,bw);
        
        wo = 2*Fo/(Fs/2);  bw = wo/60;% notch to eliminate the first harmonic
        [c,d] = iirnotch(wo,bw);
        % definitions, selections...
        
        %% multiple blocks are sotred in the same file. so avoid loading them up multiple times.
        fName1 =  BlockInfo{1,4};
        tn = 1;
        load(fName1);
        
        
        %%
        Events = [];
        for i = 1:size(BlockInfo , 1)
            clear Pall
            
            D = getrow(Dall , Dall.BN == i);
            fName = BlockInfo{i,4};
            if ~strcmp(fName , fName1)
                load(fName);
                % filter the power line noise out of the whole data and then chopp it up
                
                fName1 = fName;
            end
            %     extract the data for the block at hand
            BlockRang = [BlockInfo{i,2} : BlockInfo{i,3}];
            if ~ismember(BlockRang , -1)
                Beeg = getrow(Data , Channels);
                Beeg.values = Beeg.values(:,BlockRang);
                % get the indecies for starts ans ends of the trials
                marker = Beeg.values(find(strcmp(Beeg.label , 'TTL')) , :);
                marker = [0 diff(marker <-2*10^6)];
                for ch = 1:size(Beeg.values , 1)
                    B = Beeg.values(ch , :);
                    A = filter(b,a , B);
                    Beeg.values(ch , :) = A;
                end
                BlockPow.start_tr = find(marker == 1); % starts of trials
                % right now the end TTl pulse is being sent by the releas eof
                % the last finger, so is not aligned to the last press time. so
                % better define it this way for trials with presses
                BlockPow.end_tr = find(marker == -1);  % ends of trials
                for tn = 3:length(BlockPow.start_tr)
                    if ~isnan(D.AllPressTimes(tn , D.seqlength(tn)))
                        BlockPow.end_tr(tn) = BlockPow.start_tr(tn) + Fs*(D.AllPressTimes(tn , D.seqlength(tn))/1000)';
                    end
                end
                BlockPow.start_tr = floor(BlockPow.start_tr / DownsampleRate);
                BlockPow.end_tr = floor(BlockPow.end_tr / DownsampleRate);
                chanCatData = reshape(Beeg.values' , 1, numel(Beeg.values));
                % concatenate all the chanles to speed up analysis
                switch HilbertOrWavelet
                    case 'H'
                        [CatREG] = seqeeg_HilbertPSD(chanCatData , Fs , 'DownsampleRate' , DownsampleRate , 'PoworDec' , 'd','NumChans' , size(Beeg.values , 1)); % get the decomposition
                    case 'W'
                        [CatREG] = seqeeg_waveletPSD(chanCatData , Fs , 'DownsampleRate' , DownsampleRate , 'PoworDec' , 'd','NumChans' , size(Beeg.values , 1)); % get the decomposition
                end
                for ch = 1:size(Beeg.values , 1)
                    numsamp = length(downsample(Beeg.values(ch,:) , DownsampleRate));
                    % get the channel out of the concatenated PSD\
                    chanStart = 1+ (ch-1)*numsamp;
                    chanEnd   = ch * numsamp;
                    BlockPow.REG = CatREG(:,chanStart:chanEnd);
                    % normalize each trial to baseline : TimeDelay ms before the stim  onset
                    for tr = 1:length(BlockPow.start_tr)
                        X = nanmean(BlockPow.REG(:,BlockPow.start_tr(tr)-floor(Fs_ds*TimeDelay/2):BlockPow.start_tr(tr)-floor(Fs_ds*TimeDelay/10)) , 2);
                        eval(['D.decompBL{tr,1}(ch , :,:)  = X;']);
                        
                        X = BlockPow.REG(:,BlockPow.start_tr(tr) : BlockPow.end_tr(tr));
                        eval(['D.decompTR{tr,1}(ch , :,:)  = X;']);
                        
                        X = BlockPow.REG(:,BlockPow.start_tr(tr)-floor(Fs_ds*TimeDelay) : BlockPow.start_tr(tr)-1);
                        eval(['D.decompBefTR{tr,1}(ch , :,:)  = X;']);
                        
                        X = BlockPow.REG(:,BlockPow.end_tr(tr)+1 : BlockPow.end_tr(tr)+floor(Fs_ds*2*TimeDelay));
                        eval(['D.decompAftTR{tr,1}(ch , :,:)  = X;']);
                    end
                    disp(['PSD calculation for block ' , num2str(i) , ', channel ' , num2str(ch) , ' completed'])
                    clear REG
                end
                if length(BlockPow.start_tr)<length(D.TN)
                    missing = abs(length(BlockPow.start_tr)-length(D.TN));
                    missingTR = [length(D.TN)-(missing-1) : length(D.TN)];
                    for mtr = missingTR
                        eval(['D.decompBL{mtr,1}  = NaN']);
                        eval(['D.decompTR{mtr,1}  = NaN;']);
                        eval(['D.decompBefTR{mtr,1}= NaN;']);
                        eval(['D.decompAftTR{mtr,1} = NaN;']);
                    end
                end
                
                
                %% save the unbinned data in separate blocks for managability in size
                %             save the raw unbinned psd
                for tn = 1 :length(D.TN)
                    % prepare the individual trials to be saved as the
                    % fileds of a structre to make loading easier
                    eval(['Pall.DEC',num2str(tn),'.decompBL    = D.decompBL{', num2str(tn), '};']);
                    eval(['Pall.DEC',num2str(tn),'.decompTR    = D.decompTR{', num2str(tn), '};']);
                    eval(['Pall.DEC',num2str(tn),'.decompBefTR = D.decompBefTR{', num2str(tn), '};']);
                    eval(['Pall.DEC',num2str(tn),'.decompAftTR = D.decompAftTR{', num2str(tn), '};']);
                end
            else
                for tn = 1 :length(D.TN)
                    % prepare the individual trials to be saved as the
                    % fileds of a structre to make loading easier
                    eval(['Pall.DEC',num2str(tn),'.decompBL    = NaN;']);
                    eval(['Pall.DEC',num2str(tn),'.decompTR    = NaN;']);
                    eval(['Pall.DEC',num2str(tn),'.decompBefTR = NaN;']);
                    eval(['Pall.DEC',num2str(tn),'.decompAftTR = NaN;']);
                end
                
            end
            saveName = [saveDir,'Raw_Decomp_B',num2str(i) ,'.mat'];
            save(saveName,'-struct','Pall', '-v7.3');
            
        end
        
        %% Filter Visualization
        % h = fvtool(b,a);
        % h.Fs = Fs;
        % h.FrequencyRange='[-Fs/2, Fs/2)';
        % zplane(b,a)
        %
        % h = fvtool(c,d);
        % h.Fs = Fs;
        % h.FrequencyRange='[-Fs/2, Fs/2)';
        % zplane(c,d)
        % filt_A = Dall.EEG{tn}(10,:);
        % t = 0:(1/Fs):(length(A)/Fs) - (1/Fs);
        % figure('color' , 'white')
        % subplot(2,1,1)
        % periodogram(A(1,:),[],length(A(1,:)),Fs,'power')
        % subplot(2,1,2)
        % periodogram(filt_A(1,:),[],length(filt_A(1,:)),Fs,'power')
    case 'ParseEEG-calc_norm_PSD'
        %% preprocess EEG and filtering
        load([saveDir , 'ChanLabels.mat']);
        ChanLabels = ChanLabels(Channels);
        % Q : quality factor is the center frequency divided by the bandwidth.
        % Q = 35;
        % BW = Fo/(Fs/2);
        % [b,a] = iircomb(10,BW,'notch');
        Fo = 60;
        Fs = 1024;
        wo = Fo/(Fs/2);  bw = wo/35; % notch to eliminate 60 Hz
        [b,a] = iirnotch(wo,bw);
        Dall  = seqeeg_addEventMarker(Dall, subjNum, Fs_ds , 'addEvent' , 'NumWarpSampFast' , NumWarpSampFast, 'NumWarpSampSlow'  ,NumWarpSampSlow)';
        Events = Dall;
        save([saveDir , 'AllData_Events.mat'] , 'Events');
        wo = 2*Fo/(Fs/2);  bw = wo/60;% notch to eliminate the first harmonic
        [c,d] = iirnotch(wo,bw);
        % definitions, selections...
        
        %% multiple blocks are sotred in the same file. so avoid loading them up multiple times.
        fName1 =  BlockInfo{1,4};
        tn = 1;
        load(fName1);
        for i = 1:size(BlockInfo , 1)
            clear Pall
            
            D = getrow(Dall , Dall.BN == i);
            fName = BlockInfo{i,4};
            if ~strcmp(fName , fName1)
                load(fName);
                % filter the power line noise out of the whole data and then chopp it up
                
                fName1 = fName;
            end
            %     extract the data for the block at hand
            BlockRang = [BlockInfo{i,2} : BlockInfo{i,3}];
            if ~ismember(BlockRang , -1)
                Beeg = getrow(Data , Channels);
                Beeg.values = Beeg.values(:,BlockRang);
                % get the indecies for starts ans ends of the trials
                marker = Beeg.values(find(strcmp(Beeg.label , 'TTL')) , :);
                marker = [0 diff(marker <-2*10^6)];
                for ch = 1:size(Beeg.values , 1)
                    B = Beeg.values(ch , :);
                    A = filter(b,a , B);
                    Beeg.values(ch , :) = A;
                end
                BlockPow.start_tr = find(marker == 1); % starts of trials
                % right now the end TTl pulse is being sent by the releas eof
                % the last finger, so is not aligned to the last press time. so
                % better define it this way for trials with presses
                BlockPow.end_tr = find(marker == -1);  % ends of trials
                for tn = 3:length(BlockPow.start_tr)
                    if ~isnan(D.AllPressTimes(tn , D.seqlength(tn)))
                        BlockPow.end_tr(tn) = BlockPow.start_tr(tn) + Fs*(D.AllPressTimes(tn , D.seqlength(tn))/1000)';
                    end
                end
                BlockPow.start_tr = floor(BlockPow.start_tr / DownsampleRate);
                BlockPow.end_tr = floor(BlockPow.end_tr / DownsampleRate);
                
                chanCatData = reshape(Beeg.values' , 1, numel(Beeg.values));
                % concatenate all the chanles to speed up analysis
                switch HilbertOrWavelet
                    case 'H'
                        [CatREG] = seqeeg_HilbertPSD(chanCatData , Fs , 'DownsampleRate' , DownsampleRate , 'PoworDec' , 'p','NumChans' , size(Beeg.values , 1)); % get the decomposition
                    case 'W'
                        [CatREG] = seqeeg_waveletPSD(chanCatData , Fs , 'DownsampleRate' , DownsampleRate , 'PoworDec' , 'p','NumChans' , size(Beeg.values , 1)); % get the decomposition
                end
                
                for ch = 1:size(Beeg.values , 1)
                    numsamp = length(downsample(Beeg.values(ch,:) , DownsampleRate));
                    % get the channel out of the concatenated PSD\
                    chanStart = 1+ (ch-1)*numsamp;
                    chanEnd   = ch * numsamp;
                    BlockPow.REG = CatREG(:,chanStart:chanEnd);
                    clear baseline
                    for tr = 1:length(BlockPow.start_tr)
                        % consider the baseline to be the average of 250 ms before stim onset to 50 ms before stim onset
                        baseline(tr, :) = squeeze(nanmean(BlockPow.REG(:,BlockPow.start_tr(tr)-floor(Fs_ds*TimeDelay/2):BlockPow.start_tr(tr)-floor(Fs_ds*TimeDelay/10)) , 2));
                    end
%                      baseline = nanmean(baseline , 1)';
                    % normalize each trial to baseline : TimeDelay ms before the stim  onset
                    for tr = 1:length(BlockPow.start_tr)
                        X = squeeze(BlockPow.REG(:,BlockPow.start_tr(tr)-floor(Fs_ds*TimeDelay) : BlockPow.end_tr(tr)+floor(Fs_ds*2*TimeDelay)));
%                         X = (X - repmat( baseline(tr, :)' , 1,size(X,2)))./repmat( baseline(tr, :)' , 1,size(X,2));
                        eval(['D.PSD{tr,1}(ch , :,:)  = X;']);
                        eval(['D.BaseLine{tr,1}(ch , :,:)  = baseline(tr, :)'';']);
                    end
                    disp(['PSD calculation for block ' , num2str(i) , ', channel ' , num2str(ch) , ' completed'])
                    clear REG
                end
                % in case some trials at the end didnt get recorded on the marker
                if length(BlockPow.start_tr)<length(D.TN)
                    missing = abs(length(BlockPow.start_tr)-length(D.TN));
                    missingTR = [length(D.TN)-(missing-1) : length(D.TN)];
                    for mtr = missingTR
                        eval(['D.PSD{mtr,1}  = NaN;']);
                        eval(['D.BaseLine{mtr,1}  = NaN;']);
                    end
                end
            else
                missing = length(D.TN);
                missingTR = [length(D.TN)-(missing-1) : length(D.TN)];
                for mtr = missingTR
                    eval(['D.PSD{mtr,1}  = NaN;']);
                    eval(['D.BaseLine{mtr,1}  = NaN;']);
                end
            end
            
            %% find the event markers and normalize the power to TimeDelay ms before the stimulus came on - or press
            % complete the structure with behavior again

            %% save the unbinned data in separate blocks for managability in size
            %             save the raw unbinned psd
            for tn = 1 :length(D.TN)
                % prepare the individual trials to be saved as the
                % fileds of a structre to make loading easier
                eval(['Pall.PSD',num2str(tn),' = D.PSD{', num2str(tn), '};']);
                eval(['Pall.BaseLine',num2str(tn),' = D.BaseLine{', num2str(tn), '};']);
            end
            
            saveName = [saveDir,'Raw_PSD_B',num2str(i) ,'.mat'];
            save(saveName,'-struct','Pall', '-v7.3');
            
            % then band-average and stack up
            for tn = 1 :length(D.TN)
                if ~isnan(D.PSD{tn})
                    temp = D.PSD{tn};
                    tempbl = D.BaseLine{tn};
                    for b =1:length(BandInfo.bandid)
                        P.PSD{tn,1}(:,b, :) =  nanmean(temp(:,BandInfo.bandid{b}(1) : BandInfo.bandid{b}(end),:) , 2);
                        P.BL{tn,1}(:,b) =  nanmean(tempbl(:,BandInfo.bandid{b}(1) : BandInfo.bandid{b}(end),:) , 2);
                    end
                else
                    P.PSD{tn,1} = NaN;
                    P.BL{tn,1} = NaN;
                end
            end
            Dout = addstruct(Dout , P);
            if exist('BlockPow')
                saveName = [saveDir,'BlockPow' ,num2str(i) ,'.mat'];
                save(saveName,'-struct','BlockPow', '-v7.3')
            end
            clear D P baseline BlockPow
        end
        load([saveDir , 'AllData_Behav.mat'])
        Dall.PSD = Dout.PSD;
        Pall = Dall;
        clear Dall
        saveName = [saveDir,'AllData_PSD_StimNorm.mat'];
        save(saveName , 'Pall' , '-v7.3');        
    case 'TimeWarpPSD_Raw_Binned'
        % this case loads up it's own input so the Dall structure can be left blac
        load([saveDir,'AllData_Events.mat']);
        Fs = 1024;
        Fs_ds = floor(Fs/DownsampleRate);
        % find average event markers
        Events  = seqeeg_addEventMarker(Dall,subjNum, Fs_ds, 'addEvent' , 'NumWarpSampFast' , NumWarpSampFast, 'NumWarpSampSlow'  ,NumWarpSampSlow);
        E  = seqeeg_addEventMarker(Events, subjNum, Fs_ds , 'CalcAveragePattern' , 'NumWarpSampFast' , NumWarpSampFast, 'NumWarpSampSlow'  ,NumWarpSampSlow)';
        save([saveDir , 'AllData_AvgMarker.mat'] , 'E');
        BN = unique(Events.BN);
        
        Pall_binned = [];
        for bn = 1:length(BN)
            bg = 1;
            mem = 0;
            while mem == 0
                if ismember(BN(bn) , E.blockGroups{bg})
                    mem = 1;
                    BG = bg;
                end
                bg = bg+1;
            end
            
            clear D P
            D = getrow(Events , Events.BN == BN(bn));
            Bname = [saveDir,'Raw_PSD_B',num2str(bn) ,'.mat'];
            P = load(Bname);
            for tn = 1 :length(D.TN)
                % prepare the individual trials to be saved as the
                % fileds of a structre to make loading easier
                psd = eval(['P.PSD', num2str(tn)]);
                bl = repmat(eval(['P.BaseLine', num2str(tn)]) ,1,1, size(psd , 3)); 
                eval(['D.PSD{',num2str(tn),',1} = 100*(psd - bl)./bl ;' ] ); % percent change from baseline
            end
            % this loads up the data strutre where the raw PSD for that block is already stored
            % 'ParseEEG-calcPSD'
            % Use the average patterns of block groups to warp them
            nPSD = D.PSD;
            
            % set the sequence length for "*******" trilas to one
            D.seqlength(D.seqNumb==5) = 1;
            D.NumWarpSamp = NumWarpSampSlow* ones(size(D.TN));
            D.NumWarpSamp(logical(D.Fast)) = NumWarpSampFast;
            for tn = 1:length(D.TN)
                if length(find(D.NormEventMarker{tn})) == D.seqlength(tn) +1 & ~D.isError(tn) & ...
                        length(D.EventMarker{tn})>D.NumWarpSamp(tn) & ~isnan(nPSD{tn})
                    % check which block group this block falls into
                    
                    A = nPSD{tn};
                    % the normalized event markers in trial tn
                    idx = [0 find(D.NormEventMarker{tn}) D.NumWarpSamp(tn)];
                    
                    % find the row number corresponding to the seqNumb
                    sn = find(E.SN{BG}==D.seqNumb(tn));
                    % average normalized time stamps for the sn , BG
                    if tn>2
                        idn = [0 E.NEM{BG}(sn , ~isnan(E.NEM{BG}(sn ,:))) D.NumWarpSamp(tn)];
                        diffNEM = diff(idn);
                    else
                        % for "* * * * * * * *" trials, the average pattern is the same as the trial since there is no variability
                        idn = idx;
                        diffNEM = diff(idn);
                    end
                    
                    for evnt = 2:length(idn)
                        % make sure that each
                        idd   = linspace(1 , [idx(evnt)+1 - (idx(evnt-1) + 1)] , diffNEM(evnt-1));
                        for ch = 1:size(A,1)
                            for fi=1:size(A,2)
                                %                                 idd = floor(linspace(1, length([1:idx(e) - idx(e-1)]) , length(idd)));
                                D.PSD_stim{tn ,1}(ch,fi,idn(evnt-1)+1:idn(evnt)) = interp1([1:idx(evnt) - idx(evnt-1)] , squeeze(A(ch,fi,idx(evnt-1)+1:idx(evnt))) , idd);
                            end
                        end
                    end
                else
                    D.PSD_stim{tn ,1} = nan(size(A,1) , length(frex) , D.NumWarpSamp(tn));
                end
            end
            for tn = 1:length(D.PSD_stim)
                trialName = ['Pall.PSD',num2str(tn)];
                eval([trialName, ' = D.PSD_stim{' , num2str(tn) , '};'])
            end
            
            saveName = [saveDir,'warped_PSD_B',num2str(bn) ,'.mat'];
            save(saveName,'-struct','Pall', '-v7.3');
            % stace the binned PSDs
            
            
            D = getrow(Events , Events.BN == BN(bn));
            
            for tn = 1 :length(D.TN)
                PSD = eval(['Pall.PSD' , num2str(tn)]);
                
                for b = 1:length(BandInfo.bands)
                    D.PSD{tn,1}(:,b, :) =  nanmean(PSD(:,BandInfo.bandid{b}(1) : BandInfo.bandid{b}(2),:) , 2);
                end
                
            end
            Pall_binned = addstruct(Pall_binned, D);
            clear D
            disp(['PSD warping for block ' , num2str(bn) , ' complete'])
        end
        Pall = Pall_binned;
        saveName = [saveDir,'AllData_PSD_Warped.mat'];
        save(saveName,'Pall', '-v7.3');
        Dout = Pall;
    case 'TimeWarpPSD_Raw_Binned_seqType'
        % this case loads up it's own input so the Dall structure can be left blac
        %% the goal here is to get average time pattern for general sequence types regardless of finger
        % so all the slow single fingers, fast singel finger, random, structured, triplets, quadruples
        % so we will change the SeqNumb and average the average times of the SeqNumbs
        % within the same type
        
        
        load([saveDir,'AllData_Events.mat']);
        
        
        Fs = 1024;
        Fs_ds = floor(Fs/DownsampleRate);
        % find average event markers
        Events  = seqeeg_addEventMarker(Dall,subjNum, Fs_ds, 'addEvent' , 'NumWarpSampFast' , NumWarpSampFast, 'NumWarpSampSlow'  ,NumWarpSampSlow);
        E  = seqeeg_addEventMarker(Events, subjNum, Fs_ds , 'CalcAveragePattern_seqType' , 'NumWarpSampFast' , NumWarpSampFast, 'NumWarpSampSlow'  ,NumWarpSampSlow)';
        save([saveDir , 'AllData_AvgMarker_SeqType.mat'] , 'E');
        
        % Define sequence numbers and their transformations:
        SeqTrans = [5 11 22 33 44 55 0 1 2 3 4 103 104 203 204;...
            100 10 10 10 10 10 20 30 30 30 30 40 50 40 50];
        for sn = 1:length(SeqTrans)
            id  = Events.seqNumb == SeqTrans(1 , sn);
            Events.seqNumb(id) = SeqTrans(2 , sn);
        end
        
        BN = unique(Events.BN);
        Pall_binned = [];
        for bn = 1:length(BN)
            bg = 1;
            mem = 0;
            while mem == 0
                if ismember(BN(bn) , E.blockGroups{bg})
                    mem = 1;
                    BG = bg;
                end
                bg = bg+1;
            end
            clear D P
            D = getrow(Events , Events.BN == BN(bn));
            Bname = [saveDir,'Raw_PSD_B',num2str(bn) ,'.mat'];
            P = load(Bname);
            for tn = 1 :length(D.TN)
                % prepare the individual trials to be saved as the
                % fileds of a structre to make loading easier
                psd = eval(['P.PSD', num2str(tn)]);
                bl = repmat(eval(['P.BaseLine', num2str(tn)]) ,1,1, size(psd , 3)); 
                eval(['D.PSD{',num2str(tn),',1} = 100*(psd - bl)./bl ;' ] ); % percent change from baseline
            end
            % this loads up the data strutre where the raw PSD for that block is already stored
            % 'ParseEEG-calcPSD'
            % Use the average patterns of block groups to warp them
            nPSD = D.PSD;
            
            % set the sequence length for "*******" trilas to one
            D.seqlength(D.seqNumb==100) = 1;
            D.NumWarpSamp = NumWarpSampSlow* ones(size(D.TN));
            D.NumWarpSamp(logical(D.Fast)) = NumWarpSampFast;
            for tn = 1:length(D.TN)
                if length(find(D.NormEventMarker{tn})) == D.seqlength(tn) +1 & ~D.isError(tn) & ...
                        length(D.EventMarker{tn})>D.NumWarpSamp(tn) & ~isnan(nPSD{tn})
                    % check which block group this block falls into
                    
                    A = nPSD{tn};
                    % the normalized event markers in trial tn
                    idx = [0 find(D.NormEventMarker{tn}) D.NumWarpSamp(tn)];
                    
                    % find the row number corresponding to the seqNumb
                    sn = find(E.SN{BG}==D.seqNumb(tn));
                    % average normalized time stamps for the sn , BG
                    if tn>2
                        idn = [0 E.NEM{BG}(sn , ~isnan(E.NEM{BG}(sn ,:))) D.NumWarpSamp(tn)];
                        diffNEM = diff(idn);
                    else
                        % for "* * * * * * * *" trials, the average pattern is the same as the trial since there is no variability
                        idn = idx;
                        diffNEM = diff(idn);
                    end
                    
                    for evnt = 2:length(idn)
                        % make sure that each
                        idd   = linspace(1 , [idx(evnt)+1 - (idx(evnt-1) + 1)] , diffNEM(evnt-1));
                        for ch = 1:size(A,1)
                            for fi=1:size(A,2)
                                %                                 idd = floor(linspace(1, length([1:idx(e) - idx(e-1)]) , length(idd)));
                                D.PSD_stim{tn ,1}(ch,fi,idn(evnt-1)+1:idn(evnt)) = interp1([1:idx(evnt) - idx(evnt-1)] , squeeze(A(ch,fi,idx(evnt-1)+1:idx(evnt))) , idd);
                            end
                        end
                    end
                else
                    D.PSD_stim{tn ,1} = nan(size(A,1) , length(frex) , D.NumWarpSamp(tn));
                end
            end
            for tn = 1:length(D.PSD_stim)
                trialName = ['Pall.PSD',num2str(tn)];
                eval([trialName, ' = D.PSD_stim{' , num2str(tn) , '};'])
            end
            
            saveName = [saveDir,'warped_PSD_B_SeqType',num2str(bn) ,'.mat'];
            save(saveName,'-struct','Pall', '-v7.3');
            % stace the binned PSDs
            
            
            D = getrow(Events , Events.BN == BN(bn));
            
            for tn = 1 :length(D.TN)
                PSD = eval(['Pall.PSD' , num2str(tn)]);
                
                for b = 1:length(BandInfo.bands)
                    D.PSD{tn,1}(:,b, :) =  nanmean(PSD(:,BandInfo.bandid{b}(1) : BandInfo.bandid{b}(2),:) , 2);
                end
                
            end
            Pall_binned = addstruct(Pall_binned, D);
            clear D
            disp(['PSD warping for block ' , num2str(bn) , ' complete'])
        end
        
        Pall = Pall_binned;
        saveName = [saveDir,'AllData_PSD_Warped_SeqType.mat'];
        save(saveName,'Pall', '-v7.3');
        Dout = Pall;
    case 'AlignEvents_SeqType'
        load([saveDir,'AllData_Events.mat']);
        Fs = 1024;
        Fs_ds = floor(Fs/DownsampleRate);
        % find average event markers
        Events  = seqeeg_addEventMarker(Dall,subjNum, Fs_ds, 'addEvent' , 'NumWarpSampFast' , NumWarpSampFast, 'NumWarpSampSlow'  ,NumWarpSampSlow);
        E  = seqeeg_addEventMarker(Events, subjNum, Fs_ds , 'CalcAveragePattern' , 'NumWarpSampFast' , NumWarpSampFast, 'NumWarpSampSlow'  ,NumWarpSampSlow);
        BN = unique(Events.BN);
        
        Pall_binned = [];
        for bn = 1:length(BN)
            bg = 1;
            mem = 0;
            while mem == 0
                if ismember(BN(bn) , E.blockGroups{bg})
                    mem = 1;
                    BG = bg;
                end
                bg = bg+1;
            end
            
            clear D P
            D = getrow(Events , Events.BN == BN(bn));
            Bname = [saveDir,'Raw_PSD_B',num2str(bn) ,'.mat'];

        
            % set the sequence length for "*******" trilas to one
            D.seqlength(D.seqNumb==100) = 1;
            D.NumWarpSamp = NumWarpSampSlow* ones(size(D.TN));
            D.NumWarpSamp(logical(D.Fast)) = NumWarpSampFast;
            numEvents = max(D.seqlength) + 1;
            for tn = 1:2 % Null Trails
                for evnt = 1:length(numEvents )
                    PSD_Aligned{tn,evnt} = [];
                    PSD_betweenEvnt{tn,evnt} = [];
                end
            end
            for tn = 3:length(D.TN)
                A = load(Bname , ['PSD',num2str(tn)] , ['BaseLine',num2str(tn)]);
                psd = eval(['A.',['PSD',num2str(tn)],';']);
                bl = eval(['A.',['BaseLine',num2str(tn)],';']);
                psd = 100*(psd - bl)./bl ; % percent change
                if  ~D.isError(tn) & sum(sum(~isnan(D.EventMarker{tn}))) & sum(sum(~isnan(psd)))
                    % the normalized event markers in trial tn
                    idx = find(D.EventMarker{tn});
                    % before the first press
                    PSD_betweenEvnt{tn,1} = psd(:,:,idx(1)-10:idx(1)-5); % prior to the first press
                    for evnt = 1:length(idx)
                        % around press time
                        PSD_Aligned{tn,evnt} = psd(:,:,idx(evnt)-5:idx(evnt)+2);
                        % between presses
                        if evnt<length(idx)
                            pressmidPoint = floor(idx(evnt+1) - idx(evnt)/2);
                            PSD_betweenEvnt{tn,evnt+1} = psd(:,:,pressmidPoint-3:pressmidPoint+3);
                        end
                    end
                else
                    PSD_Aligned(tn,:) = PSD_Aligned(1,:);
                    PSD_betweenEvnt(tn,:) = PSD_betweenEvnt(1,:);
                end
            end
            for tn = 1:length(PSD_Aligned)
                trialName = ['Pall.PSD_Aligned',num2str(tn)];
                eval([trialName, ' = PSD_Aligned(' , num2str(tn) , ',:);'])
                trialName = ['Pall.PSD_betweenEvnt',num2str(tn)];
                eval([trialName, ' = PSD_betweenEvnt(' , num2str(tn) , ',:);'])
            end
            saveName = [saveDir,'EventAligned_PSD_B',num2str(bn) ,'.mat'];
            save(saveName,'-struct','Pall', '-v7.3');
            % stace the binned PSDs
            Dout = Pall;
            clear PSD_Aligned Pall            
            disp(['PSD Alignment for block ' , num2str(bn) , ' complete'])
        end        
    case 'AlignEvents_SeqType_Warped'
        load([saveDir,'AllData_Events.mat']);
        Fs = 1024;
        Fs_ds = floor(Fs/DownsampleRate);
        % find average event markers
        Events  = seqeeg_addEventMarker(Dall,subjNum, Fs_ds, 'addEvent' , 'NumWarpSampFast' , NumWarpSampFast, 'NumWarpSampSlow'  ,NumWarpSampSlow);
        E  = seqeeg_addEventMarker(Events, subjNum, Fs_ds , 'CalcAveragePattern_seqType' , 'NumWarpSampFast' , NumWarpSampFast, 'NumWarpSampSlow'  ,NumWarpSampSlow)';
        save([saveDir , 'AllData_AvgMarker_SeqType.mat'], 'E');
        BN = unique(Events.BN);
        
        
        SeqTrans = [5 11 22 33 44 55 0 1 2 3 4 103 104 203 204;...
            100 10 10 10 10 10 20 30 30 30 30 40 50 40 50];
        for sn = 1:length(SeqTrans)
            id  = Events.seqNumb == SeqTrans(1 , sn);
            Events.seqNumb(id) = SeqTrans(2 , sn);
        end
        
        BN = unique(Events.BN);
       
        
        Pall_binned = [];
        for bn = 1:length(BN)
            bg = 1;
            mem = 0;
            while mem == 0
                if ismember(BN(bn) , E.blockGroups{bg})
                    mem = 1;
                    BG = bg;
                end
                bg = bg+1;
            end
            E1 = getrow(E , BG);
            clear D P
            D = getrow(Events , Events.BN == BN(bn));
            Bname = [saveDir,'warped_PSD_B',num2str(bn) ,'.mat'];

        
            % set the sequence length for "*******" trilas to one
            D.seqlength(D.seqNumb==5) = 1;
            D.NumWarpSamp = NumWarpSampSlow* ones(size(D.TN));
            D.NumWarpSamp(logical(D.Fast)) = NumWarpSampFast;
            numEvents = max(D.seqlength) + 1;
            for tn = 1:2 % Null Trails
                for evnt = 1:numEvents
                    PSD_Aligned{tn,evnt} = [];
                end
            end
            for tn = 3:length(D.TN)
                sn = find(E1.SN{1} == D.seqNumb(tn));
                trialName = ['PSD',num2str(tn)];
                A = load(Bname , trialName);
                eval(['A = A.' , trialName , ';']);
                if   ~D.isError(tn) & sum(sum(~isnan(D.EventMarker{tn}))) & sum(sum(~isnan(A)))
                    % check which block group this block falls into
                    % the normalized event markers in trial tn
                    idx = E1.NEM{1}(sn , 1:D.seqlength(tn)+1);
                    for evnt = 1:length(idx)
                        PSD_Aligned{tn,evnt} = A(:,:,idx(evnt)-5:idx(evnt)+2);
                    end
                else
                    PSD_Aligned(tn,:) = PSD_Aligned(1,:);
                end
            end
            for tn = 1:length(PSD_Aligned)
                trialName = ['Pall.PSD_Aligned',num2str(tn)];
                eval([trialName, ' = PSD_Aligned(' , num2str(tn) , ',:);'])
            end
            saveName = [saveDir,'EventAligned_WarpedPSD_B',num2str(bn) ,'.mat'];
            save(saveName,'-struct','Pall', '-v7.3');
            % stace the binned PSDs
            Dout = Pall;
            clear PSD_Aligned Pall
            disp(['PSD Alignment for block ' , num2str(bn) , ' complete'])
        end
end

