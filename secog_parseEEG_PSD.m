function Dout  = secog_parseEEG_PSD(what , Dall , subjNum, varargin)
%% reads the all channels-packed EEG data, uses the BlockInfo file to parse the EEG into single trials and ammend the bihavioral data structure
%% It also calculates the PSD on the whole block and then parses up the PSD inot trials. This is mainly to avoid any window effect on single trials
c = 1;
%% setup the defaults and deal with the varargin
DownsampleRate = 10;
NormType = 'stim';
NumWarpSampFast = 200;
NumWarpSampSlow = 500;
TimeDelay = 0.5; % sec
FreqRange = [2 180];
numFreqBins = 90;
Channels = [1:129];

while(c<=length(varargin))
    switch(varargin{c})
        
        case {'NormType'}
            % Type of PSD normalization
            % Default : 'stim' normalize power to TimeDelay ms before the stimulation comes on
            % could be 'press' normalizez to TimeDelay ms before the stim comes on + the time before the first press
            % could be 'none' which just transforms PSD into logarithmic scale
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
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
        otherwise
            error(sprintf('Unknown option: %s',varargin{c}));
    end
end



%% HousKeeping : load up data specify which block and trial are fast and which are slow
blockGroups = {[1 2] , [3 13 26 40] , [4 14 27 41] , [5:7] , [9:11] , [8 12] , [15:17] , [19:21] , [23:25],...
            [18 22] , [28:30] , [32:34] , [36:38], [31 35 39],[42:44]}';        
blockGroupNames = {'SingleFingNat' , 'SingleFingSlow' , 'SingleFingFast' , 'Intermixed1' , 'Intermixed2' , 'ChunkDay1' , 'Intermixed3' , 'Intermixed4' , 'Intermixed5',...
            'ChunkDay2' , 'Intermixed6' , 'Intermixed7' , 'Intermixed8', 'ChunkDay3', 'Intermixed9'}';
Dall.Fast = zeros(size(Dall.TN));
fastBlock = horzcat(blockGroups{1} , blockGroups{3} , blockGroups{6}, blockGroups{10},blockGroups{14});
Dall.Fast(ismember(Dall.BN , fastBlock)) = 1;
min_freq =  FreqRange(1);
max_freq = FreqRange(2);
frex = linspace(min_freq, max_freq,numFreqBins);
BandInfo.bandsLab = {'Delta <4Hz' , 'Theta 5-8Hz' , 'Alpha 9-16Hz' , 'Beta 17-36Hz' , 'L-Gamma 37-70Hz' , 'H-Gamma 70-130' , 'NoBandLand 130-180'};
BandInfo.bands = {[0 4], [5 8] [9 16] [17 36] [37 70] [70 130] [130 180]};
for b = 1:length(BandInfo.bands)
    BandInfo.bandid{b} = [find(frex>BandInfo.bands{b}(1) ,1, 'first') , find(frex<BandInfo.bands{b}(2) ,1, 'last')];
end


subjname = {'P2'};
mainDir = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' subjname{subjNum} , '/Packed/'] ;
saveDir = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' subjname{subjNum} , '/'] ;
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
%% chop up the EEG data into trials, and filter out the power line noise
switch what
    case 'ParseEEG-freqDecomp'
       
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
            E = getrow(Dall , Dall.BN == i);
            fName = BlockInfo{i,4};
            if ~strcmp(fName , fName1)
                load(fName);
                % filter the power line noise out of the whole data and then chopp it up
                
                fName1 = fName;
            end
            %     extract the data for the block at hand
            BlockRang = [BlockInfo{i,2} : BlockInfo{i,3}];
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
            start_tr = find(marker == 1); % starts of trials
            % right now the end TTl pulse is being sent by the releas eof
            % the last finger, so is not aligned to the last press time. so
            % better define it this way for trials with presses
            end_tr = find(marker == -1);  % ends of trials
            for tn = 3:length(start_tr)
                if ~isnan(D.AllPressTimes(tn , D.seqlength(tn)))
                    end_tr(tn) = start_tr(tn) + Fs*(D.AllPressTimes(tn , D.seqlength(tn))/1000)';
                end
            end
            start_tr = floor(start_tr / DownsampleRate);
            end_tr = floor(end_tr / DownsampleRate);
            for ch = 1:size(Beeg.values , 1)
                Chname{ch} = ['RawEEGpower',num2str(ch)];
                [REG, BandInfo] = secog_waveletPSD(Beeg.values(ch , :) , Fs , 'DownsampleRate' , DownsampleRate);
                % normalize each trial to baseline : TimeDelay ms before the stim  onset
                for tr = 1:length(start_tr)
                    X = nanmean(REG(:,start_tr(tr)-floor(Fs_ds*TimeDelay):start_tr(tr)) , 2);
                    eval(['D.decompBL{tr,1}(ch , :,:)  = X;']);
                    
                    X = REG(:,start_tr(tr) : end_tr(tr));
                    eval(['D.decompTR{tr,1}(ch , :,:)  = X;']);
                    
                    X = REG(:,start_tr(tr)-floor(Fs_ds*TimeDelay) : start_tr(tr)-1);
                    eval(['D.decompBefTR{tr,1}(ch , :,:)  = X;']);
                    
                    X = REG(:,end_tr(tr)+1 : end_tr(tr)+floor(Fs_ds*2*TimeDelay));
                    eval(['D.decompAftTR{tr,1}(ch , :,:)  = X;']);
                end
                disp(['PSD calculation for block ' , num2str(i) , ', channel ' , num2str(ch) , ' completed'])
                clear REG
            end
            

            %% save the unbinned data in separate blocks for managability in size
%             save the raw unbinned psd
            for tn = 1 :length(start_tr)
                % prepare the individual trials to be saved as the
                % fileds of a structre to make loading easier
                eval(['Pall.DEC',num2str(tn),'.decompBL    = D.decompBL{', num2str(tn), '};']);
                eval(['Pall.DEC',num2str(tn),'.decompTR    = D.decompTR{', num2str(tn), '};']);
                eval(['Pall.DEC',num2str(tn),'.decompBefTR = D.decompBefTR{', num2str(tn), '};']);
                eval(['Pall.DEC',num2str(tn),'.decompAftTR = D.decompAftTR{', num2str(tn), '};']);
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
            E = getrow(Dall , Dall.BN == i);
            fName = BlockInfo{i,4};
            if ~strcmp(fName , fName1)
                load(fName);
                % filter the power line noise out of the whole data and then chopp it up
                
                fName1 = fName;
            end
            %     extract the data for the block at hand
            BlockRang = [BlockInfo{i,2} : BlockInfo{i,3}];
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
            start_tr = find(marker == 1); % starts of trials
            % right now the end TTl pulse is being sent by the releas eof
            % the last finger, so is not aligned to the last press time. so
            % better define it this way for trials with presses
            end_tr = find(marker == -1);  % ends of trials
            for tn = 3:length(start_tr)
                if ~isnan(D.AllPressTimes(tn , D.seqlength(tn)))
                    end_tr(tn) = start_tr(tn) + Fs*(D.AllPressTimes(tn , D.seqlength(tn))/1000)';
                end
            end
            start_tr = floor(start_tr / DownsampleRate);
            end_tr = floor(end_tr / DownsampleRate);
            for ch = 1:size(Beeg.values , 1)
                Chname{ch} = ['RawEEGpower',num2str(ch)];
                [REG, BandInfo] = secog_waveletPSD(Beeg.values(ch , :) , Fs , 'DownsampleRate' , DownsampleRate);
                REG = 10*log10(abs(REG));
                % normalize each trial to baseline : TimeDelay ms before the stim  onset
                for tr = 1:length(start_tr)
                    baseline = nanmean(REG(:,start_tr(tr)-floor(Fs_ds*TimeDelay):start_tr(tr)) , 2);
                    X = REG(:,start_tr(tr)-floor(Fs_ds*TimeDelay) : end_tr(tr)+floor(Fs_ds*2*TimeDelay));
                    X = (X - repmat(baseline , 1,size(X,2)))./repmat(baseline , 1,size(X,2));
                    statement1 = ['D.PSD{tr,1}(ch , :,:)  = X;'];
                    eval(statement1);
                end
                disp(['PSD calculation for block ' , num2str(i) , ', channel ' , num2str(ch) , ' completed'])
                clear REG
            end
            
          
            
            %% find the event markers and normalize the power to TimeDelay ms before the stimulus came on - or press
            % complete the structure with behavior again
            D1 = getrow(Dall , Dall.BN == i);
            D1.Pow_Norm_stim = D.PSD;
            D = D1; clear D1
            D  = secog_addEventMarker(D, subjNum, Fs_ds , 'addEvent' , 'NumWarpSampFast' , NumWarpSampFast, 'NumWarpSampSlow'  ,NumWarpSampSlow)';
            E.RealTime = D.EventMarker;
            E.NormTime = D.NormEventMarker;
            %% save the unbinned data in separate blocks for managability in size
%             save the raw unbinned psd
            for tn = 1 :length(start_tr)
                % prepare the individual trials to be saved as the
                % fileds of a structre to make loading easier
                eval(['Pall.PSD',num2str(tn),' = D.Pow_Norm_stim{', num2str(tn), '};']);
            end
            saveName = [saveDir,'Raw_PSD_B',num2str(i) ,'.mat'];
            save(saveName,'-struct','Pall', '-v7.3');

            % then band-average and stack up
            for tn = 1 :length(start_tr)
                temp = D.Pow_Norm_stim{tn};
                for b =1:length(BandInfo.bandid)
                    P.Pow_Norm_stim{tn,1}(:,b, :) =  nanmean(temp(:,BandInfo.bandid{b}(1) : BandInfo.bandid{b}(2),:) , 2);
                end
            end
            Dout = addstruct(Dout , P);
            Events  = addstruct(Events , E);
            clear D P
        end
        load([saveDir , 'AllData_Behav.mat'])
        Dall.Pow_Norm_stim = Dout.Pow_Norm_stim;
        Pall = Dall;
        clear Dall
        saveName = [saveDir,'AllData_PSD_StimNorm.mat'];
        save(saveName , 'Pall' , '-v7.3');
        
        saveName = [saveDir,'AllData_Events.mat'];
        save(saveName , 'Events' , '-v7.3');
        

    case 'NormalizePSD'
        % obtain the time normalized event stamps for Block bins
        % load up the data events to get the everage pattern
       
        Fs = 1024;
        Fs_ds = floor(Fs/DownsampleRate);
%         Dall  = secog_parseEEG_PSD('ParseEEG-calcPSD' , Dall, subjNum);
        D  = secog_addEventMarker(Dall, subjNum, Fs_ds , 'addEvent' , 'NumWarpSampFast' , NumWarpSampFast, 'NumWarpSampSlow'  ,NumWarpSampSlow)';
        
        % this takes the data strutre where the raw PSD is already stored (so the outpur of the prevous case)
        % and normalizes the PSD for each Block to the average PSD of the firs two trails with the stars.
        for tn = 1:length(D.TN)
            A = D.PSD{tn};
            if D.TN(tn)>2
                for ch = 1:size(A,1)
                    for fi=1:size(A,2)
                        switch NormType
                            case 'stim'
                                Dout.Pow_Norm_stim{tn,1}(ch,fi,:) = squeeze(10*log10(A(ch,fi,:)./mean(A(ch,fi,1:find(D.EventMarker{tn ,1} == -1)))));
                            case 'press'
                                Dout.Pow_Norm_pres{tn,1}(ch,fi,:)      = squeeze(10*log10(A(ch,fi,:)./mean(A(ch,fi,1:find(D.EventMarker{tn ,1} == 1)))));
                            case 'none'
                                Dout.Pow_non_norm{tn,1}(ch,fi,:)       = squeeze(10*log10(A(ch,fi,:)));
                        end
                    end
                end
            else
                for ch = 1:size(A,1)
                    for fi=1:size(A,2)
                        switch NormType
                            case {'stim' , 'press'}
                                Dout.Pow_Norm_stim{tn,1}(ch,fi,:) = squeeze(10*log10(A(ch,fi,:)./mean(A(ch,fi,1:find(D.EventMarker{tn ,1} == -1)))));
                            case 'none'
                                Dout.Pow_non_norm{tn,1}(ch,fi,:)       = squeeze(10*log10(A(ch,fi,:)));
                        end
                    end
                end
            end
            
        end
        
       
%%        
    case 'TimeWarpPSD_Raw_Binned'
        % this case loads up it's own input so the Dall structure can be left blac
        load([saveDir,'AllData_Events.mat']);
        Fs = 1024;
        Fs_ds = floor(Fs/DownsampleRate);
        % find average event markers
        E  = secog_addEventMarker(Events, subjNum, Fs_ds , 'CalcAveragePattern' , 'NumWarpSampFast' , NumWarpSampFast, 'NumWarpSampSlow'  ,NumWarpSampSlow)';
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
                eval(['D.Pow_Norm_stim{',num2str(tn),'} = P.PSD', num2str(tn), ';']);
            end
            % this loads up the data strutre where the raw PSD for that block is already stored 
            % 'ParseEEG-calcPSD'
            % Use the average patterns of block groups to warp them
            nPSD = D.Pow_Norm_stim;
            
            % set the sequence length for "*******" trilas to one
            D.seqlength(isnan(D.seqlength)) = 1;
            D.NumWarpSamp = NumWarpSampSlow* ones(size(D.TN));
            D.NumWarpSamp(logical(D.Fast)) = NumWarpSampFast;
            for tn = 1:length(D.TN)
                if length(find(D.NormTime{tn})) == D.seqlength(tn) +1 & ~D.isError(tn) & length(D.RealTime{tn})>D.NumWarpSamp(tn)
                    % check which block group this block falls into
                    
                    A = nPSD{tn};
                    % the normalized event markers in trial tn
                    idx = [0 find(D.NormTime{tn}) D.NumWarpSamp(tn)];
                    
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
                    
                    for e = 2:length(idn)
                        % make sure that each
                        idd   = linspace(1 , [idx(e)+1 - (idx(e-1) + 1)] , diffNEM(e-1));
                        for ch = 1:size(A,1)
                            for fi=1:size(A,2)
%                                 idd = floor(linspace(1, length([1:idx(e) - idx(e-1)]) , length(idd)));
                                D.PSD_stim{tn ,1}(ch,fi,idn(e-1)+1:idn(e)) = interp1([1:idx(e) - idx(e-1)] , squeeze(A(ch,fi,idx(e-1)+1:idx(e))) , idd);
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
                    D.Pow_Norm_stim{tn,1}(:,b, :) =  nanmean(PSD(:,BandInfo.bandid{b}(1) : BandInfo.bandid{b}(2),:) , 2);
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
        E  = secog_addEventMarker(Events, subjNum, Fs_ds , 'CalcAveragePattern_seqType' , 'NumWarpSampFast' , NumWarpSampFast, 'NumWarpSampSlow'  ,NumWarpSampSlow)';
        
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
                eval(['D.Pow_Norm_stim{',num2str(tn),'} = P.PSD', num2str(tn), ';']);
            end
            % this loads up the data strutre where the raw PSD for that block is already stored 
            % 'ParseEEG-calcPSD'
            % Use the average patterns of block groups to warp them
            nPSD = D.Pow_Norm_stim;
            
            % set the sequence length for "*******" trilas to one
            D.seqlength(isnan(D.seqlength)) = 1;
            D.NumWarpSamp = NumWarpSampSlow* ones(size(D.TN));
            D.NumWarpSamp(logical(D.Fast)) = NumWarpSampFast;
            for tn = 1:length(D.TN)
                if length(find(D.NormTime{tn})) == D.seqlength(tn) +1 & ~D.isError(tn) & length(D.RealTime{tn})>D.NumWarpSamp(tn)
                    % check which block group this block falls into
                    
                    A = nPSD{tn};
                    % the normalized event markers in trial tn
                    idx = [0 find(D.NormTime{tn}) D.NumWarpSamp(tn)];
                    
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
                    
                    for e = 2:length(idn)
                        % make sure that each
                        idd   = linspace(1 , [idx(e)+1 - (idx(e-1) + 1)] , diffNEM(e-1));
                        for ch = 1:size(A,1)
                            for fi=1:size(A,2)
%                                 idd = floor(linspace(1, length([1:idx(e) - idx(e-1)]) , length(idd)));
                                D.PSD_stim{tn ,1}(ch,fi,idn(e-1)+1:idn(e)) = interp1([1:idx(e) - idx(e-1)] , squeeze(A(ch,fi,idx(e-1)+1:idx(e))) , idd);
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
                    D.Pow_Norm_stim{tn,1}(:,b, :) =  nanmean(PSD(:,BandInfo.bandid{b}(1) : BandInfo.bandid{b}(2),:) , 2);
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
    
        
end

