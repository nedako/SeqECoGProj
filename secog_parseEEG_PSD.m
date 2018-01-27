function Dout  = secog_parseEEG_PSD(what , Dall , subjNum, varargin)
%% reads the all channels-packed EEG data, uses the BlockInfo file to parse the EEG into single trials and ammend the bihavioral data structure
%% It also calculates the PSD on the whole block and then parses up the PSD inot trials. This is mainly to avoid any window effect on single trials
c = 1;
while(c<=length(varargin))
    switch(varargin{c})
        
        case {'NormType'}
            % Type of PSD normalization
            % Default : 'stim' normalize power to 500 ms before the stimulation comes on
            % could be 'press' normalizez to 500 ms before the stim comes on + the time before the first press
            % could be 'none' which just transforms PSD into logarithmic scale
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'DownsampleRate'}
            % folds by which you wnat the PSD to be downsampled
            % default 10
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
         case {'NumWarpSamp'}
            % number of sine cycles to use in the Morlet wavelet
            % default 10
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        otherwise
            error(sprintf('Unknown option: %s',varargin{c}));
    end
end

if ~exist('DownsampleRate')
    DownsampleRate = 10;
end
if ~exist('NormType')
    NormType = 'stim';
end
if ~exist('NumWarpSamp')
    NumWarpSamp = 300;
end

%%
subjname = {'P2'};
mainDir = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' subjname{subjNum} , '/Packed/'] ;
saveDir = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' subjname{subjNum} , '/'] ;
cd(mainDir)
%% Import the BLock Info
[~, ~, BlockInfo] = xlsread([mainDir , 'BlockInfo.xlsx'],'Sheet1');
BlockInfo = BlockInfo(2:end,:);
BlockInfo(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),BlockInfo)) = {''};

idx = cellfun(@ischar, BlockInfo);
BlockInfo(idx) = cellfun(@(x) string(x), BlockInfo(idx), 'UniformOutput', false);

clearvars idx;
Dout = [];
%% chop up the EEG data into trials, and filter out the power line noise
switch what
    case 'ParseEEG-noPSD'
        Fs = 1024;
        Fo = 60;
        tn = 1;
        % multiple blocks are sotred in the same file. so avoid loading them up multiple times.
        fName1 =  BlockInfo{1,4};
        load(fName1);
        %% preprocess EEG and filtering
       
        % Q : quality factor is the center frequency divided by the bandwidth.
        % Q = 35;
        % BW = Fo/(Fs/2);
        % [b,a] = iircomb(10,BW,'notch');
        
        wo = Fo/(Fs/2);  bw = wo/35; % notch to eliminate 60 Hz
        [b,a] = iirnotch(wo,bw);
        
        wo = 2*Fo/(Fs/2);  bw = wo/60;% notch to eliminate the first harmonic
        [c,d] = iirnotch(wo,bw);
        %%
        
        for i = 1:size(BlockInfo , 1)
            D = getrow(Dall , Dall.BN == i);
            fName = BlockInfo{i,4};
            if ~strcmp(fName , fName1)
                load(fName);
                % filter the power line noise out of the whole data and then chopp it up
                for ch = 1:size(Data.values , 1)
                    A = filter(b,a , Data.values(ch , :));
                    Data.values(ch , :) = filter(c,d , A);
                end
                fName1 = fName;
            end
            %     extract the data for the block at hand
            BlockRang = [BlockInfo{i,2} : BlockInfo{i,3}];
            Beeg = Data;
            Beeg.values = Data.values(: , BlockRang);
            marker = Beeg.values(find(strcmp(Data.label , 'TTL')) , :);
            marker = [0 diff(marker <-2*10^6)];
            start_tr = find(marker == 1); % starts of trials
            end_tr = find(marker == -1);  % ends of trials
            for tn = 1:length(start_tr)
                % take 200 ms before and after the stim come on and goes off
                D.EEG{tn ,1}       = Beeg.values(: , start_tr(tn)-floor(Fs/5) : end_tr(tn)+floor(Fs/5));
                [i tn]
            end
            Dout = addstruct(Dout , D);
        end
        %% save he labels once
        saveName = [saveDir,'ChanLabels.mat'];
        ChanLabels = Data.label;
        save(saveName , 'ChanLabels');
        
        % save the data
        saveName = [saveDir,'AllData_EEG.mat'];
        save(saveName , 'Dall' , '-v7.3');
        
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
    case 'ParseEEG-calcPSD'
         %% preprocess EEG and filtering
        
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
        
        for ch = 1:size(Data.values , 1)
            A = filter(b,a , Data.values(ch , :));
            Data.values(ch , :) = filter(c,d , A);
        end
        %% 
        for i = 1:size(BlockInfo , 1)
            D = getrow(Dall , Dall.BN == i);
            fName = BlockInfo{i,4};
            if ~strcmp(fName , fName1)
                load(fName);
                % filter the power line noise out of the whole data and then chopp it up
                for ch = 1:size(Data.values , 1)
                    A = filter(b,a , Data.values(ch , :));
                    Data.values(ch , :) = filter(c,d , A);
                end
                fName1 = fName;
            end
            %     extract the data for the block at hand
            BlockRang = [BlockInfo{i,2} : BlockInfo{i,3}];
            Beeg = Data;
            Beeg.values = Data.values(: , BlockRang);
            % get the indecies for starts ans ends of the trials
            marker = Beeg.values(find(strcmp(Data.label , 'TTL')) , :);
            marker = [0 diff(marker <-2*10^6)];
            start_tr = find(marker == 1); % starts of trials
            end_tr = find(marker == -1);  % ends of trials
            for ch = 1:size(Beeg.values , 1)
                statement = ['[RawEEGpower',num2str(ch),' , BandInfo] = secog_waveletPSD(Beeg.values(ch , :) , Fs , ''DownsampleRate'' , DownsampleRate);'];
                eval(statement);
                disp(['PSD calculation for block ' , num2str(i) , ', channel ' , num2str(ch) , ' completed'])
            end
            start_tr = floor(start_tr / DownsampleRate);
            end_tr = floor(end_tr / DownsampleRate);
            Fs_ds = floor(Fs/DownsampleRate);
            for tn = 1 :length(start_tr)
                clear Pow_Norm_stim Pow_Norm_pres Pow_non_norm idx PP PEEG
                % stack up the channles
                % take 200 msec before the stim apears and half a sec after the stimulus goes away into account
                for ch = 1:size(Beeg.values , 1)
                    statement1 = ['D.PSD{tn,1}(ch , :,:)  = RawEEGpower',num2str(ch),'(: , start_tr(tn)-floor(Fs_ds/5) : end_tr(tn)+floor(Fs_ds/5));'];
                    eval(statement1);
                end
            end
            % save the un-banded version block by block
            saveName = [saveDir,'Raw_PSD_B',num2str(i) ,'.mat'];
            Pall = D;
            save(saveName , 'Pall' , '-v7.3');
             % then band-average and stack up
            temp = Pall;
            clear Pall
            for b =1:length(BandInfo.bandid)
                PSD(b, :) =  nanmean(temp(BandInfo.bandid{b}(1) : BandInfo.bandid{b}(2),:));
            end
            Dout = addstruct(Dout , D);
            clear D
        end
        Dall = Dout;
        clear Dout;
        saveName = [saveDir,'AllData_PSD.mat'];
        save(saveName , 'Dall' , '-v7.3');
    case 'NormalizePSD'
        % obtain the time normalized event stamps for Block bins
        Fs = 1024;
        Fs_ds = floor(Fs/DownsampleRate);
        Dall  = secog_parseEEG_PSD('ParseEEG-calcPSD' , Dall, subjNum);
        D  = secog_addEventMarker(Dall, subjNum, Fs_ds , 'addEvent' , 'NumWarpSamp' , NumWarpSamp)';
        
        % this takes the data strutre where the raw PSD is already stored (so the outpur of the prevous case)
        % and normalizes the PSD for each Block to the average PSD of the firs two trails with the stars.
        for tn = 1:length(D.TN)
            A = D.PSD{tn};
            if D.TN(tn)>2
                for ch = 1:size(A,1)
                    for fi=1:size(A,2)
                        switch NormType
                            case 'stim'
                                Dout.Pow_Norm_stim{tn}(ch,fi,:) = squeeze(10*log10(A(ch,fi,:)./mean(A(ch,fi,1:find(D.EventMarker{tn ,1} == -1)))));
                            case 'press'
                                Dout.Pow_Norm_pres{tn}(ch,fi,:)      = squeeze(10*log10(A(ch,fi,:)./mean(A(ch,fi,1:find(D.EventMarker{tn ,1} == 1)))));
                            case 'none'
                                Dout.Pow_non_norm{tn}(ch,fi,:)       = squeeze(10*log10(A(ch,fi,:)));
                        end
                    end
                end
            else
                for ch = 1:size(A,1)
                    for fi=1:size(A,2)
                        switch NormType
                            case {'stim' , 'press'}
                                Dout.Pow_Norm_stim{tn}(ch,fi,:) = squeeze(10*log10(A(ch,fi,:)./mean(A(ch,fi,1:find(D.EventMarker{tn ,1} == -1)))));
                            case 'none'
                                Dout.Pow_non_norm{tn}(ch,fi,:)       = squeeze(10*log10(A(ch,fi,:)));
                        end
                    end
                end
            end
            
        end
%%        
    case 'TimeWarpPSD'
        % obtain the time normalized event stamps for Block bins
        Fs = 1024;
        Fs_ds = floor(Fs/DownsampleRate);
        Pall  = secog_parseEEG_PSD('NormalizePSD' , Dall, subjNum);
        % find event markers for each trial
        D  = secog_addEventMarker(Dall, subjNum, Fs_ds , 'addEvent' , 'NumWarpSamp' , NumWarpSamp)';
        % find average event markers
        E  = secog_addEventMarker(Dall, subjNum, Fs_ds , 'CalcAveragePattern' , 'NumWarpSamp' , NumWarpSamp)';
        
        % this takes the data strutre where the raw PSD is already stored (so the outpur of the prevous case)
        % and normalizes the PSD for each Block to the average PSD of the firs two trails with the stars.
        % check for the type of normalization that has been performed on the PSD
        normName = fieldnames(Pall);
        eval(char(strcat('nPSD = Pall.' , normName , ';')))
        
        % set the sequence length for "*******" trilas to one
        D.seqlength(isnan(D.seqlength)) = 1;
        
        for tn = 1:length(D.TN)
            if length(find(D.NormEventMarker{tn})) == D.seqlength(tn) +1
                tn
                % check which block group this block falls into
                bg = 1;
                mem = 0;
                while mem == 0
                    if ismember(D.BN(tn) , E.blockGroups{bg})
                        mem = 1;
                        BG = bg;
                    end
                    bg = bg+1;
                end
                A = nPSD{tn};
                % the normalized event markers in trial tn
                idx = [0 find(D.NormEventMarker{tn})' NumWarpSamp];
                % find the row number corresponding to the seqNumb
                sn = find(E.SN{BG}==D.seqNumb(tn));
                % average normalized time stamps for the sn , BG
                if tn>2
                    idn = [0 E.NEM{BG}(sn , ~isnan(E.NEM{BG}(sn ,:))) NumWarpSamp];
                    diffNEM = diff(idn);
                else
                    % for "* * * * * * * *" trials, the average pattern is the same as the trial since there is no variability
                    idn = idx; 
                    diffNEM = diff(idn);
                end
                
                
                for e = 2:length(idn)
                    % make sure that each
                    idd   = linspace(1 , [idn(e)+1 - (idn(e-1) + 1)] , diffNEM(e-1));
                    for ch = 1:size(A,1)
                        for fi=1:size(A,2)
                            Dout.PSD_stim{tn ,1}(ch,fi,idn(e-1)+1:idn(e)) = interp1([1:idx(e) - idx(e-1)] , squeeze(nPSD{tn}(ch,fi,idx(e-1)+1:idx(e))) , idd);
                        end
                    end
                end
            else
               Dout.PSD_stim{tn ,1} = []; 
            end
        end
        Pall = Dout;
        saveName = [saveDir,'Warped_PSD.mat'];
        save(saveName , 'Pall' , '-v7.3');
end
