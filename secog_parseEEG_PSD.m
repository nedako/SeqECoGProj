function Dout  = secog_parseEEG_PSD(what , Dall , subjNum)
%% reads the all channels-packed EEG data, uses the BlockInfo file to parse the EEG into single trials and ammend the bihavioral data structure
%% It also calculates the PSD on the whole block and then parses up the PSD inot trials. This is mainly to avoid any window effect on single trials
subjname = {'P2'};
mainDir = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' subjname{subjNum} , '/Packed/'] ;
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
                % the key press trials
                if tn>2
                    D.EEG{tn ,1}       = Beeg.values(: , start_tr(tn)-floor(Fs/2) : end_tr(tn)+floor(Fs/2));
                    D.timeStim{tn ,1}  = 0:1/Fs:(length(D.EEG{tn ,1})/Fs) - 1/Fs;
                    D.timeStim{tn ,1}  = D.timeStim{tn ,1}-(500/Fs);
                    % convert the time of the forst press from hand to eeg
                    tp = floor(Fs/5) + floor(Fs*(D.pressTime0(tn)/1000));
                    D.timePress{tn ,1} = 0:1/Fs:(length(D.EEG{tn ,1})/Fs) - 1/Fs;
                    D.timePress{tn ,1}  = D.timePress{tn ,1}-(tp/Fs);
                    % the 2 idle trials at the beggining
                else
                    D.EEG{tn ,1} = Beeg.values(: , start_tr(tn) : end_tr(tn));
                    D.timeStim{tn ,1} = 0:1/Fs:(length(D.EEG{tn ,1})/Fs) - 1/Fs;
                    D.timePress{tn ,1} = 0:1/Fs:(length(D.EEG{tn ,1})/Fs) - 1/Fs;
                end
                [i tn]
            end
            Dout = addstruct(Dout , D);
        end
        %% save he labels once
        saveName = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' subjname{subjNum} '/','ChanLabels.mat'];
        ChanLabels = Data.label;
        save(saveName , 'ChanLabels');
        
        % save the data
        saveName = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' subjname{subjNum} '/','AllData.mat'];
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
        DownsampleRate = 10;
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
                statement = ['RawEEGpower',num2str(ch),'= secog_waveletPSD(Beeg.values(ch , :) , Fs , ''DownsampleRate'' , DownsampleRate);'];
                eval(statement);
                disp(['PSD calculation for block ' , num2str(i) , ', channel ' , num2str(ch) , ' completed'])
            end
            start_tr = floor(start_tr / DownsampleRate);
            end_tr = floor(end_tr / DownsampleRate);
            Fs_ds = floor(Fs/DownsampleRate);
            for tn = 1 :length(start_tr)
                clear Pow_Norm_stim Pow_Norm_pres Pow_non_norm idx PP PEEG
                % stack up the channles
                for ch = 1:size(Beeg.values , 1)
                    statement1 = ['D.PSD{tn,1}(ch , :,:)  = RawEEGpower',num2str(ch),'(: , start_tr(tn)-floor(Fs_ds/2) : end_tr(tn)+floor(Fs_ds/2));'];
                    eval(statement1);
                end
            end
            Dout = addstruct(Dout , D);
        end
    %% chop up the EEG data into trials, and filter out the power line noise
    %% also calculates the PSD on every block and and chops it up into trials
end

