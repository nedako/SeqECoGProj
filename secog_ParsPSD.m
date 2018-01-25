function secog_ParsPSD(Dall , what , GoodElec , MarkerElec)
%% %% chop up the EEG data into trials
% For this particular subject the timer was running frm the time that the START_TRIAL came on

% for patient 1 GoodElec   = [1:30,32:60,63:74,77:80];
% for patient 1 MarkerElec = 141;
alltrials = 1;
Fs = 1024;
% definitions, selections...
BN = unique(Dall.BN);
min_freq =  2;
max_freq = 150;
num_frex = 90;

% define wavelet parameters
time = -1:1/Fs:1;
frex = logspace(log10(min_freq),log10(max_freq),num_frex);
% s    = logspace(log10(3),log10(10),num_frex)./(2*pi*frex);
% s    =  3./(2*pi*frex); % with 3 cycles
s    = 10./(2*pi*frex); % with 10 cycles

% filter the power line noise
Fs = 1024;
Fo = 60;
% Q : quality factor is the center frequency divided by the bandwidth.
% Q = 35;


wo = Fo/(Fs/2);  bw = wo/60; % notch to eliminate 60 Hz
[b,a] = iirnotch(wo,bw);

wo = 2*Fo/(Fs/2);  bw = wo/100;% notch to eliminate the first harmonic
[c,d] = iirnotch(wo,bw);

% Filter Visualization
Fviz = 1;
if Fviz ==1
    h = fvtool(b,a);
    h.Fs = Fs;
    h.FrequencyRange='[-Fs/2, Fs/2)';
    zplane(b,a)
    
    h = fvtool(c,d);
    h.Fs = Fs;
    h.FrequencyRange='[-Fs/2, Fs/2)';
    zplane(c,d)
end




%% Time normalizing each trial to the average pattern
switch what
    case 'Average'
        alltrials = 251;
        % definte convolution parameters
        n_wavelet            = length(time);% sampling frequency of the wavelet
        DSR = 4;
        NormSamp = 500;
        for i = 8:length(BN )
            Fs = 1024;
            fName = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/iEEG_P1_Sep14_BlockParsed/EEG_BN' , num2str(i) , '.mat'];
            for ch  = 1:length(GoodElec)
                statement = ['clear RawEEGpower',num2str(ch)];
                eval(statement);
            end
            load(fName);
            eval(['nofilt_A = EEG_BN' , num2str(i) , '(GoodElec , :);']);
            for i2 = 1:size(nofilt_A,1)
                filt_A(i2 , :) = filter(b,a , nofilt_A(i2,:));
                A(i2 , :) = filter(c,d , filt_A(i2,:));
            end
            
            if Fviz
                t = 0:(1/Fs):(length(A)/Fs) - (1/Fs);
                figure('color' , 'white')
                subplot(2,1,1)
                periodogram(nofilt_A(1,:),[],length(nofilt_A(1,:)),Fs,'power')
                subplot(2,1,2)
                periodogram(A(1,:),[],length(A(1,:)),Fs,'power')
            end
            
            clear filt_A
            eval(['marker = EEG_BN' , num2str(i) , '(MarkerElec , :);']);
            eval(['clear ' , 'EEG_BN' , num2str(i)])
            marker = [0 diff(marker <-2*10^6)];
            start_tr = find(marker == 1); % starts of trials
            end_tr = find(marker == -1);  % ends of trials
            n_data               = size(A , 2);%...EEG.pnts*EEG.trials;
            for ch = 1:size(A,1)
                data                 = A(ch,:);
                n_convolution        = n_wavelet+n_data-1;
                n_conv_pow2          = pow2(nextpow2(n_convolution));
                half_of_wavelet_size = (n_wavelet-1)/2;
                % get FFT of data
                eegfft = fft(data,n_conv_pow2);
                
                % loop through frequencies and compute synchronization
                for fi=1:num_frex
                    wavelet = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2 );
                    
                    % convolution
                    eegconv = ifft(wavelet.*eegfft);
                    eegconv = eegconv(1:n_convolution);
                    eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
                    % assign a separate array to every chanle since for longer recordings the array size exceds the Matlab limit
                    statement = ['RawEEGpower',num2str(ch),'(fi,:) = mean(abs(reshape(eegconv,n_data,1)).^2,2);'];
                    eval(statement);
                end
                disp(['Wavelet being applied to Block ' , num2str(i) , ', Channel ' , num2str(ch)])
            end
            Fs_n = floor(Fs/DSR);
            for j = 1 :length(start_tr)
                clear Pow_Norm_stim Pow_Norm_pres Pow_non_norm idx PP PEEG
                
                if j>2
                    % stack up the channles
                    for ch = 1:size(A,1)
                        statement1 = ['PEEG(ch , :,:)  = RawEEGpower',num2str(ch),'(: , start_tr(j)-floor(Fs/2) : end_tr(j)+floor(Fs/2));'];
                        eval(statement1);
                    end
                    trialTime  = 0:1/Fs:(size(PEEG  ,3)/Fs) - 1/Fs;
                    PP.EventMarker{j ,1} = zeros(size( trialTime));
                    PP.NormEventMarker{j ,1} = zeros(NormSamp,1);
                    idx(1) = Fs/2;
                    PP.EventMarker{j ,1}(idx(1)) = -1;       % stimulus apears
                    
                    Nidx(1) = floor(idx(1)*NormSamp/length(trialTime));
                    PP.NormEventMarker{j ,1}(Nidx(1)) = -1;
                    
                    
                    idc = 2;
                    for pr = 1:Dall.seqlength(alltrials)
                        numTaskSamp = floor(( Dall.AllPressIdx(alltrials , Dall.seqlength(alltrials)) - Dall.AllPressIdx(alltrials , pr)) * (Fs/500));
                        idx(idc) = numTaskSamp + Fs/2;
                        Nidx(idc) = floor(idx(idc)*NormSamp/length(trialTime));
                        PP.EventMarker{j ,1}(end - idx(idc)) = pr;  % First press
                        PP.NormEventMarker{j ,1}(end - Nidx(idc)) = pr;
                        idc = idc + 1;
                    end
                    for ch = 1:size(A,1)
                        for fi=1:num_frex
                            Pow_Norm_stim(ch,fi,:) = squeeze(10*log10(PEEG(ch,fi,:)./mean(PEEG(ch,fi,1:find(PP.EventMarker{j ,1} == -1)))));
                            Pow_Norm_pres(ch,fi,:) = squeeze(10*log10(PEEG(ch,fi,:)./mean(PEEG(ch,fi,1:find(PP.EventMarker{j ,1} == 1)))));
                            Pow_non_norm(ch,fi,:)  = squeeze(10*log10(PEEG(ch,fi,:)));
                        end
                    end
                    
                    
                    
                    idx = [0 idx(1) size(PEEG , 3)-idx(2:end)  size(PEEG , 3)];
                    switch Dall.seqNumb(alltrials)
                        case {0}
                            lineNum = 1;
                        case {1}
                            lineNum = 2;
                        case {2}
                            lineNum = 3;
                        case {3}
                            lineNum = 4;
                        case {4}
                            lineNum = 5;
                        case {11}
                            lineNum = 6;
                        case {22}
                            lineNum = 7;
                        case {33}
                            lineNum = 8;
                        case {44}
                            lineNum = 9;
                        case {55}
                            lineNum = 10;
                        case {103}
                            lineNum = 11;
                        case {104}
                            lineNum = 12;
                        case {203}
                            lineNum = 13;
                        case {204}
                            lineNum = 14;
                    end
                    
                    idn = [0 NEM{lineNum} NormSamp];
                    for e = 2:length(idn)
                        idd   = linspace(1 , [idx(e)+1 - (idx(e-1) + 1)] , diffNEM{lineNum}(e-1));
                        for ch = 1:size(A,1)
                            for fi=1:num_frex
                                Pall.PEEG_none{j ,1}(ch,fi,idn(e-1)+1:idn(e)) = interp1([1:idx(e) - idx(e-1)] , squeeze(Pow_non_norm(ch,fi,idx(e-1)+1:idx(e))) , idd);
                                Pall.PEEG_stim{j ,1}(ch,fi,idn(e-1)+1:idn(e)) = interp1([1:idx(e) - idx(e-1)] , squeeze(Pow_Norm_stim(ch,fi,idx(e-1)+1:idx(e))) , idd);
                                Pall.PEEG_pres{j ,1}(ch,fi,idn(e-1)+1:idn(e)) = interp1([1:idx(e) - idx(e-1)] , squeeze(Pow_Norm_pres(ch,fi,idx(e-1)+1:idx(e))) , idd);
                                Pall.SeqNumb(j , 1) = Dall.seqNumb(alltrials);
                            end
                        end
                    end
                    Pall.AverageEventMarker{j ,1} = idx;
                else
                    for ch = 1:size(A,1)
                        statement1 = ['PEEG(ch , :,:)  = RawEEGpower',num2str(ch),'(: , start_tr(j)-floor(Fs/2) : end_tr(j)+floor(Fs/2));'];
                        eval(statement1);

                    end
                    
                    trialTime  = 0:1/Fs:(size(PEEG  ,3)/Fs) - 1/Fs;
                    
                    PP.EventMarker{j ,1} = zeros(size( trialTime));
                    PP.NormEventMarker{j ,1} = zeros(NormSamp,1);
                    idx = Fs/2;
                    PP.EventMarker{j ,1}(idx) = -1;       % equivalent to stimulus apears
                    Nidx = floor(idx*NormSamp/length(trialTime));
                    PP.NormEventMarker{j ,1}(Nidx) = -1;
                    
                    idx = Fs/2;
                    Nidx = floor(idx*NormSamp/length(trialTime));
                    PP.EventMarker{j ,1}(end - idx) = 7;  % equivalent to last press
                    PP.NormEventMarker{j ,1}(end - Nidx) = 7;
                    
                    
                    Pall.AverageEventMarker{j ,1} = idx;
                    Pall.SeqNumb(j , 1) = Dall.seqNumb(alltrials);
                    for ch = 1:size(A,1) 
                        for fi=1:num_frex
                            Pow_Norm_stim(ch,fi,:) = squeeze(10*log10(PEEG(ch,fi,:)./mean(PEEG(ch,fi,1:find(PP.EventMarker{j ,1} == -1)))));
                            Pow_non_norm(ch,fi,:)  = squeeze(10*log10(PEEG(ch,fi,:)));
                            
                            
                            idd   = linspace(1 , size(Pow_Norm_stim , 3) , NormSamp);
                            Pall.PEEG_stim{j ,1}(ch , fi , :) = interp1([1:size(Pow_Norm_stim  , 3)] , squeeze(Pow_Norm_stim(ch,fi,:)) , idd);
                            Pall.PEEG_none{j ,1}(ch , fi , :) = interp1([1:size(Pow_non_norm  , 3)] , squeeze(Pow_non_norm(ch,fi,:)) , idd);
                        end
                    end
                    Pall.PEEG_pres{j ,1} = NaN;
                    
                end
                disp(['PSD being normalized for Block ' , num2str(i) , ', Trial ' , num2str(j)])
                alltrials = alltrials  +1;
            end
            
            fname1 = ['AVG_Norm_PSD_b',num2str(i),'_P1_Sep14.mat'];
            fname2 = 'Pall';
            sname  = '-v7.3';
            eval('save(fname1 , fname2 , sname);')
            clear Pall A filt_A nofilt_A
        end
        
    case 'AvgPSD_all'
        %% stacking up all the
        BN = unique(Dall.BN);
        
        clear PSD
        for bn = 1:length(BN)
            fname1 = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/iEEG_P1_Sep14_BlockParsed/AVG_Norm_PSD_b',num2str(BN(bn)),'_P1_Sep14.mat'];
            load(fname1);
            D = getrow(Dall , ismember(Dall.BN , BN(bn)));
            SN = unique(D.seqNumb);
            SN = SN(SN~=5);
            for sn = 1:length(SN)
                PSD.Pow_stim{bn , sn} = zeros(size(Pall.PEEG_pres{3} , 2) , size(Pall.PEEG_stim{3} , 3));
                PSD.Pow_pres{bn , sn} = zeros(size(Pall.PEEG_pres{3} , 2) , size(Pall.PEEG_pres{3} , 3));
                PSD.Pow_non{bn , sn}  = zeros(size(Pall.PEEG_none{3} , 2) , size(Pall.PEEG_none{3} , 3));
                Pp = getrow(Pall , ~D.isError & Pall.SeqNumb == SN(sn));
                if length(Pp.SeqNumb)>0
                    for tn = 1:length(Pp.SeqNumb)
                        PSD.Pow_stim{bn , sn} = PSD.Pow_stim{bn , sn} + squeeze(mean(Pp.PEEG_stim{tn} , 1));
                        PSD.Pow_pres{bn , sn} = PSD.Pow_pres{bn , sn} + squeeze(mean(Pp.PEEG_pres{tn} , 1));
                        PSD.Pow_non{bn , sn}  = PSD.Pow_non{bn , sn} + squeeze(mean(Pp.PEEG_none{tn} , 1));
                    end
                    PSD.Pow_stim{bn , sn} = PSD.Pow_stim{bn , sn}/tn;
                    PSD.Pow_pres{bn , sn} = PSD.Pow_pres{bn , sn}/tn;
                    PSD.Pow_non{bn , sn}  = PSD.Pow_non{bn , sn}/tn;
                end
            end
        end
        PSD.SeqNumb = NaN*ones(length(BN) , 5);
        PSD.BN = NaN*ones(length(BN) , 5);
        
        for bn = 1:length(BN)
            D = getrow(Dall , ismember(Dall.BN , BN(bn)));
            SN = unique(D.seqNumb);
            SN = SN(SN~=5);
            for sn = 1:length(SN)
                PSD.SeqNumb(bn , sn) = SN(sn);
                PSD.BN(bn , sn) = BN(bn);
            end
        end    
    case 'AvgPSD_HiVar'
        %% stacking up all the
        
        BN = unique(Dall.BN);
        SN = unique(Dall.seqNumb);
        SN = SN(SN~=5);
        % sort channles based on variance based on one parsed block
        fName = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/iEEG_P1_Sep14_BlockParsed/EEG_BN1.mat'];
        load(fName);
        ChanVar = var(EEG_BN1(GoodElec , :)');
        [ChVar , ChId] = sort(ChanVar , 'descend');
        figure('color' , 'white');
        plot(ChVar , 'LineWidth' , 3);
        grid on
        hold on
        ax = gca;
        TL = cellfun(@num2str , mat2cell(ChId , 1, ones(1,length(GoodElec))),'UniformOutput',false);
        ax.XTick = 1:length(GoodElec);
        ax.XTickLabel = TL;
        
        clear PSD
        for bn = 1:length(BN)
            fname1 = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/iEEG_P1_Sep14_BlockParsed/AVG_Norm_PSD_b',num2str(BN(bn)),'_P1_Sep14.mat'];
            load(fname1);
            D = getrow(Dall , ismember(Dall.BN , BN(bn)));
            SN = unique(D.seqNumb);
            SN = SN(SN~=5);
            for sn = 1:length(SN)
                PSD.Pow_stim{bn , sn} = zeros(size(Pall.PEEG_pres{3} , 2) , size(Pall.PEEG_stim{3} , 3));
                PSD.Pow_pres{bn , sn} = zeros(size(Pall.PEEG_pres{3} , 2) , size(Pall.PEEG_pres{3} , 3));
                PSD.Pow_non{bn , sn}  = zeros(size(Pall.PEEG_none{3} , 2) , size(Pall.PEEG_none{3} , 3));
                Pp = getrow(Pall , ~D.isError & Pall.SeqNumb == SN(sn));
                if length(Pp.SeqNumb)>0
                    for tn = 1:length(Pp.SeqNumb)
                        PSD.Pow_stim{bn , sn} = PSD.Pow_stim{bn , sn} + squeeze(mean(Pp.PEEG_stim{tn}(ChId(1:7),:,:) , 1));
                        PSD.Pow_pres{bn , sn} = PSD.Pow_pres{bn , sn} + squeeze(mean(Pp.PEEG_pres{tn}(ChId(1:7),:,:) , 1));
                        PSD.Pow_non{bn , sn}  = PSD.Pow_non{bn , sn} + squeeze(mean(Pp.PEEG_none{tn}(ChId(1:7),:,:) , 1));
                    end
                    PSD.Pow_stim{bn , sn} = PSD.Pow_stim{bn , sn}/tn;
                    PSD.Pow_pres{bn , sn} = PSD.Pow_pres{bn , sn}/tn;
                    PSD.Pow_non{bn , sn}  = PSD.Pow_non{bn , sn}/tn;
                end
            end
        end
        PSD.SeqNumb = NaN*ones(length(BN) , 5);
        PSD.BN = NaN*ones(length(BN) , 5);
        
        for bn = 1:length(BN)
            D = getrow(Dall , ismember(Dall.BN , BN(bn)));
            SN = unique(D.seqNumb);
            SN = SN(SN~=5);
            for sn = 1:length(SN)
                PSD.SeqNumb(bn , sn) = SN(sn);
                PSD.BN(bn , sn) = BN(bn);
            end
        end
        
        
    case 'Visualize'
        SN = unique(Dall.seqNumb);
        SN = SN(SN~=5);
        SN = [0 1 2 3 4 11 22 33 44 55 103 203 104 204];
        % sort channles based on variance based on one parsed block
        fName = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/iEEG_P1_Sep14_BlockParsed/EEG_BN1.mat'];
        load(fName);
        ChanVar = var(EEG_BN1(GoodElec , :)');
        [ChVar , ChId] = sort(ChanVar , 'descend');
        figure('color' , 'white');
        plot(ChVar , 'LineWidth' , 3);
        grid on
        hold on
        ax = gca;
        TL = cellfun(@num2str , mat2cell(ChId , 1, ones(1,length(GoodElec))),'UniformOutput',false);
        ax.XTick = 1:length(GoodElec);
        ax.XTickLabel = TL;
        HV = 0;
        if ~HV
            fName = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/Average_TimeWarped_PSD.mat'];
        else
            fName = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/Average_TimeWarped_HighVar_PSD.mat'];
        end
        
        load(fName);
        clear EventMarker NormEventMarker EM NEM diffNEM
        for sn= 1:length(SN)
            D = getrow(Dall , ismember(Dall.seqNumb , SN(sn)) & ~Dall.isError);
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
            EM{sn, :} = floor(mean(EventMarker{sn}));
            NEM{sn, :} = floor(mean(NormEventMarker{sn}));
            diffNEM{sn, :} = diff([0 , NEM{sn} , 500]);
        end
        
        Fs = 1024;
        alltrials = 1;
        % definitions, selections...
        BN = unique(Dall.BN);
        min_freq =  2;
        max_freq = 150;
        num_frex = 90;
        % define wavelet parameters
        time = -1:1/Fs:1;
        frex = logspace(log10(min_freq),log10(max_freq),num_frex);
        SN = unique(Dall.seqNumb);
        SN = SN(SN~=5);
        bandsLab = {'Delta <4Hz' , 'Theta 4-8Hz' , 'Alpha 8-13Hz' , 'L-Beta 13-24Hz' , 'H-Beta 24-36Hz' , 'L-Gamma 36-48Hz' , 'H-Gamma >48Hz'};
        bands = {[0 4], [4 8] [8 13] [13 24] [14 36] [36 48] [48 110]};
        for b = 1:length(bands)
            bandid{b} = [find(frex>bands{b}(1) ,1, 'first') , find(frex<bands{b}(2) ,1, 'last')];
        end
        
        
        
        %%%%%%%%%%%   FINGER REPETITIONs     
        figure('color' , 'white')
        figCount = 1;
        snCount = 1;
        for sn = 6 : 10
            F = ismember(PSD.SeqNumb , SN(sn));
            figCount = snCount;
            BN = PSD.BN(F);
            Pow_non  = PSD.Pow_non(F);
            Pow_stim = PSD.Pow_stim(F);
            Pow_pres = PSD.Pow_pres(F);
            for bn = 1:length(BN)
                subplot(length(BN),5 , figCount)
                contourf([1:500],frex,Pow_stim{bn},60,'linecolor','none')
                caxis([-15 8])
                colorbar
                title (['Rep-finger ' , num2str(sn-5) ,', Block ' , num2str(BN(bn))])
                for lin = 1:length(NEM{sn})
                    line([NEM{sn}(lin) NEM{sn}(lin)] , [1 150] , 'color' , 'red' , 'LineStyle' , ':' , 'LineWidth' , 3)
                end
                xlabel('Norm Time')
                figCount = figCount + 5;
            end
            snCount = snCount+1;
        end
        
        figure('color' , 'white')
        snCount = 1;
        for sn = 6 : 10
            figCount = snCount;
            F = ismember(PSD.SeqNumb , SN(sn));
            BN = PSD.BN(F);
            Pow_non  = PSD.Pow_non(F);
            Pow_stim = PSD.Pow_stim(F);
            Pow_pres = PSD.Pow_pres(F);
            for bn = 1:length(BN)
                clear Bands T
                
                for b =1:length(bandid)
                    subplot(length(BN),5 , figCount)
                    Bands(b, :) =  nanmean(Pow_stim{bn}(bandid{b}(1) : bandid{b}(2),:));
                    plot([1:500] , Bands(b, :)+7*b , 'LineWidth' , 3)
                    T(b) = Bands(b, 1)+7*b;
                    hold on
                end
                title (['Rep-finger ' , num2str(sn-5) ,', Block ' , num2str(BN(bn))])
                for lin = 1:length(NEM{sn})
                    line([NEM{sn}(lin) NEM{sn}(lin)] , [0 60] , 'color' , 'red' , 'LineStyle' , ':' , 'LineWidth' , 3)
                end
                xlabel('Norm Time')
                if snCount == 1
                    set(gca ,'YTickLabels' , bandsLab, 'YTick' , T );
                else
                    set(gca ,'YTickLabels' , [], 'YTick' , [] );
                end
                line([100 100] , [55 60],'color' , 'black' , 'LineWidth' , 5)
                text(60,57.5,'5','FontSize' , 16 )
                set(gca , 'XLim' , [1 500] , 'YLim' , [0,60],'FontSize' , 16,'Box' , 'off')
                figCount = figCount + 5;
            end
            snCount = snCount+1;
        end
        
        %%%%%%%%%%%%%%%% CHUNKS
        figure('color' , 'white')
        
        snCount = 1;
        for sn = 11 : 14
            figCount = snCount;
            F = ismember(PSD.SeqNumb , SN(sn));
            BN = PSD.BN(F);
            Pow_non  = PSD.Pow_non(F);
            Pow_stim = PSD.Pow_stim(F);
            Pow_pres = PSD.Pow_pres(F);
            
            %     D = getrow(Dall , ismember(Dall.seqNumb , SN(sn)));
            %     BN = unique(D.BN);
            chunklab = {'Trip1' ,'Trip2' ,  'Quad1' , 'Quad2'};
            for bn = 1:length(BN)        
                subplot(length(BN),4 , figCount)
                contourf([1:500],frex,Pow_stim{bn},60,'linecolor','none')
                caxis([-15 7])
                colorbar
                ylabel('Frequency(Hz)')
                title ([chunklab{sn-10} , ' ',', Block ' , num2str(BN(bn))])
                for lin = 1:length(NEM{sn})
                    line([NEM{sn}(lin) NEM{sn}(lin)] , [1 150] , 'color' , 'red' , 'LineStyle' , ':' , 'LineWidth' , 3)
                end
                
                xlabel('Norm Time')
                figCount = figCount + 4;
            end
            snCount = snCount+1;
        end
        
        figure('color' , 'white')
        
        snCount = 1;
        for sn = 11 : 14
            figCount = snCount;
            F = ismember(PSD.SeqNumb , SN(sn));
            BN = PSD.BN(F);
            Pow_non  = PSD.Pow_non(F);
            Pow_stim = PSD.Pow_stim(F);
            Pow_pres = PSD.Pow_pres(F);
            
            %     D = getrow(Dall , ismember(Dall.seqNumb , SN(sn)));
            %     BN = unique(D.BN);
            chunklab = {'Trip1' ,'Trip2' ,  'Quad1' , 'Quad2'};
            for bn = 1:length(BN)
                subplot(length(BN),4 , figCount)
                clear Bands T
                for b =1:length(bandid)
                    Bands(b, :) =  nanmean(Pow_stim{bn}(bandid{b}(1) : bandid{b}(2),:));
                    plot([1:500] , Bands(b, :)+7*b , 'LineWidth' , 3)
                    T(b) = Bands(b, 1)+7*b;
                    hold on
                end
                set(gca , 'XLim' , [1 500] , 'YLim' , [0,60],'FontSize' , 16,'Box' , 'off')
                if snCount == 1
                    set(gca ,'YTickLabels' , bandsLab, 'YTick' , T );
                else
                    set(gca ,'YTickLabels' , [], 'YTick' , [] );
                end
                title ([chunklab{sn-10} , ' ',', Block ' , num2str(BN(bn))])
                for lin = 1:length(NEM{sn})
                    line([NEM{sn}(lin) NEM{sn}(lin)] , [0 60] , 'color' , 'red' , 'LineStyle' , ':' , 'LineWidth' , 3)
                end
                line([50 50] , [55 60],'color' , 'black' , 'LineWidth' , 5)
                text(10,57.5,'5','FontSize' , 16 )
                xlabel('Norm Time')
                figCount = figCount + 4;
            end
            snCount = snCount+1;
        end
        %% average every 3 blocks before chunk blocks for all the structures
        BG = {[1 2 3] [4 5 6] [7 8 9] [10 11 12]};
        SN = unique(Dall.seqNumb);
        SN = SN(SN~=5);
        Fs = 1024;
        % definitions, selections...
        BN = unique(Dall.BN);
        min_freq =  2;
        max_freq = 150;
        num_frex = 90;
        
        % define wavelet parameters
        time = -1:1/Fs:1;
        frex = logspace(log10(min_freq),log10(max_freq),num_frex);
        
        clear PN PS PP
        for sn = 1:5
            F = ismember(PSD.SeqNumb , SN(sn));
            BN = PSD.BN(F);
            Pow_non = PSD.Pow_non(F);
            Pow_stim = PSD.Pow_stim(F);
            Pow_pres = PSD.Pow_pres(F);
            %     D = getrow(Dall , ismember(Dall.seqNumb , SN(sn)));
            %     BN = unique(D.BN);
            for bg = 1:length(BG)
                PN{sn ,bg} = zeros(90 , 500);
                for bbg = 1:length(BG{bg})
                    PN{bg} = PN{bg} + Pow_non{BG{bg}(bbg)};
                end
                PN{sn ,bg} = PN{bg} /bbg;
                
                PS{sn ,bg} = zeros(90 , 500);
                for bbg = 1:length(BG{bg})
                    PS{bg} = PS{bg} + Pow_stim{BG{bg}(bbg)};
                end
                PS{sn , bg} = PS{bg} /bbg;
                
                PP{sn ,bg} = zeros(90 , 500);
                for bbg = 1:length(BG{bg})
                    PP{bg} = PP{bg} + Pow_pres{BG{bg}(bbg)};
                end
                PP{sn ,bg} = PP{bg} /bbg;
            end
        end
        figure('color' , 'white')
        figCount = 1;
        SeqLAb = {'Random' , 'Struct1', 'Struct2', 'Struct3', 'Struct4'};
        for bg = 1:length(BG)
            for sn = 1:5
                subplot(length(BG) , 5 , figCount)
                contourf([1:500],frex,PS{sn,bg},60,'linecolor','none')
%                 caxis([-80 70])
                colorbar
                ylabel('Frequency(Hz)')
                xlabel('Norm Time')
                title ([SeqLAb{sn} , ', Block ' , num2str(BN(bg))])
                for lin = 1:length(NEM{sn})
                    line([NEM{sn+1}(lin) NEM{sn+1}(lin)] , [2 150] , 'color' , 'red' , 'LineStyle' , ':' , 'LineWidth' , 3)
                end
                figCount = figCount + 1;
            end
        end
        figure('color' , 'white')
        figCount = 1;
        SeqLAb = {'Random' , 'Struct1', 'Struct2', 'Struct3', 'Struct4'};
        for bg = 1:length(BG)
            for sn = 1:5
                subplot(length(BG) , 5 , figCount)
                clear Bands T
                for b =1:length(bandid)
                    Bands(b, :) =  nanmean(PS{sn,bg}(bandid{b}(1) : bandid{b}(2),:));
                    plot([1:500] , Bands(b, :)+30*b , 'LineWidth' , 3)
                    T(b) = Bands(b, 1)+30*b;
                    hold on
                end
                set(gca , 'XLim' , [1 500] , 'YLim' , [0,270],'FontSize' , 16 , 'Box' , 'off')
                set(gca , 'YTickLabels' , bandsLab, 'YTick' , T);
                title ([SeqLAb{sn} , ', Block ' , num2str(BN(bg))])
                for lin = 1:length(NEM{sn})
                    line([NEM{sn+1}(lin) NEM{sn+1}(lin)] , [0 300] , 'color' , 'red' , 'LineStyle' , ':' , 'LineWidth' , 3)
                end
                line([100 100] , [220 270],'color' , 'black' , 'LineWidth' , 5)
                text(60,245,'50','FontSize' , 16 )
                xlabel('Norm Time')
%                 ylabel('Baseline-Normalized Power')
                figCount = figCount + 1;
            end
        end
        

        
        %% Time normalizing each trial separately
    case 'Self'
        
        alltrials = 1;
        
        % definte convolution parameters
        n_wavelet            = length(time);% sampling frequency of the wavelet
        DSR = 4;
        NormSamp = 500;
        for i = 1:length(BN )
            Fs = 1024;
            fName = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/iEEG_P1_Sep14_BlockParsed/EEG_BN' , num2str(i) , '.mat'];
            load(fName);
            eval(['A = EEG_BN' , num2str(i) , '(GoodElec , :);']);
            for i2 = 1:size(A,1)-1
                filt_A(i2 , :) = filter(b,a , A(i2,:));
                A(i2 , :) = filter(c,d , filt_A(i2,:));
            end
            clear filt_A
            eval(['marker = EEG_BN' , num2str(i) , '(MarkerElec , :);']);
            eval(['clear ' , 'EEG_BN' , num2str(i)])
            marker = [0 diff(marker <-2*10^6)];
            start_tr = find(marker == 1); % starts of trials
            end_tr = find(marker == -1);  % ends of trials
            n_data               = size(A , 2);%...EEG.pnts*EEG.trials;
            RawEEGpower = zeros(size(A ,1) - 1 , num_frex,n_data); % channels X frequencies X time X trials
            for ch = 1:size(A,1) - 1
                data                 = A(ch,:);
                n_convolution        = n_wavelet+n_data-1;
                n_conv_pow2          = pow2(nextpow2(n_convolution));
                half_of_wavelet_size = (n_wavelet-1)/2;
                % get FFT of data
                eegfft = fft(data,n_conv_pow2);
                
                % loop through frequencies and compute synchronization
                for fi=1:num_frex
                    
                    wavelet = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2 );
                    
                    % convolution
                    eegconv = ifft(wavelet.*eegfft);
                    eegconv = eegconv(1:n_convolution);
                    eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
                    
                    RawEEGpower(ch,fi,:) = mean(abs(reshape(eegconv,n_data,1)).^2,2);
                end
            end
            Fs_n = floor(Fs/DSR);
            for j = 1 : length(start_tr)
                
                clear Pow_Norm_stim Pow_Norm_pres Pow_non_norm
                if j>2
                    PEEG       = RawEEGpower(:,: , start_tr(j)-floor(Fs/2) : end_tr(j)+floor(Fs/2));
                    trialTime  = 0:1/Fs:(size(PEEG  ,3)/Fs) - 1/Fs;
                    Pall.EventMarker{j ,1} = zeros(size( trialTime));
                    Pall.NormEventMarker{j ,1} = zeros(NormSamp,1);
                    for pr = 1:Dall.seqlength(alltrials)
                        numTaskSamp = floor(( Dall.AllPressIdx(alltrials , Dall.seqlength(alltrials)) - Dall.AllPressIdx(alltrials , pr)) * (Fs/500));
                        idx = numTaskSamp + Fs/2;
                        Nidx = floor(idx*NormSamp/length(trialTime));
                        Pall.EventMarker{j ,1}(end - idx) = pr;  % First press
                        Pall.NormEventMarker{j ,1}(end - Nidx) = pr;
                    end
                    idx = Fs/2;
                    Pall.EventMarker{j ,1}(idx) = -1;       % stimulus apears
                    Nidx = floor(idx*NormSamp/length(trialTime));
                    Pall.NormEventMarker{j ,1}(Nidx) = -1;
                    
                    for ch = 1:size(A,1) - 1
                        for fi=1:num_frex
                            %                     Pow_Norm_stim(ch,fi,:) = squeeze(10*log10(PEEG(ch,fi,:)./mean(PEEG(ch,fi,1:find(Pall.EventMarker{j ,1} == -1)))));
                            %                     Pow_Norm_pres(ch,fi,:) = squeeze(10*log10(PEEG(ch,fi,:)./mean(PEEG(ch,fi,1:find(Pall.EventMarker{j ,1} == 1)))));
                            Pow_non_norm(ch,fi,:)  = squeeze(10*log10(PEEG(ch,fi,:)));
                            
                            idd   = linspace(1 , size(Pow_Norm_stim , 3) , NormSamp);
                            %                     Pall.PEEG_stim{j ,1}(ch , fi , :) = interp1([1:size(Pow_Norm_stim  , 3)] , squeeze(Pow_Norm_stim(ch,fi,:)) , idd);
                            %                     Pall.PEEG_pres{j ,1}(ch , fi , :) = interp1([1:size(Pow_Norm_pres  , 3)] , squeeze(Pow_Norm_pres(ch,fi,:)) , idd);
                            Pall.PEEG_none{j ,1}(ch , fi , :) = interp1([1:size(Pow_non_norm  , 3)] , squeeze(Pow_non_norm(ch,fi,:)) , idd);
                        end
                    end
                else
                    PEEG       = RawEEGpower(:,: , start_tr(j)-floor(Fs/2) : end_tr(j)+floor(Fs/2));
                    trialTime  = 0:1/Fs:(size(PEEG  ,3)/Fs) - 1/Fs;
                    
                    Pall.EventMarker{j ,1} = zeros(size( trialTime));
                    Pall.NormEventMarker{j ,1} = zeros(NormSamp,1);
                    idx = Fs/2;
                    Pall.EventMarker{j ,1}(idx) = -1;       % equivalent to stimulus apears
                    Nidx = floor(idx*NormSamp/length(trialTime));
                    Pall.NormEventMarker{j ,1}(Nidx) = -1;
                    
                    idx = Fs/2;
                    Nidx = floor(idx*NormSamp/length(trialTime));
                    Pall.EventMarker{j ,1}(end - idx) = 7;  % equivalent to last press
                    Pall.NormEventMarker{j ,1}(end - Nidx) = 7;
                    
                    
                    for ch = 1:size(A,1) - 1
                        for fi=1:num_frex
                            %                     Pow_Norm_stim(ch,fi,:) = squeeze(10*log10(PEEG(ch,fi,:)./mean(PEEG(ch,fi,1:find(Pall.EventMarker{j ,1} == -1)))));
                            Pow_non_norm(ch,fi,:)  = squeeze(10*log10(PEEG(ch,fi,:)));
                            
                            
                            idd   = linspace(1 , size(Pow_Norm_stim , 3) , NormSamp);
                            Pall.PEEG_stim{j ,1}(ch , fi , :) = interp1([1:size(Pow_Norm_stim  , 3)] , squeeze(Pow_Norm_stim(ch,fi,:)) , idd);
                            Pall.PEEG_none{j ,1}(ch , fi , :) = interp1([1:size(Pow_non_norm  , 3)] , squeeze(Pow_non_norm(ch,fi,:)) , idd);
                        end
                    end
                    %             Pall.PEEG_pres{j ,1} = NaN;
                    
                end
                [i alltrials]
                alltrials = alltrials  +1;
            end
            
            fname1 = ['RawPSD_b',num2str(i),'_P1_Sep14.mat'];
            fname1
            fname2 = 'Pall';
            sname  = '-v7.3';
            eval('save(fname1 , fname2 , sname);')
            clear Pall
            
            
            
        end
end

% %%
% figure('color' , 'white')
%
% tn  =3;
% figCount = 1;
% for sn = 1:length(SN)
%     subplot(length(SN) , 3 , figCount)
%     contourf([1:NormSamp],frex,squeeze(nanmean(Pall.PEEG_none{tn} , 1)),60,'linecolor','none')
%     set(gca,'yscale','log')
%     title (['Structure ' , num2str(SN(sn)) , '  Non-normalized absolute power'])
%     xlabel('Norm Time')
%     figCount = figCount + 1;
%
%     subplot(length(SN) , 3 , figCount)
%     contourf([1:NormSamp],frex,squeeze(nanmean(Pall.PEEG_stim{tn} , 1)),60,'linecolor','none')
%     set(gca,'yscale','log')
%     title (['Structure ' , num2str(SN(sn)) , '  Norm : 500 ms before the stimulus - stimulus'])
%     xlabel('Norm Time')
%     figCount = figCount + 1;
%
%     subplot(length(SN) , 3 , figCount)
%     contourf([1:NormSamp],frex,squeeze(nanmean(Pall.PEEG_pres{tn} , 1)),60,'linecolor','none')
%     set(gca,'yscale','log')
%     title (['Structure ' , num2str(SN(sn)) , '  Norm : 500 ms before the stimulus - the first press'])
%     xlabel('Norm Time')
%     figCount = figCount + 1;
% end
