function secog_freqAna(Dall)
%% Frequency Analysis
base = '/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/';
Fs = 1024;
% definitions, selections...

min_freq =  2;
max_freq = 100;
num_frex = 40;

% define wavelet parameters
time = -1:1/Fs:1;
frex = logspace(log10(min_freq),log10(max_freq),num_frex);
%s    = logspace(log10(3),log10(10),num_frex)./(2*pi*frex);
% s    =  3./(2*pi*frex); % with 3 cycles
s    = 10./(2*pi*frex); % with 10 cycles

% definte convolution parameters
n_wavelet            = length(time);% sampling frequency of the wavelet
% D_interest = getrow(Dall , ~ismember(Dall.TN ,[1 2]));
BN = unique(Dall.BN);
for bn = 1:length(BN)
    Dall1 = getrow(Dall , Dall.BN == BN(bn));
    for tn =1:length(Dall1.TN)
        n_data               = size(Dall1.fEEG{tn,1} , 2);%...EEG.pnts*EEG.trials;
        Pall.RawEEGpower{tn,1} = zeros(size(Dall1.fEEG{tn} ,1) - 1 , num_frex,n_data); % channels X frequencies X time X trials
        
        for ch = 1:size(Dall1.fEEG{tn} ,1) - 1
            data                 = Dall1.fEEG{tn,1}(ch,:);
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
                
                Pall.RawEEGpower{tn,1}(ch,fi,:) = mean(abs(reshape(eegconv,n_data,1)).^2,2);
            end
            disp(['Trial ' , num2str(tn) , ' , Channel ' , num2str(ch)])
        end
        
    end
    fname1 = ['PSD_b',num2str(bn),'_P1_Sep14.mat'];
    fname2 = 'Pall';
    sname  = '-v7.3';
    eval('save(fname1 , fname2 , sname);')
    
    
    %%% Time  normalize
    NormSamp = 500;
    fname1 = ['PSD_b',num2str(bn),'_P1_Sep14.mat'];
    load([base , fname1])
    for tn = 1 : length(Dall1.TN)
        clear RawEEGpower
        for ch = 1:size(Dall1.fEEG{tn} ,1) - 1
            for fi=1:num_frex
               RawEEGpower{tn,1}(ch,fi,:) = downsample(Pall.RawEEGpower{tn,1}(ch,fi,:) , 5);
            end
        end
        Fs_d = Fs/5;
        idd   = linspace(1 , size(RawEEGpower{tn,1} , 3) , NormSamp);
        AbsTime = 0:1/Fs_d:(length(RawEEGpower{tn ,1})/Fs_d)-(1/Fs_d);
        StimeTime = AbsTime - AbsTime(floor(Fs_d/5));
        tp = floor(Fs_d/5) + floor(Fs_d*(Dall1.pressTime0(tn)/1000));
        PressTime  = AbsTime-(tp/Fs_d);
        
        baseidxStim = dsearchn(StimeTime',[min(PressTime) 0]');
        baseidxPress = dsearchn(PressTime',[min(PressTime) 0]');
        
        for ch = 1:size(Dall1.fEEG{tn} ,1) - 1
            for fi=1:num_frex
                % normalize to the time before stimullus
                Pow_Norm_stim(ch,fi,:) = squeeze(10*log10(RawEEGpower{tn,1}(ch,fi,:)./mean(RawEEGpower{tn,1}(ch,fi,baseidxStim(1):baseidxStim(2)))));
                
                P.Pow_Norm_stim{tn}(ch,fi,:) = interp1([1:size(RawEEGpower{tn,1} , 3)] , squeeze(Pow_Norm_stim(ch,fi,:)) , idd);
            end
        end
        
        for ch = 1:size(Dall1.fEEG{tn} ,1) - 1
            for fi=1:num_frex
                % normalize to the time before movement
                Pow_Norm_press(ch,fi,:) = 10*log10(RawEEGpower{tn,1}(ch,fi,:)./mean(RawEEGpower{tn,1}(ch,fi,baseidxPress(1):baseidxPress(2))));
                
                P.Pow_Norm_press{tn}(ch,fi,:) = interp1([1:size(RawEEGpower{tn,1} , 3)] , squeeze(Pow_Norm_press(ch,fi,:)) , idd);
            end
        end
        
        for ch = 1:size(Dall1.fEEG{tn} ,1) - 1
            
            for fi=1:num_frex
                % No normalization
                Pow_absolute(ch,fi,:) = 10*log10(RawEEGpower{tn,1}(ch,fi,:));
                
                P.Pow_absolute{tn}(ch,fi,:) = interp1([1:size(RawEEGpower{tn,1} , 3)] , squeeze(Pow_absolute(ch,fi,:)) , idd);
            end
            
        end
        
        clear Pow_Norm_stim Pow_Norm_press Pow_absolute 
        disp(['Trial ' , num2str(tn) , ' , Block ' , num2str(bn)])
    end
    fname1 = ['PSD_Norm_b',num2str(bn),'_P1_Sep14.mat'];
    fname2 = 'P';
    sname  = '-v7.3';
    eval('save(fname1 , fname2 , sname);')
    clear  P Pall
 
end

%% power visualization for an example trial
prompt = 'Which block?';
bn = input(prompt);

prompt = 'Which trial to plot as example?';
tn = input(prompt);

load(['PSD_b' ,num2str(bn), '_P1_Sep14.mat'])
load('AllData_P1_Sep14.mat')
Dall1 = getrow(Dall , Dall.BN == bn);



% indecies of the baseline for normalization
AbsTime = 0:1/Fs:(length(Dall1.timeStim{tn ,1})/Fs)-(1/Fs);
StimeTime = AbsTime - AbsTime(floor(Fs/5));
tp = floor(Fs/5) + floor(Fs*(Dall1.pressTime0(tn)/1000));
PressTime  = AbsTime-(tp/Fs);


baseidxStim = dsearchn(StimeTime',[min(PressTime) 0]');
baseidxPress = dsearchn(PressTime',[min(PressTime) 0]');

for ch = 1:size(Dall1.fEEG{tn} ,1) - 1
    for fi=1:num_frex
        % normalize to the time before stimullus
        Pow_Norm_stim(ch,fi,:) = 10*log10(Pall.RawEEGpower{tn,1}(ch,fi,:)./mean(Pall.RawEEGpower{tn,1}(ch,fi,baseidxStim(1):baseidxStim(2))));
        
        % normalize to the time before movement
        Pow_Norm_press(ch,fi,:) = 10*log10(Pall.RawEEGpower{tn,1}(ch,fi,:)./mean(Pall.RawEEGpower{tn,1}(ch,fi,baseidxPress(1):baseidxPress(2))));
        
        % No normalization
        Pow_absolute(ch,fi,:) = 10*log10(Pall.RawEEGpower{tn,1}(ch,fi,:));
    end
end

figure('color' , 'white')
subplot(131)
contourf(AbsTime,frex,squeeze(nanmean(Pow_absolute , 1)),60,'linecolor','none')
title ('Non-normalized absolute power')
xlabel('sec')

subplot(132)
contourf(StimeTime,frex,squeeze(nanmean(Pow_Norm_stim , 1)),60,'linecolor','none')
title ('Normalized to baseline : 200 ms before the stimulus comes on till the stimulus')
xlabel('sec')

subplot(133)
contourf(PressTime,frex,squeeze(nanmean(Pow_Norm_press , 1)),60,'linecolor','none')
title ('Normalized to baseline : 200 ms before the stimulus comes on till the first press')
xlabel('sec')


%% Visualize block
prompt = 'Which block?';
bn = input(prompt);
D = getrow(Dall , ismember(Dall.BN , BN));
a = find(~D.isError);


BN   = [14 15 16];
SN = [0:4];
for bn = 1:length(BN)
    fname = ['PSD_Norm_b',num2str(BN(bn)),'_P1_Sep14.mat'];
    load([base , fname])
    Pow(bn) = P;
    clear P
end
 t = size(Pow(1).Pow_Norm_stim{1} , 3);
clear Pow_Norm_stim Pow_Norm_press Pow_absolute

min_freq =  2;
max_freq = 100;
num_frex = 40;

frex = logspace(log10(min_freq),log10(max_freq),num_frex);
for sn = 1:length(SN)
    count = 0;
    Pow_Norm_stim{sn} = zeros(size(Pow(bn).Pow_Norm_stim{1}));
    Pow_Norm_press{sn} = zeros(size(Pow(bn).Pow_Norm_press{1}));
    Pow_absolute{sn} = zeros(size(Pow(bn).Pow_absolute{1}));
    
    for bn = 1:length(BN)
        D = getrow(Dall , ismember(Dall.BN , BN(bn)));
        a = find(~D.isError & D.seqNumb == SN(sn));
        count = count + length(a);
        
        for tn = 1:length(a)
            Pow_Norm_stim{sn} = Pow_Norm_stim{sn}+Pow(bn).Pow_Norm_stim{a(tn)};
        end
        
        
        
        for tn = 1:length(a)
            Pow_Norm_press{sn} = Pow_Norm_press{sn}+Pow(bn).Pow_Norm_stim{a(tn)};
        end
        
        
        for tn = 1:length(a)
            Pow_absolute{sn} = Pow_absolute{sn}+Pow(bn).Pow_absolute{a(tn)};
        end
        
    end
    Pow_Norm_stim{sn} = Pow_Norm_stim{sn}/count;
    Pow_Norm_press{sn} = Pow_Norm_press{sn}/count;
    Pow_absolute{sn} = Pow_absolute{sn}/count;
    
end


figure('color' , 'white')
figCount = 1;
for sn = 1:length(SN)  
    subplot(length(SN) , 3 , figCount)
    contourf([1:t],frex,squeeze(nanmean(Pow_absolute{sn} , 1)),60,'linecolor','none')
    title (['Structure ' , num2str(SN(sn)) , '  Non-normalized absolute power'])
    xlabel('Norm Time')
    figCount = figCount + 1;
    
    subplot(length(SN) , 3 , figCount)
    contourf([1:t],frex,squeeze(nanmean(Pow_Norm_stim{sn} , 1)),60,'linecolor','none')
    title (['Structure ' , num2str(SN(sn)) , '  Norm : 200 ms before the stimulus - stimulus'])
    xlabel('Norm Time')
    figCount = figCount + 1;
    
    subplot(length(SN) , 3 , figCount)
    contourf([1:t],frex,squeeze(nanmean(Pow_Norm_press{sn} , 1)),60,'linecolor','none')
    title (['Structure ' , num2str(SN(sn)) , '  Norm : 200 ms before the stimulus - the first press'])
    xlabel('Norm Time')
    figCount = figCount + 1;
end


figure('color' , 'white')
figCount = 1;
for sn = 1:length(SN)
    subplot(length(SN) , 3 , figCount)
    contourf([1:t],frex,squeeze(nanmean(Pow_absolute{sn} , 1)),60,'linecolor','none')
    set(gca,'yscale','log')
    title (['Structure ' , num2str(SN(sn)) , '  Non-normalized absolute power'])
    xlabel('Norm Time')
    figCount = figCount + 1;
    
    subplot(length(SN) , 3 , figCount)
    contourf([1:t],frex,squeeze(nanmean(Pow_Norm_stim{sn} , 1)),60,'linecolor','none')
    set(gca,'yscale','log')
    title (['Structure ' , num2str(SN(sn)) , '  Norm : 200 ms before the stimulus - stimulus'])
    xlabel('Norm Time')
    figCount = figCount + 1;
    
    subplot(length(SN) , 3 , figCount)
    contourf([1:t],frex,squeeze(nanmean(Pow_Norm_press{sn} , 1)),60,'linecolor','none')
    set(gca,'yscale','log')
    title (['Structure ' , num2str(SN(sn)) , '  Norm : 200 ms before the stimulus - the first press'])
    xlabel('Norm Time')
    figCount = figCount + 1;
end

%% Time Worping
BN = unique(Dall.BN);
for bn = 12:length(BN)
   fname1 = ['PSD_Norm_b',num2str(bn),'_P1_Sep14.mat'];
   load(fname1)
   Pow.p{bn} = P.Pow_Norm_stim;
end
%% For this particular subject the timer was running frm the time that the START_TRIAL came on
for tn = 1:length(Dall.TN)
   if ~ismember(Dall.TN(tn) ,[1 2])
       behav_Len(tn) = Dall.AllPressIdx(tn , Dall.seqlength(tn)) + 1;
       eeg_Len(tn) = length(Dall.fEEG{tn}) - (1024 + 0.001*Dall.AllPressTimes(tn , 1)*1024); % length of the task minus the one second - reaction time
   else
       behav_Len(tn) = 15*500; % trialatime-out times the smapling frequency
       eeg_Len(tn) = length(Dall.fEEG{tn});
   end
end

figure
subplot(3,1,1)
plot((eeg_Len/1024))
ylabel('sec')
title('Trial length as per EEG')

subplot(3,1,2)
plot((eeg_Len/1024))
ylabel('sec')
title('Trial length as per behavior')

subplot(3,1,3)
a = (behav_Len/500)-(eeg_Len/1024) ;
plot(a)
ylabel('sec')
title('Length of behavioral - length of EEG')

