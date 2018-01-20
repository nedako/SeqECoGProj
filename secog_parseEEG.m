function Dout  = secog_parseEEG(Dall , GoodElec)


%% look for channels that are not EEG
% load('EEG_BN2.mat')
% A = EEG_BN2;
% marker = A(141,:);
% marker = [0 diff(marker <-2*10^6)];
% start_tr = find(marker == 1);
% end_tr = find(marker == -1);
% figure
% for i =20:40 ... size(EEG_BN1 , 1)
%     plot(A(i,:))
%     hold on
%     for j = 1:length(start_tr)
%         line([start_tr(j) start_tr(j)] , [min(A(i,:)) max(A(i,:))] , 'color' , 'green', 'LineWidth' , 3)
%         line([end_tr(j) end_tr(j)] , [min(A(i,:)) max(A(i,:))] , 'color' , 'red', 'LineWidth' , 3)
%     end
%     title(num2str(i))
%     hold off
%     pause()
%     drawnow
% end
% % Visually inspect the channels and identify good channels.
% GoodElec = [1:30,32:60,63:74,77:80 , 141];


%% plot to make sure the boundries are all in the right position
% figure
% plot(Data.values(141,:))
% MIN = min(Data.values(141,:));
% MAX = max(Data.values(141,:));
% hold on
% for i = 1:length(BN)
%     line([BN(i , 2) BN(i , 2)] , [MIN , MAX] , 'color' , 'red');
%     line([BN(i , 3) BN(i , 3)] , [MIN , MAX] , 'color' , 'red');
% end

%% Import the start and end ample number of each block from the manually creatd excel sheet

baseDir = '/Users/nedakordjazi/Documents/SeqECoG/';
filename = 'BlockInfo_P1_Sep-14.xlsx';

%% Import the data
[~, ~, raw] = xlsread([baseDir , filename],'Sheet1');
raw = raw(2:end,:);

%% Create output variable
BN = reshape([raw{:}],size(raw));

%% Clear temporary variables
clearvars raw;


%% chop up the EEG data into Blocks and save
% for i = 1:length(BN )
%     statement = ['EEG_BN' , num2str(i) , '= Data.values(:,BN(i,2):BN(i,3));'];
%     eval(statement)
%     save(['EEG_BN' , num2str(i) ,'.mat'] , ['EEG_BN' , num2str(i)] , '-v7.3')
%     eval(['clear EEG_BN' , num2str(i)])
% end

%% %% chop up the EEG data into trials
Fs = 1024;
alltrials = 1;
for i = 1:length(BN )
    fName = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/iEEG_P1_Sep14_BlockParsed/EEG_BN' , num2str(i) , '.mat'];
    load(fName);
    eval(['A = EEG_BN' , num2str(i) , ';']);
    marker = A(141,:);
    marker = [0 diff(marker <-2*10^6)];
    start_tr = find(marker == 1); % starts of trials
    end_tr = find(marker == -1);  % ends of trials
    for j = 1:length(start_tr)
        if j>2
            Dall.EEG{alltrials ,1}       = A(GoodElec , start_tr(j)-floor(Fs/2) : end_tr(j)+floor(Fs/2));
            Dall.timeStim{alltrials ,1}  = 0:1/Fs:(length(Dall.EEG{alltrials ,1})/Fs) - 1/Fs;
            Dall.timeStim{alltrials ,1}  = Dall.timeStim{alltrials ,1}-(500/Fs);
            % convert the time of the forst press from hand to eeg
            tp = floor(Fs/5) + floor(Fs*(Dall.pressTime0(alltrials)/1000));
            Dall.timePress{alltrials ,1} = 0:1/Fs:(length(Dall.EEG{alltrials ,1})/Fs) - 1/Fs;
            Dall.timePress{alltrials ,1}  = Dall.timePress{alltrials ,1}-(tp/Fs);
        else
            Dall.EEG{alltrials ,1} = A(GoodElec , start_tr(j) : end_tr(j));
            Dall.timeStim{alltrials ,1} = 0:1/Fs:(length(Dall.EEG{alltrials ,1})/Fs) - 1/Fs;
            Dall.timePress{alltrials ,1} = 0:1/Fs:(length(Dall.EEG{alltrials ,1})/Fs) - 1/Fs;
        end
        alltrials = alltrials  +1;
        [i alltrials]
    end
    eval(['clear EEG_BN' , num2str(i)])
end
%%
Fs = 1024;
Fo = 60;
% Q : quality factor is the center frequency divided by the bandwidth.
% Q = 35;
% BW = Fo/(Fs/2);
% [b,a] = iircomb(10,BW,'notch');

wo = Fo/(Fs/2);  bw = wo/35; % notch to eliminate 60 Hz
[b,a] = iirnotch(wo,bw);

wo = 2*Fo/(Fs/2);  bw = wo/60;% notch to eliminate the first harmonic
[c,d] = iirnotch(wo,bw);

% Filter Visualization
% h = fvtool(b,a);
% h.Fs = Fs;
% h.FrequencyRange='[-Fs/2, Fs/2)';
% zplane(b,a)
%
% h = fvtool(c,d);
% h.Fs = Fs;
% h.FrequencyRange='[-Fs/2, Fs/2)';
% zplane(c,d)

for tn = 1:length(Dall.TN)
    A = Dall.EEG{tn};
    for i = 1:size(A,1)-1
        filt_A(i , :) = filter(b,a , A(i,:));
        Dall.fEEG{tn,1}(i , :) = filter(c,d , filt_A(i,:));
    end
    clear A filt_A
end

% [filt] = ft_preproc_dftfilter(A, 1024, 60);
A      = Dall.EEG{tn}(10,:);
filt_A = Dall.fEEG{tn}(10,:);
t = 0:(1/Fs):(length(A)/Fs) - (1/Fs);
figure('color' , 'white')
subplot(2,1,1)
periodogram(A(1,:),[],length(A(1,:)),Fs,'power')
subplot(2,1,2)
periodogram(filt_A(1,:),[],length(filt_A(1,:)),Fs,'power')

