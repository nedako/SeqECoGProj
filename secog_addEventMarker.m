function Dout  = secog_addEventMarker(Dall, what , GoodElec, MarkerElec)

% for patient 1 GoodElec   = [1:30,32:60,63:74,77:80];
% for patient 1 MarkerElec = 141;

switch what
    case 'ParsEEG'
        
        
        % Import the start and end sample number of each block from the manually creatd excel sheet
        
        baseDir = '/Users/nedakordjazi/Documents/SeqECoG/';
        filename = 'BlockInfo_P1_Sep-14.xlsx';
        
        % Import the data
        [~, ~, raw] = xlsread([baseDir , filename],'Sheet1');
        raw = raw(2:end,:);
        
        % Create output variable
        BN = reshape([raw{:}],size(raw));
        
        % Clear temporary variables
        clearvars raw;
        
        
        chop up the EEG data into Blocks and save
        for i = 1:length(BN )
            statement = ['EEG_BN' , num2str(i) , '= Data.values(:,BN(i,2):BN(i,3));'];
            eval(statement)
            save(['EEG_BN' , num2str(i) ,'.mat'] , ['EEG_BN' , num2str(i)] , '-v7.3')
            eval(['clear EEG_BN' , num2str(i)])
        end
        
        % Load one of the parsed blocks and look for channels that are not EEG to define GoodElec and MarkerElec
        load('EEG_BN2.mat')
        A = EEG_BN2;
        marker = A(141,:);
        marker = [0 diff(marker <-2*10^6)];
        start_tr = find(marker == 1);
        end_tr = find(marker == -1);
        figure
        for i =20:40 ... size(EEG_BN1 , 1)
                plot(A(i,:))
            hold on
            for j = 1:length(start_tr)
                line([start_tr(j) start_tr(j)] , [min(A(i,:)) max(A(i,:))] , 'color' , 'green', 'LineWidth' , 3)
                line([end_tr(j) end_tr(j)] , [min(A(i,:)) max(A(i,:))] , 'color' , 'red', 'LineWidth' , 3)
            end
            title(num2str(i))
            hold off
            pause()
            drawnow
        end
        Visually inspect the channels and identify good channels.
        plot to make sure the boundries are all in the right position
        figure
        plot(Data.values(141,:))
        MIN = min(Data.values(141,:));
        MAX = max(Data.values(141,:));
        hold on
        for i = 1:length(BN)
            line([BN(i , 2) BN(i , 2)] , [MIN , MAX] , 'color' , 'red');
            line([BN(i , 3) BN(i , 3)] , [MIN , MAX] , 'color' , 'red');
        end
        Dout = [];
        
    case 'addEvent'
        %% find the normal event marker and averages
        % definte convolution parameters
        alltrials = 1;
        NormSamp = 500;
        for i = 1:length(BN)
            Fs = 1024;
            fName = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/iEEG_P1_Sep14_BlockParsed/EEG_BN' , num2str(i) , '.mat'];
            load(fName);
            eval(['marker = EEG_BN' , num2str(i) , '(MarkerElec , :);']);
            eval(['clear ' , 'EEG_BN' , num2str(i)])
            marker = [0 diff(marker <-2*10^6)];
            start_tr = find(marker == 1); % starts of trials
            end_tr   = find(marker == -1);  % ends of trials
            for j = 1 : length(start_tr)
                if j>2
                    PEEG       =  start_tr(j)-floor(Fs/2) : end_tr(j)+floor(Fs/2);
                    trialTime  = 0:1/Fs:(length(PEEG)/Fs) - 1/Fs;
                    Tall.EventMarker{alltrials , 1} = zeros(size( trialTime));
                    Tall.NormEventMarker{alltrials , 1} = zeros(NormSamp,1);
                    for pr = 1:Dall.seqlength(alltrials)
                        numTaskSamp = floor(( Dall.AllPressIdx(alltrials , Dall.seqlength(alltrials)) - Dall.AllPressIdx(alltrials , pr)) * (Fs/500));
                        idx = numTaskSamp + Fs/2;
                        Nidx = floor(idx*NormSamp/length(trialTime));
                        Tall.EventMarker{alltrials , 1}(end - idx) = pr;  % First press
                        Tall.NormEventMarker{alltrials , 1}(end - Nidx) = pr;
                    end
                    idx = Fs/2;
                    Tall.EventMarker{alltrials , 1}(idx) = -1;       % stimulus apears
                    Nidx = floor(idx*NormSamp/length(trialTime));
                    Tall.NormEventMarker{alltrials , 1}(Nidx) = -1;
                else
                    PEEG       = start_tr(j)-floor(Fs/2) : end_tr(j)+floor(Fs/2);
                    trialTime  = 0:1/Fs:(length(PEEG)/Fs) - 1/Fs;
                    
                    Tall.EventMarker{alltrials , 1} = zeros(size( trialTime));
                    Tall.NormEventMarker{alltrials , 1} = zeros(NormSamp,1);
                    idx = Fs/2;
                    Tall.EventMarker{alltrials , 1}(idx) = -1;       % equivalent to stimulus apears
                    Nidx = floor(idx*NormSamp/length(trialTime));
                    Tall.NormEventMarker{alltrials , 1}(Nidx) = -1;
                    
                    idx = Fs/2;
                    Nidx = floor(idx*NormSamp/length(trialTime));
                    Tall.EventMarker{alltrials , 1}(end - idx) = 7;  % equivalent to last press
                    Tall.NormEventMarker{alltrials , 1}(end - Nidx) = 7;
                end
                [i alltrials]
                alltrials = alltrials  +1;
            end
        end
        Dall.EventMarker = Tall.EventMarker;
        Dall.NormEventMarker = Tall.NormEventMarker;
        Dout = Dall;
end