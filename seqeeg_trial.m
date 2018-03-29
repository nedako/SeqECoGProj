function [D]=seqeeg_trial(MOV,D,fig,varargin)
% D is the data structure including the forces

close all

% ------------------------------------------------
% Defaults
sample=2; % what proportion of 1000 is the sampling frequency
forceThres=1;
minDist = 100/sample;
PH = {};
PT = {};

% ------------------------------------------------
% extract data
if (isempty(MOV))
    return;
end;
state    = MOV(:,1);
realTime = MOV(:,2);
time     = MOV(:,3);
sampfreq = 1000/sample;

% The force recordings

D.F{1}        = (MOV(:,4:8));
D.F{1}        = smooth_kernel(D.F{1},4);


fingers = 5;
peakHeight = [];
peakTime   = [];

% ------------------------ Force peak detection
for i = 1:fingers
    if max(D.F{1}(:,i)) > forceThres
        [PH,PT]=findpeaks(D.F{1}(:,i),'minpeakdistance',minDist,'minpeakheight',forceThres);
        peakHeight          = [peakHeight ; PH];
        peakTime            = [peakTime   ; PT];
        D.numPeaks(i)       = length(PH);
        D.avrgPeakHeight(i) = nanmean(PH);
        i = i+1;
    else
       D.numPeaks(i)      = NaN;
       D.avrgPeakHeight(i)= NaN; 
    end
end

% -------- Find the number of finger presses by counting "responsexx" fields
A = fieldnames(D);
D.AllResponse   = [];
D.AllPress   = []; 
presscnt = 0; 
for i = 1:length(A)
    if length(A{i}) > 8
        if strcmp(A{i}(1:8) , 'response')
            presscnt = presscnt + 1;
            eval(['D.AllResponse = [D.AllResponse D.',A{i} , ']']);
        end
    end
end

presscnt = 0; 
for i = 1:length(A)
    if length(A{i}) > 5
        if strcmp(A{i}(1:5) , 'press') & ~strcmp(A{i}(6) , 'T')
            presscnt = presscnt + 1;
            eval(['D.AllPress = [D.AllPress D.',A{i} , ']']);
        end
    end
end



disp(['Maximum number of presses detected: ' , num2str(presscnt)])
disp('Putting together the press times...')

for prs = 0:presscnt - 1 
    pressName                    = ['D.pressTime' , num2str(prs) ];    
    D.AllPressIdx(:,prs + 1)     = (sampfreq/1000)* eval(pressName);
    D.AllPressTimes(:,prs + 1)   = eval(pressName);   
end

D.AllResponse(~D.AllResponse) = NaN;
D.AllPress(~D.AllPress) = NaN;
D.AllPressTimes(~D.AllPressTimes) = NaN;
D.AllPressIdx(~D.AllPressIdx) = NaN;


% -----------------Display trial
if fig
    figure('color' , 'white')
    subplot(2,1,1)
    % ---------------- Plot forces only during state 5 (WAIT_PRESS)
    %idx = (state == 5);
    idx = 1:length(time);
    plot(time(idx), D.F{1}(idx , :),'LineWidth' , 3)
    hold on
    legend('Finger1','Finger2','Finger3','Finger4','Finger5')
    xlabel('Time(ms)')
    ylabel('Force')
    for i = 1:length(find(~isnan(D.AllPressIdx)))
        line([time(D.AllPressIdx(i)) , time(D.AllPressIdx(i))],[-.1 5],'LineStyle' , ':','LineWidth' , 3 , 'color' , 'r');
    end
    hold on
    ss = unique(state);
    for j = 1:length(ss)
        a = find(state==ss(j) , 1, 'first');
        line([a*sample a*sample] , [-.1 5] , 'color' , [0 0 0],'LineStyle' , ':','LineWidth' , 3 );
    end
    axis([min(time(idx)) max(time(idx)) -.1 5]);
    err = D.isError;
    title(['Block' , num2str(D.BN) , ' - Trial ' , num2str(D.TN) , ' - Error ' , num2str(err)]);
   
    
    % --------------- Plot state
    subplot(2,1,2);
    temp = state(idx);    
    plot(time(idx), temp,'LineWidth' , 3) 
    hold on
    ss = unique(state);
    for j = 1:length(ss)
        a = find(state==ss(j) , 1, 'first');
        line([a*sample a*sample] , [min(temp) max(temp)] , 'color' , [0 0 0],'LineStyle' , ':','LineWidth' , 3 );
    end
    axis([min(time(idx)) max(time(idx)) min(temp) max(temp)]);
    xlabel('Time(ms)')
    ylabel('State')
    
    figure('color' , 'white')
    subplot(2,1,1)
    plot(diff(realTime))
    title('realTime')
    subplot(2,1,2)
    plot(diff(time))
    title('time')
    keyboard; 
end


