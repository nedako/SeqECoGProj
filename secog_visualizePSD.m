function secog_visualizePSD(Pall , subjNum, what)
% Pall needs to be the structure containing time normalized, average
% patterned PSDs so the output of Pall  = secog_parseEEG_PSD('TimeWarpPSD' , Dall, subjNum);
while(c<=length(varargin))
    switch(varargin{c})
        case {'NumWarpSamp'}
            % number of sine cycles to use in the Morlet wavelet
            % default 10
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
         case {'Chan2Plot'}
            % Channel(s) to plot
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        otherwise
            error(sprintf('Unknown option: %s',varargin{c}));
    end
end

if ~exist('NumWarpSamp')
    NumWarpSamp = 300;
end


subjname = {'P2'};
mainDir = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' subjname{subjNum} , '/'] ;
cd(mainDir)
load('AveragePatternMarker.mat')
load('secog_all.mat')
load('ChanLabels.mat')

bandsLab = {'Delta <4Hz' , 'Theta 4-8Hz' , 'Alpha 8-13Hz' , 'L-Beta 13-24Hz' , 'H-Beta 24-36Hz' , 'L-Gamma 36-48Hz' , 'H-Gamma >48Hz'};
E = getrow(E , strcmp(E.blockGroupNames , what));
if isempty(E.EM)
    error('Nothing to plot!')
end
figure('color' , 'white')
figCount = 1;
for sn = 1:length(E.SN{1})
    NEM = E.NEM{1}(sn , :);
    id = ismember(Dall.BN , E.blockGroups{1}) & ismember(Dall.seqNumb , E.SN{1}(sn));
    F = getrow(Pall , id);
    % sum the warped PSDs insife the F structure
    tcount = 1;
    for tn = 1:length(F.PSD_stim)
        if isequal(size(F.PSD_stim{tn}) ,[length(ChanLabels) , length(bandsLab) , NumWarpSamp])
            tempPow(tcount , :,:,:) = F.PSD_stim{tn};
            tcount = tcount +1;
        end
    end
    AvgPow{sn} = squeeze(nanmean(tempPow , 1));
    
    for ch = 1:length(Chan2Plot)
        subplot(length(E.SN{1}) ,length(Chan2Plot), figCount)
        for b =1:length(bandsLab)
            plot([1:NumWarpSamp] , squeeze(AvgPow{sn}(ch , b , :))+7*b , 'LineWidth' , 3)
            hold on
        end
        title (['Average Time-Warped PSD for ' ,E.blockGroupNames{1}, ' Blocks'])
        for lin = 1:length(NEM)
            line([NEM(lin) NEM(lin)] , [0 60] , 'color' , 'red' , 'LineStyle' , ':' , 'LineWidth' , 3)
        end
        xlabel('Norm Time')
%         if figCount == 1
%             set(gca ,'YTickLabels' , bandsLab, 'YTick' , T );
%         else
%             set(gca ,'YTickLabels' , [], 'YTick' , [] );
%         end
        line([100 100] , [55 60],'color' , 'black' , 'LineWidth' , 5)
        text(60,57.5,'5','FontSize' , 16 )
        set(gca , 'XLim' , [1 500] , 'YLim' , [0,60],'FontSize' , 16,'Box' , 'off')
        
    end
    figCount = figCount + 1;
end


