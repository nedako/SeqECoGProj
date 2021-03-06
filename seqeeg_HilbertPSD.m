function [PSD] = seqeeg_HilbertPSD(EEG , Fs , varargin)
%% psd = seqeeg_waveletPSD(EEG , Fs , varargin)
% EEG is a single channle EEG recording
% Fs is the sampling frequency
% PSD is a numFreqBins by length(EEG) tensor
% Neda Kordjazi


FreqRange = [4 184];
numFreqBins = 45;
DownsampleRate = 10;
AvgOverBand = 0;


c = 1;
while(c<=length(varargin))
    switch(varargin{c})
        
        case {'numCycles'}
            % number of sine cycles to use in the Morlet wavelet
            % default 10
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
        case {'DownsampleRate'}
            % folds by which you wnat the PSD to be downsampld
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'PoworDec'}
            % 'p' 'P' for power (10log10) 
            % anything else for decompositions
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'NumChans'}
            % number of channels that have been concatenated 
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        otherwise
            error(sprintf('Unknown option: %s',varargin{c}));
    end
end






nyqst = Fs/2;
min_freq =  FreqRange(1);
max_freq = FreqRange(2);
frex = round(linspace(min_freq, max_freq,numFreqBins));
time = -1:1/Fs:1;
n_data = length(EEG); %...EEG.pnts*EEG.trials;
EEG = double(EEG);
for fi=1:numFreqBins
    tempPSD = [];
    freqspread  =  2; % Hz +/- the center frequency
    center_freq = frex(fi);
    transwid    = .15;
    % construct filter kernels
    ffrequencies  = [ 0 (1-transwid)*(center_freq-freqspread) (center_freq-freqspread) (center_freq+freqspread) (1+transwid)*(center_freq+freqspread) nyqst ]/nyqst;
    idealresponse = [ 0 0 1 1 0 0 ];
    filterweights = firls(3*round(Fs/(center_freq-freqspread)),ffrequencies,idealresponse);
    filtdat= filtfilt(filterweights,1,EEG);
    if  sum(strcmp(PoworDec , {'p' , 'P'}))
        temp = 10*log10(abs(nanmean(reshape(filtdat,n_data,1).^2,2)').^2);
    else
        temp = nanmean(reshape(filtdat,n_data,1).^2,2)';
    end
    dataLen = length(temp)/NumChans;
    temp = reshape(temp' ,dataLen,NumChans)';
    % downsample each channel separately, cuz done together the number of samples wont match
    for ch = 1:NumChans
        tempPSD = [tempPSD , downsample(temp(ch , :) , DownsampleRate)];
    end
    PSD(fi,:) = tempPSD;
end




