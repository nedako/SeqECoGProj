function [PSD , bandsLab] = secog_waveletPSD(EEG , Fs , varargin)
%% psd = secog_waveletPSD(EEG , Fs , varargin)
% EEG is a single channle EEG recording
% Fs is the sampling frequency
% PSD is a numFreqBins by length(EEG) tensor
% Neda Kordjazi
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
        case {'AvgOverBand'}
            % 0 or 1 meaning you want every indudual band or averaged over each band respectively Default = 1
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        otherwise
            error(sprintf('Unknown option: %s',varargin{c}));
    end
end

if ~exist('numCycles')
    numCycles = 10;
end
if ~exist('FreqRange')
    FreqRange = [2 150];
end
if ~exist('numFreqBins')
    numFreqBins = 90;
end
if ~exist('numTrials')
    numTrials = 1000;
end
if ~exist('DownsampleRate')
    DownsampleRate = 1;
end
if ~exist('AvgOverBand')
    AvgOverBand = 1;
end


min_freq =  FreqRange(1);
max_freq = FreqRange(2);
frex = logspace(log10(min_freq),log10(max_freq),numFreqBins);
time = -1:1/Fs:1;
s    = numCycles./(2*pi*frex); % with 10 cycles
n_wavelet            = length(time);% sampling frequency of the wavelet
n_data = length(EEG); %...EEG.pnts*EEG.trials;
n_convolution        = n_wavelet+n_data-1;
n_conv_pow2          = pow2(nextpow2(n_convolution));
half_of_wavelet_size = (n_wavelet-1)/2;
% get FFT of data
eegfft = fft(EEG,n_conv_pow2);
bandsLab = {'Delta <4Hz' , 'Theta 4-8Hz' , 'Alpha 8-13Hz' , 'L-Beta 13-24Hz' , 'H-Beta 24-36Hz' , 'L-Gamma 36-48Hz' , 'H-Gamma >48Hz'};
bands = {[0 4], [4 8] [8 13] [13 24] [14 36] [36 48] [48 110]};
for b = 1:length(bands)
    bandid{b} = [find(frex>bands{b}(1) ,1, 'first') , find(frex<bands{b}(2) ,1, 'last')];
end

for fi=1:numFreqBins
    wavelet = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2 );
    
    % convolution
    eegconv = ifft(wavelet.*eegfft);
    eegconv = eegconv(1:n_convolution);
    eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
    % assign a separate array to every chanle since for longer recordings the array size exceds the Matlab limit
    PSD(fi,:) = downsample(mean(abs(reshape(eegconv,n_data,1)).^2,2)' , DownsampleRate);
end

if AvgOverBand
    temp = PSD;
    clear PSD
    for b =1:length(bandid)
        PSD(b, :) =  nanmean(temp(bandid{b}(1) : bandid{b}(2),:));
    end
else
    bandsLab = [];
end