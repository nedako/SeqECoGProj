function [PSD] = seqeeg_waveletPSD(EEG , Fs , varargin)
%% psd = seqeeg_waveletPSD(EEG , Fs , varargin)
% EEG is a single channle EEG recording
% Fs is the sampling frequency
% PSD is a numFreqBins by length(EEG) tensor
% Neda Kordjazi
c = 1;

numCycles = 7;
FreqRange = [2 180];
numFreqBins = 90;
DownsampleRate = 10;
AvgOverBand = 0;



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




min_freq =  FreqRange(1);
max_freq = FreqRange(2);
frex = linspace(min_freq, max_freq,numFreqBins);
time = -1:1/Fs:1;
s    = numCycles./(2*pi*frex); % with 10 cycles
n_wavelet            = length(time);% sampling frequency of the wavelet
n_data = length(EEG); %...EEG.pnts*EEG.trials;
n_convolution        = n_wavelet+n_data-1;
n_conv_pow2          = pow2(nextpow2(n_convolution));
half_of_wavelet_size = (n_wavelet-1)/2;
% get FFT of data
eegfft = fft(EEG,n_conv_pow2);

%% visualize wavelet
numCycles1 = 10;
s1    = numCycles1./(2*pi*frex); % with 10 cycles
plotWave = 0;

if plotWave
    figure('color' , 'white')
    for fi=1:numFreqBins
        subplot(211)
        wave = sqrt(1/(s1(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s1(fi)^2)));
        plot(real(wave) , 'LineWidth' , 2);
        hold on
        plot(imag(wave) , 'LineWidth' , 2);
        hold off
        legend({'Real' , 'Imaginary'})
        set(gca , 'FontSize' , 20 , 'Box' , 'off')
        title('Wavelet')
        
        subplot(212)
        wavelet = fft(wave);
        L = length(wave);
        P2 = wavelet/L;
        f = Fs*(0:L-1)/L;
        plot(f(1:FreqRange(2)) , real(P2(1:FreqRange(2))) , 'LineWidth' , 2);
        hold on
        plot(f(1:FreqRange(2)) , imag(P2(1:FreqRange(2))) , 'LineWidth' , 2);
        plot(f(1:FreqRange(2)) , abs(P2(1:FreqRange(2))), 'LineWidth' , 3 , 'color' , 'k')
        set(gca , 'XTick' , [5:5:90] , 'XTickLabel' , [10:10:180] , 'FontSize' , 20 , 'Box' , 'off')
        hold off
        title('FFT of the wavelet')
        pause()
        drawnow
    end
end

%%

for fi=1:numFreqBins
    wavelet = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2 );
    
    % convolution
    eegconv = ifft(wavelet.*eegfft);
    eegconv = eegconv(1:n_convolution);
    eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
    % assign a separate array to every chanle since for longer recordings the array size exceds the Matlab limit
    if  sum(strcmp(PoworDec , {'p' , 'P'}))
        temp = 10*log10(abs(nanmean(reshape(eegconv,n_data,1).^2,2)').^2);
    else
        temp = nanmean(reshape(eegconv,n_data,1).^2,2)';
    end
    dataLen = length(temp)/NumChans;
    temp = reshape(temp' ,dataLen,NumChans)';
    % downsample each channel separately, cuz done together the number of samples wont match
    for ch = 1:NumChans
        tempPSD = [tempPSD , downsample(temp(ch , :) , DownsampleRate)];
    end
    PSD(fi,:) = tempPSD;
  
end
