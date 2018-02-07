% complex Morlet Wavelet
N = 1000;

% this is the Gaussian tapering horizon -->  effective support
Lb = -8;
Ub = 8;

% The higher the fc, the more cycles you get. this defines the width of the Gaussian Kernel
% Therefore, increasing Fb results in a narrower concentration of energy around the center frequency.


fb = 8;
fc = 1;
[psi,x] = cmorwavf(Lb,Ub,N,fb,fc);
subplot(311)
plot(x,real(psi)); title('Real Part');
subplot(312)
plot(x,imag(psi)); title('Imaginary Part');

subplot(313)
f = Lb:.01:Ub;
psihat = exp(-pi^2*fb*(f-fc).^2);  
plot(f,psihat)   
title(['fb = ' , num2str(fb)])

%% create and add a wavelet to the library 

% Create the pattern signal
t    = linspace(0, 1, 1024);
y    = (1./(1+(4*(t-0.5)).^4)).*cos(20*(t-0.5));
% Create a polynomial fit of the pattern signal
[psi,tval,nc] = pat2cwav(y, 'polynomial', 10, 'continuous') ;  
% Plot the pattern compared to the polynomial fit
plot(t,y,'-',tval,nc*psi,'--'),  
title('Original Pattern and Adapted Wavelet (dashed line)')
save('mother', 'nc', 'psi', 'tval');
wavemngr('del','CosWave');
wavemngr('add','CosWave','cosw',4,'','mother'); 

%% 