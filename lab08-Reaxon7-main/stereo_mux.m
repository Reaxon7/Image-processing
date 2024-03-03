% EE 438 Digital Signal Processing with Applications - Spring 2000
% Digital stereo multiplexer
%---------------------------------------------------------------
% Initialize
clear all;
close all;
%Read and play back data sampled at 8192Hz
Fs = 8192;
%---------------------------------------------------------------
%Read data for left and right channels

data=load('erf1s1t0.txt')';
%plot(data)
%end1=input('end1?');
%end2=input('end2?');
xl=data;
Nl = length(xl)
input('Left channel');
soundsc(xl,Fs);

data=load('ysf1s1t0.txt')';
%plot(data)
%end1=input('end1?');
%end2=input('end2?');
%xr=data(end1:end2);
xr=data;
Nr = length(xr)
input('Right channel');
soundsc(xr,Fs);

%---------------------------------------------------------------
%Bandlimit both channels to half Nyquist rate
N_spect = 128;
del_mu = 1./N_spect;
mu = -0.5:del_mu:0.5-del_mu;
N_filt = 64;
delta = 0.02;
bands = [0.0, 0.5-delta, 0.5+delta, 1.0];
H_ideal = [1.0, 1.0, 0.0, 0.0];
h = remez(N_filt, bands, H_ideal);
%---------------------------------------------------------------
%Prepare plots
figure(1);
H = fft(h,N_spect);
H = fftshift(H);
Hmag = abs(H);
Hpha = angle(H);
subplot(3,1,1), stem(h);
ylabel('Filter');
subplot(3,1,2), plot(mu,Hmag);
ylabel('Magnitude');
subplot(3,1,3), plot(mu,Hpha);
ylabel('Phase');
pause;
%---------------------------------------------------------------
%Filter left and right channels
xlbl = conv(h, xl);
xrbl = conv(h, xr);
%---------------------------------------------------------------
%Prepare plots and playback results
figure(2);
N_spect = pow2(ceil(log2(Nl)))
Xl = fft(xl,N_spect);
Xl = fftshift(Xl);
Xlmag = abs(Xl);
Xlbl = fft(xlbl,N_spect);
Xlbl = fftshift(Xlbl);
Xlblmag = abs(Xlbl);
del_mu = 1./N_spect;
mu = -0.5:del_mu:0.5-del_mu;
subplot(4,1,1), plot(xl);
ylabel('Original');
subplot(4,1,2), plot(mu,Xlmag);
ylabel('Spectrum');
subplot(4,1,3), plot(xlbl);
ylabel('Bandlimited');
subplot(4,1,4), plot(mu,Xlblmag);
ylabel('Spectrum');
input('Left channel');
soundsc(xl,Fs);
input('Left channel bandlimited to 1/2 Nyquist rate');
soundsc(xlbl,Fs);
xrbl = conv(h, xr);
input('Right channel');
soundsc(xr,Fs);
input('Right channel bandlimited to 1/2 Nyquist rate');
soundsc(xrbl,Fs);
%---------------------------------------------------------------
%Pad signals to same length with zeros
Nlbl = length(xlbl);
Nrbl = length(xrbl);
Nlr = max(Nlbl,Nrbl);
xl(1:Nlbl) = xlbl;
xl(Nlbl+1:Nlr) = zeros(1,Nlr-Nlbl);
xr(1:Nrbl) = xrbl;
xr(Nrbl+1:Nlr) = zeros(1,Nlr-Nrbl);
%Get rid of unneeded variables to free some memory
clear xlbl; clear Xlbl; clear Xlblmag; clear mu;
clear xrbl;
clear H; clear Hmag; clear Hpha;
%---------------------------------------------------------------
%Matrix left and right channels to form sum and difference
%components.
xsum = xl+xr;
xdif = xl-xr;
%---------------------------------------------------------------
%Modulate difference channel and add to sum channel
omega_0 = pi;
Nplot = 256;
n = 0:1:Nlr-1;
xmod = xdif.*cos(omega_0.*n);
xmux = xsum+xmod;
%---------------------------------------------------------------
%Prepare plots and playback results
figure(3);
N_spect = pow2(ceil(log2(Nlr)))
Xdif = fft(xdif,N_spect);
Xdif = fftshift(Xdif);
Xdifmag = abs(Xdif);
Xmod = fft(xmod,N_spect);
Xmod = fftshift(Xmod);
Xmodmag = abs(Xmod);
del_mu = 1./N_spect;
mu = -0.5:del_mu:0.5-del_mu;
xdifclip = xdif(round(Nlr./2):round(Nlr./2+Nplot-1));
subplot(4,1,1), plot(xdifclip);
ylabel('Difference');
subplot(4,1,2), plot(mu,Xdifmag);
ylabel('Spectrum');
xmodclip = xmod(round(Nlr./2):round(Nlr./2+Nplot-1));
subplot(4,1,3), plot(xmodclip);
ylabel('Modulated');
subplot(4,1,4), plot(mu,Xmodmag);
ylabel('Spectrum');
input('Difference channel');
soundsc(xdif,Fs);
input('Modulated difference channel');
soundsc(xmod,Fs);
%Get rid of unneeded variables to free some memory
clear xdif; clear xdifclip; clear Xdif; clear Xdifmag;
clear xmod; clear xmodclip; clear Xmod; clear Xmodmag;
Xmux = fft(xmux,N_spect);
Xmux = fftshift(Xmux);
Xmuxmag = abs(Xmux);
xmuxclip = xmux(round(Nlr./2):round(Nlr./2+Nplot-1));
figure(4);
subplot(2,1,1), plot(xmuxclip);
title('Multiplexed signal');
subplot(2,1,2), plot(mu,Xmuxmag);
title('Spectrum of multiplexed signal');
input('Final multiplexed signal');
soundsc(xmux,Fs);
