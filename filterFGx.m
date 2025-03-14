function [filtdat,empVals,fx] = filterFGx(data,srate,f,fwhm,showplot)
% filterFGx   Narrow-band filter via frequency-domain Gaussian
%  [filtdat,empVals] = filterFGx(data,srate,f,fwhm,showplot)
% 
% 
%    INPUTS
%       data : 1 X time or chans X time
%      srate : sampling rate in Hz
%          f : peak frequency of filter
%       fhwm : standard deviation of filter, 
%              defined as full-width at half-maximum in Hz
%   showplot : set to true to show the frequency-domain filter shape
% 
%    OUTPUTS
%    filtdat : filtered data
%    empVals : the empirical frequency and FWHM (in Hz and in ms)
% 
% Empirical frequency and FWHM depend on the sampling rate and the
% number of time points, and may thus be slightly different from
% the requested values.
% 
% mikexcohen@gmail.com

%% input check

if size(data,1)>size(data,2)
%     help filterFGx
%     error('Check data size')
end

if (f-fwhm)<0
%     help filterFGx
%     error('increase frequency or decrease FWHM')
end

if nargin<4
    help filterFGx
    error('Not enough inputs')
end

if fwhm<=0
    error('FWHM must be greater than 0')
end

if nargin<5
    showplot=false;
end

%% compute filter

% frequencies
hz = linspace(0,srate,size(data,2));

% create Gaussian
s  = fwhm*(2*pi-1)/(4*pi); % normalized width
x  = hz-f;                 % shifted frequencies
fx = exp(-.5*(x/s).^2);    % gaussian
fx = fx./max(fx);          % gain-normalized

%% filter

filtdat = 2*real( ifft( bsxfun(@times,fft(data,[],2),fx) ,[],2) );

%% compute empirical frequency and standard deviation

idx = dsearchn(hz',f);
empVals(1) = hz(idx);

% find values closest to .5 after MINUS before the peak
empVals(2) = hz(idx-1+dsearchn(fx(idx:end)',.5)) - hz(dsearchn(fx(1:idx)',.5));

% also temporal FWHM
tmp = abs(hilbert(real(fftshift(ifft(fx)))));
tmp = tmp./max(tmp);
tx = (0:length(data)-1)/srate;
[~,idxt] = max(tmp);
empVals(3) = (tx(idxt-1+dsearchn(tmp(idxt:end)',.5)) - tx(dsearchn(tmp(1:idxt)',.5)))*1000;

%% inspect the Gaussian (turned off by default)

if showplot
    figure(10001+showplot),clf
    %subplot(211)
    plot(hz,fx,'k', 'LineWidth' , 1.3)
    hold on
    %plot([hz(dsearchn(fx(1:idx)',.5)) hz(idx-1+dsearchn(fx(idx:end)',.5))],[fx(dsearchn(fx(1:idx)',.5)) fx(idx-1+dsearchn(fx(idx:end)',.5))],'k--')
    xline(1,'r--');
    set(gca,'xlim',[max(f-10,0) f+10]);
    legend({'Filter kernel','High-pass cut-off'})
    title(['Empirical fitler: ' num2str(empVals(1)) ', ' num2str(empVals(2)) ' Hz' ])
    xlabel('Frequency (Hz)'), ylabel('Filter gain')
%     
%     subplot(212)
%     tmp1 = real(fftshift(ifft(fx))); tmp1 = tmp1./max(tmp1);
%     tmp2 = abs(hilbert(tmp1));
%     plot(tx,tmp1, tx,tmp2), zoom on
%     xlabel('Time (s)'), ylabel('Amplitude gain')
end

%% done.
