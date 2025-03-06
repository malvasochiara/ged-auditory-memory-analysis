function P = morlet_transform_singlesub(ts,fc, fwhm,srate,baseline)

% INPUTS:
%   ts: time series for each trial (time x trial)
%   fc: central frequency of the wavelet
%   fwhm: full width half maximum of the wavelet 
%   srate : sampling rate of the data
%   baseline: array that contains the time indices for normalization
% OUTPUT:
%   P: (ntimeseries) power time serie averaged over trials
%
    P_db = zeros(size(ts));
    for tt = 1:(size(ts,2)) %over trials
        clear filtdat
        % applying a gaussian filter
        [filtdat,~] = filterFGx(ts(:,tt)',srate,fc,fwhm);
        P_raw = abs(hilbert(filtdat)).^2;
        reference = mean(P_raw(baseline));
        P_db(:,tt) = 10*log10(P_raw./reference);
        
    end
    P = mean(P_db(:,:),2);
    

end