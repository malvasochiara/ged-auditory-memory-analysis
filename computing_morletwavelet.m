function [P_singlesub] = computing_morletwavelet(freqq, delta_f, path, indexes, compp, srate,baseline, norm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calls the function morlet_transform_singlesub and compute the transform
% for each subject.
%
% INPUT:
%      -freqq    = (array) all the frequencies for which you want to compute the amplitude
%      -delta_f  = (array) frequency amplitude corresponding to each frequency
%      -path     = (string) path where the input data are stored
%      -indexeses  = (array) all the indexeses of subjects 
%      -compp    = (array) components you want to consider (time seires you are loading have dim (compp, time, trials) )
%      -srate    = (double) sampling rate of the data
%      -baseline = (array) time indices for normalization
%      -norm     = (double) set to 1 to normalize the power time serie for each subject
% OUTPUT:
%      -P_singlesub = power time serie foer each subject with dimensions (frequency, time, subjects )

ntime = 876; % number of time points

freq_list = dir([path '/freq*']);
P_singlesub = zeros(length(freq_list),ntime,length(indexes));
P_singlesub_temp = zeros(length(freq_list),ntime,length(indexes));
for ff = 1:length(freq_list)
    disp(ff)
    sub_list = dir([path '/' freq_list(ff).name '/SUB*']);
    P_subj = zeros(ntime, length(indexes));
    fc = freqq(ff);
    fwhm = delta_f(ff);
    for sub = 1:length(indexes)
        clear GEDts
        
        load([path '/' freq_list(ff).name '/' sub_list(indexes(sub)).name], 'GEDts')
        P_subj(:,sub) = morlet_transform_singlesub(squeeze(GEDts(compp,:,:)),fc, fwhm,srate,baseline);
        P_singlesub_temp(ff,:,sub) = P_subj(:,sub); 
        
    end
end

if norm == 1
    for sub = 1:length(indexes)
        P_singlesub(:,:,sub) = P_singlesub_temp(:,:,sub)./(max(max(abs(P_singlesub_temp(:,:,sub)))));
    end
else
    P_singlesub(:,:,:) = P_singlesub_temp(:,:,:);
end