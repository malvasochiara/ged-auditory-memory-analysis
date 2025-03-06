%%%%%%%%%
% This code contains various plot and analysis that we tried after GED
%%
addpath('/projects/MINDLAB2023_MEG-AuditMemDement/scripts/chiaramalvaso/GED_TSA2021');
addpath('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021');
%% *** START UP FUNCTIONS.. (LBPD_startup_D) ***
 
%starting up some functions for LBPD toolbox.

pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions
addpath(pathl);
LBPD_startup_D(pathl);
%% VAREXP PLOT FOR EACH FREQUENCY- GED ON SINGLE TRIALS

%%%%%%%%%%%%%% loading explained variance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the evals dividing in young and old, average over participants, plot
% the variance and perform the ANOVA test
load('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/groups_Age.mat');
index_old = cat(2,S.subjs{1,1},S.subjs{1,2});
index_young = S.subjs{1,3};

result_path = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/GED_results';
freq_path_M = [result_path '/Memory'];
freq_path_R = [result_path '/Resting'];

%%%%%%%%%%%%%%%%%% discarding subject with less than 14 correct answers %%%%%%%%%%%%%%%%%%

thresh = 14;
disp('Memory')
[evals_M_O, evals_M_Y, ~,~] = young_or_old_evals(freq_path_M, index_old, index_young,thresh);

disp('Resting')
[evals_R_O, evals_R_Y, index_old, index_young] = young_or_old_evals(freq_path_R, index_old, index_young,thresh);

evals_M_av_O = squeeze(mean(evals_M_O(:,:,:),3));
evals_M_av_Y = squeeze(mean(evals_M_Y(:,:,:),3));

evals_R_av_O = squeeze(mean(evals_R_O(:,:,:),3));
evals_R_av_Y = squeeze(mean(evals_R_Y(:,:,:),3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%% getting the frequency array for the plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_ref = 2.857;  %1000/350 
%%%%%%% 0.5 < f < 2.857  %%%%%%%%
freq_1 = [(1/4)*f_ref, (3/8)*f_ref, (1/2)*f_ref, (3/4)*f_ref, f_ref];

%%%%%% 2.857 < f < 22.856 %%%%%%%% 
freq_2 = zeros(1,14);
for ii = 1:14
    freq_2(1,ii) = ((ii+2)/2)*f_ref;
end

%%%%%% 22.856 < f < 100 %%%%%%%%%%%
freq_3 = zeros(1,9);
for jj = 1:9
    freq_3(1,jj) = (3*jj +8)*f_ref;
end
%%%%%%%% complete array %%%%%%%%%%
central_freq = [freq_1, freq_2, freq_3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% actual plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
compp = 1;
figure;
plot(central_freq, evals_M_av_O(compp,:),'Linewidth', 2, 'Marker', '*', 'DisplayName', 'Listening')
hold on;

plot(central_freq, evals_R_av_O(compp,:),'Linewidth', 2, 'Marker', '*','DisplayName', 'Resting')
hold on;
title('Old participants-cov on concatenated trials')
legend('show')

grid minor
if isempty(thresh)
    saveas(gcf,[result_path '/Varexp_OLD_conctrial.jpg']);
else
    saveas(gcf,[result_path '/Varexp_OLD_conctrial_THRESH_' num2str(thresh) '.jpg']);
end

figure;
plot(central_freq, evals_M_av_Y(compp,:),'Linewidth', 2, 'Marker', '*', 'DisplayName', 'Listening')
hold on;

plot(central_freq, evals_R_av_Y(compp,:),'Linewidth', 2, 'Marker', '*','DisplayName', 'Resting')
hold on;
title('Young participants-cov on concatenated trials')
legend('show')

grid minor

if isempty(thresh)
    saveas(gcf,[result_path '/Varexp_YOUNG_conctrial.jpg']);
else
    saveas(gcf,[result_path '/Varexp_YOUNG_conctrial_THRESH_' num2str(thresh) '.jpg']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ANOVA test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% should we test all frequencies or just the one we are interested in? 
cool_freqs = [2 5 10]; % index of the 10 Hz frequency, change later
for ii = 1:length(cool_freqs)
    clear ff
    ff = cool_freqs(ii);
    numyoung = length(index_young);
    numold = length(index_old);
    % data order : YOUNG RESTING, YOUNG MEMORY, OLD RESTING, OLD MEMORY
    varexp= [squeeze(evals_R_Y(compp,ff,:));squeeze(evals_M_Y(compp,ff,:));squeeze(evals_R_O(compp,ff,:));squeeze(evals_M_O(compp,ff,:))];
    age = [repmat({'Y'},numyoung,1);repmat({'Y'},numyoung,1);repmat({'O'},numold,1);repmat({'O'},numold,1)];
    condition = [repmat({'resting'},numyoung,1);repmat({'memory'},numyoung,1);repmat({'resting'},numold,1);repmat({'memory'},numold,1)];
    dataTable = table(varexp,age,condition);
    [p,tbl,stats] = anovan(dataTable.varexp,{dataTable.age,dataTable.condition}, 'model', 'interaction', 'varname', {'age', 'condition'});
    headers= tbl(1,:);
    headers = strrep(headers, ' ', '_');
    headers = strrep(headers, '.', '');
    headers = strrep(headers, '?', '');
    headers = strrep(headers, '>', '_biggerthan_');
    data = tbl(2:end,:);
    T = cell2table(data, 'VariableNames', headers);
    if isempty(thresh)
        save([result_path '/ANOVA_freq_' num2str(central_freq(ff)) '.mat'], 'T');
    else
        save([result_path '/ANOVA_freq_' num2str(central_freq(ff)) '_THRESH_' num2str(thresh) '.mat'], 'T');
        
    end
end
%% T test on ratio 
% statistic? Freq by freq?
compp = 1;


P = zeros(length(central_freq),1);
T = zeros(length(central_freq),1);
for ff = 1:length(central_freq)
    ratio_Y = squeeze(evals_M_Y(compp, ff,:)./evals_R_Y(compp, ff,:));
    ratio_O = squeeze(evals_M_O(compp, ff,:)./evals_R_O(compp, ff,:));
    [~,p,~,stats] = ttest(squeeze(evals_M_O(compp,ff,:)),squeeze(evals_R_O(compp,ff,:)));
    P(ff,1) = p;
    T(ff,1) = stats.tstat;
    clear ratio_Y
    clear ratio_O
    clear p
    clear stats
end

[fdr_siglevel, ~,~] = fdr(P);
ii = 1; % counting significant frequencies
for ff = 1:length(central_freq)
    if P(ff) < 0.05 % no correction
        disp(central_freq(ff))
        sig_freq_index(ii) = ff;
        ii = ii+1;
    end
end


%% VAREXP PLOT FOR EACH FREQUENCY- GED ON CONCATENATED TRIALS

%%%%%%%%%%%%%% loading explained variance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the evals dividing in young and old, average over participants, plot
% the variance and perform the ANOVA test
load('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/groups_Age.mat');
index_old = cat(2,S.subjs{1,1},S.subjs{1,2});
index_young = S.subjs{1,3};

result_path = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/GED_results';
freq_path_M = [result_path '/Memory'];
freq_path_R = [result_path '/Resting'];
%%%%%%%%%%%%%%%%%% discarding subject with less than 14 correct answers %%%%%%%%%%%%%%%%%%
thresh = 14;

disp('Memory')
[evals_M_O, evals_M_Y, ~, ~] = young_or_old_evals(freq_path_M, index_old, index_young,thresh);

disp('Resting')
[evals_R_O, evals_R_Y, index_old, index_young] = young_or_old_evals(freq_path_R, index_old, index_young,thresh);

evals_M_av_O = squeeze(mean(evals_M_O(:,:,:),3));
evals_M_av_Y = squeeze(mean(evals_M_Y(:,:,:),3));

evals_R_av_O = squeeze(mean(evals_R_O(:,:,:),3));
evals_R_av_Y = squeeze(mean(evals_R_Y(:,:,:),3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%% getting the frequency array for the plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_ref = 2.857;  %1000/350 
%%%%%%% 0.5 < f < 2.857  %%%%%%%%
freq_1 = [(1/4)*f_ref, (3/8)*f_ref, (1/2)*f_ref, (3/4)*f_ref, f_ref];

%%%%%% 2.857 < f < 22.856 %%%%%%%% 
freq_2 = zeros(1,14);
for ii = 1:14
    freq_2(1,ii) = ((ii+2)/2)*f_ref;
end

%%%%%% 22.856 < f < 100 %%%%%%%%%%%
freq_3 = zeros(1,9);
for jj = 1:9
    freq_3(1,jj) = (3*jj +8)*f_ref;
end
%%%%%%%% complete array %%%%%%%%%%
central_freq = [freq_1, freq_2, freq_3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% actual plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
compp = 1;
figure;
plot(central_freq, evals_M_av_O(compp,:),'Linewidth', 2, 'Marker', '*', 'DisplayName', 'Listening')
hold on;

plot(central_freq, evals_R_av_O(compp,:),'Linewidth', 2, 'Marker', '*','DisplayName', 'Resting')
hold on;
title('Old participants-cov on concatenated trials trials')
legend('show')

grid minor
% if isempty(thresh)
%     saveas(gcf,'/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/GED_results/varexp_old.fig');
% else
%     saveas(gcf,['/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/GED_results/varexp_old_THRESH_' num2str(thresh) '.fig']);
% end

figure;
plot(central_freq, evals_M_av_Y(compp,:),'Linewidth', 2, 'Marker', '*', 'DisplayName', 'Listening')
hold on;

plot(central_freq, evals_R_av_Y(compp,:),'Linewidth', 2, 'Marker', '*','DisplayName', 'Resting')
hold on;
title('Young participants-cov on concatenated trials')
legend('show')

grid minor

% if isempty(thresh)
%     saveas(gcf,'/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/GED_results/varexp_young.fig');
% else
%     saveas(gcf,['/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/GED_results/varexp_young_THRESH_' num2str(thresh) '.fig']);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ANOVA test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% should we test all frequencies or just the one we are interested in? 
cool_freqs = [2 5 10]; % index of the 10 Hz frequency, change later
for ii = 1:length(cool_freqs)
    clear ff
    ff = cool_freqs(ii);
    numyoung = length(index_young);
    numold = length(index_old);
    % data order : YOUNG RESTING, YOUNG MEMORY, OLD RESTING, OLD MEMORY
    varexp= [squeeze(evals_R_Y(compp,ff,:));squeeze(evals_M_Y(compp,ff,:));squeeze(evals_R_O(compp,ff,:));squeeze(evals_M_O(compp,ff,:))];
    age = [repmat({'Y'},numyoung,1);repmat({'Y'},numyoung,1);repmat({'O'},numold,1);repmat({'O'},numold,1)];
    condition = [repmat({'resting'},numyoung,1);repmat({'memory'},numyoung,1);repmat({'resting'},numold,1);repmat({'memory'},numold,1)];
    dataTable = table(varexp,age,condition);
    [p,tbl,stats] = anovan(dataTable.varexp,{dataTable.age,dataTable.condition}, 'model', 'interaction', 'varname', {'age', 'condition'});
    headers= tbl(1,:);
    headers = strrep(headers, ' ', '_');
    headers = strrep(headers, '.', '');
    headers = strrep(headers, '?', '');
    headers = strrep(headers, '>', '_biggerthan_');
    data = tbl(2:end,:);
    T = cell2table(data, 'VariableNames', headers);
%     if isempty(thresh)
%         save([result_path '/ANOVA_freq_' num2str(central_freq(ff)) '.mat'], 'T');
%     else
%         save([result_path '/ANOVA_freq_' num2str(central_freq(ff)) '_THRESH_' num2str(thresh) '.mat'], 'T');
%     end
end
%% OBTAINING THE INDEXES OF THOSE WHO GAVE MORE THAN 50 % CORRECT ANSWERS

addpath('/projects/MINDLAB2023_MEG-AuditMemDement/scripts/chiaramalvaso/GED_TSA2021/Final_scripts')
thresh = 14;
[index_old, index_young] = young_or_old(thresh);


%%%%%%%%%%%%%%%%%%%%%%%%%%% getting the frequency array %%%%%%%%%%%%%%%%%%%%%%%%%%%
f_ref = 2.857;  %1000/350 
%%%%%%% 0.5 < f < 2.857  %%%%%%%%
freq_1 = [(1/4)*f_ref, (3/8)*f_ref, (1/2)*f_ref, (3/4)*f_ref, f_ref];

%%%%%% 2.857 < f < 22.856 %%%%%%%% 
freq_2 = zeros(1,14);
for ii = 1:14
    freq_2(1,ii) = ((ii+2)/2)*f_ref;
end

%%%%%% 22.856 < f < 100 %%%%%%%%%%%
freq_3 = zeros(1,9);
for jj = 1:9
    freq_3(1,jj) = (3*jj +8)*f_ref;
end
%%%%%%%% complete array %%%%%%%%%%
central_freq = [freq_1, freq_2, freq_3];
clear freq_1
clear freq_2
clear freq_3
clear f_ref
% define the coordinates of the reference points computed by checking the filter with filter_properties.m
x_r = 2.857; %1000/350 = target frequency
y_r = 0.6;
x_low = 1.429;
y_low = 0.6;
x_high = 32;
y_high = 1.429; %half of the target frequency
m_high = (y_high - y_r)/ (x_high-x_r);
q_high = y_high - m_high*x_high;
m_low = (y_r-y_low)/(x_r-x_low);
q_low = y_low - m_low*x_low;

ampl = @(freq) (freq <= 2.875).* (m_low.*freq + q_low) + (freq > 2.875 ).*(m_high.*freq + q_high);
delta_freq = ampl(central_freq); % array containing the interval amplitude corresponding to each central frequency
clear x_r
clear y_r
clear x_low
clear y_low
clear x_high
clear y_high
clear m_high
clear m_low
clear q_high
clear q_low
clear ampl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% actual wavelet computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% setting up variablese that must be the same in all the computations %%%%%%%%%%%%%%% 
srate = 250; % sampling rate in Hz
baseline = 1:25; % index for baseline correction 
compp = 1; % GED component you want to use for the wavelet transform
norm = []; % normalizing independently for each subject
memory_path = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/GED_covsingletrial_results/Memory';

disp('Computing memory young')
P_M_Y_singlesub_norm = computing_morletwavelet(central_freq, delta_freq, memory_path, index_young, compp, srate,baseline, norm);
disp('Computing memory old')
P_M_O_singlesub_norm = computing_morletwavelet(central_freq, delta_freq, memory_path, index_old, compp, srate,baseline, norm);

rest_path = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/GED_covsingletrial_results/Resting';
disp('Computing resting young')
P_R_Y_singlesub_norm = computing_morletwavelet(central_freq, delta_freq, rest_path, index_young, compp, srate,baseline, norm);
disp('Computing resting old')
P_R_O_singlesub_norm = computing_morletwavelet(central_freq, delta_freq, rest_path, index_old, compp, srate,baseline, norm);

%%%%%%%%%%%%%%%%%%%%%%%%%%% ANOVA %%%%%%%%%%%%%%%%%%%%%%%%%%%
% should we test all frequencies or just the one we are interested in? 
numyoung = length(index_young);
numold = length(index_old);
time = linspace(0,3.5,876);
anova_results = cell(length(central_freq),1);
for ff = 1:length(central_freq)
    temp_results = {};
    count = 0;
    disp(ff)
    for tt = 1:length(time)
        
        % data order : YOUNG RESTING, YOUNG MEMORY, OLD RESTING, OLD MEMORY
        power_ts= [squeeze(P_R_Y_singlesub_norm(ff,tt,:));squeeze(P_M_Y_singlesub_norm(ff,tt,:));squeeze(P_R_O_singlesub_norm(ff,tt,:));squeeze(P_M_O_singlesub_norm(ff,tt,:))];
        age = [repmat({'Y'},numyoung,1);repmat({'Y'},numyoung,1);repmat({'O'},numold,1);repmat({'O'},numold,1)];
        condition = [repmat({'resting'},numyoung,1);repmat({'memory'},numyoung,1);repmat({'resting'},numold,1);repmat({'memory'},numold,1)];
        dataTable = table(power_ts,age,condition);
        clear p
        [p,tbl,stats] = anovan(dataTable.power_ts,{dataTable.age,dataTable.condition}, 'model', 'interaction', 'varname', {'age', 'condition'}, 'display', 'off');
        if any(p <= 0.005)
            count = count +1;
            temp_results{count, 1} = time(tt);
            headers= tbl(1,:);
            headers = strrep(headers, ' ', '_');
            headers = strrep(headers, '.', '');
            headers = strrep(headers, '?', '');
            headers = strrep(headers, '>', '_biggerthan_');
            data = tbl(2:end,:);
            T = cell2table(data, 'VariableNames', headers);
            temp_results{count, 2} = T;
        end
        clear headers
        clear data
        clear T
    end
    anova_results{ff,1} = temp_results;
end

result_path = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021';
if isempty(thresh)
    if norm == 1
        save([result_path '/ANOVA_normalized_power_ts.mat'], 'anova_results');
    else
        save([result_path '/ANOVA_power_ts.mat'], 'anova_results');
    end
else
    if norm == 1
        save([result_path '/ANOVA_normalized_power_ts_THRESH_' num2str(thresh) '.mat'], 'anova_results');
    else
        save([result_path '/ANOVA_power_ts_THRESH_' num2str(thresh) '.mat'], 'anova_results');
    end
end
%% PLOTTING THE POWER TIME SERIES
dest_path = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/wavelet_results';
%%%%% MEMORY YOUNG %%%%%%%
figure;
imagesc(time, 1:length(central_freq), squeeze(mean(P_M_Y_singlesub_norm(:,:,:),3))./max(max(abs(mean(P_M_Y_singlesub_norm(:,:,:),3)))));
yticks(1:length(central_freq));
yticklabels(arrayfun(@num2str, round(central_freq,1), 'UniformOutput', false))
hold on;
title('Memory young - normalized')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
colorbar;
caxis([-1 0.2])
set(gcf, 'color', 'w')
saveas(gcf,[dest_path '/wavelet_MY.jpg']);


%%%%% MEMORY OLD %%%%%%%
figure;
imagesc(time,1:length(central_freq), squeeze(mean(P_M_O_singlesub_norm(:,:,:),3))./max(max(abs(mean(P_M_O_singlesub_norm(:,:,:),3)))));
yticks(1:length(central_freq));
yticklabels(arrayfun(@num2str, round(central_freq,1), 'UniformOutput', false))
hold on;
title('Memory old - normalized')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
colorbar;
caxis([-1 0.2])
set(gcf, 'color', 'w')
saveas(gcf,[dest_path '/wavelet_MO.jpg']);

%%%%% RESTING YOUNG %%%%%%%
figure;
imagesc(time,1:length(central_freq), squeeze(mean(P_R_Y_singlesub_norm(:,:,:),3))./max(max(abs(mean(P_R_Y_singlesub_norm(:,:,:),3)))));
yticks(1:length(central_freq));
yticklabels(arrayfun(@num2str, round(central_freq,1), 'UniformOutput', false))
hold on;
title('Resting young - normalized')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
colorbar;
caxis([-1 0.2])
set(gcf, 'color', 'w')
saveas(gcf,[dest_path '/wavelet_RY.jpg']);

%%%%% RESTING OLD %%%%%%%
figure;
imagesc(time,1:length(central_freq), squeeze(mean(P_R_O_singlesub_norm(:,:,:),3))./max(max(abs(mean(P_R_O_singlesub_norm(:,:,:),3)))));
yticks(1:length(central_freq));
yticklabels(arrayfun(@num2str, round(central_freq,1), 'UniformOutput', false))
hold on;
title('Resting old - normalized')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
colorbar;
caxis([-1 0.2])
set(gcf, 'color', 'w')
saveas(gcf,[dest_path '/wavelet_RO.jpg']);

%%%%%%%%%%%%% trying to plot the difference of the mean to see which group has the highest values %%%%%%%%%%%%%
% MEMORY - YOUNG VS OLD
diff_M = squeeze(mean(P_M_Y_singlesub_norm(:,:,:),3))- squeeze(mean(P_M_O_singlesub_norm(:,:,:),3));
figure;
imagesc(time,1:length(central_freq), diff_M./(max(max(abs(diff_M)))) );
yticks(1:length(central_freq));
yticklabels(arrayfun(@num2str, round(central_freq,1), 'UniformOutput', false))
hold on;
title('MEMORY: young - old')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
colorbar;
caxis([-1 1])
set(gcf, 'color', 'w')
saveas(gcf,[dest_path '/wavelet_DIFF_memory.jpg']);

% RESTING - YOUNG VS OLD
diff_R = squeeze(mean(P_R_Y_singlesub_norm(:,:,:),3))- squeeze(mean(P_R_O_singlesub_norm(:,:,:),3));
figure;
imagesc(time,1:length(central_freq), diff_R./(max(max(abs(diff_R)))) );
yticks(1:length(central_freq));
yticklabels(arrayfun(@num2str, round(central_freq,1), 'UniformOutput', false))
hold on;
title('RESTING: young - old')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
colorbar;
caxis([-1 1])
set(gcf, 'color', 'w')
saveas(gcf,[dest_path '/wavelet_DIFF_resting.jpg']);

% YOUNG - Memory - Resting
diff_Y = squeeze(mean(P_M_Y_singlesub_norm(:,:,:),3))- squeeze(mean(P_R_Y_singlesub_norm(:,:,:),3));
figure;
imagesc(time,1:length(central_freq), diff_Y./(max(max(abs(diff_Y)))) );
yticks(1:length(central_freq));
yticklabels(arrayfun(@num2str, round(central_freq,1), 'UniformOutput', false))
hold on;
title('YOUNG: memory - resting')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
colorbar;
caxis([-1 1])
set(gcf, 'color', 'w')
saveas(gcf,[dest_path '/wavelet_DIFF_young.jpg']);

% OLD - Memory - Resting
diff_O = squeeze(mean(P_M_O_singlesub_norm(:,:,:),3))- squeeze(mean(P_R_O_singlesub_norm(:,:,:),3));
figure;
imagesc(time,1:length(central_freq), diff_O./(max(max(abs(diff_O)))) );
yticks(1:length(central_freq));
yticklabels(arrayfun(@num2str, round(central_freq,1), 'UniformOutput', false))
hold on;
title('OLD: memory - resting')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
colorbar;
caxis([-1 1])
set(gcf, 'color', 'w')
saveas(gcf,[dest_path '/wavelet_DIFF_old.jpg']);

%% Plotting F values on signifcance
load('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/ANOVA_power_ts_THRESH_14.mat')

F_age = zeros(length(central_freq));
F_cond = zeros(length(central_freq));
F_agecond = zeros(length(central_freq));

for ff = 1:length(central_freq)
    disp(ff)
    clear temp
    if ~isempty(anova_results{ff,1})
        temp = anova_results{ff,1};
        for jj = 1:size(temp,1)
            clear tbl
            clear index
            index = find(time == temp{jj,1}); % finding the time index
            tbl = temp{jj,2};
            if cell2mat(tbl{1,7}) <= 0.05 % if the age is significant
                F_age(ff,index) = cell2mat(tbl{1,6});
            end
            
            if cell2mat(tbl{2,7}) <= 0.05 % if the age is significant
                F_cond(ff,index) = cell2mat(tbl{2,6});
            end
            
            if cell2mat(tbl{3,7}) <= 0.05 % if the age is significant
                F_agecond(ff,index) = cell2mat(tbl{3,6});
            end
        
        end
    end
end

freq_lim = 28;

F_age_norm = F_age;%./(max(max(abs(F_age))));
h =figure;
imagesc(time, 1:length(central_freq), F_age_norm(1:freq_lim,:));
hold on;
yticks(1:length(central_freq));
yticklabels(arrayfun(@num2str, round(central_freq,1), 'UniformOutput', false))
title('F values - age significance')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
colorbar;
caxis([0 25])
set(gcf, 'color', 'w')
saveas(gcf,[dest_path '/wavelet_ANOVA_age.jpg']);

F_cond_norm = F_cond;%./(max(max(abs(F_cond))));
h =figure;
imagesc(time, 1:length(central_freq), F_cond_norm(1:freq_lim,:));
hold on;
yticks(1:length(central_freq));
yticklabels(arrayfun(@num2str, round(central_freq,1), 'UniformOutput', false))
title('F values - condition significance')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
colorbar;
caxis([0 25])
set(gcf, 'color', 'w')
saveas(gcf,[dest_path '/wavelet_ANOVA_condition.jpg']);


F_agecond_norm = F_agecond;%./(max(max(abs(F_agecond))));
h =figure;
imagesc(time, 1:length(central_freq), F_agecond_norm(1:freq_lim,:));
hold on;
yticks(1:length(central_freq));
yticklabels(arrayfun(@num2str, round(central_freq,1), 'UniformOutput', false))
title('F values - age*condition significance')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
colorbar;
caxis([0 25])
set(gcf, 'color', 'w')
saveas(gcf,[dest_path '/wavelet_ANOVA_agecondition.jpg']);

%%