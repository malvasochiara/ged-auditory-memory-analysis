%% 
%% NOTES FOR FUTURE CHIARA
% this code include all the steps, from preprocessing to GED. You didn't
% actually perform all these steps, you did everything starting from
% epoching. Most of these steps need to be repeated for memory and
% resting, some of the are different among the two conditions (for example epoching is different since you don't have any trigger in resting) 
%% TEMP SEQ AGES 2021 - PREPROCESSING

%% Maxfilter


%OBS! before running maxfilter you need to close matlab, open the terminal and write: 'use anaconda', then open matlab and run maxfilter script

maxfilter_path = '/neuro/bin/util/maxfilter';
project = 'MINDLAB2021_MEG-TempSeqAges';
maxDir = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter'; %output path
movement_comp = 1; %1 = yes; 0 = no


path = '/raw/sorted/MINDLAB2021_MEG-TempSeqAges'; %path with all the subjects folders
jj = dir([path '/0*']); %list all the folders starting with '0' in order to avoid hidden files
for ii = 75%:length(jj) %over subjects
    if ii ~= 49
        cart = [jj(ii).folder '/' jj(ii).name]; %create a path that combines the folder name and the file name
        pnana = dir([cart '/2*']); %search for folders starting with '2'
        for pp = 1:length(pnana) %loop to explore ad analyze all the folders inside the path above
            cart2 = [pnana(1).folder '/' pnana(pp).name];
            pr = dir([cart2 '/ME*']); %looks for meg folder
            if ~isempty(pr) %if pr is not empty, proceed with subfolders inside the meg path
                if ii ~= 48 %fixing issue of SUBJ0049 archived as additional 5 blocks of SBUJ0048
                    pnunu = dir([pr(1).folder '/' pr(1).name '/00*']);
                else
                    pnunu = dir([pr(1).folder '/' pr(1).name '/0*']); %only one '0' here since the folder has 10 files and so '009-010'
                    movement_comp = 0; %we cannot do movement compensation for SUBJ0049 since we have movement data of SUBJ0048
                end
                for dd = 1:length(pnunu)
                    fpath = dir([pnunu(1).folder '/' pnunu(dd).name '/files/*.fif']); % looks for .fif file
                    rawName = ([fpath.folder '/' fpath.name]); %assigns the final path of the .fif file to the rawName path used in the maxfilter command
                    if ii ~= 48 %SUBJ0049 was not archived properly (not loaded its preparation..), so we needed a manual operation to deal with it, as reported in the next section
                        maxfName = ['SUBJ' jj(ii).name '_' fpath.name(1:end-4)]; %define the output name of the maxfilter processing
                    else
                        if dd < 6 %first 5 block are of SUBJ0048
                            maxfName = ['SUBJ' jj(ii).name '_' fpath.name(1:end-4)]; %define the output name of the maxfilter processing
                        else %blocks from 6 to 10 are of SUBJ0049
                            maxfName = ['SUBJ0049_' fpath.name(1:end-4)]; %define the output name of the maxfilter processing
                        end
                    end
                    if movement_comp == 1
                        %movement compensation
                        cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ' ,project, ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsssdsm.fif'] ' -st 4 -corr 0.98 -movecomp -ds 4 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsssdsm.log"']];
                    else %no movement compensation (to be used if HPI coils did not work properly)
                        cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ' ,project, ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsssdsm.fif'] ' -st 4 -corr 0.98 -ds 4 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsssdsm.log"']];
                    end
                    system(cmd);
                end
            end
        end
    end
end


%% Starting up OSL

%OBS! run this before converting the .fif files into SPM objects

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/osl-core'); %adds the path to OSL functions
osl_startup %starts the osl package


%% Converting the .fif files into SPM objects

%OBS! remember to run 'starting up OSL' first

%setting up the cluster
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path to the function that submits the jobs to the cluster
clusterconfig('scheduler', 'cluster');
clusterconfig('long_running', 1); %there are different queues for the cluster depending on the number and length of the jobs you want to submit 
clusterconfig('slot', 1); %slot in the queu

%% conversion to SPM objects

fif_list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/*.fif'); %path to SPM objects

for ii = 372:376%length(fif_list) %over the .fif files
    S = []; %structure 'S'                   
    S.dataset = [fif_list(ii).folder '/' fif_list(ii).name];
    D = spm_eeg_convert(S);
%     D = job2cluster(@cluster_spmobject, S); %actual function for conversion
end

%% Removing bad segments using OSLVIEW

%checks data for potential bad segments (periods)
%marking is done by right-clicking in the proximity of the event and click on 'mark event'
%a first click (green dashed label) marks the beginning of a bad period
%a second click indicates the end of a bad period (red)
%this will mean that we are not using about half of the data, but with such bad artefacts this is the best we can do
%we can still obtain good results with what remains
%NB: Push the disk button to save to disk (no prefix will be added, same name is kept)

%OBS! remember to check for bad segments of the signal both at 'megplanar' and 'megmag' channels (you can change the channels in the OSLVIEW interface)

%OBS! remember to mark the trial within the bad segments as 'badtrials' and use the label for removing them from the Averaging (after Epoching) 

spm_list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/spmeeg*.mat'); %path to SPM objects

for ii = 31%372:length(spm_list) %over experimental blocks %OBS!
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]);
    D = oslview(D);
    D.save(); %save the selected bad segments and/or channels in OSLVIEW
    disp(ii)
end

%% UPDATE FIDUCIALS FOR SUBJECT0049

%loading empty room with SUBJ0049 fiducials
Dfid = spm_eeg_load('/scratch7/MINDLAB2021_MEG-TempSeqAges/spmeeg_emptyroom.mat');
fid = Dfid.fiducials; %extracting fiducials
spm_list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/spmeeg_SUBJ0049*.mat'); %path to SPM objects of SUBJ0049
%loading files preprocessed for SUBJ0049

for ii = 1:length(spm_list)
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]);
    %pasting fiducials (here you need also to adapt the function "subsasgn.m" (you find comments there to help you..))
    D.fiducials = fid;
    D.save();
    disp(ii)
end


%% AFRICA denoising (part I)

%setting up the cluster
clusterconfig('scheduler', 'cluster');
clusterconfig('long_running', 1); %there are different queues for the cluster depending on the number and length of the jobs you want to submit 
clusterconfig('slot', 1); %slot in the queu
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path to the function that submits the jobs to the cluster


%%

%ICA calculation
spm_list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/spmeeg*.mat');

for ii = 372:length(spm_list) %OBS!
    S = [];
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]);
    S.D = D;
    
    jobid = job2cluster(@cluster_africa,S);
%   D = osl_africa(D,'do_ica',true,'do_ident',false,'do_remove',false,'used_maxfilter',true); 
%   D.save();
end

%% AFRICA denoising (part II)

% v = [11 12 19 32];
%visual inspection and removal of artifacted components
%look for EOG and ECG channels (usually the most correlated ones, but check a few more just in case)
spm_list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/spmeeg*.mat');

for ii = 372:length(spm_list) %OBS!%38:41
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]);
    D = osl_africa(D,'do_ident','manual','do_remove',false,'artefact_channels',{'EOG','ECG'});
    %hacking the function to manage to get around the OUT OF MEMORY problem..
    S = [];
    S.D = D;
    jobid = job2cluster(@cluster_rembadcomp,S);
%   D.save();
    disp(ii)
end

%%

%% NOTHING UNITL HERE

%%

%%

%% CHIARA, FROM HERE..

%% COPYING MEMORY TASK DATA FILES
% never run this again, it's just copying data
spm_list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/spmeeg*minor*.mat');
spm_list2 = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/spmeeg*minor*.dat');
outdir = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/Continuous/Memory';
mkdir(outdir);

for ii = 1:length(spm_list)
   copyfile([spm_list(ii).folder '/' spm_list(ii).name],[outdir '/' spm_list(ii).name]) 
   copyfile([spm_list2(ii).folder '/' spm_list2(ii).name],[outdir '/' spm_list2(ii).name])     
end


%% COPYING RESTING STATE DATA FILES
% never run this again, it's just copying data
spm_list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/spmeeg*res*.mat');
spm_list2 = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/spmeeg*res*.dat');
outdir = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/Continuous/Resting';
mkdir(outdir);

for ii = 1:length(spm_list)
   copyfile([spm_list(ii).folder '/' spm_list(ii).name],[outdir '/' spm_list(ii).name]) 
   copyfile([spm_list2(ii).folder '/' spm_list2(ii).name],[outdir '/' spm_list2(ii).name])     
end


%%

%% Epoching: one epoch per old/new excerpt (baseline = (-)100ms)
 %%%%% RUN SECTION Starting up OSL BEFORE RUNNING THIS ONE TO SET UP FUNCTIONS %%%%%
% Epoching the memory data. This is automatic since you have the
% information of the triggers
prefix_tobeadded = 'e'; %adds this prefix to epoched files
epoch_end = 3.4; %time in seconds of epoch length !The actual length of the epoch is 3.5  s bc you use the 0.1 s BEFORE the trigger for beseline correction

spm_list = dir('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/Continuous/Memory/spmeeg*.mat');

for ii = 1:length(spm_list) %over .mat files
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]); %load spm_list .mat files
    D = D.montage('switch',0);
    dummy = D.fname; %OBS! D.fname does not work, so we need to use a 'dummy' variable instead
    %if strcmp(dummy(22:26), 'speed') %checks whether characters 22 to 26 are equal to 'speed'; the loop continues if this is true (1) and it stops if this is false (0)
    events = D.events; %look for triggers
    
    %takes the correct triggers sent during the recording
    clear trigcor
    count_evval = 0; %???
    for ieve = 1:length(events) %over triggers
        if strcmp(events(ieve).type,'STI101_up') %only triggers at the beginning of each stimuli
            if strcmp('rec',spm_list(ii).name(17:19)) || strcmp('rac',spm_list(ii).name(17:19))
                if events(ieve).value == 10 || events(ieve).value == 50 %10 and 50 are old and new in recogminor (block 3), while 11 and 21 are old and new in blocks 4 and 5 (aud and vis)
                    count_evval = count_evval + 1;
                    trigcor(count_evval,1) = events(ieve).time; %+ 0.010; %this takes the correct triggers and add 10ms of delay of the sound travelling into the tubes
                    %variable with all the triggers we need
                end
            end
        end
    end
    % matrix with the information about the beginning of the epoch in seconds
    trl_sam = zeros(length(trigcor),3); %prepare the samples matrix with 0's in all its cells
    % matrix with the information about the beginning of the epoch in points (use sample frequency for the conversion
    trl_sec = zeros(length(trigcor),3); %prepare the seconds matrix with 0's in all its cells
    %deftrig = zeros(length(trigcor),1); %this is not useful
    for k = 1:length(trigcor) %over selected triggers
        %deftrig(k,1) = 0.012 + trigcor(k,1); %adding a 0.012 seconds delay to the triggers sent during the experiment (this delay was due to technical reasons related to the stimuli)
        trl_sec(k,1) = trigcor(k,1) - 0.1000; %beginning time-window epoch in s (please note that we computed this operation two times, obtaining two slightly different pre-stimulus times.
        %this was done because for some computations was convenient to have a slightly longer pre-stimulus time
        %remove 1000ms of baseline
        trl_sec(k,2) = trigcor(k,1) + epoch_end; %end time-window epoch in seconds
        trl_sec(k,3) = trl_sec(k,2) - trl_sec(k,1); %range time-windows in seconds
        trl_sam(k,1) = round(trl_sec(k,1) * 250) + 1; %beginning time-window epoch in samples %250Hz per second
        trl_sam(k,2) = round(trl_sec(k,2) * 250) + 1; %end time-window epoch in samples
        trl_sam(k,3) = -25; %sample before the onset of the stimulus (corresponds to 0.100ms)
    end
    dif = trl_sam(:,2) - trl_sam(:, 1); %difference between the end and the beginning of each sample (just to make sure that everything is fine)
    if ~all(dif == dif(1)) %checking if every element of the vector are the same (i.e. the length of the trials is the same; we may have 1 sample of difference sometimes because of different rounding operations..)
        trl_sam(:,2) = trl_sam(:,1) + dif(1);
    end
    %creates the epochinfo structure that is required for the source reconstruction later
    epochinfo.trl = trl_sam;
    epochinfo.time_continuous = D.time;
    %switch the montage to 0 because for some reason OSL people prefer to do the epoching with the not denoised data
    D = D.montage('switch',0);
    %build structure for spm_eeg_epochs
    S = [];
    S.D = D;
    S.trl = trl_sam;
    S.prefix = prefix_tobeadded;
    % function that actually perform the epoching
    D = spm_eeg_epochs(S);
    %store the epochinfo structure inside the D object
    D.epochinfo = epochinfo;
    D.save();
    
    %take bad segments registered in OSLVIEW and check if they overlap with the trials. if so, it gives the number of overlapped trials that will be removed later
    count = 0;
    Bad_trials = zeros(length(trl_sec),1);
    for kkk = 1:length(events) %over events
        if strcmp(events(kkk).type,'artefact_OSL')
            for k = 1:length(trl_sec) %over trials
                if events(kkk).time - trl_sec(k,2) < 0 %if end of trial is > than beginning of artifact
                    if trl_sec(k,1) < (events(kkk).time + epoch_end) %if beginning of trial is < than end of artifact
                        Bad_trials(k,1) = 1; %it is a bad trial (stored here)
                        count = count + 1;
                    end
                end
            end
        end
    end
    %if bad trials were detected, their indices are stored within D.badtrials field
    disp(spm_list(ii).name);
    if count == 0
        disp('there are no bad trials marked in oslview');
    else
        D = badtrials(D,find(Bad_trials),1); %get the indices of the badtrials marked as '1' (that means bad)
        %         D = conditions(D,find(Bad_trials),1); %get the indices of the badtrials marked as '1' (that means bad)
        epochinfo = D.epochinfo;
        xcv = find(Bad_trials == 1);
        %this should be done only later.. in any case.. not a problem..
        for jhk = 1:length(xcv)
            D = D.conditions(xcv(jhk),'Bad');
            epochinfo.conditionlabels(xcv(jhk)) = {'Bad'};
            disp([num2str(ii) ' - ' num2str(jhk) ' / ' num2str(length(xcv))])
        end
        D.epochinfo = epochinfo;
        D.save(); %saving on disk
        disp('bad trials are ')
        length(D.badtrials)
    end
    D.save();
    disp(ii)
end


%% 

%%

%% Defining the conditions - All blocks
% This section gives the correct label to the trials (what we used for now
% is only Old Correct)
%%% CHANGE PATHS, ETC. %%%
%is this the right path? are they going to be saved here?
epoch_list = dir('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/Continuous/Memory/espmeeg*.mat'); %dir to epoched files


xlsx_dir_behav = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/MEG_behavioural_TSA2021'; %dir to MEG behavioral results (.xlsx files)
for ii = 1:length(epoch_list) %over epoched data
       %WHEN YOU HAVE THE DATA CHECK IF 18:20 IS CORRECT...
    if ~strcmp('mmn',epoch_list(ii).name(18:20)) %if the block is not the MMN block..
        %loading SPM object (epoched data)
        D = spm_eeg_load([epoch_list(ii).folder '/' epoch_list(ii).name]);
        dummy = D.fname; %%% OBS!! (SUBJ0041 IS NOT AUD1 AS IN THE NAME BUT IT IS AUD2!!) HERE YOU DID MANUAL ADJUSTMENT OF VARIABLE dummy TO FIX IT..
        %barbaric solution.. to build the name to be read for the excel files with the MEG behavioral tasks performance
        
%         MAYBE YOU CAN OPTIMIZE THIS PART IF YOU ARE SURE THAT YOU CARE ONLY ABOUT BLOCK 3

        if strcmp(dummy(18:20),'rec') || strcmp(dummy(18:20),'rac')
            dumbloc = 'Block_3.xlsx';
            bl = 3;
        elseif strcmp(dummy(18:20),'aud')
            dumbloc = 'Block_4_Auditory.xlsx';
            bl = 4;
        elseif strcmp(dummy(18:20),'vis')
            dumbloc = 'Block_5_Visual.xlsx';
            bl = 5;
        end
        dumls = ['Subj_' dummy(13:16) '_' dumbloc];
        [~,~,raw_recog] = xlsread([xlsx_dir_behav '/' dumls]); %excel files
        %picking the current block
        if bl == 3 %block 3
            for k = 1:length(D.trialonset)
                if raw_recog{(k + 1),3} == 0 %if there was no response
                    D = D.conditions(k,'No_response');
                elseif strcmp(raw_recog{(k + 1),2}(8:9),'ol') && raw_recog{(k + 1),3} == 1 %old correct
                    D = D.conditions(k,'Old_Correct'); %assign old correct
                elseif strcmp(raw_recog{(k + 1),2}(8:9),'ol') && raw_recog{(k + 1),3} == 2 %old incorrect
                    D = D.conditions(k,'Old_Incorrect'); %otherwise assign new correct
                elseif strcmp(raw_recog{(k + 1),2}(14:15),'t1') && raw_recog{(k + 1),3} == 2 %new t1 correct
                    D = D.conditions(k,'New_T1_Correct');
                elseif strcmp(raw_recog{(k + 1),2}(14:15),'t1') && raw_recog{(k + 1),3} == 1 %new t1 incorrect
                    D = D.conditions(k,'New_T1_Incorrect');
                elseif strcmp(raw_recog{(k + 1),2}(14:15),'t3') && raw_recog{(k + 1),3} == 2 %new t3 correct
                    D = D.conditions(k,'New_T3_Correct');
                elseif strcmp(raw_recog{(k + 1),2}(14:15),'t3') && raw_recog{(k + 1),3} == 1 %new t3 incorrect
                    D = D.conditions(k,'New_T3_Incorrect');
                end
            end
        end
        %this is for every block
        if ~isempty(D.badtrials) %overwriting badtrials (if any) on condition labels
            BadTrials = D.badtrials;
            for badcount = 1:length(BadTrials) %over bad trials
                D = D.conditions(BadTrials(badcount),'Bad_trial');
            end
        end
        D = D.montage('switch',1);
        D.epochinfo.conditionlabels = D.conditions; %to add for later use in the source reconstruction
        D.save(); %saving data on disk
    end
    disp(num2str(ii))
end

%%


%% RESTING STATE - EPOCHING - DONE 
% This section makes tyhe epoch in resting state data. Since there is no
% external trigger, startin time for the epochs are randomly chose
% following these criteria: - for each subjects, the number of trials is the same as in memory
%                           - avoid the first 60 s of recordings because the initial data is 0
%                           - trials should not overlap
%                           - epoch length is the same as in memory (3.5 s considering baseline correction)
%                           - trials must be ordered (trial 1 happens before trial 2 and so on)

% %% Epoching: one epoch per old/new excerpt (baseline = (-)100ms)
 %%%%% RUN SECTION Starting up OSL BEFORE RUNNING THIS ONE TO SET UP FUNCTIONS %%%%%

prefix_tobeadded = 'e'; %adds this prefix to epoched files

epoch_end = 3.5; %time in seconds of epoch length
epoch_end_sam = round(epoch_end*250); %length of the epoch in sample
% this is the list from where you should get the number of trials
xlsx_dir_behav = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/MEG_behavioural_TSA2021'; %dir to MEG behavioral results (.xlsx files)

spm_list = dir('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/Continuous/Resting/spmeeg*.mat');

% test = zeros(length(spm_list),1);
for ii = 1:length(spm_list) %over .mat files
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]); %load spm_list .mat files
    D = D.montage('switch',0);
    dummy = D.fname; %OBS! D.fname does not work, so we need to use a 'dummy' variable instead
    
    %%%%%%%% counting the number of old correct for this subject %%%%%%%%%%
    clear raw_recog
    dumls = ['Subj_' dummy(12:15) '_Block_3.xlsx'];
    [~,~,raw_recog] = xlsread([xlsx_dir_behav '/' dumls]); %excel files
    number_oldC = 0;
    for tt = 2:size(raw_recog,1)
        if strcmp(raw_recog{tt,2}(8:10),'old')
            if raw_recog{tt,3} == 1
                number_oldC = number_oldC + 1;
            end
        end
    end
%     test(ii,1) = number_oldC;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % total_time are all the possible time indexes from where you can
    % extract the trigger for the epoch. I'm removing the first minutes
    % since at the beginning of acquisition data may be always 0
    total_time = D.time(15001:end);
    time_window = floor(length(total_time)/number_oldC); % in sample points
    
    % matrix with the information about the beginning of the epoch in seconds
    trl_sam = zeros(number_oldC,3); %prepare the samples matrix with 0's in all its cells
    % matrix with the information about the beginning of the epoch in points (use sample frequency for the conversion
    trl_sec = zeros(number_oldC,3); %prepare the seconds matrix with 0's in all its cells
    
    for k = 1:number_oldC %over selected triggers
        
        clear temp_time
        %defining the time interval from which the trigger will be extracted
        start_index = 1 +(k-1)*time_window; 
        end_index = start_index + time_window - epoch_end_sam; % remove the length of 1 epoch to make sure that there are no overlaps
        temp_time = total_time(start_index:end_index);
        
        %deftrig(k,1) = 0.012 + trigcor(k,1); %adding a 0.012 seconds delay to the triggers sent during the experiment (this delay was due to technical reasons related to the stimuli)
        trl_sec(k,1) = temp_time(randi(length(temp_time))); % extracting a random trigger time
        %remove 1000ms of baseline
        trl_sec(k,2) = trl_sec(k,1) + epoch_end ; %end time-window epoch in seconds
        trl_sec(k,3) = trl_sec(k,2) - trl_sec(k,1); %range time-windows in seconds
        trl_sam(k,1) = round(trl_sec(k,1) * 250) + 1; %beginning time-window epoch in samples %250Hz per second
        trl_sam(k,2) = round(trl_sec(k,2) * 250) + 1; %end time-window epoch in samples
        trl_sam(k,3) = -25; %sample before the onset of the stimulus (corresponds to 0.100ms)
    end
    dif = trl_sam(:,2) - trl_sam(:, 1); %difference between the end and the beginning of each sample (just to make sure that everything is fine)
    if ~all(dif == dif(1)) %checking if every element of the vector are the same (i.e. the length of the trials is the same; we may have 1 sample of difference sometimes because of different rounding operations..)
        trl_sam(:,2) = trl_sam(:,1) + dif(1);
    end
    %creates the epochinfo structure that is required for the source reconstruction later
    epochinfo.trl = trl_sam;
    epochinfo.time_continuous = D.time;
    %switch the montage to 0 because for some reason OSL people prefer to do the epoching with the not denoised data
    D = D.montage('switch',0);
    %build structure for spm_eeg_epochs
    S = [];
    S.D = D;
    S.trl = trl_sam;
    S.prefix = prefix_tobeadded;
    D = spm_eeg_epochs(S);
    %store the epochinfo structure inside the D object
    D.epochinfo = epochinfo;
    D = D.montage('switch',1);
    D.epochinfo.conditionlabels = D.conditions; %to add for later use in the source reconstruction
    D.save();
    disp(['Subject ' num2str(ii) ' is done'])
end


%% NO TO BE USED
% this section was made bc we believed that we had to write the fake
% triggers in the spm object for epoching the resting state. Actually it's
% not necessary and the section is still here for future reference



% %%% HACKERAGGIO %%%
% 
% %%% HERE YOU NEED TO MATCH THE NUMBER OF TRIALS BETWEEN MEMORY AND RESTING %%%
% %%% YOU NEED TO LOAD THE EXCEL FILES AND CHECK HOW MANY CORRECT TRIALS WERE PRESENT %%%
% 
% 
% %list of resting state data
% spm_list = dir('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/Continuous/Resting/spmeeg*.mat');
% %list of memory data
% spm_list2 = dir('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/Continuous/Memory/spmeeg*.mat');
% workdir = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/Continuous/RestingH'; %new working directory for resting data
% 
% for ii = 1:length(spm_list) %over .mat files
%     D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]); %load spm_list .mat files (resting data)
%     D2 = spm_eeg_load([spm_list2(ii).folder '/' spm_list2(ii).name]); %load spm_list .mat files (memory data)
%     %cloning SPM object
%     [workdiroriginal, fname, ext] = fileparts(D2.fullfile); %getting information about the path of the memory data
%     Dnew = D2.clone( fullfile(workdir,['RestH' spm_list(ii).name]) ); %cloning memory data
%     data = D(:,:,:); %extracting data from resting state
%     Dnew(:,:,:) = 0; %zeroing cloning data to be sure that the "extra" data from memory data does not remain there after you paste the resting state data
%     Dnew(:,1:size(data,2),:) = data; %pasting the resting state data into the cloned file
%     Dnew.save(); %saving the clone file
% end
% 
% %%% AFTER THIS PROCEDURE THE FINAL RESTING STATE DATA 
% 
% %%% REMEMBER THAT YOU NEED TO USE WHEN EPOCHING THE SAME NUMBER OF TRIALS
% %%% USED FOR THE CORRESPONIDNG MEMORY DATA (FOR EACH SUBJECT) AND THAT THIS
% %%% DATA MUST NOT BE TAKEN FOR THE "0s" !! %%%



%%


%%

%% *** SOURCE RECONSTRUCTION ***

%%
%% LBPD_startup_D

pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions
addpath(pathl);
LBPD_startup_D(pathl);

%% CREATING 8mm PARCELLATION FOR EASIER INSPECTION IN FSLEYES

%%% *** (NOT TO BE RUN!!) *** %%%

%OBS!! This section is done only for better handling of some visualization purposes, but it does not affect any of the beamforming algorithm;
% it is just important not to mix up the MNI coordinates, thus I would recommend to use the following lines

%1) USE load_nii TO LOAD A PREVIOUS NIFTI IMAGE
imag_8mm = load_nii('/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/MNI152_T1_8mm_brain.nii.gz');
Minfo = size(imag_8mm.img); %get info about the size of the original image
M8 = zeros(Minfo(1), Minfo(2), Minfo(3)); %Initialize an empty matrix with the same dimensions as the original .nii image
cc = 0; %set a counter
M1 = imag_8mm.img;
for ii = 1:Minfo(1) %loop across each voxel of every dimension
    for jj = 1:Minfo(2)
        for zz = 1:Minfo(3)
            if M1(ii,jj,zz) ~= 0 %if we have an actual brain voxel
                cc = cc+1;
                M8(ii,jj,zz) = cc;
            end
        end
    end
end
%2) PUT YOUR MATRIX IN THE FIELD ".img"
imag_8mm.img = M8; %assign values to new matrix 
%3) SAVE NIFTI IMAGE USING save_nii
save_nii(imag_8mm,'/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/MNI152_8mm_brain_diy.nii.gz');
%4) USE FSLEYES TO LOOK AT THE FIGURE
%Create parcellation on the 8mm template
for ii = 1:3559 %for each 8mm voxel
    cmd = ['fslmaths /scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/MNI152_8mm_brain_diy.nii.nii.gz -thr ' num2str(ii) ' -uthr ' num2str(ii) ' -bin /scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/AAL_80mm_3559ROIs/' num2str(ii) '.nii.gz'];
    system(cmd)
    disp(ii)
end
%5) GET MNI COORDINATES OF THE NEW FIGURE AND SAVE THEM ON DISK
MNI8 = zeros(3559,3);
for mm = 1:3559 %over brain voxel
    path_8mm = ['/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/parcel_80mm_3559ROIs/' num2str(mm) '.nii.gz']; %path for each of the 3559 parcels
    [mni_coord,pkfo] = osl_mnimask2mnicoords(path_8mm);  %getting MNI coordinates
    MNI8(mm,:) = mni_coord; %storing MNI coordinates
end
%saving on disk
save('/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/MNI152_8mm_coord_dyi.mat', 'MNI8');

%% CONVERSION T1 - DICOM TO NIFTI - DONE

% THIS HAS BEEN ALREADY RUN, YOU NOW HAVE THIS DATA. 
% DON'T DO IT  AGAIN

%%% BLOCK 3 = MEMORY TASK
%%% BLOCK 4 = RESTING STATE

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/dicm2nii'); %adds path to the dcm2nii folder in osl
MRIsubj = dir('/projects/MINDLAB2021_MEG-TempSeqAges/raw/0*');
MRIsubj(79:end) = [];
MRIoutput = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/Source_LBPD/MRI_nifti';
mkdir(MRIoutput);
MRIout_block{1} = 'Block_1'; MRIout_block{2} = 'Block_2'; MRIout_block{3} = 'Block_3'; MRIout_block{4} = 'Block_4'; MRIout_block{5} = 'Block_5'; MRIout_block{6} = 'Block_6';

for bb = 3%:4 %:length(MRIout_block) %over experimental blocks
    for ii = 38 %1:length(MRIsubj) %over subjects
        asd = [MRIoutput '/Block_' num2str(bb) '/' MRIsubj(ii).name];
        if ~exist(asd,'dir') %checking whether the directory exists
            mkdir(asd); %if not, creating it
        end
        if strcmp(MRIsubj(ii).name,'0072') %we previously obtained subject 0072's MRI so we convert it from a different folder
            dcmSource = '/raw/sorted/MINDLAB2020_MEG-AuditoryPatternRecognition/0011/20201208_173601/MR/005.t1_mprage_3D_sag_fatsat/files/';
            niiFolder = asd;
            dicm2nii(dcmSource, niiFolder, '.nii');
        else %actual subjects (only) from the new data collection
            if isempty(dir([asd '/*.nii'])) %if there are no nifti images.. I need to convert them
                flagg = 0;
                MRIMEGdate = dir([MRIsubj(ii).folder '/' MRIsubj(ii).name '/20*']);
                niiFolder = [MRIoutput '/' MRIout_block{bb} '/' MRIsubj(ii).name];
                for jj = 1:length(MRIMEGdate) %over dates of recording
                    %                     if ~isempty(dir([MRIMEGdate(jj).folder '/' MRIMEGdate(jj).name '/MR*'])) %if we get an MRI recording
                    MRI2 = dir([MRIMEGdate(jj).folder '/' MRIMEGdate(jj).name '/MR/*fatsat']); %looking for T1
                    if ~isempty(MRI2) %if we have one T1
                        flagg = 1; %determining that I could convert MRI T1
                        dcmSource = [MRI2(1).folder '/' MRI2(1).name '/files/'];
                        dicm2nii(dcmSource, niiFolder, '.nii');
%                     elseif length(MRI2) ~= 1 && jj == length(MRIMEGdate)
%                         warning(['subject ' MRIsubj(ii).name ' has no MRI T1 or has more than 1 MRI T1']);
%                         warning('copying brain template..')
%                         copyfile('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_2mm.nii',[niiFolder '/MNI152_T1_2mm.nii'])
                    end
                    %                     else
                    %                         warning(['subject ' MRIsubj(ii).name ' has no MRI T1']);
                    %                         warning('copying brain template..')
                    %                         copyfile('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_2mm.nii',[niiFolder '/MNI152_T1_2mm.nii'])
                    %                     end
                end
%                 if isempty(dir([niiFolder '/*.nii'])) %if something goes wrong with the conversion of the MRI file, I copy-paste a template
%                     copyfile('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_2mm.nii',[niiFolder '/MNI152_T1_2mm.nii'])
%                 end
            end
        end
        disp(ii)
    end
end

for bb = 3:4
    subfolder = dir([MRIoutput '/Block_' num2str(bb) '/00*' ]);
    for ii = 1:length(subfolder) %over subjects
        niiFolder = [subfolder(ii).folder '/' subfolder(ii).name];
        if isempty(dir([niiFolder '/*.nii'])) %if something goes wrong with the conversion of the MRI file, I copy-paste a template
            copyfile('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_2mm.nii',[niiFolder '/MNI152_T1_2mm.nii'])
            disp([ 'Copying in block ' num2str(bb) ' subject ' num2str(ii)])
        end
    %     bum = dir([niiFolder '/*.txt']);
    %     for pp = 1:length(bum)
    %         delete([bum(pp).folder '/' bum(pp).name]);
    %     end
        
    end
end

%% SETTING FOR CLUSTER (PARALLEL COMPUTING)

clusterconfig('scheduler', 'none'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
% clusterconfig('scheduler', 'cluster'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('long_running', 1); % This is the cue we want to use for the clsuter. There are 3 different cues. Cue 0 is the short one, which should be enough for us
clusterconfig('slot', 1); %slot is memory, and 1 memory slot is about 8 GB. Hence, set to 2 = 16 GB
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing')

%% RHINO coregistration - DONE 
% No need to run this again even if you have to perform again source reconstruction, the information remains in the data in the
% field D.inv
%%% FROM HERE, CHAIRA.. CHANG DIRECTORIES AND CHECK.. SIMPLIFY.. ETC.. %%%
%%% RUN THS FOR BOTH BLOCK 3 (MEMORY) AND 4 (RESTING STATE)

%change this directory once you have the epoched files
% list = dir ('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/e*cogminor*_tsssdsm.mat'); %dir to epoched files (encoding)
extr2 = 13:16;

%running rhino
%OBS! check that all MEG data are in the same order and number as MRI nifti files!
for block = 3:4 % over the two blocks I'm interested in (just 3 = memory and 4 = resting)
    clear list
    a = ['/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/Source_LBPD/MRI_nifti/Block_' num2str(block)]; %set path to MRI subjects' folders
    if block == 3 % MEMORY
        list = dir('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/Continuous/Memory/espmeeg*.mat');
    else % the nested if is not necessary but I want to be sure and avoid mistakes
        if block == 4 % RESTING
            list = dir('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/Continuous/Resting/espmeeg*.mat');
        end
    end
    for ii = 1:length(list) %OBS! change this depending on atonal vs. major
        S = [];
        S.ii = ii;
        S.D = [list(ii).folder '/' list(ii).name]; %path to major files
        D = spm_eeg_load(S.D);

        if isfield(D,'inv') %checking if the coregistration was already run
            D.inv = [];
            D.save();
            if isempty(D.inv)
                dummyname = D.fname;
                % exist output is 7 if NAME is a folder
                if exist([a '/' dummyname(extr2)],'dir') == 7 %if you have the MRI folder
                    dummymri = dir([a '/' dummyname(extr2) '/*.nii']); %path to nifti files (ending with .nii)
                    if ~isempty(dummymri)
                        S.mri = [dummymri(1).folder '/' dummymri(1).name];
                        %standard parameters
                        S.useheadshape = 1;
                        S.use_rhino = 1; %set 1 for rhino, 0 for no rhino
                        %         S.forward_meg = 'MEG Local Spheres';
                        S.forward_meg = 'Single Shell'; %CHECK WHY IT SEEMS TO WORK ONLY WITH SINGLE SHELL!!
                        S.fid.label.nasion = 'Nasion';
                        S.fid.label.lpa = 'LPA';
                        S.fid.label.rpa = 'RPA';
                        jobid = job2cluster(@coregfunc,S); %running with parallel computing
                    else
                        warning(['subject ' dummyname(extr2) ' does not have the MRI'])
                    end
                end
            end
%         else
%             if isempty(D.inv{1}) %checking whether the coregistration was run but now it is empty..
%                 warning(['subject ' D.fname ' has an empty rhino..']);
%             end
        end
        disp(ii)
   
    end
end

%% checking (or copying) RHINO - DONE
% No need to run this again even if you have to perform again source reconstruction, the information remains in the data in the
% field D.inv
%%% RUN THS FOR BOTH BLOCK 3 (MEMORY) AND 4 (RESTING STATE)

copy_label = 1; % 1 = pasting inv RHINO from epoched data (where it was computed) to continuous data (not supported for block 6); 0 = simply showing RHINO coregistration
folders = {'Memory','Resting'};

for ff = 1:length(folders)
    path = ['/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/Continuous/' folders{1,ff}];
    list = dir ([path '/espmeeg*.mat']); %dir to epoched files (encoding)
    extr = 13:16;
    for ii = 1:length(list)
        D = spm_eeg_load([list(ii).folder '/' list(ii).name]);
        if isfield(D,'inv')
            if copy_label == 0 %simply displaying RHINO coregistration
                if isfield(D,'inv') %checking if the coregistration was already run
                    rhino_display(D)
                end
            else %pasting inv RHINO from epoched data (where it was computed) to continuous data
                if isfield(D,'inv') %checking if the coregistration was already run
                    inv_rhino = D.inv;
                end
                D2 = spm_eeg_load([list(ii).folder '/' list(ii).name(2:end)]); %loading spm files by leaving out the 'e'
                D2.inv = inv_rhino;
                D2.save();
            end
        end
        disp([ folders{1,ff} ' - Subject ' num2str(ii)])
    end
end
%block 4 subj 1

%%
%% SETTING UP FOR SOURCE RECONSTRUCTION (LBPD functions and cluster settings)
% RUN THIS SECTION BEFORE PERFORMING SOURCE RECONSTRUCTION
%%%%%% *** START UP FUNCTIONS.. (LBPD_startup_D) *** %%%%%%%%%%%%%
%starting up some functions for LBPD toolbox.

pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions
addpath(pathl);
LBPD_startup_D(pathl);

%%%%%%%%%%%%% SETTING FOR CLUSTER (PARALLEL COMPUTING) %%%%%%%%%%%%%%%%%%%%%%%%%%%

% clusterconfig('scheduler', 'none'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('scheduler', 'cluster'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('long_running', 1); % This is the cue we want to use for the clsuter. There are 3 different cues. Cue 0 is the short one, which should be enough for us
clusterconfig('slot', 1); %slot is memory, and 1 memory slot is about 8 GB. Hence, set to 2 = 16 GB
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing');

%% NARROWBAND SOURCE RECONSTRUCTION -BOTH MEMORY AND RESTING
% The cluster sometimes lose some jobs. This code can be run as many times
% as needed to complete source reconstruction for each frequency and each
% subject. If you run this section more than once, subjects that have
% already been completed will be skipped.


%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFINE THE FREQUENCY ARRAY %%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT AND DEFINE THE FUNCTION FOR THE AMPLITUDE CORRESPONDING TO EACH FREQUENCY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
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
figure;
delta_freq = ampl(central_freq); % array containing the interval amplitude corresponding to each central frequency
plot(central_freq, delta_freq, 'LineWidth', 1,'marker', '*', 'color', 'r');
grid minor;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILTERING AND SOURCE RECONSTRUCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


folders = {'Memory','Resting'};

for jj = 1:length(folders) % over memory ad resting
    disp(folders{1,jj})
    % user settings
    clust_l = 1; %1 = using cluster of computers (CFIN-MIB, Aarhus University); 0 = running locally
    timek = []; %time-points, empty bc we are using the whole epoch  
    sensl = 1; %1 = magnetometers only; 2 = gradiometers only; 3 = both magnetometers and gradiometers (SUGGESTED 1!)
    workingdir2 = ['/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/Source_rec_data/' folders{1,jj}]; %high-order working directory (a subfolder for each analysis with information about frequency, time and absolute value will be created)
    %block = 3; % 3 = recogminor; 4 = recogspeed; 5 = samemel; 6 = visualpat; 7 = pd
    invers = 5; %5; %1-4 = different ways (e.g. mean, t-values, etc.) to aggregate trials and then source reconstruct only one trial; 5 for single trial independent source reconstruction

    absl = 0; % 1 = absolute value of sources; 0 = not

    %actual computation 
    
    %import continuous data
    list_c = dir(['/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/Continuous/' folders{1,jj} '/spmeeg*.mat' ]); %dir to continuos files
    
    if strcmp(folders{1,jj}, 'Memory')
        condss = {'Old_Correct'}; %selecting only the old correct condition
    else 
        condss = {'Undefined'};
    end
    %load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat');
    load('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/time_normal.mat');

    % import epoched data
    list_e = dir (['/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/Continuous/' folders{1,jj} '/espmeeg*.mat' ]); %dir to new epoched files



    for ff = 1:size(central_freq,2) %iterate over frequencies
        freqq = [central_freq(1,ff)-(delta_freq(1,ff)/2), central_freq(1,ff)+(delta_freq(1,ff)/2)];
        if ff < 10
            workingdir = [workingdir2 '/Beam_abs_' num2str(absl) '_sens_' num2str(sensl) '_freq_000' num2str(ff) '_' num2str(freqq(1)) '_' num2str(freqq(2)) '_invers_' num2str(invers)];
        else
            workingdir = [workingdir2 '/Beam_abs_' num2str(absl) '_sens_' num2str(sensl) '_freq_00' num2str(ff) '_' num2str(freqq(1)) '_' num2str(freqq(2)) '_invers_' num2str(invers)];
        end

        if ~exist(workingdir,'dir') %creating working folder if it does not exist
            mkdir(workingdir)
        end

        computed_subjs = dir([workingdir '/SUBJ*.mat']);
        if length(computed_subjs) == length(list_e) % if all the subjects are already there skip this folder
            disp(['Skipping folder ' num2str(ff) ' because it is already complete'])
            continue
        else

            allsubj_IDs = cell(length(list_e),1);
            for s = 1:length(list_e)
                allsubj_IDs{s,1} = list_e(s).name(13:16);
            end
            somesubj_IDs = cell(length(computed_subjs),1);
            for s = 1:length(computed_subjs)
                somesubj_IDs{s,1} = computed_subjs(s).name(6:9);
            end
            clear missing_subj
            clear missing_index
            [missing_subj, missing_index ] = setdiff(allsubj_IDs, somesubj_IDs, 'stable');
        end

        for ii = 1:length(missing_index) %over missing subjects
            S = [];
            S.norm_megsensors.MEGdata_c = [list_c(missing_index(ii)).folder '/' list_c(missing_index(ii)).name]; %continuous file %perchè non il nome intero?

            S.Aarhus_cluster = clust_l; %1 for parallel computing; 0 for local computation

            S.norm_megsensors.zscorel_cov = 1; % 1 for zscore normalization; 0 otherwise
            S.norm_megsensors.workdir = workingdir;
            S.norm_megsensors.MEGdata_e = [list_e(missing_index(ii)).folder '/' list_e(missing_index(ii)).name]; %epoched files
            S.norm_megsensors.freq = freqq; %frequency range
            S.norm_megsensors.forward = 'Single Shell'; %forward solution (for now better to stick to 'Single Shell')

            S.beamfilters.sensl = sensl; %1 = magnetometers; 2 = gradiometers; 3 = both MEG sensors (mag and grad) (SUGGESTED 3!)
            S.beamfilters.maskfname = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz'; % path to brain mask: (e.g. 8mm MNI152-T1: '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz')

            S.inversion.znorml = 0; % 1 for inverting MEG data using the zscored normalized one; (SUGGESTED 0 IN BOTH CASES!)
            %                                 0 to normalize the original data with respect to maximum and minimum of the experimental conditions if you have both magnetometers and gradiometers.
            %                                 0 to use original data in the inversion if you have only mag or grad (while e.g. you may have used zscored-data for covariance matrix)
            %
            S.inversion.timef = timek; %data-points to be extracted (e.g. 1:300); leave it empty [] for working on the full length of the epoch
            S.inversion.conditions = condss; %cell with characters for the labels of the experimental conditions (e.g. {'Old_Correct','New_Correct'})
            S.inversion.bc = [1 26]; %extreme time-samples for baseline correction (leave empty [] if you do not want to apply it)
            S.inversion.abs = absl; %1 for absolute values of sources time-series (recommendnded 1!)
            S.inversion.effects = invers;

            S.smoothing.spatsmootl = 0; %1 for spatial smoothing; 0 otherwise
            S.smoothing.spat_fwhm = 100; %spatial smoothing fwhm (suggested = 100)
            S.smoothing.tempsmootl = 0; %1 for temporal smoothing; 0 otherwise
            S.smoothing.temp_param = 0.01; %temporal smoothing parameter (suggested = 0.01)
            S.smoothing.tempplot = [1 2030 3269]; %vector with sources indices to be plotted (original vs temporally smoothed timeseries; e.g. [1 2030 3269]). Leave empty [] for not having any plot.

            S.nifti = 1; %1 for plotting nifti images of the reconstructed sources of the experimental conditions
            S.out_name = ['SUBJ_' list_e(missing_index(ii)).name(13:16)]; %name (character) for output nifti images (conditions name is automatically detected and added)

            if clust_l ~= 1 %useful  mainly for begugging purposes
                MEG_SR_Beam_LBPD(S);
            else
                jobid = job2cluster(@MEG_SR_Beam_LBPD,S); %running with parallel computing
            end
        end
    end

end

%% BROADBAND SOURCE RECONSTRUCTION -BOTH MEMORY AND RESTING
% The cluster sometimes loose some jobs. This code can be run as many times
% as needed to complete source reconstruction for each subject. 
%Subjects that have already been completed will be skipped.

folders = {'Memory','Resting'};

for jj = 1:length(folders) % over memory ad resting
    disp(folders{1,jj})
    % user settings
    clust_l = 1; %1 = using cluster of computers (CFIN-MIB, Aarhus University); 0 = running locally
    timek = []; %time-points, empty bc we are using the whole epoch  
    sensl = 1; %1 = magnetometers only; 2 = gradiometers only; 3 = both magnetometers and gradiometers (SUGGESTED 1!)
    workingdir2 = ['/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/Source_rec_data/' folders{1,jj}]; %high-order working directory (a subfolder for each analysis with information about frequency, time and absolute value will be created)
    %block = 3; % 3 = recogminor; 4 = recogspeed; 5 = samemel; 6 = visualpat; 7 = pd
    invers = 5; %5; %1-4 = different ways (e.g. mean, t-values, etc.) to aggregate trials and then source reconstruct only one trial; 5 for single trial independent source reconstruction

    absl = 0; % 1 = absolute value of sources; 0 = not

    %actual computation 
    %import continuous data

    list_c = dir(['/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/Continuous/' folders{1,jj} '/spmeeg*.mat' ]); %dir to continuos files
    
    if strcmp(folders{1,jj}, 'Memory')
        condss = {'Old_Correct'}; %selecting only the old correct condition
    else 
        condss = {'Undefined'};
    end
    %load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat');
    load('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/time_normal.mat');

    % import epoched data
    list_e = dir (['/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/Continuous/' folders{1,jj} '/espmeeg*.mat' ]); %dir to new epoched files

    freqq = []; % considering broadband data

    
        
    workingdir = [workingdir2 '/Beam_abs_' num2str(absl) '_sens_' num2str(sensl) '_BROADBAND_invers_' num2str(invers)];
        

    if ~exist(workingdir,'dir') %creating working folder if it does not exist
        mkdir(workingdir)
    end

    computed_subjs = dir([workingdir '/SUBJ*.mat']);
    if length(computed_subjs) == length(list_e) % if all the subjects are already there skip this folder
        disp(['Skipping ' folders{1,jj} ' BROAD BAND folder because it is already complete'])
        continue
    else

        allsubj_IDs = cell(length(list_e),1);
        for s = 1:length(list_e)
            allsubj_IDs{s,1} = list_e(s).name(13:16);
        end
        somesubj_IDs = cell(length(computed_subjs),1);
        for s = 1:length(computed_subjs)
            somesubj_IDs{s,1} = computed_subjs(s).name(6:9);
        end
        clear missing_subj
        clear missing_index
        [missing_subj, missing_index ] = setdiff(allsubj_IDs, somesubj_IDs, 'stable');
    end

    for ii = 1:length(missing_index) %over missing subjects
        S = [];
        S.norm_megsensors.MEGdata_c = [list_c(missing_index(ii)).folder '/' list_c(missing_index(ii)).name]; %continuous file %perchè non il nome intero?

        S.Aarhus_cluster = clust_l; %1 for parallel computing; 0 for local computation

        S.norm_megsensors.zscorel_cov = 1; % 1 for zscore normalization; 0 otherwise
        S.norm_megsensors.workdir = workingdir;
        S.norm_megsensors.MEGdata_e = [list_e(missing_index(ii)).folder '/' list_e(missing_index(ii)).name]; %epoched files
        S.norm_megsensors.freq = freqq; %frequency range
        S.norm_megsensors.forward = 'Single Shell'; %forward solution (for now better to stick to 'Single Shell')

        S.beamfilters.sensl = sensl; %1 = magnetometers; 2 = gradiometers; 3 = both MEG sensors (mag and grad) (SUGGESTED 3!)
        S.beamfilters.maskfname = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz'; % path to brain mask: (e.g. 8mm MNI152-T1: '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz')

        S.inversion.znorml = 0; % 1 for inverting MEG data using the zscored normalized one; (SUGGESTED 0 IN BOTH CASES!)
        %                                 0 to normalize the original data with respect to maximum and minimum of the experimental conditions if you have both magnetometers and gradiometers.
        %                                 0 to use original data in the inversion if you have only mag or grad (while e.g. you may have used zscored-data for covariance matrix)
        %
        S.inversion.timef = timek; %data-points to be extracted (e.g. 1:300); leave it empty [] for working on the full length of the epoch
        S.inversion.conditions = condss; %cell with characters for the labels of the experimental conditions (e.g. {'Old_Correct','New_Correct'})
        S.inversion.bc = [1 26]; %extreme time-samples for baseline correction (leave empty [] if you do not want to apply it)
        S.inversion.abs = absl; %1 for absolute values of sources time-series (recommendnded 1!)
        S.inversion.effects = invers;

        S.smoothing.spatsmootl = 0; %1 for spatial smoothing; 0 otherwise
        S.smoothing.spat_fwhm = 100; %spatial smoothing fwhm (suggested = 100)
        S.smoothing.tempsmootl = 0; %1 for temporal smoothing; 0 otherwise
        S.smoothing.temp_param = 0.01; %temporal smoothing parameter (suggested = 0.01)
        S.smoothing.tempplot = [1 2030 3269]; %vector with sources indices to be plotted (original vs temporally smoothed timeseries; e.g. [1 2030 3269]). Leave empty [] for not having any plot.

        S.nifti = 1; %1 for plotting nifti images of the reconstructed sources of the experimental conditions
        S.out_name = ['SUBJ_' list_e(missing_index(ii)).name(13:16)]; % name (character) for output nifti images (conditions name is automatically detected and added)

        if clust_l ~= 1 %useful  mainly for begugging purposes
            MEG_SR_Beam_LBPD(S);
        else
            jobid = job2cluster(@MEG_SR_Beam_LBPD,S); %running with parallel computing
        end
    end
    

end
%% CONFIGURATION FOR LBPD AND CLUSTER - TO BE RUN BEFORE GED

addpath('/projects/MINDLAB2023_MEG-AuditMemDement/scripts/chiaramalvaso/GED_scripts');
addpath('/projects/MINDLAB2023_MEG-AuditMemDement/scripts/chiaramalvaso/GED_TSA2021');
addpath('/projects/MINDLAB2023_MEG-AuditMemDement/scripts/chiaramalvaso/GED_scripts2')

%%%%%%%%%%%%%%%% *** START UP FUNCTIONS.. (LBPD_startup_D) *** %%%%%%%%%%%%%%%% 
 
%starting up some functions for LBPD toolbox.

pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions
addpath(pathl);
LBPD_startup_D(pathl);
%%%%%%%%%%%%%%%%%%%% SETTING FOR CLUSTER (PARALLEL COMPUTING) %%%%%%%%%%%%%%%%%%%% 

clusterconfig('scheduler', 'none'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
% clusterconfig('scheduler', 'cluster'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('long_running', 1); % This is the cue we want to use for the clsuter. There are 3 different cues. Cue 0 is the short one, which should be enough for us
clusterconfig('slot', 2); %slot is memory, and 1 memory slot is about 8 GB. Hence, set to 2 = 16 GB
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing');

%% GED...
% The following sections include the code for three different ways of
% performing GED, each one performed for each single subject: 
%       - covariance matrix computed on concatenated trials
%       - covariance matrix computed on single trials and then averaged
%       - covariance matrix on single trials removing the outlier: 
%           the covariance matrix is computed independently for each trial,
%           averaged and then cov matrix that are more distant than three
%           standard deviation from the mean are discarded. The new
%           averaged covariance matrix with no outliers is computed and
%           that is the one used for calculations. The procedure is done
%           both for Reference and Signal independently.
% Since the cluster may loose some jobs, the following code for GED can be
% run as many times as needed to compute GED for each subject and
% frequency. If you have to run it more than once, subjects that have
% already been computed are skipped.
% 
%% GED ON CONCATENATED TRIALS - BOTH MEMORY AND LISTENING
folders = {'Memory','Resting'};

for jj = 1:length(folders) % over memory ad resting
    disp(folders{1,jj})
    destpath = ['/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/GED_results/' folders{1,jj}];
    freq_list = dir(['/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/Source_rec_data/' folders{1,jj} '/Beam_abs*freq*']);
    for ff=1:length(freq_list)
        freq_index = strfind(freq_list(ff).name, '_invers') -1;
        destdir= [destpath '/' freq_list(ff).name(19:freq_index)];
        if ~exist(destdir,'dir') %creating working folder if it does not exist
            mkdir(destdir)
        end
        % List of broad-band subjects
        sub_list_B = dir(['/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/Source_rec_data/' folders{1,jj} '/Beam_abs*BROAD*/SUBJ*.mat']);
        
        % List of narrow-band subjects
        sub_list_N = dir([freq_list(ff).folder '/' freq_list(ff).name '/SUBJ*.mat']); % complete list of subjects from the input directory
        
        computed_subjs = dir([destdir '/SUBJ*.mat']); % subjects already present in the folder with the results
        
        if length(computed_subjs) == length(sub_list_N) % if all the subjects are already there skip this folder
            disp(['Skipping folder ' num2str(ff) ' because it is already complete'])
            continue
        else

            allsubj_IDs = cell(length(sub_list_N),1);
            for s = 1:length(sub_list_N)
                allsubj_IDs{s,1} = sub_list_N(s).name(6:9); % check if the index is correct after source rec (it should be)
            end
            somesubj_IDs = cell(length(computed_subjs),1);
            for s = 1:length(computed_subjs)
                somesubj_IDs{s,1} = computed_subjs(s).name(6:9);
            end
            clear missing_subj
            clear missing_index
            [missing_subj, missing_index ] = setdiff(allsubj_IDs, somesubj_IDs, 'stable');
        end
        
        for sub=1:length(missing_subj) % over subjects that still needs to be computed
            S = struct(); % initializing the input structure
            S.filenameNarrow = sub_list_N(missing_index(sub)).name;
            S.pathNarrow = sub_list_N(missing_index(sub)).folder;
            S.filenameBroad = sub_list_B(missing_index(sub)).name;
            S.pathBroad = sub_list_B(missing_index(sub)).folder;
            S.destpath = destdir;
            S.ncomps = 10;
            jobid = job2cluster(@GED_single_subject_LandR_separately,S); %running with parallel computing
        end
        
    end   
end
%% GED WITH COVARIANCE ON SINGLE TRIALS - BOTH MEMORY AND LISTENING
folders = {'Memory','Resting'};

for jj = 1:length(folders) % over memory ad resting
    disp(folders{1,jj})
    destpath = ['/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/GED_covsingletrial_results/' folders{1,jj}];
    freq_list = dir(['/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/Source_rec_data/' folders{1,jj} '/Beam_abs*freq*']);
    for ff=1:length(freq_list)
        freq_index = strfind(freq_list(ff).name, '_invers') -1;
        destdir= [destpath '/' freq_list(ff).name(19:freq_index)];
        if ~exist(destdir,'dir') %creating working folder if it does not exist
            mkdir(destdir)
        end
        % List of broad-band subjects
        sub_list_B = dir(['/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/Source_rec_data/' folders{1,jj} '/Beam_abs*BROAD*/SUBJ*.mat']);
        
        % List of narrow-band subjects
        sub_list_N = dir([freq_list(ff).folder '/' freq_list(ff).name '/SUBJ*.mat']); % complete list of subjects from the input directory
        
        computed_subjs = dir([destdir '/SUBJ*.mat']); % subjects already present in the folder with the results
        
        if length(computed_subjs) == length(sub_list_N) % if all the subjects are already there skip this folder
            disp(['Skipping folder ' num2str(ff) ' because it is already complete'])
             continue
        else

            allsubj_IDs = cell(length(sub_list_N),1);
            for s = 1:length(sub_list_N)
                allsubj_IDs{s,1} = sub_list_N(s).name(6:9); % check if the index is correct after source rec (it should be)
            end
            somesubj_IDs = cell(length(computed_subjs),1);
            for s = 1:length(computed_subjs)
                somesubj_IDs{s,1} = computed_subjs(s).name(6:9);
            end
            clear missing_subj
            clear missing_index
            [missing_subj, missing_index ] = setdiff(allsubj_IDs, somesubj_IDs, 'stable');
        end
        
        for sub=1:length(missing_subj) % over subjects that still needs to be computed
            S = struct(); % initializing the input structure
            S.filenameNarrow = sub_list_N(missing_index(sub)).name;
            S.pathNarrow = sub_list_N(missing_index(sub)).folder;
            S.filenameBroad = sub_list_B(missing_index(sub)).name;
            S.pathBroad = sub_list_B(missing_index(sub)).folder;
            S.destpath = destdir;
            S.ncomps = 10;
            jobid = job2cluster(@GED_single_subject_singletrial,S); %running with parallel computing
        end
        
    end

end
%% GED WITH COVARIANCE ON SINGLE TRIALS AND OUTLIERS REMOVAL - BOTH MEMORY AND LISTENING
folders = {'Memory','Resting'};

for jj = 1:length(folders) % over memory ad resting
    disp(folders{1,jj})
    destpath = ['/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/GED_covsingletrial_outlier_results/' folders{1,jj}];
    if ~exist(destpath,'dir') %creating working folder if it does not exist
        mkdir(destpath)
    end
    freq_list = dir(['/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/Source_rec_data/' folders{1,jj} '/Beam_abs*freq*']);
    for ff=1:length(freq_list)
        freq_index = strfind(freq_list(ff).name, '_invers') -1;
        destdir= [destpath '/' freq_list(ff).name(19:freq_index)];
        if ~exist(destdir,'dir') %creating working folder if it does not exist
            mkdir(destdir)
        end
        % List of broad-band subjects
        sub_list_B = dir(['/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/Source_rec_data/' folders{1,jj} '/Beam_abs*BROAD*/SUBJ*.mat']);
        
        % List of narrow-band subjects
        sub_list_N = dir([freq_list(ff).folder '/' freq_list(ff).name '/SUBJ*.mat']); % complete list of subjects from the input directory
        
        computed_subjs = dir([destdir '/SUBJ*.mat']); % subjects already present in the folder with the results
        
        if length(computed_subjs) == length(sub_list_N) % if all the subjects are already there skip this folder
            disp(['Skipping folder ' num2str(ff) ' because it is already complete'])
             continue
        else

            allsubj_IDs = cell(length(sub_list_N),1);
            for s = 1:length(sub_list_N)
                allsubj_IDs{s,1} = sub_list_N(s).name(6:9); % check if the index is correct after source rec (it should be)
            end
            somesubj_IDs = cell(length(computed_subjs),1);
            for s = 1:length(computed_subjs)
                somesubj_IDs{s,1} = computed_subjs(s).name(6:9);
            end
            clear missing_subj
            clear missing_index
            [missing_subj, missing_index ] = setdiff(allsubj_IDs, somesubj_IDs, 'stable');
        end
        
        for sub=1:length(missing_subj) % over subjects that still needs to be computed
            S = struct(); % initializing the input structure
            S.filenameNarrow = sub_list_N(missing_index(sub)).name;
            S.pathNarrow = sub_list_N(missing_index(sub)).folder;
            S.filenameBroad = sub_list_B(missing_index(sub)).name;
            S.pathBroad = sub_list_B(missing_index(sub)).folder;
            S.destpath = destdir;
            S.ncomps = 10;
            jobid = job2cluster(@GED_single_subject_singletrial_outlier,S); %running with parallel computing
        end
        
    end
 
end


