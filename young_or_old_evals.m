function [old, young, index_old, index_young] = young_or_old_evals(freq_path, index_old, index_young, thresh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the explained variance for each frequency and divide them between
% young and old. If you want to discard subjects that gave the correct
% answer less than a certain amount of times you can add the minimum number
% of correct answer needed.
% INPUT:
%       -freq_path    = path to the folder containing all the frequencies folder where you stored the output from GED
%       -index_old    = array with the indices of the old participants
%       -index_young  = array with the indices of the young participants
%       -thresh       = (OPTIONAL). Minimum number of correct answers
%                       If a subject ha give less than thresh correct answer, the subject is discarded.
%                       Don't insert a thresh is you want to consider all
%                       subjects.
% OUTPUT:
%       -old   = variance explained for the old with dim (voxel, frequency, subjects)
%       -young = variance explained for the young with dim (voxel, frequency, subjects)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (nargin == 3) || (nargin == 4 && (isempty(thresh)))
        disp('No threshold applied, considering all subjects')
        thresh = [];
    else
        if (nargin == 4) && (~isempty(thresh))
            disp(['Discarding subjects with less than ' num2str(thresh) ' correct answers'])
        else
            error('Incorrect number of inputs')
        end
    end
    if isempty(thresh)
        freq_list = dir([freq_path '/freq*']);
        old = zeros(3559,length(freq_list), length(index_old));
        young = zeros(3559, length(freq_list), length(index_young));
        for ff = 1:length(freq_list)

            sub_list = dir([freq_list(ff).folder '/' freq_list(ff).name '/SUBJ*.mat']);

            disp(['Loading frequency ' num2str(ff) ' - young'])
            for yy = 1:length(index_young)
                load([sub_list(index_young(yy)).folder '/' sub_list(index_young(yy)).name], 'evals');
                young(:,ff,yy) = evals;
                clear evals
            end

            disp(['Loading frequency ' num2str(ff) ' - old'])
            for oo = 1:length(index_old)
                load([sub_list(index_old(oo)).folder '/' sub_list(index_old(oo)).name], 'evals');
                old(:,ff,oo) = evals;
                clear evals
            end
        end
    else
        %getting information about the number of trials
        list = dir('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/MEG_behavioural_TSA2021/Subj_*_Block_3.xlsx');
        sub_count = 0;
        for sub = 1:length(list)

            clear data;
            clear headers;
            clear count;
            clear D;
            clear D_old;
            [data, headers, raw]= xlsread([list(sub).folder '/' list(sub).name]);
            count = 0;
            for ii = 2:size(headers,1) % skipping the first row with the labels
                if strcmp(headers{ii,2}(8:10), 'old') == 1 % check if it's condition 'old'
                    if data(ii-1,3) == 1 % if the answer is correct

                        count = count +1;
                    end
                end
            end
            if count <= 14
                sub_count = sub_count +1;
                badsubj(sub_count) = sub;
            end
        end

        % check if there are bad subjects among old subjects
        for ii = 1:length(badsubj)
            clear index
            index = find(index_old == badsubj(ii));
            if ~isempty(index)
                index_old(index) = [];
            end
        end
        % check if there are bad subjects among young subjects
        for ii = 1:length(badsubj)
            clear index
            index = find(index_young == badsubj(ii));
            if ~isempty(index)
                index_young(index) = [];
            end
        end
        % loading only good subjects
        freq_list = dir([freq_path '/freq*']);
        old = zeros(3559,length(freq_list), length(index_old));
        young = zeros(3559, length(freq_list), length(index_young));
        for ff = 1:length(freq_list)

            sub_list = dir([freq_list(ff).folder '/' freq_list(ff).name '/SUBJ*.mat']);

            disp(['Loading frequency ' num2str(ff) ' - young'])
            for yy = 1:length(index_young)
                load([sub_list(index_young(yy)).folder '/' sub_list(index_young(yy)).name], 'evals');
                young(:,ff,yy) = evals;
                clear evals
            end

            disp(['Loading frequency ' num2str(ff) ' - old'])
            for oo = 1:length(index_old)
                load([sub_list(index_old(oo)).folder '/' sub_list(index_old(oo)).name], 'evals');
                old(:,ff,oo) = evals;
                clear evals
            end
        end
        
        
        
        
        
    end

end