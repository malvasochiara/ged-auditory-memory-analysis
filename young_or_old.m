function [index_old, index_young] = young_or_old(thresh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the indexes of young and old subjects that have given more than
% thresh correct answers. If thresh is empty, it just get the informations
% about old and young indexes.
% INPUT:
%       -thresh       = (OPTIONAL). Minimum number of correct answers
%                       If a subject ha give less than thresh correct answer, the subject is discarded.
%                       Don't insert a thresh is you want to consider all
%                       subjects.
% OUTPUT:
%       -index_old    = array with the indices of the old participants
%       -index_young  = array with the indices of the young participants

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/groups_Age.mat');
    index_old = cat(2,S.subjs{1,1},S.subjs{1,2});
    index_young = S.subjs{1,3};
    if isempty(thresh)
        disp('No threshold applied, considering all subjects')
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
    end