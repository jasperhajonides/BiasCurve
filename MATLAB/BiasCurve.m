% Bias curve fitting
% Computes circular bias in responses towards for secondary orientation
% here called 'distractors'
%
%   key variables:
%     vector of target orientations [-pi, pi]
%     vector of responses by participants [-pi, pi]
%     vector of distractor orientations [-pi, pi]
%
% Jasper Hajonides, Dec. 2020, and Nicholas Myers, 2018

% read in data from .mat file
dir = pwd;
load('../data_raw/data.mat');

% read in circular stat functions by Philipp Berens, 2009
circ_stats   = [pwd '/CircStat_2012a'];
addpath(circ_stats);

%% Compute biases

% loop over subjects to fit the bias curve for a single participant.
for sb = [1:20];
    
    %select trials 
    incl = data.subject == sb;
    itrl = find(mean(incl,2)==1);

    % now calculate
    config = [];
    config.pbin = 1/4; % window size
    config.nbin = 64; % amount of bins for target-distractor distance
    [biases(sb,:), angmus(sb,:), precs(sb,:)]= calculate_biasCurve(data.response,...
        data.target,...
        data.distractor,...
        itrl,...
        config);

end

config.same_diff_subplot = [1]; % a vector of which the length indicates the amount of conditions, put same numbers in the same plot e.g. [1 1 2 2] will make 2 subplots with 2 lines in each
config.psubs = [1:20]; %amount of subjects
config.rootdir = cd;
config.do_print = 0; % save figure = 1, not save figure = 0;
[y,theta]= Plot_BiasCurve(biases,angmus,config);