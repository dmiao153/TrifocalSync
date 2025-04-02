%%% Need to be in the main folder
%% Main Script for running EPFL experiments
addpath(genpath(pwd))
clear all;
 
rng(1)
s = rng;

dataset = "FountainP11";
% dataset = "CastleP19";
% dataset = "EntryP10";
% dataset = "HerzP25";
% dataset = "HerzP8";
% dataset = "CastleP30";

[K, R_true, t_true, im_names] = EPFL_get_ground_truth(dataset);
t_scale = sqrt(size(t_true,2))/norm(t_true);

%% BLOCK A % % %%% comment Block B and C. 
% test_main_EPFL;
% save("EPFL_code/"+dataset+"_sift_Z.mat")
% % load("EPFL_code/"+dataset+"_sift_Z.mat")
% EPFL_FCC;

%% BLOCK B %%% Need to go to the 2024a version to use one of the functions. Comment both Block A and C. 

% EPFL_two_view_data;

%% BLOCK C Need to come back to the 2021a version to use the graphconncomp. Comment both Block A and B. 
%%% function, after running BLOCK B, just need these 
%  
% EPFL_experiments_all_time;
EPFL_GlueStick_Experiments;
% EPFL_run2view; 

% EPFL_cheirality_check;

