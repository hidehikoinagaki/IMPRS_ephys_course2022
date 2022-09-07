%
% A compilation of analysis scripts for the 2022 IMPRS Ephys course
%
%
% 

%% Cell-based analyses

%%
% 1. Plot rasters
%   each spike is a dot
%   trials are arrayed in the vertical dimension, time in the horizontal dimension
plot_raster

%%
% 2. Estimate and plot mean spike rate for each trial type as sr(R), sr(L)
% 3. Compute selectivity as sr(R) - sr(L)
plot_PSTH_with_selectivity


%%
% 4. Error trial
% Change the name of variable in each function and explore the error trials

%% Session-based analysis

%%
% 5. Grand average (average across trials and neurons) ...
% 6. Grand average selectivity and its absolute value
plot_session_PSTH_with_selectivity


%% Dimensionality reduction

%%
% 7. Perform PCA
% 8. Compute the coding direction (CD) 
plot_PCA_CD
