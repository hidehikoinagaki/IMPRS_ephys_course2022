% This is a script to plot PSTH and selectivity of individual cell
% 
% cellId: determines which cell to plot
%
% Plot
% blue: lick right
% red: lick left


%% load data
cellId    = 1; % cell to plot
timeBin   = 0.001; % time bin for PSTH (sec) 
tAxis     = -3.5:timeBin:2; % tAxis for PSTH
smoothBin = 100; % bin size to smooth PSTH

load('ephysDataset.mat') % load data


%% Calculate the mean spike rate & selectivity
% lick R trials
st_right        = ephysDataset(cellId).st_right; % spike timing in lick R trials
numTrials_right = numel(st_right); % numbef of trials
spkTime_right   = cell2mat(st_right); % pool all trials together

PSTH_tmp   = hist(spkTime_right,tAxis)/numTrials_right/timeBin; % use histogram to calucalte spike rate & convert to spikes per s
PSTH_right = smooth(PSTH_tmp,smoothBin); % smooth PSTH 


% lick L trials
st_left        = ephysDataset(cellId).st_left;
numTrials_left = numel(st_left);
spkTime_left   = cell2mat(st_left);

PSTH_tmp  = hist(spkTime_left,tAxis)/numTrials_left/timeBin;  
PSTH_left = smooth(PSTH_tmp,smoothBin); 

% selectivity
selectivity = PSTH_right - PSTH_left; % contra selectivity: difference in spike rate between two trial types (R - L)

%% plot the PSTH
figure;set(gcf,'Color','w')
hold on
plot(tAxis,PSTH_right,'b')
plot(tAxis,PSTH_left,'r')
xlim([-3  1.5]); % range of X axis
gridxy([-2.6 -1.3 0],'Color','k','Linestyle','--') ; % plot timing of each epoch
xlabel('Time (s)')
ylabel('Spikes per s')
hold off
set(gca,'box','off','tickdir','out','fontsize',16)

%% plot the contra selectivity
figure;set(gcf,'Color','w')
hold on
plot(tAxis,selectivity,'k')
xlim([-3  1.5]); % range of X axis
gridxy([-2.6 -1.3 0],'Color','k','Linestyle','--') ;  % plot timing of each epoch
xlabel('Time (s)')
ylabel('Contra selectivity (Spikes per s)')
set(gca,'box','off','tickdir','out','fontsize',16)



