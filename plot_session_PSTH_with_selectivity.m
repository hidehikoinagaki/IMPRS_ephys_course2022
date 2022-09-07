% This script plots the grand average PSTH and selectivity (across units
% and trials wtihin one session)
% We exclude fast spiking cells for this analysis
% ephysDataset.cell_type: cell type. 1 is regular spiking cell, 0 is fast
% spiking cell
%
% Plot
% blue: lick right
% red: lick left


%% load data
sessionId = 1; %  ID of session to analyze
timeBin   = 0.001; % time bin for PSTH (sec) 
tAxis     = -3.5:timeBin:2; % tAxis for PSTH
smoothBin = 100; % bin size to smooth PSTH

load('ephysDataset.mat') % load data


%% find the regular spiking units from the session
% sessionIndex should be the assigned sessionId & cell_type needs to be 1 (regular spiking cell)
sessionData = ephysDataset([ephysDataset.sessionIndex]==sessionId & [ephysDataset.cell_type]==1);
numUnit     = length(sessionData); % number of units
numTime     = length(tAxis);       % number of time bins


%% Calculate the mean spike rate of each unit
% predefine matrix of nan's (i.e. 'not a number')
PSTH_right  = nan(numUnit, numTime);
PSTH_left   = nan(numUnit, numTime);
selectivity = nan(numUnit, numTime);

% run a loop to extract data from each unit
for cellID = 1:numUnit

    % calcualte mean PSTH of each unit
    st_right = sessionData(cellID).st_right;
    sr_right = acquireSpikeRatePerTrial(st_right,timeBin,tAxis);
    
    st_left  = sessionData(cellID).st_left;
    sr_left  = acquireSpikeRatePerTrial(st_left,timeBin,tAxis);
    
    PSTH_right(cellID,:) = smooth(mean(sr_right,1),smoothBin);
    PSTH_left(cellID,:)  = smooth(mean(sr_left,1),smoothBin);
    
    % calculate selectivity of each cell
    % Only cells with siginificant selectivity during the delay epoch is included 
    % negative selectivity during the delay epoch is flipped to be positive 
        
    % extarct spike rate during the delay epoch
    delayTimbin = tAxis>-1.3 & tAxis<0; % timbin of delay epoch
    srDelayR    = mean(sr_right(:,delayTimbin),2); % spike rate during delay lick R trial
    srDelayL    = mean(sr_left(:,delayTimbin),2);  % spike rate during delay lick L trial
    
    % ranksum test to check if spike rates are significantly different between two trial types
    p = ranksum(srDelayR,srDelayL);
    
    % if spike rates are significantly different calcualte selectivity
    % reverse the direction if slectivity is negative
    % non selective trials are kept as nan
    if     p < 0.05 && mean(srDelayR) > mean(srDelayL)
        
        selectivity(cellID,:) = PSTH_right(cellID,:) - PSTH_left(cellID,:);
        
    elseif p < 0.05 && mean(srDelayR) < mean(srDelayL)
        
        selectivity(cellID,:) = PSTH_left(cellID,:) - PSTH_right(cellID,:); 
        
    end
end


%% plot the PSTH
figure;set(gcf,'Color','w')
hold on
plot(tAxis,mean(PSTH_right),'b') % mean spike rate among cells
plot(tAxis,mean(PSTH_left),'r')
xlim([-3  1.5]);
gridxy([-2.6 -1.3 0],'Color','k','Linestyle','--') ;
xlabel('Time (s)')
ylabel('Spikes per s')
set(gca,'box','off','tickdir','out','fontsize',16)
hold off



%% plot the selectivity
figure;set(gcf,'Color','w')
hold on
plot(tAxis,nanmean(selectivity),'k') % mean selectivty  among cells
xlim([-3  1.5]);
gridxy([-2.6 -1.3 0],'Color','k','Linestyle','--') ;
xlabel('Time (s)')
ylabel('Selectivity (Spikes per s)')
set(gca,'box','off','tickdir','out','fontsize',16)
hold off





