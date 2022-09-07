% plot rasters 
% 
% This code will create raster plot for a selected single neuron at
% for each behavioral conditions.
%
% blue: lick right trials
% red: lick left trials


% load data
load('ephysDataset.mat')

% cell number (ranged from 1 to length of datasets)
cellId = 1;


% plot
figure;set(gcf,'Color','w');hold on;
title(['Raster plot for cell #' num2str(cellId)]);

numTrial = 0;

% lick R trials
st_right = ephysDataset(cellId).st_right;

for ntrial = 1:length(st_right)
    spkTime = st_right{ntrial};
    numTrial = numTrial + 1;
    plot(spkTime, numTrial * ones(length(spkTime), 1), '.b');
end

% lick L trials
st_left = ephysDataset(cellId).st_left;

for ntrial = 1:length(st_left)
    spkTime = st_left{ntrial};
    numTrial = numTrial + 1;
    plot(spkTime, numTrial * ones(length(spkTime), 1), '.r');
end

gridxy([-2.6 -1.3 0],'Color','k','Linestyle','--') ;
xlim([-3.0  1.5]);
ylim([0.5 numTrial+0.5]);
xlabel('Time from movement (sec)')
ylabel('Trial index')
set(gca,'box','off','tickdir','out','fontsize',16)
