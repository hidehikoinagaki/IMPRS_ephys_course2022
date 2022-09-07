% plot coding direction
% 
% This code will go through coding direction plot for a selected single
% session for two behavioral conditions.

%% load data
sessionId = 1; %  ID of session to analyze
timeBin   = 0.001; % time bin for PSTH (sec) 
tAxis     = -3.5:timeBin:2; % tAxis for PSTH
smoothBin = 100; % bin size to smooth PSTH
timeMask  = tAxis>-3 & tAxis<1.5; % time range to analyze (and remove hist artifact)
tAxisToPlot = tAxis(timeMask); % tAxis to be used for plot

load('ephysDataset.mat') % load data


%% find the regular spiking units from the session
% sessionIndex should be the assigned sessionId & cell_type needs to be 1 (regular spiking cell)
sessionData = ephysDataset([ephysDataset.sessionIndex]==sessionId & [ephysDataset.cell_type]==1);
numUnit     = length(sessionData); % number of units
numTime     = sum(timeMask);     % number of time bins



%% coding direction

% coding direction is defined as a vector spanned by neurons;
% each element is the difference of the mean activity of neurons in two
% trial type conditions.


PSTH_right_fit  = zeros(numUnit, numTime);
PSTH_left_fit   = zeros(numUnit, numTime);
PSTH_right_test = zeros(numUnit, numTime);
PSTH_left_test  = zeros(numUnit, numTime);

for cellID = 1:numUnit
    % calcualte mean PSTH of each unit
    st_right = sessionData(cellID).st_right;
    sr_right = acquireSpikeRatePerTrial(st_right,timeBin,tAxis);
    numTrials_right = numel(st_right);
    
    st_left  = sessionData(cellID).st_left;
    sr_left  = acquireSpikeRatePerTrial(st_left,timeBin,tAxis);
    numTrials_left  = numel(st_left);
    
    % select randomely half trials to estimate mode and use others for plot
    % (cross validation)
    
    tr_fit_right  = randsample(1:numTrials_right,floor(numTrials_right/2));
    tr_test_right = setdiff(1:numTrials_right,tr_fit_right);
    
    tr_fit_left   = randsample(1:numTrials_left,floor(numTrials_left/2));
    tr_test_left  = setdiff(1:numTrials_left,tr_fit_left);

    PSTH_right_fit(cellID,:) = smooth(mean(sr_right(tr_fit_right,timeMask),1),smoothBin);
    PSTH_left_fit(cellID,:)  = smooth(mean(sr_left(tr_fit_left,timeMask),1),smoothBin);
    
    PSTH_right_test(cellID,:) = smooth(mean(sr_right(tr_test_right,timeMask),1),smoothBin);
    PSTH_left_test(cellID,:)  = smooth(mean(sr_left(tr_test_left,timeMask),1),smoothBin);

end

%% PCA

% remove the mean variable-wise (row-wise)
X = (PSTH_right_fit+PSTH_left_fit)/2;
X = X-repmat(mean(X,2),1,size(X,2));

% calculate eigenvectors (loadings) W, and eigenvalues of the covariance matrix
[W, EvalueMatrix] = eig(cov(X'));
Evalues = diag(EvalueMatrix);

% order by largest eigenvalue
Evalues = Evalues(end:-1:1);
W = W(:,end:-1:1); W=W';  

% Eigen spectrum
figure;set(gcf,'Color','w')
plot(1:numUnit,Evalues/sum(Evalues));
xlim([0 10])
xlabel('PCs');title('Eigen spectrum');ylabel('Variance explained')
set(gca,'box','off','tickdir','out','fontsize',16)

% Principle components

figure;set(gcf,'Color','w','Position',[300 400 1000 350])
for i = 1:3
    subplot(1,3,i);hold on
    plot(tAxisToPlot,W(i,:)*PSTH_right_test,'b')
    plot(tAxisToPlot,W(i,:)*PSTH_left_test,'r')
    title(['PC ',num2str(i)])
    gridxy([-2.6 -1.3 0],'Color','k','Linestyle','--') ;
    xlim([-3.0  1.5]);
    xlabel('Time from movement (sec)')
    ylabel('Activity projected to PC')
    set(gca,'box','off','tickdir','out','fontsize',16)
end


%% CD

diffRL_fit  = PSTH_right_fit - PSTH_left_fit;
diffRL_test = PSTH_right_test - PSTH_left_test;

% We now explore the similarity of coding direction across time using corr
% function in matlab.
figure;set(gcf,'Color','w')
title(['Coding direction correlation across time for Session #' num2str(sessionId)])
hold on
imagesc(tAxisToPlot, tAxisToPlot, corr(diffRL_fit,diffRL_test));
gridxy([-2.6 -1.3 0],[-2.6 -1.3 0],'Color','k','Linestyle','--') ;
xlim([-3.0  1.5]);
ylim([-3.0  1.5]);
xlabel('Time from movement (sec)')
ylabel('Time from movement (sec)')
set(gca,'box','off','tickdir','out','fontsize',16)


%% projection of data to delay-epoch coding direction 

% delay-epoch coding direction is defined as an average of coding direction
% over whole delay period -1.3 to 0 sec in time tag.

cdDelay = mean(diffRL_fit(:, tAxisToPlot > -1.3 & tAxisToPlot < 0), 2);
cdDelay = cdDelay/norm(cdDelay);
CDprojR    = cdDelay' * PSTH_right_test;
CDprojL    = cdDelay' * PSTH_left_test;

% We now project the neuronal activity in different trial type conditions
% to the coding direction.

figure;set(gcf,'Color','w')
title(['Coding direction projection for Session #' num2str(sessionId)])
hold on
plot(tAxisToPlot, CDprojR, '-b')
plot(tAxisToPlot, CDprojL, '-r')
gridxy([-2.6 -1.3 0],'Color','k','Linestyle','--') ;
xlim([-3.0  1.5]);
xlabel('Time from movement (sec)')
ylabel('Activity projected coding direction')
set(gca,'box','off','tickdir','out','fontsize',16)

% how much of selectivity does it explain?=
varSel   = sum(diffRL_test.^2,1);

% square of projection to CD
varCD   = (CDprojR-CDprojL).^2;

figure;set(gcf,'Color','w')
title(['Coding direction projection for Session #' num2str(sessionId)])
hold on
plot(tAxisToPlot, varCD./varSel*100, '-k')
gridxy([-2.6 -1.3 0],'Color','k','Linestyle','--') ;
xlim([-3.0  1.5]);ylim([0 100])
xlabel('Time from movement (sec)')
ylabel('Selectivity explained by CD (%)')
set(gca,'box','off','tickdir','out','fontsize',16)



