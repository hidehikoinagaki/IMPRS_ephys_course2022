%% PCA tutorial

clear;close all

numCell = 1000;

%% create test data
% two baiss
mode(1,:) = sin(pi*(1:201)/60)*1;
mode(2,:) = [zeros(1,100),[0:0.01:1]]*1;
tAxis = 1:size(mode,2);

% plot basis
figure;set(gcf,'Color','w')
plot(tAxis,mode)
xlabel('Time');title('Two basis');ylabel('Activity (a.u.)');xlim([0 max(tAxis)])
set(gca,'box','off','tickdir','out','fontsize',16)

% weight of basis
w = randn(numCell,2);

% creata data set using the weight
X = w*mode;

% show data
% all data
figure;set(gcf,'Color','w')
plot(tAxis,X)
xlabel('Time');title('Activity of neurons');ylabel('Activity (a.u.)');xlim([0 max(tAxis)])
set(gca,'box','off','tickdir','out','fontsize',16)

% % plot example individual cells
% for i=1:20
% figure;set(gcf,'Color','w')
% plot(tAxis,X(i,:))
% xlabel('Time');title(['w = ',num2str(w(i,1)),' , ',num2str(w(i,2))]);ylabel('Activity (a.u.)');xlim([0 max(tAxis)])
% set(gca,'box','off','tickdir','out','fontsize',16)
% end

%% PCA
% remove the mean variable-wise (row-wise)
X=X-repmat(mean(X,2),1,size(X,2));

% calculate eigenvectors (loadings) W, and eigenvalues of the covariance matrix
[W, EvalueMatrix] = eig(cov(X'));
Evalues = diag(EvalueMatrix);

% order by largest eigenvalue
Evalues = Evalues(end:-1:1);
W = W(:,end:-1:1); W=W';  

% Eigen spectrum
figure;set(gcf,'Color','w')
plot(1:numCell,Evalues/sum(Evalues));
xlim([0 10])
xlabel('PCs');title('Eigen spectrum');ylabel('Variance explained')
set(gca,'box','off','tickdir','out','fontsize',16)

% Principle components
pc = W(1:2,:) * X;
figure;set(gcf,'Color','w')
plot(tAxis,pc)
xlabel('Time');title('Two PCs');ylabel('Activity (a.u.)');xlim([0 max(tAxis)])
set(gca,'box','off','tickdir','out','fontsize',16)

% comapre the original weight w, and PCA weight W
w(:,1)=w(:,1)/norm(w(:,1));
w(:,2)=w(:,2)/norm(w(:,2));
W(1,:)=W(1,:)/norm(W(1,:));
W(2,:)=W(2,:)/norm(W(2,:));
figure;set(gcf,'Color','w')
subplot(1,2,1);plot(W(1,:),w(:,1),'ko');title('Mode1')
xlabel('Original weight');title('PC 1');ylabel('PC weight')
set(gca,'box','off','tickdir','out','fontsize',16)
subplot(1,2,2);plot(W(2,:),w(:,2),'ko');title('Mode2')
xlabel('Original weight');title('PC 2');ylabel('PC weight')
set(gca,'box','off','tickdir','out','fontsize',16)

%% we can do the same using matlab PCA function
[coeff,score,latent] = pca(X);
pc2 = score'*X;

figure;set(gcf,'Color','w')
plot(tAxis,pc2(1:2,:)')
xlabel('Time');title('Two PCs');ylabel('Activity (a.u.)');xlim([0 max(tAxis)])
set(gca,'box','off','tickdir','out','fontsize',16)
