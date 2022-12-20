clear all, close all, clc
load sampleEEGdata

%% Task 1
times2save     = 200:1:600; % in ms


% convert from ms to index
times2saveidx = zeros(size(times2save));
for i=1:length(times2save)
    [junk,times2saveidx(i)]=min(abs(EEG.times-times2save(i)));
end

times = EEG.times(min(times2saveidx) + 1: max(times2saveidx) - 1);
data = EEG.data(:, min(times2saveidx) + 1: max(times2saveidx) - 1, :);

erp = squeeze(mean(data,3));

erp = bsxfun(@minus,erp,mean(erp,2));
covar = (erp*erp')./(length(times)-1);

figure(100)
subplot(121)
imagesc(covar)
axis square
set(gca,'xticklabel',{EEG.chanlocs(get(gca,'xtick')).labels},'yticklabel',{EEG.chanlocs(get(gca,'ytick')).labels},'clim',[-1 5])
title('Covariance of ERP')

subplot(122)
covar = zeros(EEG.nbchan);

for tr=1:EEG.trials
    eeg = bsxfun(@minus,squeeze(data(:,:,tr)),squeeze(mean(data(:,:,tr),2)));
    covar = covar + (eeg*eeg')./(length(times)-1);
end
covar = covar./tr;
imagesc(covar)
axis square
set(gca,'xticklabel',{EEG.chanlocs(get(gca,'xtick')).labels},'yticklabel',{EEG.chanlocs(get(gca,'ytick')).labels},'clim',[20 150])
title('Average covariance of single-trial EEG')

% Here we can see a difference: for the single-trial EEG, covariance matrix
% reflects mostly negative relationship, although for the ERP data, covariance
% matrix reflects positive relationship for some electrodes

%% Task 2
times2save1     = -500:1:0; % in ms


% convert from ms to index
times2saveidx1 = zeros(size(times2save1));
for i=1:length(times2save1)
    [junk,times2saveidx1(i)]=min(abs(EEG.times-times2save1(i)));
end

times1 = EEG.times(min(times2saveidx1): max(times2saveidx1));
data1 = EEG.data(:, min(times2saveidx1): max(times2saveidx1), :);

times2save2     = 100:1:600; % in ms


% convert from ms to index
times2saveidx2 = zeros(size(times2save2));
for i=1:length(times2save2)
    [junk,times2saveidx2(i)]=min(abs(EEG.times-times2save2(i)));
end

times2 = EEG.times(min(times2saveidx2): max(times2saveidx2));
data2 = EEG.data(:, min(times2saveidx2): max(times2saveidx2), :);

erp1 = squeeze(mean(data1,3));
erp1 = bsxfun(@minus,erp1,mean(erp1,2));
covar1 = (erp1*erp1')./(length(times1)-1);

% principle components analysis via eigenvalue decomposition on ERP covariance
[pc1,eigvals1] = eig(covar1);

% components are listed in increasing order, and converted here to descending order for convenience
pc1      = pc1(:,end:-1:1);
eigvals1 = diag(eigvals1);
eigvals1 = 100*eigvals1(end:-1:1)./sum(eigvals1); % convert to percent change

erp2 = squeeze(mean(data2,3));
erp2 = bsxfun(@minus,erp2,mean(erp2,2));
covar2 = (erp2*erp2')./(length(times2)-1);

% principle components analysis via eigenvalue decomposition on ERP covariance
[pc2,eigvals2] = eig(covar2);

% components are listed in increasing order, and converted here to descending order for convenience
pc2      = pc2(:,end:-1:1);
eigvals2 = diag(eigvals2);
eigvals2 = 100*eigvals2(end:-1:1)./sum(eigvals2); % convert to percent change

eeglabpath = '/Volumes/MY_DRIVE/MA_CS/RSII/eeglab2019_1'
addpath(genpath(eeglabpath))

for i=1:4 % only first 6 are shown in the real figure
    figure(1)
    subplot(2,2,i)
    topoplot(double(pc1(:,i)),EEG.chanlocs,'electrodes','off','plotrad',.53);
    title([ 'PC #' num2str(i) ', eigval=' num2str(eigvals1(i)) ])

    figure(2)
    subplot(4,1,i)
    plot(times1,pc1(:,i)'*erp1)
    hold on
    plot(get(gca,'xlim'),[0 0],'k')
    set(gca,'xlim',[-200 0])
    title([ 'PC #' num2str(i) ', eigval=' num2str(eigvals1(i)) ])
end

for i=1:4 % only first 6 are shown in the real figure
    figure(3)
    subplot(2,2,i)
    topoplot(double(pc2(:,i)),EEG.chanlocs,'electrodes','off','plotrad',.53);
    title([ 'PC #' num2str(i) ', eigval=' num2str(eigvals2(i)) ])

    figure(4)
    subplot(4,1,i)
    plot(times2,pc1(:,i)'*erp2)
    hold on
    plot(get(gca,'xlim'),[0 0],'k')
    set(gca,'xlim',[100 600])
    title([ 'PC #' num2str(i) ', eigval=' num2str(eigvals2(i)) ])
end;

% There is a difference in the topographical maps before and after stimuli:
% the activity patters look similar for PC 1,2, and 4 components, but
% differ for the PC 3 component. The magnitude of the EEG components differs
% before and after stimuli. This can be explained that the stimuli changes
% activity pattern, 