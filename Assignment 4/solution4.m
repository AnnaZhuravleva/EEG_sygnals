%% solution 4

clear all, close all, clc
load sampleEEGdata

%% Figure 1

timewin           = 300; % in ms, for stFFT
times2save     = -300:25:1000; % in ms


% convert from ms to index
times2saveidx = zeros(size(times2save));
for i=1:length(times2save)
    [junk,times2saveidx(i)]=min(abs(EEG.times-times2save(i)));
end
timewinidx = round(timewin/(1000/EEG.srate)); % window length in samples
chan2useidx = strcmpi(channel2plot,{EEG.chanlocs.labels});


% create hann taper
hann_win = .5*(1-cos(2*pi*(0:timewinidx-1)/(timewinidx-1)));

% define frequencies
frex = linspace(0,EEG.srate/2,floor(timewinidx/2)+1);

% initialize power output matrix
tf = zeros(length(frex),length(times2save));

% loop over time points and perform FFT
for timepointi=1:length(times2save)
    
    % extract time series data for this center time point
    time_interval = times2saveidx(timepointi)-floor(timewinidx/2):times2saveidx(timepointi)+floor(timewinidx/2)-mod(timewinidx+1,2) % make sure it is an even number
    tempdat = squeeze(EEG.data(chan2useidx,time_interval,:)); % note: the 'mod' function here corrects for even or odd number of points
    
    % taper data (using bsxfun instead of repmat... note sizes of tempdat
    % and hann_win)
    taperdat = bsxfun(@times,tempdat,hann_win'); %[time x trials]
    % with repmat would be:
    % taperdat =  tempdat.*repmat(hann_win',1,size(tempdat,2));

    fdat = fft(taperdat,[],1)/timewinidx; % FFT(X,[],DIM), 3rd input is to make sure fft is over time
    tf(:,timepointi) = mean(abs(fdat(1:floor(timewinidx/2)+1,:)).^2,2); % average over trials
end

% plot
figure(1)
subplot(2,2,1) % plot time trend for a specific frequency
[junk,freq2plotidx]=min(abs(frex-frequency2plot));
plot(times2save,mean(log10(tf(freq2plotidx-2:freq2plotidx+2,:)),1))
title([ 'Sensor ' channel2plot ', ' num2str(frequency2plot) ' Hz' ])
axis square
set(gca,'xlim',[times2save(1) times2save(end)])

subplot(2,2,2) % plot FFT for a specific time 
[junk,time2plotidx]=min(abs(times2save-timepoint2plot));
plot(frex,log10(tf(:,time2plotidx)))
title([ 'Sensor ' channel2plot ', ' num2str(timepoint2plot) ' ms' ])
axis square
set(gca,'xlim',[frex(1) 40])

subplot(2,2,[3 4])
contourf(times2save,frex,log10(tf),40,'linecolor','none')
set(gca,'clim',[-2 1])
title([ 'Sensor ' channel2plot ', power plot (no baseline correction)' ])

disp([ 'Overlap of ' num2str(100*(1-mean(diff(times2save))/timewin)) '%' ])


%% Figure 2

%%%%%%%%%% CWT - RSII_211012_a

num_frex = 1
frex = 6

% definitions, selections...
chan2use = 'fcz'; 

% define wavelet parameters
time = -1:1/EEG.srate:1;
 % s    = logspace(log10(3),log10(10),num_frex)./(2*pi*frex); % number of cycles for the wavelets
s    =  ones(size(frex))*3./(2*pi*frex); % this line is for figure 13.14
 
n_wavelet            = length(time);
n_data               = EEG.pnts*EEG.trials;
n_convolution        = n_wavelet+n_data-1;             % length of fft
n_conv_pow2          = pow2(nextpow2(n_convolution));
half_of_wavelet_size = (n_wavelet-1)/2;

% eegdata = reshape(EEG.data(strcmpi(chan2use,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials);
eegfft = fft(reshape(EEG.data(strcmpi(chan2use,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials),n_conv_pow2);
 

% initialize
eegpower_cwt = zeros(num_frex,EEG.pnts); % frequencies X time X trials

baseidx = dsearchn(EEG.times',[-500 -200]'); % interval for baseline normalization

% loop through frequencies and compute synchronization
for fi=1 
    
    wavelet = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2 );
    
    % convolution
    eegconv = ifft(wavelet.*eegfft); % some zeros at the end. why?
    eegconv = eegconv(1:n_convolution); % take away zeros
    eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size); % length of original signal
    
    % Average power over trials (this code performs baseline transform,
    % which you will learn about in chapter 18)
    temppower = mean(abs(reshape(eegconv,EEG.pnts,EEG.trials)).^2,2);
        eegpower_cwt(fi,:) =  temppower ;

%     eegpower_cwt(fi,:) = 10*log10(temppower./mean(temppower(baseidx(1):baseidx(2))));
end
 


%%%%%%%%%% FIR + Hilbert

for fi=1 
    
nyquist            = EEG.srate/2;
lower_filter_bound = frex(fi)-2; % Hz
upper_filter_bound = frex(fi)+2; % Hz
transition_width   = 0.2; % defines the steepness of the slopes of the filter
filter_order       = round(3*(EEG.srate/lower_filter_bound));

% create the filter shape (this is explained more in the text around figure 14.4)
ffrequencies  = [ 0 (1-transition_width)*lower_filter_bound lower_filter_bound upper_filter_bound (1+transition_width)*upper_filter_bound nyquist ]/nyquist;
idealresponse = [ 0 0 1 1 0 0 ];
filterweights = double(firls(filter_order,ffrequencies,idealresponse));
input2filter = double(squeeze(EEG.data(strcmpi(chan2use,{EEG.chanlocs.labels}),:,:)));
% filtfilt(filterweights,1,input2filter)

temppower = mean(abs(hilbert(filtfilt(filterweights,1,input2filter))).^2,2);
 eegpower_hilbert(fi,:) = temppower;
 
end

%%%%%%%%%% STFT
 
timewin           = 300; % in ms, for stFFT
times2save     = -800:25:1300; % in ms
 frequency2plot = 6;  % in Hz
% timepoint2plot = 100; % ms

% convert from ms to index
times2saveidx = zeros(size(times2save));
for i=1:length(times2save)
    [junk,times2saveidx(i)]=min(abs(EEG.times-times2save(i)));
end
timewinidx = round(timewin/(1000/EEG.srate)); % window length in samples
chan2useidx = strcmpi(chan2use,{EEG.chanlocs.labels});


% create hann taper
hann_win = .5*(1-cos(2*pi*(0:timewinidx-1)/(timewinidx-1)));

% define frequencies
 f = linspace(0,EEG.srate/2,floor(timewinidx/2)+1);

% initialize power output matrix
tf = zeros(length(f),length(times2save));

% loop over time points and perform FFT
for timepointi=1:length(times2save)
    
    % extract time series data for this center time point
    time_interval = times2saveidx(timepointi)-floor(timewinidx/2):times2saveidx(timepointi)+floor(timewinidx/2)-mod(timewinidx+1,2) % make sure it is an even number
    tempdat = squeeze(EEG.data(chan2useidx,time_interval,:)); % note: the 'mod' function here corrects for even or odd number of points
    
    % taper data (using bsxfun instead of repmat... note sizes of tempdat
    % and hann_win)
    taperdat = bsxfun(@times,tempdat,hann_win'); %[time x trials]
    % with repmat would be:
    % taperdat =  tempdat.*repmat(hann_win',1,size(tempdat,2));

    fdat = fft(taperdat,[],1)/timewinidx; % FFT(X,[],DIM), 3rd input is to make sure fft is over time
    tf(:,timepointi) = mean(abs(fdat(1:floor(timewinidx/2)+1,:)).^2,2); % average over trials
end

[junk,freq2plotidx]=min(abs(f-frequency2plot));

eegpower_stft = tf(freq2plotidx,:);


%%%%%%%%%%% Tapers

nw_product        = 3;  
tapers                 = dpss(timewinidx,nw_product); % note that in practice, you'll want to set the temporal resolution to be a function of frequency

 
% initialize output matrix
multitaper_tf = zeros(floor(timewinidx/2)+1,length(times2save));

% loop through time bins
for ti=1:length(times2saveidx)
    
    % initialize power vector (over tapers)
    taperpow = zeros(floor(timewinidx/2)+1,1);
    
    % loop through tapers
    for tapi = 1:size(tapers,2)-1
        
        % window and taper data, and get power spectrum
        sampleswindow = times2saveidx(ti)+[-floor(timewinidx/2)+1:ceil(timewinidx/2)];
        
        data      = bsxfun(@times,squeeze(EEG.data(chan2useidx,sampleswindow,:)),tapers(:,tapi));
        pow       = fft(data,timewinidx)/timewinidx;
        pow       = pow(1:floor(timewinidx/2)+1,:);
        taperpow  = taperpow + mean(abs(pow).^2,2); % power accumulated by each taper convolution
    end
    
    % finally, get power from closest frequency
    multitaper_tf(:,ti) = taperpow/tapi;
end

[junk,freq2plotidx]=min(abs(f-frequency2plot));

eegpower_mt = multitaper_tf(freq2plotidx,:);

figure(2),

subplot(411), plot(EEG.times, eegpower_cwt)
subplot(412), plot(EEG.times, eegpower_hilbert)
subplot(413), plot(times2save, eegpower_stft)
subplot(414), plot(times2save, eegpower_mt)


%% Figure 3
close all
clc
clear eegpower*
%%%%%%%%%% CWT - RSII_211012_a

 min_freq =  5;
max_freq = 80;
num_frex = 20;

frex = logspace(log10(min_freq),log10(max_freq),num_frex); % frequencies for the sines


% definitions, selections...
chan2use = 'fcz'; 

% define wavelet parameters
time = -1:1/EEG.srate:1;
 % s    = logspace(log10(3),log10(10),num_frex)./(2*pi*frex); % number of cycles for the wavelets
s    =  ones(size(frex))*3./(2*pi*frex); % this line is for figure 13.14
 
n_wavelet            = length(time);
n_data               = EEG.pnts*EEG.trials;
n_convolution        = n_wavelet+n_data-1;             % length of fft
n_conv_pow2          = pow2(nextpow2(n_convolution));
half_of_wavelet_size = (n_wavelet-1)/2;

% eegdata = reshape(EEG.data(strcmpi(chan2use,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials);
eegfft = fft(reshape(EEG.data(strcmpi(chan2use,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials),n_conv_pow2);
 

% initialize
eegpower_cwt = zeros(num_frex,EEG.pnts); % frequencies X time X trials

baseidx = dsearchn(EEG.times',[-500 -200]'); % interval for baseline normalization

% loop through frequencies and compute synchronization
for fi=1:length(frex)
    
    wavelet = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2 );
    
    % convolution
    eegconv = ifft(wavelet.*eegfft); % some zeros at the end. why?
    eegconv = eegconv(1:n_convolution); % take away zeros
    eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size); % length of original signal
    
    % Average power over trials (this code performs baseline transform,
    % which you will learn about in chapter 18)
    temppower = mean(abs(reshape(eegconv,EEG.pnts,EEG.trials)).^2,2);
%         eegpower_cwt(fi,:) =  temppower ;

 eegpower_cwt(fi,:) = (temppower./mean(temppower(baseidx(1):baseidx(2))));
end
 


%%%%%%%%%% FIR + Hilbert

 

frex = logspace(log10(min_freq),log10(max_freq),num_frex); % frequencies for the sines


for fi=1 :length(frex)
    
nyquist            = EEG.srate/2;
lower_filter_bound = frex(fi)-frex(fi)/3; % Hz
upper_filter_bound = frex(fi)+frex(fi)/3; % Hz
transition_width   = .1; % defines the steepness of the slopes of the filter
filter_order       = 100;round(2*(EEG.srate/lower_filter_bound));

% create the filter shape (this is explained more in the text around figure 14.4)
ffrequencies  = [ 0 (1-transition_width)*lower_filter_bound lower_filter_bound upper_filter_bound (1+transition_width)*upper_filter_bound nyquist ]/nyquist;
idealresponse = [ 0 0 1 1 0 0 ];
filterweights = double(firls(filter_order,ffrequencies,idealresponse));

temppower = mean(abs(hilbert(filtfilt(filterweights,1,double(squeeze(EEG.data(strcmpi(chan2use,{EEG.chanlocs.labels}),:,:)))))).^2,2);
eegpower_hilbert(fi,:) = (temppower./mean(temppower(baseidx(1):baseidx(2))));

% % to check the quality of the filter
% filterweightsW = (firls(filter_order,ffrequencies,idealresponse));
% hz_filtkern   = linspace(0,nyquist,filter_order/2+1); % list of frequencies in Hz corresponding to filter kernel
%         
%         fft_filtkern   = abs(fft(filterweightsW));
%         fft_filtkern   = fft_filtkern./max(fft_filtkern); % normalized to 1.0 for visual comparison ease
%         freqsidx = dsearchn(hz_filtkern',ffrequencies'*nyquist);
%         
%         sse(fi) = sum( (idealresponse-fft_filtkern(freqsidx)).^2 );
%         figure
%         plot(ffrequencies*nyquist,idealresponse)
%         hold on
%         plot(hz_filtkern,fft_filtkern(1:filter_order/2+1))
 
end

%%%%%%%%%% STFT
 
timewin           = 50; % in ms, for stFFT
times2save     = -800:5:1300; % in ms
 frequency2plot = 6;  % in Hz
% timepoint2plot = 100; % ms

baseidx = dsearchn(times2save',[-500 -200]'); 

% convert from ms to index
times2saveidx = zeros(size(times2save));
for i=1:length(times2save)
    [junk,times2saveidx(i)]=min(abs(EEG.times-times2save(i)));
end
timewinidx = round(timewin/(1000/EEG.srate)); % window length in samples
chan2useidx = strcmpi(chan2use,{EEG.chanlocs.labels});


% create hann taper
hann_win = .5*(1-cos(2*pi*(0:timewinidx-1)/(timewinidx-1)));

% define frequencies
 f = linspace(0,EEG.srate/2,floor(timewinidx/2)+1);

% initialize power output matrix
tf = zeros(length(f),length(times2save));

% loop over time points and perform FFT
for timepointi=1:length(times2save)
    
    % extract time series data for this center time point
    time_interval = times2saveidx(timepointi)-floor(timewinidx/2):times2saveidx(timepointi)+floor(timewinidx/2)-mod(timewinidx+1,2) ;% make sure it is an even number
    tempdat = squeeze(EEG.data(chan2useidx,time_interval,:)); % note: the 'mod' function here corrects for even or odd number of points
    
    % taper data (using bsxfun instead of repmat... note sizes of tempdat
    % and hann_win)
    taperdat = bsxfun(@times,tempdat,hann_win'); %[time x trials]
    % with repmat would be:
    % taperdat =  tempdat.*repmat(hann_win',1,size(tempdat,2));

    fdat = fft(taperdat,[],1)/timewinidx; % FFT(X,[],DIM), 3rd input is to make sure fft is over time
    eegpower_stft(:,timepointi) = mean(abs(fdat(1:floor(timewinidx/2)+1,:)).^2,2); % average over trials
end
 
baseline = repmat(mean(eegpower_stft(:,baseidx(1):baseidx(2)),2),1,size(eegpower_stft,2));

eegpower_stft = eegpower_stft./baseline;

%%%%%%%%%%% Tapers

nw_product        = 3;  
tapers                 = dpss(timewinidx,nw_product); % note that in practice, you'll want to set the temporal resolution to be a function of frequency

 
% initialize output matrix
multitaper_tf = zeros(floor(timewinidx/2)+1,length(times2save));

% loop through time bins
for ti=1:length(times2saveidx)
    
    % initialize power vector (over tapers)
    taperpow = zeros(floor(timewinidx/2)+1,1);
    
    % loop through tapers
    for tapi = 1:size(tapers,2)-1
        
        % window and taper data, and get power spectrum
        sampleswindow = times2saveidx(ti)+[-floor(timewinidx/2)+1:ceil(timewinidx/2)];
        
        data      = bsxfun(@times,squeeze(EEG.data(chan2useidx,sampleswindow,:)),tapers(:,tapi));
        pow       = fft(data,timewinidx)/timewinidx;
        pow       = pow(1:floor(timewinidx/2)+1,:);
        taperpow  = taperpow + mean(abs(pow).^2,2); % power accumulated by each taper convolution
    end
    
    % finally, get power from closest frequency
    eegpower_mt(:,ti) = taperpow/tapi;
end

baseline = repmat(mean(eegpower_mt(:,baseidx(1):baseidx(2)),2),1,size(eegpower_mt,2));

eegpower_mt = eegpower_mt./baseline;
 
figure(3),

subplot(411), contourf(EEG.times, frex, 10*log10(eegpower_cwt),40,'linecolor','none'),colorbar
subplot(412), contourf(EEG.times, frex, 10*log10(eegpower_hilbert),40,'linecolor','none'),colorbar
subplot(413), contourf(times2save, f, 10*log10(eegpower_stft),40,'linecolor','none'),colorbar
subplot(414), contourf(times2save, f, 10*log10(eegpower_mt),40,'linecolor','none'),colorbar
