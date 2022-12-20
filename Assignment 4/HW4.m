load sampleEEGdata;
%% Task 1

timewin           = 400; % in ms, for stFFT
times2save     = -300:50:1000; % in ms
channel2plot   = 'p7';
frequency2plot = 10;  % in Hz
timepoint2plot = 100; % ms

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

for timepointi=1:length(times2save)
   
    time_interval = times2saveidx(timepointi)-floor(timewinidx/2):times2saveidx(timepointi)+floor(timewinidx/2)-mod(timewinidx+1,2) 
    tempdat = squeeze(EEG.data(chan2useidx,time_interval,:));
    taperdat = bsxfun(@times,tempdat,hann_win'); 
    fdat = fft(taperdat,[],1)/timewinidx; 
    ifdat = ifft(fdat(1:floor(timewinidx/2)+1,:), [], 1);
    tf(:,timepointi) = mean(abs(fdat(1:floor(timewinidx/2)+1,:)).^2,2); 
    
end

% plot
figure
subplot(131) % plot time trend for a specific frequency
[junk,freq2plotidx]=min(abs(frex-frequency2plot));
plot(times2save,mean(log10(tf(freq2plotidx-2:freq2plotidx+2,:)),1))
title([ 'Sensor ' channel2plot ', ' num2str(frequency2plot) ' Hz' ])
axis square
set(gca,'xlim',[times2save(1) times2save(end)])

subplot(132) % plot FFT for a specific time 
[junk,time2plotidx]=min(abs(times2save-timepoint2plot));
plot(frex,log10(tf(:,time2plotidx)))
title([ 'Sensor ' channel2plot ', ' num2str(timepoint2plot) ' ms' ])
axis square
set(gca,'xlim',[frex(1) 40])

subplot(133) % plot power 
f           = linspace(0,EEG.srate/2,floor(length(hann_win)/2)+1); % frequencies of FFT
power  = abs(tf(1:floor(length(hann_win)/2)+1)).^2;

plot(f(2:end),power(2:end),'.-');
title('power spectrum from that time window')
%set(gca,'xlim',[1 128],'ylim',[-1000 25000],'xtick',0:10:EEG.srate/2)



%% Task 2
% other wavelet parameters
time             = -1:1/EEG.srate:1;
half_of_wavelet_size = (length(time)-1)/2;

% FFT parameters (use next-power-of-2)
n_wavelet     = length(time);
n_data        = EEG.pnts;
n_convolution = n_wavelet+n_data-1;
n_conv_pow2   = pow2(nextpow2(n_convolution));
wavelet_cycles= 4; 
% FFT of data (note: this doesn't change on frequency iteration)
fft_data = fft(squeeze(EEG.data(46,:,1)),n_conv_pow2);


A = (pi*2*sqrt(pi))^-.5;
Gaussian_fi = exp(-time.^2./(2*( wavelet_cycles /(2*pi*2))^2));
wavelet = A * exp(2*1i*pi*2.*time) .* Gaussian_fi/2;
fft_wavelet = fft(wavelet,n_conv_pow2);
convolution_result_fft = ifft(fft_wavelet.*fft_data,n_conv_pow2);
convolution_result_fft = convolution_result_fft(1:n_convolution); 
convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
tf_data = abs(convolution_result_fft).^2;
figure
subplot(131)
plot(EEG.times, tf_data)
xlabel('Time (ms)'), ylabel('Amplitude')
title('Wavelet')

%%
timewin           = 400; % in ms, for stFFT
times2save     = -300:50:1000; % in ms
frequency2plot = 2;  % in Hz

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

for timepointi=1:length(times2save)
   
    time_interval = times2saveidx(timepointi)-floor(timewinidx/2):times2saveidx(timepointi)+floor(timewinidx/2)-mod(timewinidx+1,2) 
    tempdat = squeeze(EEG.data(chan2useidx,time_interval,:));
    taperdat = bsxfun(@times,tempdat,hann_win'); 
    fdat = fft(taperdat,[],1)/timewinidx; 
    ifdat = ifft(fdat(1:floor(timewinidx/2)+1,:), [], 1);
    tf(:,timepointi) = mean(abs(fdat(1:floor(timewinidx/2)+1,:)).^2,2); 
    
end

% plot
subplot(132) % plot time trend for a specific frequency
[junk,freq2plotidx]=min(abs(frex-frequency2plot)); % plot power for a specific time 
f           = linspace(0,EEG.srate/2,floor(length(hann_win)/2)+1); % frequencies of FFT
power  = abs(tf(1:floor(length(hann_win)/2)+1)).^2;

plot(f(2:end),power(2:end),'.-');
title('STFT')
%%
nyquist = EEG.srate/2;
lower_filter_bound = 4; % Hz
upper_filter_bound = 10; % Hz
transition_width   = 0.2;
filter_order       = round(3*(EEG.srate/lower_filter_bound));

% create the filter shape (this is explained more in the text around figure 14.4)
ffrequencies  = [ 0 (1-transition_width)*lower_filter_bound lower_filter_bound upper_filter_bound (1+transition_width)*upper_filter_bound nyquist ]/nyquist;
idealresponse = [ 0 0 1 1 0 0 ];
filterweights = firls(filter_order,ffrequencies,idealresponse);

% apply the filter kernal to the data to obtain the band-pass filtered signal
filtered_data = filtfilt(filterweights,1,double(EEG.data(46,:,1)));
filtered_data = reshape(filtered_data, [640, 1]);


hilbert  = hilbert(filtered_data); % time should be in the first dimension.
subplot(133)
plot(EEG.times,abs(hilbert).^2,'b');
title('Hilbert')
xlabel('Time (ms)'), ylabel('Amplitude')
set(gca,'xlim',[-1000 1500])