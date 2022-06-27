clc;
close all;
clear all;



filename = 'eric.wav';
[message,fs]= audioread(filename);

%Plotting the original signal in TIME domain
time = linspace(0, length(message)/fs, length(message));
figure()
subplot(2,1,1);
plot(time,message);
title('Message Signal in time domain');
xlabel('Time');
ylabel('Amplitude');

%Plotting of original signal in FREQUENCY domain
message_in_freq = fftshift(fft(message));
frequency_of_Message = linspace(-fs/2, fs/2, length(message_in_freq));
subplot(2,1,2);
plot(frequency_of_Message, abs(message_in_freq));
title('Message Signal in frequency domain');
xlabel('Frequency');
ylabel('Amplitude');

%filter
filter_freq = 4000;
lp = ones(length(message_in_freq),1);
lp(abs(frequency_of_Message)>filter_freq) = 0;
yf_filtered = lp.*message_in_freq;

%plot the filter
f2 = linspace(-fs/2,fs/2,length(lp));
figure()
plot(f2,abs(lp));
title('Filter');
xlabel('Frequency');
ylabel('Amplitude;');

%plot the filtered signal in FREQUENCY
f2 = linspace(-fs/2,fs/2,length(yf_filtered));
figure()
subplot(2,1,1);
plot(f2,abs(yf_filtered));
title('Filtered Signal in frequency domain');
xlabel('Frequency');
ylabel('Amplitude;');

%plot the filtered signal in TIME domain
yt_filtered= ifft(ifftshift(yf_filtered));
time2 = linspace(0, length(yt_filtered)/fs, length(yt_filtered));
%figure()
subplot(2,1,2);
plot(time2,yt_filtered); %imaginary nnumber warning without abs
title('Filtered signal in time domain');
xlabel('Time');
ylabel('Amplitude');

sound(real(yt_filtered),fs);

%NBFM
Fc_NBFM = 100000;
%sampling frequency
Fs_NBFM = 5*Fc_NBFM;
%Amplitude of the message
Am = max(abs(message));
%equation:
%cos(2*pi*fc*t + 2*pi*kf*integ_message);
%expression : 2*pi*kf*integ_message < 1 rad

Resampled_NBFM = resample(yt_filtered,Fs_NBFM,fs);
integ_message = cumsum(Resampled_NBFM)/Fs_NBFM;
Kf = 0.2/(2*pi*max(integ_message));
%carrier amplitude
Ac_NBFM = max(abs(Resampled_NBFM));

time_NBFM = linspace(0, (length(Resampled_NBFM))/Fs_NBFM, length(Resampled_NBFM));
time_NBFM = transpose(time_NBFM);
freq_NBFM = linspace((-Fs_NBFM)/2, Fs_NBFM/2, length(Resampled_NBFM));
NBFM_Signal = Ac_NBFM*cos(2*pi*Fc_NBFM*time_NBFM + 2*pi*Kf*integ_message);
NBFM_signal_spectrum = fftshift (fft(NBFM_Signal)./Fs_NBFM);

figure()
subplot(2,1,1);
plot(time_NBFM, real(NBFM_Signal));
title('NBFM modulated signal in time domain');
xlabel('Time');
ylabel('Amplitude');
subplot(2,1,2);
plot(freq_NBFM, abs(NBFM_signal_spectrum));
title('NBFM modulated signal in frequency domain');
xlabel('Frequency');
ylabel('Amplitude');

%Demodulation
Output_Signal_NBFM_diff = diff(NBFM_Signal)*Fs_NBFM;
Output_Signal_demod_time = abs(hilbert(Output_Signal_NBFM_diff));
Output_Signal_demod_time = Output_Signal_demod_time - abs(mean(Output_Signal_demod_time));
Output_Signal_demod_time_resampled = resample(Output_Signal_demod_time,fs,Fs_NBFM);
Output_Signal_demod_freq = fftshift(fft(Output_Signal_demod_time_resampled)./Fs_NBFM);

figure()
subplot(2,1,1);
time_demod =linspace(0, (length(Output_Signal_demod_time_resampled))/fs, length(Output_Signal_demod_time_resampled));
plot(time_demod,0.01*Output_Signal_demod_time_resampled);
ylim([-0.2 0.2]);
title('Demodulated signal in time domain');
xlabel('Time');
ylabel('Amplitude');

sound(Output_Signal_demod_time_resampled,fs);
subplot(2,1,2);
freq_demod = linspace(-fs/2, fs/2, length(Output_Signal_demod_freq));
plot(freq_demod,abs(Output_Signal_demod_freq));
title('Demodulated signal in frequency domain');
xlabel('Frequency');
ylabel('Amplitude');

