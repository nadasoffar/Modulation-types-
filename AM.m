clc;
close all;
clear all;
filename = 'eric.wav';
[y,fs]= audioread(filename);

%Plotting of original signal in TIME domain
t = linspace(0, length(y)/fs, length(y));
figure();
subplot(2,1,1);
plot(t,y);
title('Signal in time domain');
xlabel('Time');
ylabel('Amplitude');

%Plotting of original signal in FREQUENCY domain
yf = fftshift(fft(y));
f = linspace(-fs/2, fs/2, length(yf));
subplot(2,1,2);
plot(f, abs(yf)./fs);
title('Signal in frequency domain');
xlabel('Frequency');
ylabel('Amplitude');

%filter
filter_freq = 4000;
lp = ones(length(yf),1);  
lp(abs(f)>filter_freq) = 0;
yf_filtered = (lp.*yf);

%plot the filter
f2 = linspace(-fs/2,fs/2,length(lp));
figure();
plot(f2,abs(lp)./fs);
title('Filter');
xlabel('Frequency');
ylabel('Amplitude');

%plot the filtered signal
f2 = linspace(-fs/2,fs/2,length(yf_filtered));
figure();
subplot(2,1,2);
plot(f2,abs(yf_filtered)./fs);
title('Filtered Signal in frequency domain');
xlabel('Frequency');
ylabel('Amplitude');

yt_filtered = ifft(ifftshift(yf_filtered));
subplot(2,1,1);
plot(t,real(yt_filtered));
title('Filtered Signal in time domain');
xlabel('Time');
ylabel('Amplitude');
%sound(real(yt_filtered), fs);

%DSBSC mod
fc=100000;
fs_new=5*fc;
yt_new = resample(yt_filtered,fs_new,fs);
t_new=linspace(0, length(yt_new)/fs_new, length(yt_new));
carrier=transpose(cos(2*pi*fc*t_new));
yt_dsbsc=yt_new.*carrier;
yf_dsbsc=real(fftshift(fft(yt_dsbsc)));
f_new=linspace(-fs_new/2,fs_new/2,length(yt_dsbsc)); 

%plot
figure();
subplot(2,1,1);
plot(t_new,real(yt_dsbsc));
title('DSBSC modulated Signal in time domain');
xlabel('Time');
ylabel('Amplitude');
subplot(2,1,2);
plot(f_new,abs(yf_dsbsc)./fs_new);
title('DSBSC modulated Signal in frequency domain');
xlabel('Frequency');
ylabel('Amplitude');

%DSBSC env demod
yt_sc_env = resample(abs(hilbert(real(yt_dsbsc))),fs,fs_new);
yf_sc_env = real(fftshift(fft(yt_sc_env)));

%plot
figure();
subplot(2,1,1);
t1=linspace(0,length(yt_sc_env)/fs,length(yt_sc_env));
plot(t1,real(yt_sc_env));
title('DSBSC envelope demodulated Signal in time domain');
xlabel('Time');
ylabel('Amplitude');
subplot(2,1,2);
f1 = linspace(-fs/2, fs/2, length(yf_sc_env));
plot(f1,abs(yf_sc_env)./fs);
title('DSBSC envelope demodulated Signal in frequency domain');
xlabel('Frequency');
ylabel('Amplitude');
%sound(real(yt_sc_env),fs); %sound of DSBSC env demod

%DSBSC coherent demodulation and SNR
t_coh=linspace(0,length(yt_dsbsc)/fs_new,length(yt_dsbsc));
carrier_coh_no_snr=transpose(cos(2*pi*fc*t_coh));
yt_sc_coh=yt_dsbsc.*carrier_coh_no_snr;
yf_sc_coh=fftshift(fft(yt_sc_coh));
lp_sc_coh= 2*ones(length(yt_dsbsc),1);  
lp_sc_coh(abs(f)>4000) = 0;
yf_sc_no_snr_filtered = lp_sc_coh.*yf_sc_coh;
temp_after_filter=ifft(ifftshift(yf_sc_no_snr_filtered));
yt_sc_coh_no_snr=resample(temp_after_filter,fs,fs_new);
yf_sc_coh_no_snr=real(fftshift(fft(yt_sc_coh_no_snr)));

%plot
figure();
subplot(2,1,1);
t_after_resample=linspace(0,length(yt_sc_coh_no_snr)/fs,length(yt_sc_coh_no_snr));
plot(t_after_resample,real(yt_sc_coh_no_snr));
title('DSBSC coherent demodulated Signal no snr in time domain');
xlabel('Time');
ylabel('Amplitude');
subplot(2,1,2);
f_after_resample=linspace(-fs/2,fs/2,length(yf_sc_coh_no_snr));
plot(f_after_resample,abs(yf_sc_coh_no_snr)./fs);
title('DSBSC coherent demodulated Signal no snr in frequency domain');
xlabel('Frequency');
ylabel('Amplitude');
%sound(real(yt_sc_coh_no_snr),fs); 

%SNR=0
yt_coh_0_snr=awgn(yt_dsbsc,0);
t_coh=linspace(0,length(yt_coh_0_snr)/fs_new,length(yt_coh_0_snr));
carrier_coh_0_snr=transpose(cos(2*pi*fc*t_coh));
temp_before_filter=yt_coh_0_snr.*carrier_coh_0_snr;
yf_coh_0_snr=fftshift(fft(temp_before_filter));
f_SNR0=linspace(-(fs_new)/2,(fs_new)/2,length(yf_coh_0_snr));
lp_sc_coh= 2*ones(length(yt_coh_0_snr),1);  
lp_sc_coh(abs(f_SNR0)>4000) = 0;
yf_sc_0_snr_filtered = lp_sc_coh.*yf_coh_0_snr;
temp_after_filter=ifft(ifftshift(yf_sc_0_snr_filtered));
yt_sc_coh_0_snr=resample(temp_after_filter,fs,fs_new);
yf_sc_coh_0_snr=real(fftshift(fft(yt_sc_coh_0_snr)));

%plot
figure();
subplot(2,1,1);
t_after_resample=linspace(0,length(yt_sc_coh_0_snr)/fs,length(yt_sc_coh_0_snr));
plot(t_after_resample,real(yt_sc_coh_0_snr));
title('DSBSC coherent demodulated Signal SNR=0 in time domain');
xlabel('Time');
ylabel('Amplitude');
subplot(2,1,2);
f_after_resample=linspace(-fs/2,fs/2,length(yf_sc_coh_0_snr));
plot(f_after_resample,abs(yf_sc_coh_0_snr)./fs);
title('DSBSC coherent demodulated Signal SNR=0 in frequency domain');
xlabel('Frequency');
ylabel('Amplitude');
%sound(real(yt_sc_coh_0_snr),fs);

%SNR=10
yt_coh_10_snr=awgn(yt_dsbsc,10);
t_coh=linspace(0,length(yt_coh_10_snr)/fs_new,length(yt_coh_10_snr));
carrier_coh_10_snr=transpose(cos(2*pi*fc*t_coh));
temp_before_filter=yt_coh_10_snr.*carrier_coh_10_snr;
yf_coh_10_snr=fftshift(fft(temp_before_filter));
f_SNR10=linspace(-(fs_new)/2,(fs_new)/2,length(yf_coh_10_snr));
lp_sc_coh= 2*ones(length(yt_coh_10_snr),1);  
lp_sc_coh(abs(f_SNR10)>4000) = 0;
yf_sc_10_snr_filtered = lp_sc_coh.*yf_coh_10_snr;
temp_after_filter=ifft(ifftshift(yf_sc_10_snr_filtered));
yt_sc_coh_10_snr=resample(temp_after_filter,fs,fs_new);
yf_sc_coh_10_snr=real(fftshift(fft(yt_sc_coh_10_snr)));

%plot
figure();
subplot(2,1,1);
t_after_resample=linspace(0,length(yt_sc_coh_10_snr)/fs,length(yt_sc_coh_10_snr));
plot(t_after_resample,real(yt_sc_coh_10_snr));
title('DSBSC coherent demodulated Signal SNR=10 in time domain');
xlabel('Time');
ylabel('Amplitude');
subplot(2,1,2);
f_after_resample=linspace(-fs/2,fs/2,length(yf_sc_coh_10_snr));
plot(f_after_resample,abs(yf_sc_coh_10_snr)./fs);
title('DSBSC coherent demodulated Signal SNR=10 in frequency domain');
xlabel('Frequency');
ylabel('Amplitude');
%sound(real(yt_sc_coh_10_snr),fs);

%SNR=30
yt_coh_30_snr=awgn(yt_dsbsc,30);
t_coh=linspace(0,length(yt_coh_30_snr)/fs_new,length(yt_coh_30_snr));
carrier_coh_30_snr=transpose(cos(2*pi*fc*t_coh));
temp_before_filter=yt_coh_30_snr.*carrier_coh_30_snr;
yf_coh_30_snr=fftshift(fft(temp_before_filter));
f_SNR30=linspace(-(fs_new)/2,(fs_new)/2,length(yf_coh_30_snr));
lp_sc_coh= 2*ones(length(yt_coh_30_snr),1);  
lp_sc_coh(abs(f_SNR30)>4000) = 0;
yf_sc_30_snr_filtered = lp_sc_coh.*yf_coh_30_snr;
temp_after_filter=ifft(ifftshift(yf_sc_30_snr_filtered));
yt_sc_coh_30_snr=resample(temp_after_filter,fs,fs_new);
yf_sc_coh_30_snr=real(fftshift(fft(yt_sc_coh_30_snr)));

%plot
figure();
subplot(2,1,1);
t_after_resample=linspace(0,length(yt_sc_coh_30_snr)/fs,length(yt_sc_coh_30_snr));
plot(t_after_resample,real(yt_sc_coh_30_snr));
title('DSBSC coherent demodulated Signal SNR=30 in time domain');
xlabel('Time');
ylabel('Amplitude');
subplot(2,1,2);
f_after_resample=linspace(-fs/2,fs/2,length(yf_sc_coh_30_snr));
plot(f_after_resample,abs(yf_sc_coh_30_snr)./fs);
title('DSBSC coherent demodulated Signal SNR=30 in frequency domain');
xlabel('Frequency');
ylabel('Amplitude');
%sound(real(yt_sc_coh_30_snr),fs);

%phase error
t_coh=linspace(0,length(yt_dsbsc)/fs_new,length(yt_dsbsc));
carrier_coh_phase_error=transpose(cos(2*pi*fc*t_coh+(20*(pi/180))));
temp_before_filter=yt_dsbsc.*carrier_coh_phase_error;
yf_sc_coh_phase_error=fftshift(fft(temp_before_filter));
lp_sc_coh= 2*ones(length(yt_dsbsc),1);  
lp_sc_coh(abs(f)>fs_new) = 0;
yf_sc_phase_error_filtered = lp_sc_coh.*yf_sc_coh_phase_error;
temp_after_filter=ifft(ifftshift(yf_sc_phase_error_filtered));
yt_sc_coh_phase_error=resample(temp_after_filter,fs,fs_new);
yf_sc_coh_phase_error=real(fftshift(fft(yt_sc_coh_phase_error)));

%plot
figure();
subplot(2,1,1);
t_after_resample=linspace(0,length(yt_sc_coh_phase_error)/fs,length(yt_sc_coh_phase_error));
plot(t_after_resample,real(yt_sc_coh_phase_error));
title('DSBSC coherent demodulated Signal phase error in time domain');
xlabel('Time');
ylabel('Amplitude');
subplot(2,1,2);
f_after_resample=linspace(-fs/2,fs/2,length(yf_sc_coh_phase_error));
plot(f_after_resample,abs(yf_sc_coh_phase_error)./fs);
title('DSBSC coherent demodulated Signal phase error in frequency domain');
xlabel('Frequency');
ylabel('Amplitude');
%sound(real(yt_sc_coh_phase_error),fs); 

%frequency error
t_coh=linspace(0,length(yt_dsbsc)/fs_new,length(yt_dsbsc));
carrier_coh_freq_error=transpose(cos(2*pi*100100*t_coh));
temp_before_filter=yt_dsbsc.*carrier_coh_freq_error;
yf_sc_coh_freq_error=fftshift(fft(temp_before_filter));
lp_sc_coh= 2*ones(length(yt_dsbsc),1);  
lp_sc_coh(abs(f)>fs_new) = 0;
yf_sc_freq_error_filtered = lp_sc_coh.*yf_sc_coh_freq_error;
temp_after_filter=ifft(ifftshift(yf_sc_freq_error_filtered));
yt_sc_coh_freq_error=resample(temp_after_filter,fs,fs_new);
yf_sc_coh_freq_error=real(fftshift(fft(yt_sc_coh_freq_error)));

%plot
figure();
subplot(2,1,1);
t_after_resample=linspace(0,length(yt_sc_coh_freq_error)/fs,length(yt_sc_coh_freq_error));
plot(t_after_resample,real(yt_sc_coh_freq_error));
title('DSBSC coherent demodulated Signal frequency error in time domain');
xlabel('Time');
ylabel('Amplitude');
subplot(2,1,2);
f_after_resample=linspace(-fs/2,fs/2,length(yf_sc_coh_freq_error));
plot(f_after_resample,abs(yf_sc_coh_freq_error)./fs);
title('DSBSC coherent demodulated Signal frquency error in frequency domain');
xlabel('Frequency');
ylabel('Amplitude');
%sound(real(yt_sc_coh_freq_error),fs); 

%DSBTC mod
Am=max(abs(real(yt_new))); %max gain of message
yt_dsbtc=(2*Am+real(yt_new)).*carrier;
yf_dsbtc=real(fftshift(fft(yt_dsbtc)));

%plot
figure();
subplot(2,1,1);
plot(t_new,real(yt_dsbtc));
title('DSBTC modulated Signal in time domain');
xlabel('Time');
ylabel('Amplitude');
subplot(2,1,2);
plot(f_new,abs(yf_dsbtc));
xlim([-10500 10500])
title('DSBTC modulated Signal in frequency domain');
xlabel('Frequency');
ylabel('Amplitude');

%DSBTC envelope demodulation
yt_tc_env = resample(abs(hilbert((yt_dsbtc))),fs,fs_new)-(2*Am);
yf_tc_env = real(fftshift(fft(yt_tc_env)));

%plot
figure();
subplot(2,1,1);
t1=linspace(0,length(yt_tc_env)/fs,length(yt_tc_env));
plot(t1,yt_tc_env);
title('DSBTC envelope demod Signal in time domain');
xlabel('Time');
ylabel('Amplitude');
subplot(2,1,2);
f1 = linspace(-fs/2, fs/2, length(yf_tc_env));
plot(f1,abs(yf_tc_env)./fs);
title('DSBTC envelope demod Signal in frequency domain');
xlabel('Frequency');
ylabel('Amplitude');
sound(real(yt_tc_env),fs);
%%%%%%%%%%%%%%%%%%%SSB_SC%%%%%%%%%%%%%%%%%%%%%%%%%
%(5)Obtain the SSB by filtering out the USB (we need to get LSB only) of the DSB-SC modulated 
%signal using an ideal filter then Plot the spectrum again.
f_ssb=linspace(-(fs_new/2),(fs_new/2),length(yf_dsbsc));
lpfc = 2.*ones(length(yf_dsbsc),1);   
lpfc(abs(f_ssb)>fc) = 0;
yf_lssb = lpfc.*yf_dsbsc;
yt_ssb = ifft(ifftshift(yf_lssb));
t_ssb=linspace(0, length(yt_ssb)/fs_new, length(yt_ssb));
figure();
subplot(2,1,1);
plot(f_ssb,abs(yf_lssb)./fs_new);
title('Single SideBand Modulation (LSB) frequency domain');
xlabel('Frequency');
ylabel('Amplitude');

subplot(2,1,2);
plot(t_ssb,real(yt_ssb));
title('Single SideBand Modulation (LSB) time domain');
xlabel('Time');
ylabel('Amplitude');

%(6)-Use coherent detection with no noise interference to get the received signal (to demodulate the 
%SSB-SC) and play the file back also sketch the received waveform and spectrum
%SSB coherent demod time domain
carrier1=cos(2*pi*fc*t_ssb).';
yt_ssbdemod = carrier1.*yt_ssb;
%SSB coherent demod freq domain
yf_ssbdemod = fftshift(fft(yt_ssbdemod));
%filter ssb
yf_ssbdemod(f_ssb>4000) = 0;
yf_ssbdemod(f_ssb<-4000) = 0;
yf_ssbdemod2 = yf_ssbdemod.*2;
%filtered ssb demodulated

f_ssbf=linspace(-(fs_new/2),(fs_new/2),length(yf_ssbdemod2));
yt_ssbdemodf = ifft(ifftshift(yf_ssbdemod2));
yt_ssbdemodf = resample(yt_ssbdemodf,fs,fs_new);
t_ssbf=linspace(0, length(yt_ssbdemodf)/fs_new, length(yt_ssbdemodf));
figure();
subplot(2,1,1);
plot(f_ssbf,abs(yf_ssbdemod2)./fs);
title('SSB cohDemodulated (Frequency)');
xlabel('Frequency');
ylabel('Amplitude');


subplot(2,1,2);
plot(t_ssbf,real(yt_ssbdemodf));
title('SSB cohDemodulated (Time)');
xlabel('Time');
ylabel('Amplitude');

%sound(real(yt_ssbdemodf),fs);


%(7)Repeat steps 5 and 6, only this time. Use a practical 4th order Butterworth filter.
%ssb butterworth
flow=(fc-4e3)/(fs_new/2);
fhigh=fc/(fs_new/2);
fcut=[flow,fhigh];
[a,b] = butter(4,fcut);
yt_ssb2 = filter(a,b,yt_dsbsc).*2;
yf_ssb2 = fftshift(fft(yt_ssb2));
f_ssb2=linspace(-(fs_new/2),(fs_new/2),length(yf_ssb2));
t_ssb2 = linspace(0, length(yt_ssb2)/fs_new, length(yt_ssb2));

figure();
subplot(2,1,1);
plot(f_ssb2,abs(yf_ssb2)./fs_new);
title('Single SideBand butter (frequency)');
xlabel('Frequency');
ylabel('Amplitude');
subplot(2,1,2);
plot(t_ssb2,real(yt_ssb2));
title('Single SideBand butter (time)');
xlabel('Time');
ylabel('Amplitude');

%ssb butter cohdemodulation
carrier5=cos(2*pi*fc*t_ssb2).';
yt_ssbbutterdemod = carrier5.*yt_ssb2;
%SSB cohbutter demod freq domain
yf_ssbbutterdemod = fftshift(fft(yt_ssbbutterdemod));
%filter ssb
filter_freq = 4000;
lp4 = 2.*ones(length(yf_ssbbutterdemod),1);
lp4(abs(f_ssb2)>filter_freq) = 0;
%filtered ssb demodulated
yf_ssbbutterdemodf = lp4.*yf_ssbbutterdemod;
f_ssbbutterf=linspace(-(fs_new/2),(fs_new/2),length(yf_ssbbutterdemodf));
yt_ssbbutterdemodf = ifft(ifftshift(yf_ssbbutterdemodf));
yt_ssbbutterdemodf = resample(yt_ssbbutterdemodf,fs,fs_new);
t_ssbbutterf = linspace(0, length(yt_ssbbutterdemodf)/fs_new, length(yt_ssbbutterdemodf));
figure();
subplot(2,1,1);
plot(f_ssbbutterf,abs(yf_ssbbutterdemodf)./fs);
title('SSB butter cohdemod(frequency)');
xlabel('Frequency');
ylabel('ssb');

subplot(2,1,2);
plot(t_ssbbutterf,real(yt_ssbbutterdemodf));
title('SSB butter cohdemod(time)');
xlabel('Time');
ylabel('ssb');

%sound(real(yt_ssbbutterdemodf),fs);

%(8) For the ideal filter case, get the received signal again but when noise is added to SSB-SC with 
%SNR = 0, 10, and 30 also play the received sound and sketch the received waveform/spectrum in 
%each case. 

%SNR 0
yt_SNR0 = awgn(yt_ssb, 0);
t_SNR0=linspace(0, length(yt_SNR0)/fs_new, length(yt_SNR0));
carrier2=cos(2*pi*fc*t_SNR0).';
yt_SNR0demod = carrier2.*yt_SNR0;
yf_SNR0demod = fftshift(fft(yt_SNR0demod));
f_SNR0=linspace(-(fs_new/2),(fs_new/2),length(yf_SNR0demod));

%LPF SNR0
filter_freq = 4000;
lp3 = 2.*ones(length(yf_SNR0demod),1);
lp3(abs(f_SNR0)>filter_freq) = 0;

%filtered ssb demodulated
yf_SNR0demodf = lp3.*yf_SNR0demod;
f_SNR0f=linspace(-(fs_new/2),(fs_new/2),length(yf_SNR0demodf));
yt_SNR0demodf= ifft(ifftshift(yf_SNR0demodf));
yt_SNR0demodf = resample(yt_SNR0demodf,fs,fs_new);
t_SNR0f=linspace(0, length(yt_SNR0demodf)/fs_new, length(yt_SNR0demodf));
figure();
subplot(2,1,1);
plot(f_SNR0f,abs(yf_SNR0demodf)./fs);
title('Single SideBand SNR0 demod (frequency)');
xlabel('Frequency');
ylabel('ssb');

subplot(2,1,2);
plot(t_SNR0f,real(yt_SNR0demodf));
title('Single SideBand SNR0 demod (time)');
xlabel('Time');
ylabel('ssb');


%sound(real(yt_SNR0demodf),fs);

%SNR 10
yt_SNR10 = awgn(yt_ssb, 10);
t_SNR10=linspace(0, length(yt_SNR10)/fs_new, length(yt_SNR10));
carrier3=cos(2*pi*fc*t_SNR10).';
yt_SNR10demod = carrier3.*yt_SNR10;

yf_SNR10demod = fftshift(fft(yt_SNR10demod));
f_SNR10=linspace(-(fs_new/2),(fs_new/2),length(yf_SNR10demod));
%lpf 
filter_freq = 4000;
lp10 = 2.*ones(length(yf_SNR10demod),1);
lp10(abs(f_SNR10)>filter_freq) = 0;
%filtered ssb demodulated
yf_SNR10demodf = lp10.*yf_SNR10demod;
yt_SNR10demodf= ifft(ifftshift(yf_SNR10demodf));
f_SNR10f=linspace(-(fs_new/2),(fs_new/2),length(yf_SNR10demodf));
yt_SNR10demodf = resample(yt_SNR10demodf,fs,fs_new);
t_SNR10f=linspace(0, length(yt_SNR10demodf)/fs_new, length(yt_SNR10demodf));
figure();
subplot(2,1,1);
plot(f_SNR10f,abs(yf_SNR10demodf)./fs);
title('Single SideBand SNR10 demod (frequency)');
xlabel('Frequency');
ylabel('ssb');

subplot(2,1,2);
plot(t_SNR10f,real(yt_SNR10demodf));
title('Single SideBand SNR10 demod (time)');
xlabel('Time');
ylabel('ssb');

%sound(real(yt_SNR10demodf),fs);

%SNR 30
yt_SNR30 = awgn(yt_ssb, 30);
t_SNR30=linspace(0, length(yt_SNR30)/fs_new, length(yt_SNR30));
carrier4=cos(2*pi*fc*t_SNR30).';
yt_SNR30demod = carrier3.*yt_SNR30;
yf_SNR30demod = fftshift(fft(yt_SNR30demod));
f_SNR30=linspace(-(fs_new/2),(fs_new/2),length(yf_SNR30demod));
%LPF SNR30
filter_freq = 4000;
lp30 = 2.*ones(length(yf_SNR30demod),1);
lp30(abs(f_SNR30)>filter_freq) = 0;
%filtered ssb demodulated
yf_SNR30demodf = lp30.*yf_SNR30demod;
yt_SNR30demodf = ifft(ifftshift(yf_SNR30demodf));
f_SNR30f=linspace(-(fs_new/2),(fs_new/2),length(yf_SNR30demodf));
yt_SNR30demodf = resample(yt_SNR30demodf,fs,fs_new);
t_SNR30f=linspace(0, length(yt_SNR30demodf)/fs_new, length(yt_SNR30demodf));
figure();
subplot(2,1,1);
plot(f_SNR30f,abs(yf_SNR30demodf)./fs_new);
title('Single SideBand SNR30 demod (frequency)');
xlabel('Frequency');
ylabel('ssb');

subplot(2,1,2);
plot(t_SNR30f,real(yt_SNR30demodf));
title('Single SideBand SNR30 demod (time)');
xlabel('Time');
ylabel('ssb');

%sound(real(yt_SNR30demodf),fs);


%%%%%%%%%%%%%%%%%%%SSB-TC%%%%%%%%%%%%%%%%%%%%%%%%%

f_ssbtc=linspace(-(fs_new/2),(fs_new/2),length(yf_dsbtc));
yf_dsbtc(f_ssbtc>100001) = 0;
yf_dsbtc(f_ssbtc<-100001) = 0;
yt_dsbtcf=ifft(ifftshift(yf_dsbtc));

t_ssbtc3 = linspace(0, length(yt_dsbtcf)/fs_new, length(yt_dsbtcf));
figure();
subplot(2,1,1);
plot(t_ssbtc3,real(yt_dsbtc));
title('SSBTC Signal modulated');
xlabel('Time');
ylabel('Amplitude');

subplot(2,1,2);
plot(f_new,abs(yf_dsbtc)./fs_new);
xlim([-10500 10500])
title('SSBTC Signal modulated');
xlabel('Frequency');
ylabel('Amplitude');



yt_ssbtc_env = resample(abs(hilbert(yt_dsbtc)),fs,fs_new)-(2*Am);
yf_ssbtc_env = real(fftshift(fft(yt_ssbtc_env)));

%yf_tc_env = real(fftshift(fft(yt_tc_env)));
t_ssbtc = linspace(0, length(yt_ssbtc_env)/fs, length(yt_ssbtc_env));
figure();
plot(t_ssbtc,yt_ssbtc_env);
title('SSBTC Signal demodulated');
xlabel('Time');
ylabel('Amplitude');  

sound(real(yt_ssbtc_env),fs);