%getting the signal in the time domain and it's power spectrum
load('\Users\user\Desktop\ecg_421n');
Mdata = M(20000:22500,2); %taking the last two columns containing the data
Mtime = M(20000:22500,1); %taking the first column containing the time
Mperiod = Mtime(2)-Mtime(1);
Mfreq = 1/Mperiod;
figure;
plot(Mtime,Mdata)
title('Signal in time domain')
xlabel('Time')
ylabel('Magnitude')
%% 
N = 10*Mfreq; %N-point DFT
X = fft(Mdata);
n = length(X);
f = (0:n-1)*Mfreq/n;
power = abs(X).^2;
figure;
plot(f,power);
title('Frequency plot of unfiltered signal')
xlabel('Frequency')
ylabel('Power')

%% plot in dB
figure;
plot(f,20*log10(power))%spectral power magnitude\\
title('Data Power Spectrum')
xlabel('Frequency')
ylabel('Power (dB)')
%% IIR Filter Butterworth
freqc = [(2*atan(2*pi*Mperiod))/pi (2*atan(45*pi*Mperiod))/pi];
[b1,a1] = butter(6,freqc,'bandpass');
[h,w] = freqz(b1,a1);

figure;
subplot(2,1,1)
plot(w,abs(h));
title('IIR Magnitude Filter Response')
xlabel('Frequency')
ylabel('Magnitude')
subplot(2,1,2)
plot(w,unwrap(angle(h)*180/pi))
title('IIR Phase Filter Response')
xlabel('Frequency')
ylabel('Phase')
%% plot in dB
figure;
plot(w,20*log10(abs(h)));
title('IIR Filter Response (dB)')
xlabel('Frequency')
ylabel('Magnitude (dB)')
%% impulse response
[h,w] = impz(b1,a1);
figure;
plot(w,abs(h));
title('IIR Filter impulse response')
xlabel('Frequency (Hz)')
ylabel('Impulse Response')
%% filtering the data in M
Y = filter(b1,a1,Mdata);
figure
plot(Mtime,Mdata)
xlabel('Time')
ylabel('Magnitude')
hold on
plot(Mtime,Y,'r')
legend('Original Signal','Filtered Signal IIR Butter')
%% IIR Filter Chebychev
n = 6; 
Rp = 0.1; %peak to peak bandpass ripple
W = [2*atan(2*pi*Mperiod) 2*atan(45*pi*Mperiod)];
ftype = 'bandpass';
[z,p,k] = cheby1(n,Rp,W/pi,ftype);
sos = zp2sos(z,p,k);
hfvt = fvtool(sos,'FrequencyScale','log');

%% filtering the data in M
Y = filter(b1,a1,Mdata);
R = filter(p/4000,k,Mdata);
figure
plot(Mtime,Mdata,'y')
xlabel('Time')
ylabel('Magnitude')
hold on
plot(Mtime,Y,'r')
hold on
plot(Mtime,R,'b');
legend('Original Signal','Filtered Signal IIR Butter','Filtered Signal IIR Cheb')
%%
N1 = 10*Mfreq; %N-point DFT
X1 = fft(Y);
X2 = fft(R);
n = length(X1);
f1 = (0:n-1)*Mfreq/n;
power1 = abs(X1).^2;
power2 = abs(X2).^2;
figure;
plot(f1,power1);
hold on
plot(f1,power2,'y');
title('Frequency plot of IIR filtered signal')
xlabel('Frequency')
ylabel('Power')
legend('Filtered by Butter','Filtered by Cheby')

%% plot in dB
figure;
plot(f1,20*log10(power1))%spectral power magnitude\\
hold on
plot(f1,20*log10(power2),'y');
title('Data Power Spectrum of filtered signal')
xlabel('Frequency')
ylabel('Power (dB)')
legend('Filtered by Butter','Filtered by Cheby')
%% FIR filter
f3 = [2 45]/Mfreq;
c = fir1(100,f3); %100 isthe length of the signal
out = filter(c,1,Mdata);
freqz(c,1);
figure;
out = filter(FIR2,1,Mdata);

%% 
figure;
subplot(4,1,1)
plot(Mtime,Mdata)
title('Original Signal')
xlabel('Time')
ylabel('Magnitude')
subplot(4,1,2)
plot(Mtime,R)
title('Filtered Signal IIR Cheby')
xlabel('Time')
ylabel('Magnitude')
subplot(4,1,3)
plot(Mtime,Y)
title('Filtered Signal IIR Butter')
xlabel('Time')
ylabel('Magnitude')
subplot(4,1,4)
plot(Mtime,out)
title('Filtered Signal FIR')
xlabel('Time')
ylabel('Magnitude')