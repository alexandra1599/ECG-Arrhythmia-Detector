%getting the signal in the time domain and it's power spectrum
load('C:\Users\user\Downloads\ecg_421n');
Mdata = M(:,2); %taking the last two columns containing the data
Mtime = M(:,1); %taking the first column containing the time
Mperiod = Mtime(2)-Mtime(1);
Mfreq = 1/Mperiod;
figure;
plot(Mtime,Mdata)
title('Signal in time domain')
xlabel('Time')
ylabel('Magnitude')

%% plotting the frequency content of the signal
Mdata1 = M(1:2500,2); %first 10 seconds of the data
Mdata3 = M(27500:30000,2); %duplets in heart beat
Mdata2 = M(17500:20000,2); %VF period of the data
N = 10*Mfreq; %N-point DFT
X = fft(Mdata1);
Y = fft(Mdata2);
Z = fft(Mdata3);
n = length(X);
n1 = length(Y);
n2 = length(Z);
f = (0:n-1)*Mfreq/n;
f1 = (0:n1-1)*Mfreq/n1;
f2 = (0:n2-1)*Mfreq/n2;
power1 = abs(X).^2;
power2 = abs(Y).^2;
power3 = abs(Z).^2;
figure;
subplot(3,1,1)
plot(f,power1);
title('Frequency plot of normal ECG')
xlabel('Frequency')
ylabel('Power')
subplot(3,1,3)
plot(f1,power2);
title('Frequency plot of VF')
xlabel('Frequency')
ylabel('Power')
subplot(3,1,2)
plot(f2,power3);
title('Frequency plot of Duplets')
xlabel('Frequency')
ylabel('Power')

%% VF detector
%continuously create 10 sec segments of the data
%plot the frequency content
%if desired frequency found go to next 10 sec segment
%if not then signal VF
i=1;
while i < 30000
    for j = i:1:i+2499
    time = Mtime(i:j,1);
    segment = M(i:j,2);
    end
    N5 = 10*Mfreq; %N-point DFT
    X5 = fft(segment);
    n = length(X5);
    f = (0:n-1)*Mfreq/n;
    power11 = abs(X5).^2;
    b = interp1(f,power11,1.7);
    if b<20000
        Y = -10:0.1:10;
        X = Mtime(i,1) * ones(size(Y));
        figure;
        plot(Mtime,Mdata);
        hold on
        plot(X, Y);
        xlabel('Time')
        ylabel('Magnitude')
        title('ECG signal')
        amp=10; 
        fs=250;  % sampling frequency
        duration=10;
        freq=100;
        values=0:1/fs:duration;
        a=amp*sin(2*pi* freq*values);
        sound(a)
        break
    end
    i=i+2500;
end    