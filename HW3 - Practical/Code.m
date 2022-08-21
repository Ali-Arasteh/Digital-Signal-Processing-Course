%%
% Q1)
%%
% 7.5)
delta1 = 0.01;
delta2 = 0.05;
delta3 = 0.01;
delta = min([delta1,delta2,delta3]);
A = -20*log10(delta);
if (A < 21)
    beta = 0;
elseif (A >= 21 && A <= 50)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
elseif (A > 50)
    beta = 0.1102*(A-8.7);
end
deltaw1 = 0.1*pi;
deltaw2 = 0.05*pi;
deltaw = min(deltaw1,deltaw2);
M = ceil((A-8)/(2.285*deltaw));
n = 0:M;
hd = sin(0.3*pi*(n-M/2))./(pi*(n-M/2)) - sin(0.625*pi*(n-M/2))./(pi*(n-M/2));
if M/2 == floor(M/2)
    hd(M/2+1) = 0.3 - 0.625;
end
KaiserFilter = kaiser(M+1,beta)'; 
h = hd.*KaiserFilter;
H = fftshift(fft(h));
figure();
stem(n,h)
xlabel('n')
ylabel('h[n]')
title('impulse response in time domain')
figure();
f = (0:length(H)-1)*2*pi/(length(H)-1) - pi;
plot(f,abs(H))
xlabel('\omega')
ylabel('H(e^{j\omega})')
title('impulse response in frequency domain')
%%
% 7.6)
delta1 = 0.1;
delta2 = 0.06;
delta3 = 0.05;
delta = min([delta1,delta2,delta3]);
A = -20*log10(delta);
if (A < 21)
    beta = 0;
elseif (A >= 21 && A <= 50)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
elseif (A > 50)
    beta = 0.1102*(A-8.7);
end
deltaw1 = 0.1*pi;
deltaw2 = 0.05*pi;
deltaw = min(deltaw1,deltaw2);
M = ceil((A-8)/(2.285*deltaw));
if M/2 ~= floor(M/2)
    M = M+1;
end
KaiserFilter = kaiser(M+1,beta)'; 
n = -M/2:M/2;
hd = sin(0.25*pi*(n))./(pi*(n)) + 2*sin(pi*(n))./(pi*(n)) - 2*sin(0.5*pi*(n))./(pi*(n));
if M/2 == floor(M/2)
    hd(M/2+1) = 0.25 - 2*0.5 + 2*1;
end
h = hd.*KaiserFilter;
H = fftshift(fft(h));
figure();
stem(n,h)
xlabel('n')
ylabel('h[n]')
title('impulse response in time domain')
figure();
f = (0:length(H)-1)*2*pi/(length(H)-1) - pi;
plot(f,abs(H))
xlabel('\omega')
ylabel('H(e^{j\omega})')
title('impulse response in frequency domain')
%%
% Q2)
%%
% Q3)
%%
[signal,fs] = audioread('MultiFreq_Sig.wav');
signal = (signal(:,1)+signal(:,2))'/2;
frequency = -fs/2:fs/length(signal):fs/2-fs/length(signal);
signalF = fftshift(fft(signal));
figure();
plot(frequency,abs(signalF),'b','LineWidth',1)
xlabel('frequency')
ylabel('signal(e^{j2\pif})')
title('signal in fequency domain')
%%
% a)
M = 40;
[signal,fs] = audioread('MultiFreq_Sig.wav');
signal = (signal(:,1)+signal(:,2))'/2;
h = ones(1,M)/M;
H = fftshift(fft(h,length(signal)));
frequency = -fs/2:fs/length(signal):fs/2-fs/length(signal);
figure();
subplot(2,1,1);
plot(frequency,abs(H),'b','LineWidth',1)
xlabel('frequency')
ylabel('|H(e^{j2\pif})|')
title('magnitude of frequency response')
subplot(2,1,2);
plot(frequency,angle(H),'b','LineWidth',1)
xlabel('frequency')
ylabel('<H(e^{j2\pif})')
title('phase of frequency response')
%%
% b)
%%
% N = 3
N = 3;
fc = 2000;
[~,fs] = audioread('MultiFreq_Sig.wav');
[bmax,amax] = butter(N,fc/fs);
H = freqz(bmax,amax)';
frequency = 0:fs/(2*length(H)):fs/2-fs/(2*length(H));
figure();
subplot(2,1,1);
plot(frequency,abs(H),'b','LineWidth',1)
xlabel('frequency')
ylabel('|H(e^{j2\pif})|')
title('magnitude of frequency response')
subplot(2,1,2);
plot(frequency,angle(H),'b','LineWidth',1)
xlabel('frequency')
ylabel('<H(e^{j2\pif})')
title('phase of frequency response')
%%
% N = 10
N = 10;
fc = 2000;
[~,fs] = audioread('MultiFreq_Sig.wav');
[bmax,amax] = butter(N,fc/fs);
H = freqz(bmax,amax)';
frequency = 0:fs/(2*length(H)):fs/2-fs/(2*length(H));
figure();
subplot(2,1,1);
plot(frequency,abs(H),'b','LineWidth',1)
xlabel('frequency')
ylabel('|H(e^{j2\pif})|')
title('magnitude of frequency response')
subplot(2,1,2);
plot(frequency,angle(H),'b','LineWidth',1)
xlabel('frequency')
ylabel('<H(e^{j2\pif})')
title('phase of frequency response')
%%
% c)
N = 5;
fc = 2000;
[~,fs] = audioread('MultiFreq_Sig.wav');
[bmax,amax] = cheby1(N,0.01,fc/fs);
H = freqz(bmax,amax)';
frequency = 0:fs/(2*length(H)):fs/2-fs/(2*length(H));
figure();
subplot(2,1,1);
plot(frequency,abs(H),'b','LineWidth',1)
xlabel('frequency')
ylabel('|H(e^{j2\pif})|')
title('magnitude of frequency response')
subplot(2,1,2);
plot(frequency,angle(H),'b','LineWidth',1)
xlabel('frequency')
ylabel('<H(e^{j2\pif})')
title('phase of frequency response')
%%
% d)
N = 8;
fc = 2000;
[~,fs] = audioread('MultiFreq_Sig.wav');
[bmax,amax] = ellip(N,0.01,50,fc/fs);
H = freqz(bmax,amax)';
frequency = 0:fs/(2*length(H)):fs/2-fs/(2*length(H));
figure();
subplot(2,1,1);
plot(frequency,abs(H),'b','LineWidth',1)
xlabel('frequency')
ylabel('|H(e^{j2\pif})|')
title('magnitude of frequency response')
subplot(2,1,2);
plot(frequency,angle(H),'b','LineWidth',1)
xlabel('frequency')
ylabel('<H(e^{j2\pif})')
title('phase of frequency response')
%%
fc = 17000;
M = 40;
[signal,fs] = audioread('MultiFreq_Sig.wav');
signal = (signal(:,1)+signal(:,2))'/2;
n = 0:length(signal)-1;
signal = signal.*cos(2*pi*fc.*n/fs);
h = ones(1,M)/M;
fsignal = filter(h,1,signal);
fsignal = 2*fsignal.*cos(2*pi*fc.*n/fs);
H = fftshift(fft(fsignal));
frequency = -fs/2:fs/length(fsignal):fs/2-fs/length(fsignal);
figure();
plot(frequency,abs(H),'b','LineWidth',1)
xlabel('frequency')
ylabel('|filtered signal(e^{j2\pif})|')
title('filtered signal in frequency domain')
%%
fc = 17000;
N = 3;
fp = 2000;
[signal,fs] = audioread('MultiFreq_Sig.wav');
signal = (signal(:,1)+signal(:,2))'/2;
n = 0:length(signal)-1;
signal = signal.*cos(2*pi*fc.*n/fs);
[bmax,amax] = butter(N,fp/fs);
fsignal = filter(bmax,amax,signal);
fsignal = 2*fsignal.*cos(2*pi*fc.*n/fs);
H = fftshift(fft(fsignal));
frequency = -fs/2:fs/length(fsignal):fs/2-fs/length(fsignal);
figure();
plot(frequency,abs(H),'b','LineWidth',1)
xlabel('frequency')
ylabel('|filtered signal(e^{j2\pif})|')
title('filtered signal in frequency domain')
%%
fc = 17000;
N = 10;
fp = 2000;
[signal,fs] = audioread('MultiFreq_Sig.wav');
signal = (signal(:,1)+signal(:,2))'/2;
n = 0:length(signal)-1;
signal = signal.*cos(2*pi*fc.*n/fs);
[bmax,amax] = butter(N,fp/fs);
fsignal = filter(bmax,amax,signal);
fsignal = 2*fsignal.*cos(2*pi*fc.*n/fs);
H = fftshift(fft(signal));
frequency = -fs/2:fs/length(fsignal):fs/2-fs/length(fsignal);
figure();
plot(frequency,abs(H),'b','LineWidth',1)
xlabel('frequency')
ylabel('|filtered signal(e^{j2\pif})|')
title('filtered signal in frequency domain')
%%
fc = 17000;
N = 5;
fp = 2000;
[signal,fs] = audioread('MultiFreq_Sig.wav');
signal = (signal(:,1)+signal(:,2))'/2;
n = 0:length(signal)-1;
signal = signal.*cos(2*pi*fc.*n/fs);
[bmax,amax] = cheby1(N,0.01,fp/fs);
fsignal = filter(bmax,amax,signal);
fsignal = 2*fsignal.*cos(2*pi*fc.*n/fs);
H = fftshift(fft(fsignal));
frequency = -fs/2:fs/length(fsignal):fs/2-fs/length(fsignal);
figure();
plot(frequency,abs(H),'b','LineWidth',1)
xlabel('frequency')
ylabel('|filtered signal(e^{j2\pif})|')
title('filtered signal in frequency domain')
%%
fc = 17000;
N = 8;
fp = 2000;
[signal,fs] = audioread('MultiFreq_Sig.wav');
signal = (signal(:,1)+signal(:,2))'/2;
n = 0:length(signal)-1;
signal = signal.*cos(2*pi*fc.*n/fs);
[bmax,amax] = ellip(N,0.01,50,fp/fs);
fsignal = filter(bmax,amax,signal);
fsignal = 2*fsignal.*cos(2*pi*fc.*n/fs);
H = fftshift(fft(fsignal));
frequency = -fs/2:fs/length(fsignal):fs/2-fs/length(fsignal);
figure();
plot(frequency,abs(H),'b','LineWidth',1)
xlabel('frequency')
ylabel('|filtered signal(e^{j2\pif})|')
title('filtered signal in frequency domain')
%%
% Q4)
%%
% a)
bmax = [1 -4 8 -64 201 -450 650];
amax = 1;
zeros = roots(bmax)';
zplane(bmax,amax)
title('zero-poles of maximum phase system')
%%
% b)
bmax = [1 -4 8 -64 201 -450 650];
amax = 1;
[Zmax,Pmax,Kmax] = tf2zpk(bmax,amax);
Zap = conj(1./Zmax);
Pap = Zmax;
Kap = Kmax*prod(conj(Zmax));
[bap,aap] = zp2tf(Zap,Pap,Kap);
Hap = freqz(bap,aap);
figure()
zplane(bap,aap);
title('zero-poles of allpass system')
f = 0:pi/length(Hap):(length(Hap)-1)*pi/length(Hap);
figure();
plot(f,abs(Hap))
xlim([0 (length(Hap)-1)*pi/length(Hap)])
ylim([0 2])
xlabel('frequency')
ylabel('|Hap(e^{j2\pif})|')
title('magnitude of frequency response of allpass system')
%%
% c)
bHBmin = firhalfband(25,0.45,'minphase');
aHBmin = 1;
HHBmin = freqz(bHBmin,aHBmin);
figure();
zplane(bHBmin,aHBmin)
title('zero-poles of halfband minimum phase filter')
f = 0:pi/length(HHBmin):(length(HHBmin)-1)*pi/length(HHBmin);
figure();
subplot(2,1,1);
plot(f,abs(HHBmin))
xlabel('frequency')
ylabel('|HHBmin(e^{j2\pif})|')
title('magnitude of frequency response halfband minimum phase filter')
subplot(2,1,2);
plot(f,angle(HHBmin))
xlabel('frequency')
ylabel('<HHBmin(e^{j2\pif})')
title('phase of frequency response halfband minimum phase filter')
figure();
impz(bHBmin,aHBmin);
xlabel('time')
ylabel('HHBmin(t)')
title('impulse response of halfband minimum phase filter')
[ZHBmin,PHBmin,KHBmin] = tf2zpk(bHBmin,aHBmin);
ZHBmax = conj(1./ZHBmin);
PHBmax = PHBmin;
KHBmax = KHBmin;
[bHBmax,aHBmax] = zp2tf(ZHBmax,PHBmax,KHBmax);
HHBmax = freqz(bHBmax,aHBmax);
figure();
zplane(bHBmax,aHBmax)
title('zero-poles of halfband maximum phase filter')
f = 0:pi/length(HHBmax):(length(HHBmax)-1)*pi/length(HHBmax);
figure();
subplot(2,1,1);
plot(f,abs(HHBmax))
xlabel('frequency')
ylabel('|HHBmax(e^{j2\pif})|')
title('magnitude of frequency response halfband maximum phase filter')
subplot(2,1,2);
plot(f,angle(HHBmax))
xlabel('frequency')
ylabel('<HHBmax(e^{j2\pif})')
title('phase of frequency response halfband maximum phase filter')
figure();
impz(bHBmax,aHBmax);
xlabel('time')
ylabel('HHBmax(t)')
title('impulse response of halfband maximum phase filter')