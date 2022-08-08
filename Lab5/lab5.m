%%Student #:500913092
%           ABCDEFGHI
% H = 9, I = 2
%% A) Discrete Fourier Transform and Zero Padding
%%% Part 1 - Discrete Fourier Transform
n = 0:9;
N = length(n);
fr = linspace(-0.5, 0.49, 10);
x1 = exp(1i*2*pi*(1)*n)+exp(1i*2*pi*(33/100)*n);
x2 = cos(2*pi*(1)*n)+0.5*cos(2*pi*(0.3*n));

X1 = fftshift(fft(x1));
X2 = fftshift(fft(x2));

figure();
stem(fr, abs(X1));
title("Magnitude of X1 with Respect to fr");
xlabel('|X1[fr]|')
ylabel('fr')

figure();
stem(fr, abs(X2));
title("Magnitude of X2 with Respect to fr");
xlabel('|X2[fr]|')
ylabel('fr')

%%
% Part 2 - 490-zero padded DFT with 500 samples
n = 0:9;
N = 500;
fr = linspace(-0.5, 0.49, N);
x1 = exp(1i*2*pi*(1)*n)+exp(1i*2*pi*(33/100)*n);
x3 = [zeros(1, 245), x1, zeros(1, 245)];
x2 = cos(2*pi*(1)*n)+0.5*cos(2*pi*(0.3*n));
x4 = [zeros(1, 245), x2, zeros(1, 245)];

X3 = fftshift(fft(x3));
X4 = fftshift(fft(x4));

figure();
stem(fr, abs(X3));
title("Magnitude of X1 with Respect to fr");
xlabel('|X1[fr]|')
ylabel('fr')

figure();
stem(fr, abs(X4));
title("Magnitude of X2 with Respect to fr");
xlabel('|X2[fr]|')
ylabel('fr')

%%
% Part 3 - Part 1 DFT with 100 samples
n = 0:99;
N = 100;
fr = linspace(-0.5, 0.49, N);
x1 = exp(1i*2*pi*(1)*n)+exp(1i*2*pi*(33/100)*n);
x2 = cos(2*pi*(1)*n)+0.5*cos(2*pi*(0.3*n));

X3 = fftshift(fft(x1));
X4 = fftshift(fft(x2));

figure();
stem(fr, abs(X3));
title("Magnitude of X1 with Respect to fr");
xlabel('|X1[fr]|')
ylabel('fr')

figure();
stem(fr, abs(X4));
title("Magnitude of X2 with Respect to fr");
xlabel('|X2[fr]|')
ylabel('fr')

%%
% Part 4 - 400-zero padded DFT with 500 samples
n = 0:100;
N = 501;
fr = linspace(-0.5, 0.49, N);
x1 = exp(1i*2*pi*(1)*n)+exp(1i*2*pi*(33/100)*n);
x3 = [zeros(1, 200), x1, zeros(1, 200)];
x2 = cos(2*pi*(1)*n)+0.5*cos(2*pi*(0.3*n));
x4 = [zeros(1, 200), x2, zeros(1, 200)];

X3 = fftshift(fft(x3));
X4 = fftshift(fft(x4));

figure();
stem(fr, abs(X3));
title("Magnitude of X1 with Respect to fr");
xlabel('|X1[fr]|')
ylabel('fr')

figure();
stem(fr, abs(X4));
title("Magnitude of X2 with Respect to fr");
xlabel('|X2[fr]|')
ylabel('fr')

%% B) Sampling
%%% Part 1 - Signal Properties
clc 
clear
load laughter.mat
filename = 'laughter.wav';
audiowrite(filename, y, Fs);
clear y Fs
[y, Fs] = audioread("laughter.wav");

No = length(y);
To = length(y)/Fs;
Ti = 1/Fs;

%Number of Samples
No

%Duration of Singal
To

%Sampling Interval
Ti

%%% Part 2 - Plotted Signal
t = 0:(No-1);
figure();
plot(t, y);
title("Signal y vs Time");


%%% Part 3 - Plot of Discrete Fourier Transformed Signal
Y = fftshift(fft(y));
fr = (-No/2) : ((No/2)-1);

figure();
plot(fr, abs(Y));
title("DFT of audio signal");

%%% Part 4 - Rate 2 Subsampled Signal Properties
rate = 2;
y1 = y(1:rate:end);
No2 = length(y1);
To2 = 2*length(y1)/Fs;
Ti2 = 2/Fs;

%Number of Samples
No2

%Duration of Singal
To2

%Sampling Interval
Ti2

%%% Part 5 - Plot of Rate 2 Subsampled Signal 
figure();
plot(y1);
title("Signal y vs Time");

%%% Part 6 - Plot of Discrete Fourier Transformed Subsampled Signal
Y1 = fftshift(fft(y1));
fr = (-No2/2) : ((No2/2)-1);

figure();
plot(fr, abs(Y1));
title("DFT of audio signal");


%%% Part 7 - Listening to the two signals
sound(y, Fs);
sound(y1, Fs);

%%% Part 8 - Subsample of signal with rate 5
rate = 5;
y2 = y(1:rate:end);
No3 = length(y2);
To3 = 2*length(y2)/Fs;
Ti3 = 2/Fs;

No3

To3

Ti3

figure();
plot(y2);
title("Signal y vs Time");

Y2 = fftshift(fft(y2));
fr = (-No3/2) : ((No3/2)-1);

figure();
plot(fr, abs(Y2));
title("DFT of audio signal");

%% C) Filter Design
%%% Part 1 - Rect Filter
clc
clear
load laughter.mat; 
filename = "laughter.wav"; 
audiowrite (filename, y, Fs); 
[y, Fs] = audioread ("laughter.wav")
; n = length(y);
range = -n/2:((n/2)-1);
period = 1/Fs; p = n*period;
t = (0 : (n-1)); j = 1/p;
f = range*j;
filter_2000 = abs(f) < 2000; Y = fftshift(fft(y.'));
Yfiltered = Y.*filter_2000;
ytime = ifft(fftshift(Yfiltered));
figure(); plot(t,real(ytime));
title ("Plot of filtered sound (2000kHz) in time domain");
xlabel('t'); ylabel('amplitude');
figure(); plot(f, abs(Yfiltered));
title ("Plot of filtered sound (2000kHz) in frequency domain");
xlabel('frequency'); ylabel('|Y(Ω)|');

%%% Part 2 - Playing the Rect Filtered Signal
sound(real(ytime),Fs);

%%% Part 3 - Bass Frequency (16Hz - 256Hz) Filter
[y, Fs] = audioread ("laughter.wav");
n = length(y);
range = -n/2:((n/2)-1);
period = 1/Fs; p = n*period;
t = (0 : (n-1)); j = 1/p;
f = range*j;
filter_16_256 = abs(f) < 16 | abs(f) > 256;
Y = fftshift(fft(y.'));
Yfiltered = Y.*filter_16_256;
ytime = ifft(fftshift(Yfiltered));
figure(); plot(t,real(ytime));
title ("Plot of filtered sound (16-256kHz) in time domain");
xlabel('t'); ylabel('amplitude');
figure(); plot(f, abs(Yfiltered));
title ("Plot of filtered sound (16-256kHz) in frequency domain");
xlabel('frequency');
ylabel('|Y(Ω)|');
sound(real(ytime),Fs);

%%% Part 4 - Treble Frequency (2048Hz - 16384Hz) Amplifying Filter
[y, Fs] = audioread ("laughter.wav");
n = length(y); range = -n/2:((n/2)-1);
period = 1/Fs; p = n*period;
t = (0 : (n-1)); j = 1/p;
f = range*j;
filter_2048_16384 = abs(f) > 2048 | abs(f) < 16384;
Y = fftshift(fft(y.'));
Yfiltered = Y.*filter_2048_16384.*1.25; ytime = ifft(fftshift(Yfiltered));
figure(); plot(t,real(ytime));
title ("Plot of filtered and amplified sound (2048-16384kHz) in time domain");
xlabel('t');
ylabel('amplitude');

figure();
plot(f, abs(Yfiltered));
title ("Plot of filtered and amplified sound (2048-16384kHz) in frequency domain");
xlabel('frequency');
ylabel('|Y(Ω)|');
sound(real(ytime),Fs);

%Part 5 - Question Response