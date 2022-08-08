%% A) Discrete-Time Fourier Series (DTFT)
%%% Part 1 - DTFT
clc
clear
N0 = 128;
n = 0:N0-1;
ohm = (2*pi/128).*(-64:63);

x = [1 (6/7) (5/7) (4/7) (3/7) (2/7) (1/7) zeros(1, 121)];
X = fft(x);

figure;
plot(ohm, fftshift(abs(X)));
title('Magnitude of X[Ω]')
xlabel('Ω')
ylabel('|X[Ω]|')

figure;
plot(ohm, fftshift(angle(X)));
title('Phase of X[Ω]')
xlabel('Ω')
ylabel('<X[Ω]')

% figure;
% plot(ohm, abs(X));
% title('Magnitude of X[Ω]')
% xlabel('Ω')
% ylabel('|X[Ω]|')
% 
% figure;
% plot(ohm, angle(X));
% title('Phase of X[Ω]')
% xlabel('Ω')
% ylabel('<X[Ω]')

%% Part 2 - FTDT by Hand
clc
clear
ohm = (2*pi/128).*(-64:63);

X = 1 + (6/7).*exp(-1i.*ohm) ...
    + (5/7).*exp(-1i.*2.*ohm) ...
    + (4/7).*exp(-1i.*3.*ohm) ...
    + (3/7).*exp(-1i.*4.*ohm) ...
    + (2/7).*exp(-1i.*5.*ohm) ...
    + (1/7).*exp(-1i.*6.*ohm);

figure;
plot(ohm, abs(X));
title('Magnitude of X[Ω]')
xlabel('Ω')
ylabel('|X[Ω]|')



figure;
plot(ohm, angle(X));
title('Phase of X[Ω]')
xlabel('Ω')
ylabel('<X[Ω]')


%%
% 
%  yes the results of part 2 are consistent with part 1.
% 

%% Part 3 - IFFT of Hand Calculated DTFT
clc
clear
N0 = 128;
n = 0:N0-1;
ohm = (2*pi/128).*(-64:63);

X = 1 + (6/7).*exp(-1i.*ohm) ...
    + (5/7).*exp(-1i.*2.*ohm) ...
    + (4/7).*exp(-1i.*3.*ohm) ...
    + (3/7).*exp(-1i.*4.*ohm) ...
    + (2/7).*exp(-1i.*5.*ohm) ...
    + (1/7).*exp(-1i.*6.*ohm);

% X = fftshift(X);

x = ifft(X);

figure;
stem(n, x);
title('Plot of x[n]')
xlabel('n')
ylabel('x[n]')

%%
% 
%  No the result isn't the same as x[n] since the definition of the IFFT
%  command in MATLAB requires an additional shift using the FFTSHIFT
%  command.
% 

%% B) Time Convolution
%%% Part 1 - DTFT plot of x[n]
clc
clear
n = (0:15);
u_c = @(t) 1.0.*(t>=0);
u = @(n) u_c(n).*(mod(n,1) == 0);
x = sin(2*pi*n/10).*(u(n) - u(n-10));

omega= linspace(-pi,pi,1001);
W_omega = exp(-1i).^((0:length(x)-1)'*omega);
X = (x*W_omega);

figure;
stem(n, x);
title('Plot of x[n]')
xlabel('n')
ylabel('x[n]')

figure;
plot(omega, abs(X));
title('Magnitude of X[Ω]')
xlabel('Ω')
ylabel('|X[Ω]|')

figure;
plot(omega, angle(X));
title('Phase of X[Ω]')
xlabel('Ω')
ylabel('<X[Ω]')


% Part 2 - DTFT plot of h[n]
n = (0:9);
x = u(n)-u(n-9);


omega= linspace(-pi,pi,1001);
W_omega = exp(-1i).^((0:length(x)-1)'*omega);
H = (x*W_omega);

figure;
stem(n, x);
title('Plot of h[n]')
xlabel('n')
ylabel('h[n]')

figure;
plot(omega, abs(H));
title('Magnitude of H[Ω]')
xlabel('Ω')
ylabel('|H[Ω]|')

figure;
plot(omega, angle(H));
title('Phase of H[Ω]')
xlabel('Ω')
ylabel('<H[Ω]')

% Part 3 - Convolution Plot of X[Ω] and H[Ω]

Y = X.*H;

figure;
plot(omega, abs(Y));
title('Magnitude of Y[Ω]')
xlabel('Ω')
ylabel('|Y[Ω]|')

figure;
plot(omega, angle(Y));
title('Phase of Y[Ω]')
xlabel('Ω')
ylabel('<Y[Ω]')



%% Part 4 - Convolution Plot of x[n] and h[n] by conv command
clc
clear
n = (0:15);
u_c = @(t) 1.0.*(t>=0);
u = @(n) u_c(n).*(mod(n,1) == 0);

h = u(0:9);
x = sin(2*pi*n/10).*(u(n) - u(n-10));
n=0:24;

y = conv(x, h);

figure;
stem(n, y);
title('plot of y[n] using conv command')
xlabel('n')
ylabel('y[n]')

% Part 5 - DTFT Plot of y[n] from Part 4

omega= linspace(-pi,pi,1001);
W_omega = exp(-1i).^((0:length(y)-1)'*omega);
Y = (y*W_omega);

figure;
plot(omega, abs(Y));
title('Magnitude of Y[Ω] from part 4')
xlabel('Ω')
ylabel('|Y[Ω]|')

figure;
plot(omega, angle(Y));
title('Phase of Y[Ω]')
xlabel('Ω')
ylabel('<Y[Ω]')

%%
% 
%  Yes the same results were achieved in part 3 and 5, this is due to a
%  property of the fourier transform; The covolution of two functions in 
%  the time domain is equivelent to the product in the frequency domain.
% 

%% C) FIR Filter Design by Frequency Sampling
%%% Part 1 - High Pass FIR Filter
ohm0 = 2*pi/3;
N = 35;
n = 0:N-1;
Omega = linspace(0,2*pi*(1-1/N),N);
H_d = @(Omega) (mod(Omega,2*pi)>ohm0).*(mod(Omega,2*pi)<2*pi-ohm0);

H = H_d(Omega).*exp(-1i.*Omega.*((N-1)/2)); 
h = ifft(H);

figure;
plot(Omega, abs(H));
title('Original Highpass Filter Magnitude |Y[Ω]|')
xlabel('Ω')
ylabel('|H[Ω]|')

figure;
plot(Omega, angle(H));
title('Original Highpass Filter Phase <H[Ω]')
xlabel('Ω')
ylabel('<H[Ω]')

figure;
stem(n, h);
title('Impulse Response (h[n]) of Filter H[Ω]')
xlabel('n')
ylabel('h[n]')


%% Part 2 - Frequency Response from h[n] by freqz Command
H = freqz(h,1, 0:2*pi/1001:2*pi);
figure;
plot(0:2*pi/1001:2*pi, abs(H));
title('Highpass Filter Magnitude by freqz command')
xlabel('Ω')
ylabel('|Y[Ω]|')

%%
% 
%  Part 3 The result in part 2 is different from the ideal filter we 
%  started with in terms of the frequencies it permits. The ideal filter
%  has a strict cutoff frequency, whereas the result in part 2 has ripples
%  beyond the desired cutoff frequency. This is due to the impulse response
%  of an ideal highpass filter being noncausal and therefore physically
%  unrealizable 
% 

%% Part 4 - 71 points
clc
clear
ohm0 = 2*pi/3;
N = 71;
n = 0:N-1;
Omega = linspace(0,2*pi*(1-1/N),N);
H_d = @(Omega) (mod(Omega,2*pi)>ohm0).*(mod(Omega,2*pi)<2*pi-ohm0);

H = H_d(Omega).*exp(-1i.*Omega.*((N-1)/2)); 
h = ifft(H);

figure;
plot(Omega, abs(H));
title('Original Highpass Filter')
xlabel('Ω')
ylabel('|Y[Ω]|')

figure;
plot(Omega, angle(H));
title('Original Highpass Filter Phase <H[Ω]')
xlabel('Ω')
ylabel('<H[Ω]')

figure;
stem(n, h);
title('Impulse Response (h[n]) of Filter Y[Ω]')
xlabel('n')
ylabel('h[n]')


%%
% 
%  Part 5
%  Increasing N increases the 'resolution' since we get a more accurate
%  representation of the impulse response due to the existance of more points 
%  and a more unique signal representation.
% 

