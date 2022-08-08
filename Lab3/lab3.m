%% A) Discrete-Time Fourier Series (DTFS)
%%% Part 1 - Fundamental Period & Fundamental Frequency
clc
clear
n = 0:20;
x = 4*cos(2.4*pi*n)+2*sin(3.2*pi*n);
mN1 = (2.4*pi)/(2*pi);
mN2 = (3.2*pi)/(2*pi);
%m/N1 = 6/5
%m/N2 = 8/5
N0 = 5;
Om0 = 2*pi/N0;

%%% Part 2 - DTFS
clc
clear
N0 = 5;
Om0 = 2*pi/N0;
n = 0:N0-1;
x = 4*cos(2.4*pi*n)+2*sin(3.2*pi*n);
for(r = 0:N0-1)
    Dr(r+1)=(sum(x.*exp(-1i.*r.*n.*Om0)))/N0;
end

figure;
stem(n, x);
title('Plot of x(n) DTFS')
xlabel('n')
ylabel('x[n]')

figure;
stem(abs(Dr));
title('Magnitude of x(n) DTFS')
xlabel('r')
ylabel('|Dr|')

figure;
stem(angle(Dr));
title('Phase of x(n) DTFS')
xlabel('r')
ylabel('<Dr')

%% Part 3 - P1 and P2 with Lab3 Manual Figure 1
clc
clear
N0 = 6;
Om0 = 2*pi/N0;

n = 0:N0-1;
y = [3 2 1 0 1 2];
for(r = 0:N0-1)
    Dr(r+1)=(sum(y.*exp(-1i.*r.*n.*Om0)))/N0;
end

figure;
stem(n, y);
title('Plot of y(n) DTFS')
xlabel('n')
ylabel('y[n]')

figure;
stem(abs(Dr));
title('Magnitude of y(n) DTFS')
xlabel('r')
ylabel('|Dr|')
grid on

figure;
stem(angle(Dr));
title('Phase of y(n) DTFS')
xlabel('r')
ylabel('<Dr')
grid on


%% B) Inverse DTFS and Time Shifting Property
%%% Part 1 - Inverse DTFS
clc
clear
N0 = 32;
Om0 = 2*pi/N0;
n = 0:N0-1;

xr = [ones(1,5) zeros(1,23) ones(1,4)];
% x = real(ifft(xr).*N0);
x = ifft(xr).*N0;

figure;
stem(n, x);
title('Plot of x(n) IDTFS')
xlabel('n')
ylabel('x[n]')


%% Part 2 - Time Shifting
clc
clear
N0 = 32;
Om0 = 2*pi/N0;
n = 0:N0-1;

xr = [ones(1,5) zeros(1,23) ones(1,4)];
X = xr.*exp(-1i.*5.*n.*Om0);
% X = real(ifft(X)*N0);
X = ifft(X)*N0;

figure;
stem(n, X);
title('Plot of Shifted X(r) IDTFS')
xlabel('n')
ylabel('X[r]')


%% C) System Response
%%% Part 1 - Plot of H[r] for Lab3 Figure 3 with Respect to Omega
clc
clear
N0 = 32;
Om0 = 2*pi/N0;
r = 0:N0-1;

Hr = [ones(1,5) zeros(1, 23) ones(1,4)];

figure;
stem(Om0.*r, Hr);
title('H[r] for Lab3 Figure 3 with Respect to Omega')
xlabel('r')
ylabel('H[r]')

%% Part 2 - Plot of System Response for Lab3 Figure 3 with input X1[n]
clc
clear
N0 = 32;
Om0 = 2*pi/N0;
r = 0:N0-1;

Hr = [ones(1,5) zeros(1, 23) ones(1,4)];
x = 4*cos(pi*r/8);
X = fft(x);

figure;
stem(r, X);
title('Plot of X[r]')
xlabel('n')
ylabel('X[r]')

Y = X.*Hr;
y = ifft(Y);

figure;
stem(r, y);
title('Plot of y[n]')
xlabel('n')
ylabel('y[n]')

%% Part 3 - Plot of System Response for Lab3 Figure 3 with input X2[n]
clc
clear
N0 = 32;
Om0 = 2*pi/N0;
r = 0:N0-1;

Hr = [ones(1,5) zeros(1, 23) ones(1,4)];
x = 4*cos(pi*r/2);
X = fft(x);




Y = X.*Hr;
y = ifft(Y);

figure;
stem(r, y);
title('Plot of y[n]')
xlabel('n')
ylabel('y[n]')