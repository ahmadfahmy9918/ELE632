%%
%% A) Unit Impulse Response
%%% Part 1 - Filter to Recieve Unit Impulse Response
% * Function I
n =  0:10;

b = [1/3 0 0];
a = [1 1/6 -1/6];
delta = (n)==0;
% figure;
% stem(n, delta);
% title('Part A.1 - Unit Impulse')
% axis([-0.5 0.5 -0.1 1.2])
% xlabel('n')


h = filter(b, a, delta);

figure;
stem(n, h, 'k');
title('Part A.1 - Function I Unit Impulse Response')
axis([-0.5 10.5 -0.1 0.5])
xlabel('n')
ylabel('h[n]')

%%
%%% Function II
n =  0:10;

b = [1 0 0];
a = [1 0 1/4];
delta = (n)==0;

h = filter(b, a, delta);

figure;
stem(n, h, 'k');
title('Part A.1 - Function II Unit Impulse Response')
axis([-0.5 10.5 -0.3 1.1])
xlabel('n')
ylabel('h[n]')

%% B) Zero Input Response
n =  0:50;

b = [2 0 0];
a = [1 -3/10 -1/10];


z_i = filtic(b, a, [1 2]);
y_0 = filter(b, a, zeros(size(n)), z_i);

figure
stem(n, y_0, 'k');
title('Part B - Zero Input Response')
axis([-2 10.5 -0.3 1.1])
xlabel('n')
ylabel('y0[n]')

%% C) Zero State Response
n =  0:50;
u = (mod(n,1)==0)*1.0.*(n>=0);
u2 = (mod(n,1)==0)*1.0.*(n>=10);

b = [2 0 0];
a = [1 -3/10 -1/10]; 
x = 2.*cos((2.*pi.*n)./6).*(u-u2);

z_i = filtic(b, a, 0);
y_0 = filter(b, a, x, z_i);

figure
stem(n, y_0, 'k');
title('Part C - Zero State Response')
axis([-2 18 -4 4.1])
xlabel('n')
ylabel('y0[n]')

%% D) Total Response
%%% Part I - Total response
n =  -50:50;
u = (mod(n,1)==0)*1.0.*(n>=0);
u2 = (mod(n,1)==0)*1.0.*(n>=10);

b = [2 0 0];
a = [1 -3/10 -1/10]; 
x = 2.*cos((2.*pi.*n)./6).*(u-u2);

z_i = filtic(b, a, [1 2]);
y_0 = filter(b, a, x, z_i);

figure
stem(n, y_0, 'k');
title('Part D.I - Total Response')
axis([-1 18 -4 4.1])
xlabel('n')
ylabel('y0[n]')

%%
%%% Part II - Total response by Summation of Zero Input and Zero state
%%% Responses
n =  -50:50;
u = (mod(n,1)==0)*1.0.*(n>=0);
u2 = (mod(n,1)==0)*1.0.*(n>=10);
b = [2 0 0];
a = [1 -3/10 -1/10];
x = 2.*cos((2.*pi.*n)./6).*(u-u2);

z_i = filtic(b, a, [1 2]);
y_0i = filter(b, a, zeros(size(n)), z_i);

z_i2 = filtic(b, a, 0);
y_0 = filter(b, a, x, z_i2);

yt = y_0 + y_0i;

figure
stem(n, y_0, 'k');
title('Part D.II - Total Response by Summation')
axis([-1 18 -4 4.1])
xlabel('n')
ylabel('y0[n]')

%% E) Convolution and System Stability
%%% Part I - conv command
n =  0:50;
u = (mod(n,1)==0)*1.0.*(n>=0);
u2 = (mod(n,1)==0)*1.0.*(n>=10);
b = [2 0 0];
a = [1 -3/10 -1/10];
x = 2.*cos((2.*pi.*n)./6).*(u-u2);
% delta = @(n) 1.0.*(n==0);
delta = (n)==0;

h = filter(b, a, delta);

y = conv(x, h);

n =  0:100;

figure
stem(n, y, 'k');
title('Part E - Zero State Response with Conv Function')
axis([-1 18 -4 4.1])
xlabel('n')
ylabel('y0[n]')

%%
%%% Part II
%Yes the result is the same

%%
%%% Part III
%Yes the system is asymptomatically stable, this is because the system's
% Characteristic equation has two real  roots within the unit circle in the
% Complex Plane 

%% F) Moving Average Filter
%%% Part I - Constant Coefficient Difference Equation with h[n] Impulse
%%% Response

%%
%%% Part II & III - MATLAB Function for N-point moving-average filter
n = 0:1:45;
d = (n-30)==0;
d2 = (n-35)==0;
a = 1;

x = cos((pi.*n)./5)+d-d2;

%original function
figure
stem(n, x);
title('Part F - x[n]')
axis([0 45 -2.2 2.2])
xlabel('n')
ylabel('h[n]')

%N=4
filterN = 4;

b = (1/filterN)*ones(1, filterN);
h = filter(b, a, x);
figure
stem(n, h);
title('Part F - x[n] after Moving Average Filter of Window Size 4')
axis([0 45 -1.05 1.05])
xlabel('n')
ylabel('h[n]')

%N=8
filterN = 8;
b = (1/filterN)*ones(1, filterN);
h = filter(b, a, x);
figure
stem(n, h);
title('Part F - x[n] after Moving Average Filter of Window Size 8')
axis([0 45 -0.5 0.5])
xlabel('n')
ylabel('h[n]')

%N=12
filterN = 12;
b = (1/filterN)*ones(1, filterN);
h = filter(b, a, x);
figure
stem(n, h);
title('Part F - x[n] after Moving Average Filter of Window Size 12')
axis([0 45 -0.25 0.25])
xlabel('n')
ylabel('h[n]')

%As the window size (N) of the moving-average filter increases the
%amplitude of the resulting sinusoidal function decreases. This is
%because the filter takes an average over a larger number of values while
%the number of repeating y-values and the number of y-values that 
% cancel out increases 