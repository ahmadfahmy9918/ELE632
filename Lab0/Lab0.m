%% Intro to Matlab
% Matlab is an interpreted language. Variables do not need to be declared,
% but Matalb assume the variable type. Be careful not to name variables
% with the same name as functions. If you have trouble with any Matlab
% functions, you can visit Matlab's help website as it usually comes with
% tutorials with examples.
%% Things that the TA should show you
% * Workspace (how to view individual variables)
% * The command window and the editor windows
% * Difference between a script and a function
% * How to publish a document
% * How to run, run section, and run with default inputs (for functions)
%% How to plot
% You should already know how to plot from previous courses with Matlab,
% but in this course, we want to make a strong distinction between plot and
% stem commands since we are differentiating between discrete and
% continuous domains.

%%% 
% Here is an example for the plot command. Make sure to get used to the
% version where you give both x and y arguments instead of the version
% where you only give y arguments. In this example, note how I am using
% plot on a "continuous" signal where I just specify the x to have
% non-integer increments

x = -1:0.01:1;
y = sin(pi*x);
figure;
plot(x,y);
title('Continuous sinusoid')


%%
% In this next example, x is going to be a discrete signal where I specify
% x to be in integer increments. Therefore, I will be using the stem
% command. Note that when I define x, I use the notation where I don't
% specify the increment amount, since the default increment is 1.
% Additionally, to publish in such a way that 

x = -10:10;
y = sin(0.1*pi*x);
figure;
stem(x,y);
title('Discrete sinusoid')


%% close, clear, and clc commands
% * close: closes the currently selected figure. close all  closes all
% figures
% * clear: removes variables from the workspace. clear all removes all
% variables from the workspace while clear <variable> removes only the
% variable specified by <variable>
% * clc: clears the command window. This one is mostly user preference
%
% When coding your labs, get into the habit of using close all and clear
% all at the beginning of each section of code. The most common mistake
% students will make is to accidentally use variables from previous
% sections while thinking they represent something else. The best example
% is accidentally reusing time variables, especially if the time regions
% are different in each section. Additionally, students sometimes wonder
% what happened to their figures, when in the previous section they had a
% figure with subfigures, and then they call the plot command in their
% current section.

% create some variables
x = 1:10;
y = 11:20;

% clear only x

clear x

%% save and load
% save and load are pretty self explanitory as they save and load variables
% to the workspace. Their main usage would be to save results of scripts
% that take a long time to process, or if you want to save a specific
% variariable before clearing the workspace, and then reload those
% variables after you clear it.

% create variables x and y (unsupressed so we can see output)
x = 1:10
y = 11:20

%%

% save only variable x
save('save_example', 'x');

% call clear to remove the variables
clear;


% create a new variable x with different values
x = 21:30

%%

% load your save
load('save_example');

% your x variable should now contain the old values
x

%% Loops and enumeration
% The syntax for for loops is: 
% for index = values
%   statements
% end
%
% where index is the iterator for the for loop which you will reference in
% your statements, and values is an array that you can initialize inside of
% the statement for before hand. Typically the normal use case is:

for i = 1:3
    disp(sprintf('this is iteration %d', i)); 
end

%%
% you can also do this:

x = [5 2 10 -1];

for i = x
    disp(sprintf('the current value of i is %d', i)); 
end

%% Matrix indexing
% In Matlab, we will mostly be dealing with matricies and matrix operations
% as they are generally more computationally efficient that loops (this is
% specific to Matlab). If you have a matrix, you can index into it like so:

x = [101:105; 111:115; 121:125]
indices = [10, 5, 3, 14, 7]
x(indices)

%%%
% note that we are using matrix indices, therefore we go from top to bottom
% before going left to right.
%%
% Additionally, if you have appropriately sized arrays, you can concatenate
% then like this:

y = [131:135]
xy_vertical = [x;y]

%%
% or like this:

z = [106:110;116:120; 126:130]
xz_horizontal = [x,z]

%% Displaying Results
% Matlab has multiple different ways of printing the results to the
% standard output. The standard out (stdout) will depend on if you are
% running the code or publishing the code. When running the code the stdout
% will be the command window and when publishing the code stdout with be
% the document you are publishing to. The easiest ways to outputing to
% stdout are unsupressing the output and using the disp command.

%%%
% Unsuppressing an output will cause both the variable name and it's value
% to be outputted to std out. You can unsuppress a code snippit by removing
% the semicolon at the end of the statement.
X1 = 53; % This is a suppressed output and will not show up in stdout
X2 = 21  % This is now unsuppressed and will now output to stdout


%%
% The disp command will only show values of the variable and so needs to be
% combined with text for context. This can be done using formatting from
% sprintf.

text = sprintf('X1 = \n\n\t\t\t\t%d', X1);
disp(text);
%% Booleans
% Boolean logic can be used in multiple different ways in matlab to access
% arrays or supress values. In matlab boolean are stored as 1s and 0s,
% however it does have it's own type that behaves differently.

DoubleMatrix = [1.0 1.0 0; 0 1.0 1.0; 0 0 0]
BooleanMatrix = (DoubleMatrix == 1);
A = [11 12 13; 14 15 16; 17 18 19]
% B(DoubleMatrix); This will error.
C = A(BooleanMatrix)
%%
% When setting a matrix like this the value of the new array is a 1D array
% with the values of the original matrix which had a true statement from
% the boolean matrix.
%
% Alternatively we can adjust the matrix using a boolean matrix on the left
% hand side of the equation

D = A;
D(BooleanMatrix) = D(BooleanMatrix) - 10

%%

E = A;
E(E < 14) = E(E < 14) * -1  % This is equal to E(A < 14) = E(A < 14) * -1;

%%
% When you multiply with a boolean matrix it will first be typecast into a
% double, behaving as a mask.

F = A .* BooleanMatrix

%%SAMPLING -------------------
% In multiple parts of your labs, you will be dealing with discrete
% functions which need to have the behaviour where the output needs to be
% zero if the input is non-integer. This behaviour can be replicated by
% using logicals:

n = [-5:0.1:5];
u = (mod(n,1)==0)*1.0.*(n>=0);
figure;
stem(n,u);

%%%
% In the discrete unit step function we just created, the logic statement
% on the left makes any output 0 if the input is non-integer. The logic
% statement on the right makes the output 0 if the input is negative. The
% middle *1.0 is to just have the output be saved as a double instead of a
% boolean.

%%
% Alternatively, if we want a continuous unit step function, we can skip
% the first logic statement:

n = [-5:0.1:5];
u = 1.0.*(n>=0);
figure;
plot(n,u);

%% Versioning
% When using different versions of Matlab certain functions change name or
% functionallity. To get around this issue you can check the version and
% use an if statement to get around the issue or throw an error.

%%%
% Here is an example in practice:
if verLessThan('matlab','8.4')
    error('Matlab version outdated');
elseif verLessThan('simulink','8.4')
    error('Simulink 4.0 is required');
else
    % myFunction();
end

%%
% A specific example is when you want a title over a plot with multiple
% subplots in older versions you use suptitle, however in version 9.5 the
% function was removed in favour of sgtitle.

y1 = 1:2:20;
y2 = 10:-1:1;
x = 1:10;
figure;
subplot(1,2,1);
plot(x,y1);
title('Plot 1');
subplot(1,2,2);
plot(x,y2);
title('Plot 2');

if verLessThan('matlab','9.5')
    suptitle('Plot 1 vs Plot 2');
else
    sgtitle('Plot 1 vs Plot 2');
end

%% Anonymous functions
% There are three points to the function: 
%
%     Name - This is what is used to call the function
%     Variables - These are your inputs
%     Calculation - This is what you want your function to do
%
% In anonymous functions the name takes the space of the variable to the
% left of the equals sign.
%
% Name = ...
%
%
% The variables are then placed in paranthesis after the @ symbol and are
% separated by commas.
%
% Name = @(Variable1,Variable2,...) ...
%
%
% The calculation behaves the same as if you were just writing the function
% normally.
%
% Name = @(Variable1,Variable2,...) 5*Variable1 + 2*Variable2;
%
%
% Anonymous functions are called the same way regular functions are called,
% so in this example:
%
% Name(x,y)

%%%
% Here is an example in practice:
a = 1;
b = 5;
c = 6;

root1 = @(x1,x2,x3) (x2 + sqrt(x2.^2-4*x1*x3))/(2*x1);
root2 = @(x1,x2,x3) (x2 - sqrt(x2.^2-4*x1*x3))/(2*x1);

ans1 = root1(a,b,c)
ans2 = root2(a,b,c)

