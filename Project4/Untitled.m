%% Stochastics Project #5 - Omar Thenmalai
% Part 2 - MMSE estimation for filters length N=4,6,10
N = [4,6,10];
num_values = 100;
S = randi([0,1],1,num_values);
S(S == 0) = -1;
