%% Stochastics Project #5 - Omar Thenmalai
% Part 2 - MMSE estimation for filters length N=4,6,10
N = [4;6;10];
num_values = 100000;
s = randi([0,1],1,num_values);
s(s == 0) = -1;
sigma_squared = 1; % assume variance of noise is 1, zero-mean
c = [1, 0.2, 0.4]; % impulse response of filter c
r = filter(c, 1, s); % s after passing through filter c 
d = sqrt(sigma_squared)*randn(1,length(r)); % additive gaussian white noise
r = r + d; % signal after noise has been added
MSE = zeros(length(N),1);
for i=1:length(N)
    Rss = 1;
    Rsr = conv(Rss,c);
    C_n = [zeros(1,N(i)-2), c];
    C_negn = flip(C_n);
    Rrr = conv(conv(Rss,C_n),C_negn); % Get Rrr
    Rrr(N(i)+1) = Rrr(N(i)+1) + sigma_squared; % Add sigma_squared to Rrr(0)
    
    temp = [c, zeros(1,N(i)-length(c))]';
    Rxx_mat = zeros(N(i),N(i));
    for j=N(i)+1:2*N(i);
        Rxx_mat(j-N(i),:) = flip(Rrr(1, j-N(i)+1:j));
    end
    h = (Rxx_mat\temp)'; % impulse response of filter h
    s_hat = filter(h, 1, r); % MMSE estimation of s_hat
    MSE(i) = mean((s_hat-s).^2); % Mean squared error of input signal s and MMSE prediction s_hat
end

table(N, MSE) % Table of Filter lengths N and their respective MSEs
% MSE does not seem to significantly change as the length of the filter, H,
% increases.



