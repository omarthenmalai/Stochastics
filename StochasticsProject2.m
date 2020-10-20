%% Stochastics Project #2
% Omar Thenmalai

%% Scenario 1

% Bayes MMSE
clc;
clear;
num_simulations = 10000; % number of simulations
w_max = 2; % Max value of W
w_min = -2; % Min value of W
y_max = 1; % Max value of Y
y_min = -1; % Min value of Y
Y = unifrnd(y_min, y_max, 1, num_simulations); % Uniform distribution from y_min to y_max
W = unifrnd(w_min, w_max, 1, num_simulations); % Uniform distribution from w_min to w_max
X = Y + W; 
Y_est_bayes = zeros(1, num_simulations);
for i=1:num_simulations
    % The estimate of y is the mean of all of the possible values of Y
    % given a specific value of X.
    % Since Y is bounded by (-1,1), make sure that the value being taken
    % for the mean falls within these bounds using the max and min
    % functions
    Y_est_bayes(i) = (max(X(i)+w_min, y_min) + min(X(i)+w_max, y_max))/2;
end
bayes_empirical_mse = immse(Y, Y_est_bayes); % MSE from observed predictions for Y
bayes_theoretical_mse = 1/4; % Theoretical MSE found in textbook

% Linear MMSE
mu_y = mean(Y); % Mean of Y
mu_x = mean(X); % Mean of X
var_y = var(Y); % variance of Y
var_w = var(W); % Variance of W 
var_x = var(X); % Variance of X
covar_yx = var_y; % The covariance of Y equals the variance of Y 
Y_est_linear = mu_y + (covar_yx/var_x)*(X-mu_x); % Formula for the linear estimate
linear_empirical_mse = immse(Y, Y_est_linear); % The MSE based on the observed predictions for Y
linear_theoretical_mse = var_y*(1-covar_yx/var_x); % The theoretical MSE

% Table Output
bayes_table = table(bayes_empirical_mse, bayes_theoretical_mse);
linear_table = table(linear_empirical_mse, linear_theoretical_mse);
display(bayes_table);
display(linear_table);

%% Scenario 2
num_points = 100000;
variances = linspace(1, 10, 50);
empirical_mse = zeros(1, length(variances));
theoretical_mse = zeros(1, length(variances));
num_obs = 2;
for i=1:length(variances)
    [empirical_mse(i), theoretical_mse(i)] = multiple_noisy_observations(variances(i), 1, num_obs, num_points);
end
figure;
hold on;
title('Empirical and Theoretical MSE for 2 noisy observations with respect to the Variance of Y')
plot(variances, empirical_mse);
plot(variances, theoretical_mse);
legend('Empirical MSE', 'Theoretical MMSE')
ylabel("MSE")
xlabel("Variance of Y (VAR(R) = 1)")
hold off


for i=1:length(variances)
    [empirical_mse(i), theoretical_mse(i)] = multiple_noisy_observations(1, variances(i), num_obs, num_points);
end
figure;
hold on;
title('Empirical and Theoretical MSE for 2 noisy observations with respect to the Variance of R')
plot(variances, empirical_mse);
plot(variances, theoretical_mse);
legend('Empirical MSE', 'Theoretical MMSE')
ylabel("MSE")
xlabel("Variance of R (VAR(Y) = 1)")
hold off


function [empirical_mse, theoretical_mse]= multiple_noisy_observations(var_y, var_r, num_observations, num_points)
    mu_y = 1; % Assume that the mean of Y is 1
    mu_r = 0; % Assume that the mean of all of R's is 1
    Y = normrnd(mu_y, var_y, [num_points, 1]); % Generate a set of points from a gaussian distribution with mean mu_y, and variance var_y
    X_sum = var_r * mu_y; % The estimated values for Y are determined by some coefficient based on the variances of Y and R times a sum of the variance of Y and the values of X
    for i=1:num_observations 
        X = Y + normrnd(mu_r, var_r, [num_points, 1]); % Xi = Y + Ri
        X_sum = X_sum + var_y*X; % Add var_y*Xi to the sum
    end
    Y_est = (1/(2*var_y + var_r))*X_sum; % the estimate for Y equals the sum of these X values and VAR(Y) times a coefficient
    empirical_mse = immse(Y, Y_est); % Get the observed MSE
    theoretical_mse = var_y*var_r/(2*var_y + var_r);  % Use the formula from the notes to find the theoretical MSE for the given variances
end


