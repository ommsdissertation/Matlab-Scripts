%% Read first

% This script finds optimal investment allocation solving
% Optimizaiton Program 6 (Investment Only) and generates figure 2.4.

% Additional requirements: 1. Optimization Toolbox must be downloaded
%                          2. Functions required: optimal_inv.m, invs.m, f_male_f.m, co.m 

% This script is part of the project The Evolution of Sex Ratios and Sex Allocation: a
% Theoretical Study, submitted on April 2021 by candidate of the Oxford
% Masters in Mathematical Sciences.


%% RunOptimalInv.m script

% Choose number of females considered
N=50;
% Choose array of maternal condition
cms=linspace(0.01,1,N);
% Choose array of newborn population sex ratio
alphas=linspace(0.01,0.99,N);

% Initialize matrix to store optimal investemnt ratio for every 
% combination of maternal condition and newborn population sex ratio
inv_ratio=zeros(N,N);

% Run loops over all values of maternal condition and newborn population sex ratio
for j=1:N
   a= alphas(j);
for i=1:N
   cm=cms(i);

   [i1,i2]=optimal_inv(cm,a); % Find optimal investemnt to male using function optimal_inv.m, 
                              % which implements using the Optimization
                              % Toolbox function fmincon.

   inv_ratio(j,i)=i1/invs(cm); % Find investment ratio, defined as the proportion of total
                               % investment resources invested in the male
                               % of a mixed brood.


end
end

% Plot results as a Heatmap
figure(1)
heatmap(inv_ratio)
colormap parula
xlabel('Maternal Condition');
ylabel('Newborn Sex-Ratio \tilde{\alpha}');
titel('Optimal Investment Ratio')
grid off