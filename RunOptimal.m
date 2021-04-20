%% Read first

% This script finds optimal investment allocation solving
% Optimizaiton Program 1 (Basic) and generates figure 2.1.

% Additional requirements: 1. Optimization Toolbox must be downloaded
%                          2. Functions required: optimal_inv.m, inv.m, 

% This script is part of the project The Evolution of Sex Ratios and Sex Allocation: a
% Theoretical Study, submitted on April 2021 by candidate of the Oxford
% Masters in Mathematical Sciences.


%% RunOptimal.m script

% Choose number of females considered
N=50;
% Choose array of maternal condition
cms=linspace(0.001,1,N);

% Initialize variables to store optimal progeny sex and investemnt ratios

% For linear fitness
sex_ratio_l=zeros(1,N);
inv_ratio_l=zeros(1,N);
im_l=zeros(2,N);
if_l=zeros(2,N);

% For male-convex fitness
sex_ratio_x=zeros(1,N);
inv_ratio_x=zeros(1,N);
im_x=zeros(2,N);
if_x=zeros(2,N);

% For logistic fitness
sex_ratio_g=zeros(1,N);
inv_ratio_g=zeros(1,N);
im_g=zeros(2,N);
if_g=zeros(2,N);


% Find optimal sex ratio and investment ratio 

for i=1:N
cm=cms(i);
% For linear fitness
Func=1;
[s,i1,i2,i3,fval]=optimal(cm,Func);
sex_ratio_l(i)=s;
inv_ratio_l(i)=i2(1)/invs(cm);
im_l(:,i)=i1;
if_l(:,i)=i3;


% For linear fitness
Func=2;
[s,i1,i2,i3,fval]=optimal(cm,Func);
sex_ratio_x(i)=s;
inv_ratio_x(i)=i2(1)/invs(cm);
im_x(:,i)=i1;
if_x(:,i)=i3;


% For logistic fitness
Func=3;
[s,i1,i2,i3,fval]=optimal(cm,Func);
sex_ratio_g(i)=s;
inv_ratio_g(i)=i2(1)/invs(cm);
im_g(:,i)=i1;
if_g(:,i)=i3;

end

% Plot results optimal ratios and fitness functions (figure 2.1)

figure(1)
y = tiledlayout(3,2);

% plot linear
nexttile
plot(cms,f_male_l(cms),'k','LineWidth',2)
hold on
plot(cms,f_female_l(cms),'g--','LineWidth',2)
legend('Male Fitness','Female Fitness');
xlabel('Offspring Condition c');
title('Linear fitness')

nexttile
plot(cms,sex_ratio_l,'b.','LineWidth',2)
hold on
plot(cms,inv_ratio_l,'r.','LineWidth',2)
legend('Sex-ratio','Investment ratio');
xlabel('Maternal Condition c_m');

% plot male-convex
nexttile
plot(cms,f_male_x(cms),'k','LineWidth',2)
hold on
plot(cms,f_female_x(cms),'g--','LineWidth',2)
legend('Male Fitness','Female Fitness');
xlabel('Offspring Condition c');
title('Male-convex fitness')

nexttile
plot(cms,sex_ratio_x,'b.','LineWidth',2)
hold on
plot(cms,inv_ratio_x,'r.','LineWidth',2)
legend('Sex-ratio','Investment ratio');
xlabel('Maternal Condition c_m');

% plot logistic
nexttile
plot(cms,f_male_g(cms),'k','LineWidth',2)
hold on
plot(cms,f_female_g(cms),'g--','LineWidth',2)
legend('Male Fitness','Female Fitness');
xlabel('Offspring Condition c');
title('Male-convex fitness')

nexttile
plot(cms,sex_ratio_g,'b.','LineWidth',2)
hold on
plot(cms,inv_ratio_g,'r.','LineWidth',2)
legend('Sex-ratio','Investment ratio');
xlabel('Maternal Condition c_m');




% Check that investment in all-male and all-female brood is equitable

figure(2)
y = tiledlayout(1,3);

half=invs(cms)/2;

nexttile
plot(cms,im_l(1,:),cms,if_l(1,:),cms, half,'LineWidth',2)

legend('inv in male of all-male','inv in female of all-female','1/2 of max inv');
xlabel('Offspring Condition c');
title('investment in all-male and all-female brood, linear fitness')

nexttile
plot(cms,im_x(1,:),cms,if_x(1,:),cms, half,'LineWidth',2)
legend('inv in male of all-male','inv in female of all-female','1/2 of max inv');
xlabel('Offspring Condition c');
title('investment in all-male and all-female brood, male convex fitness')

nexttile
plot(cms,im_g(1,:),cms,if_g(1,:),cms, half,'LineWidth',2)
legend('inv in male of all-male','inv in female of all-female','1/2 of max inv');
xlabel('Offspring Condition c');
title('investment in all-male and all-female brood, logistic fitness')
