%% Read first

% This script finds optimal investment allocation solving
% Optimizaiton Program 5 (with Uncertainty) and generates figures 2.2 and
% 2.3, which appear in Section 2.3, Considering uncertainty (chapter 2). 

% Additional requirements: 1. Optimization Toolbox must be downloaded
%                          2. Functions required:
%                          optimal_uncertainty.m,invs.m, co.m, f_female_g.m,
%                          f_male_g.m,  
%                          limitdist.m  (this was downloaded from )
%                          3. Data: files deltasp1.mat, deltasp2.mat,
%                          deltasp3.mat. deltasp4.mat (these can also be
%                          generated with the script) and cmap.mat (just to
%                          a personalized colormap)

% This script is part of the project The Evolution of Sex Ratios and Sex Allocation: a
% Theoretical Study, submitted on April 2021 by candidate of the Oxford
% Masters in Mathematical Sciences.


%% RunOptimalUncertainty.m script


% Choose number of females considered
N=500;
% Choose array of maternal condition normally distributed around threshold
% value
cms=normrnd(0.6,0.2,1,N);


newdeltas=0; % Change to 1 to generate new data for Figure 2.3. 
             % Otherwise data is imported to speed computation.
             
% Parameters to generate Figure 2.2 (only one set of uncertainty
% parameters)
k=1;
h=1;
kappas=0.05; % value of kappa (determines uncertainty in mothers prediction)
phis=0.15; % value of phi (determines dimension)
ps=0.2;


if newdeltas==1
% If the user chooses to generate new data for Figure 2.3, 
% choose ranges parameter values:
k=10; % number of phi and kappa values consdiered
h=4; % number of ps values consdiered
kappas=linspace(0,0.4,k);
phis=linspace(0,0.5,k);
ps=linspace(0.1,0.9,h);
% And initialize matrix to store values of Delta (see Report)
Delta_kph=zeros(k,k);
end



% Run loop over all values of uncertainty parameters

for m=1:h
    
   p=ps(m);
   
   % Generate probability transition matrix
   PM=[0.5,0.4,0.1;(p)*0.5,1-p,(p)*0.5;0.1,0.4,0.5];
   
   % Find limiting distribution 
   y=limitdist(PM);

for j=1:k
for t =1:k
        
 
    phi=phis(t);
     
        kappa=kappas(j);

    % Initialize vectors to store predicted fitness of mothers optimizing
    fv1=zeros(1,N); % for s_t=good
    fv2=zeros(1,N); % for s_t=neutral
    fv3=zeros(1,N); % for s_t=bad


    % Initialize vectors to store realized fitness of mothers optimizing
    fr11=zeros(1,N); % for s_t=good, s_t+1=good
    fr12=zeros(1,N); % for s_t=good, s_t+1=neutral
    fr13=zeros(1,N); % for s_t=good, s_t+1=bad

    fr21=zeros(1,N);
    fr22=zeros(1,N);
    fr23=zeros(1,N);

    fr31=zeros(1,N);
    fr32=zeros(1,N);
    fr33=zeros(1,N);

    % Initialize vectors to store realized fitness of random mothers (FR)
    fn1=zeros(1,N); % for s_t+1=good
    fn2=zeros(1,N); % for s_t+1=neutral
    fn3=zeros(1,N); % for s_t+1=bad


    % Run loop over all values of maternal condition:

    for i=1:N
    cm=cms(i);
        % Find optimal sex ratio and investment ratio, and use that to
        % determine predicted and realized fitness.

        % for s_t=good:
        % find optimal strategy
        [sb,i2b,fvalb]=optimal_uncertainty(cm,PM,1,phi);  
        % store predicted fitness
        fv1(i)=fvalb;

        % find realized fitness of mothers optimizing and randomizing sex
        % allocation 
        % for s_t=good, s_t+1=good
        [x]=realized_fitness(cm,sb,i2b,fvalb,1,phi,kappa) ;
        fr11(i)=x(1); 
        fn1(i)=x(2);
        % for s_t=good, s_t+1=neutral
        [x]=realized_fitness(cm,sb,i2b,fvalb,2,phi,kappa) ;
        fr12(i)=x(1);
        fn2(i)=x(2);
        % for s_t=good, s_t+1=bad
        [x]=realized_fitness(cm,sb,i2b,fvalb,3,phi,kappa) ;
        fr13(i)=x(1);
        fn3(i)=x(2);


        % for s_t=neutral:
        [sb,i2b,fvalb]=optimal_uncertainty(cm,PM,2,phi);
        % store predicted fitness
        fv2(i)=fvalb;

        % find realized fitness of mothers optimizing
        % for s_t=neutral, s_t+1=good
        [x]=realized_fitness(cm,sb,i2b,fvalb,1,phi,kappa) ;
        fr21(i)=x(1);
        % for s_t=neutral, s_t+1=neutral
        [x]=realized_fitness(cm,sb,i2b,fvalb,2,phi,kappa) ;
        fr22(i)=x(1);
        % for s_t=neutral, s_t+1=bad
        [x]=realized_fitness(cm,sb,i2b,fvalb,3,phi,kappa) ;
        fr23(i)=x(1);


        % for s_t=bad:
        [sb,i2b,fvalb]=optimal_uncertainty(cm,PM,3,phi);
        % store predicted fitness
        fv3(i)=fvalb;

        % find realized fitness of mothers optimizing
        % for s_t=bad, s_t+1=good
        [x]=realized_fitness(cm,sb,i2b,fvalb,1,phi,kappa) ;
        fr31(i)=x(1);
        % for s_t=bad, s_t+1=neutral
        [x]=realized_fitness(cm,sb,i2b,fvalb,2,phi,kappa) ;
        fr32(i)=x(1);
        % for s_t=bad, s_t+1=bad
        [x]=realized_fitness(cm,sb,i2b,fvalb,3,phi,kappa) ;
        fr33(i)=x(1);

    end

    
    

% find expected realized fitness of mothers optimizing (expectations taken
% over limiting distribution for present and future states)
fM=y(1)*(y(1)*fr11+y(2)*fr21+y(3)*fr31) +y(2)*(y(1)*fr12+y(2)*fr22+y(3)*fr32) +y(3)*(y(1)*fr13+y(2)*fr23+y(3)*fr33) ;

% find expected realized fitness of mothers randomizing(expectations taken
% over limiting distribution for future states)
fR=y(1)*fn1+y(2)*fn2 + y(3) *fn3; 


if newdeltas==1  % if new data is generated, save value of delta for each combination
Delta_kph(j,t)=mean(fM)-mean(fR);
end

end
end


if newdeltas==1 % if new data is generated, save matrix at each
titl=sprintf('p%f.mat', m);
save(titl,'Delta_kph')
end


end


%% Make Figure 2.2

figure(1)
x = tiledlayout(1,4);

nexttile
plot(cms,fv1,'b.',cms,fv2,'k.',cms,fv3,'r.','LineWidth',2)
xlim([0 1])
ylim([0 2])
title('Predicted fitness')
xlabel('Maternal Condition c_m')

nexttile
plot(cms,fr11,'b.',cms,fr21,'k.',cms,fr31,'r.','LineWidth',2)
title('s_{t+1} good')
xlabel('Maternal Condition c_m')
xlim([0 1])
ylim([0 2])

nexttile
plot(cms,fr12,'b.',cms,fr22,'k.',cms,fr32,'r.','LineWidth',2)
xlim([0 1])
ylim([0 2])
title('s_{t+1} neutral')
xlabel('Maternal Condition c_m')

nexttile
plot(cms,fr13,'b.',cms,fr23,'k.',cms,fr33,'r.','LineWidth',2)
xlim([0 1])
ylim([0 2])
title('s_{t+1} bad')
xlabel('Maternal Condition c_m')



% plot figure to see the difference between random
% and optimizing fitness for the chosen uncertainty parameters

figure(2)
plot(cms,fM,'co',cms,fR,'gx')
xlim([0 1])
ylim([0 2])
legend('Optimizing Mothers','Random Mothers')
title('Expected Mother Fitness')
xlabel('Maternal Condition c_m')




%% Make Figure 2.3

% Load data for heatmaps 

x1=load('deltasp1.mat').Delta_kph;
x2=load('deltasp2.mat').Delta_kph;
x3=load('deltasp3.mat').Delta_kph;
x4=load('deltasp4.mat').Delta_kph;

% If new data has been generated, load that.
if newdeltas==1 
x1=load('p1.mat').Delta_kph;
x2=load('p2.mat').Delta_kph;
x3=load('p3.mat').Delta_kph;
x4=load('p4.mat').Delta_kph ;  
end


% plot heatmaps with loaded or generated data 
figure(3)
y = tiledlayout(1,4);
mycolormap=load('cmap.mat').cmap; % import colormap
nexttile
heatmap(x1,'ColorLimits',[-0.3 0.3],'CellLabelColor','none');
grid off
colormap(mycolormap);

nexttile
heatmap(x2,'ColorLimits',[-0.3 0.3],'CellLabelColor','none');
grid off
colormap(mycolormap);

nexttile
heatmap(x3,'ColorLimits',[-0.3 0.3],'CellLabelColor','none');
colormap(mycolormap);
grid off

nexttile
heatmap(x4,'ColorLimits',[-0.3 0.3],'CellLabelColor','none');
grid off
colormap(mycolormap);


