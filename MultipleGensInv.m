%% Read first

% This script simulates the evolution of a population in which mothers
% determine investment allocation according to Optimizaiton Program 6 (Investment Only).
% The pseudo-code of the alogrithm can be found in the Appendix.

% Additional requirements: 1. Optimization Toolbox must be downloaded
%                          2. Functions called:
%                          optimal_inv.m,invs.m, co.m, f_female_x.m,
%                          f_male_x.m,  
%                          limitdist.m  (this was adapted from  http://www.math.wustl.edu/~feres/Math450Lect04.pdf )
% 

% This script is part of the project The Evolution of Sex Ratios and Sex Allocation: a
% Theoretical Study, submitted on April 2021 by candidate 1048690 of the Oxford
% Masters in Mathematical Sciences.



%% MultipleGensInv.m script

addpath('supp_functions')

% choose max number of generations simulated
gens=100;

% choose uncertainty parameters
phi=0.15;
p=0.2;
% choose number of values of parameter phi considered (total k1*k2)
k=4;

kappas=0.05;
phis=linspace(0,0.8,k1*k2);

% choose number of repetitions of the algorithm
itermax=50;

iter=1;
while  iter<itermax
% make probability transition matrix
PM=[0.5,0.4,0.1;(p)*0.5,1-p,(p)*0.5;0.1,0.4,0.5];
y=limitdist(PM); % find limiting distribution
mc=dtmc(PM); % generate markov chain object

% initialize figure 
figure(iter)
u = tiledlayout(1,k);

% contribution from mother vs father 
v=0.3;
l=0.7;

for ii=1:k
    
phi=phis(ii);

% choose initial population size (and maximal)
P=100;

% choose initial population sex ratio
alpha=0.5;


% Generate a markov chain with probability transition matrix PM and as
% many states as the max number of generations considered
X = simulate(mc,gens);

% Initial number of females and males 
F=P*(1-alpha);
S=P*(alpha);

% Initial array of maternal condition 
cms=normrnd(0.5,0.2,1,F);

% Initital array of fathers condition
cfs=normrnd(0.5,0.2,1,S);

% scaled fitness of fathers
fitf=f_male_g(cfs);
fitf=fitf/sum(fitf);


% initialize vector to store results : generation number, number of
% females, number of sons, population size,

results=zeros(8,gens);



% generate array of mothers       
mothers=[cms; arrayfun(@(x) invs(x),cms);zeros(1,F);zeros(1,F);zeros(1,F)];


n=1;

results(1,n)=n;
results(2,n)=F;
results(6,n)=mean(cms);
results(7,n)=var(cms);
results(8,n)=S/(S+F);


while  F>0 && S>0 && n<(gens)
    
        results(1,n)=n;
        results(2,n)=F;
        results(6,n)=mean(cms);
        results(7,n)=var(cms);
        results(8,n)=S/(S+F);
 
    % initialize array of sons and daughters, where we'll store their
    % condition, fitness and sex (0=male, 1=female)
    
    sons=zeros(3,100*P);
    daughters=zeros(3,100*P);
    
    % find adult condition depending on future state of the markov chain
    if X(n+1)==1 %if good
        cf = randsrc(1,1,[cfs,fitf]); 
        ca= @(cm,i) max(0,co(cm,i)*l+ cf*v + normrnd( phi, kappa));
    elseif X(n+1)==2 %if neutral
        cf = randsrc(1,1,[cfs,fitf]); 
        ca= @(cm,i) max(0,co(cm,i)*l+ cf*v  + normrnd( 0, kappa));
    elseif X(n+1)==3 %if bad
        cf = randsrc(1,1,[cfs,fitf]); 
        ca= @(cm,i) max(0,co(cm,i)*l+ cf*v  + normrnd( -phi, kappa));
    end
    
    % For every mother, generate brood randomly (probability of male 0.5)

    for i=1:F  
        mothers(3,i)=binornd(2,0.5); % number of sons
    end
    
    % Find newborns sex ratio
    alpha=sum(mothers(3,:))/(length(mothers(3,:))*2);
    
    % For every mother, find her optimal investment strategy and generate
    % brood
    for i=1:F
        mothers(4,i)=2-mothers(3,i); % number of daughters
       
        [i1,i2]=optimal_inv(mothers(1,i),alpha) ;

        mothers(5,i)=i1; % store investment to sons in mized brood
        I=invs(mothers(1,i)); % find max investment
    
        
        if mothers(3,i)==2 % if all-male
            sons(1,i)=ca(mothers(1,i),I/2); % son 1 condition
            sons(2,i)=f_male_x(sons(2,i)); % son 1 fitness
            sons(3,i)=0; % son sex
            sons(1,100*P-i)=ca(mothers(1,i),I/2); % son 2 condition
            sons(2,100*P-i)=f_male_x(sons(2,i)); % son 2 fitness
            sons(3,100*P-i)=0; % son sex
            
        elseif mothers(3,i)==1 % if mixed-brood
            sons(1,i)=ca(mothers(1,i), mothers(5,i)); % son 1 condition
            sons(2,i)=f_male_x(sons(2,i)); % son 1 fitness
            sons(3,i)=0;
            
            daughters(1,i)=ca(mothers(1,i),I- mothers(5,i)); % daughter 1 condition
            daughters(2,i)=f_female_x(daughters(2,i)); % daughter 1 fitness
            daughters(3,i)=1;
            
        elseif mothers(3,i)==0
            daughters(1,i)=ca(mothers(1,i), I/2); % daughter 1condition
            daughters(2,i)=f_female_x(daughters(2,i)); %% daughter 1 fitness
            daughters(3,i)=1;
            daughters(1,100*P-i)=ca(mothers(1,i),I/2); % daughter 2 condition
            daughters(2,100*P-i)=f_female_x(daughters(2,i)); % daughter 2 fitness
            daughters(3,100*P-i)=1;
        
        end
        
        
         
    end
    
      % Remove zero entries
        daughters(:,all(daughters==0))=[];
        sons(:,all(sons==0))=[];
  
       % Find and save number of sons and daughters born
       S=length(sons(1,:)); % number of sons born
       D=length(daughters(1,:)); % number of daughters born
       results(3,n)=S; 
       results(4,n)=D; 
       
       % Find and save total number of newborns
       Pn=S+D; % Newborns population size
       results(5,n)=Pn;

      
        
        % Generate offspring array storing their condition, fitness and sex
        offspring=[sons(1,:),daughters(1,:);sons(2,:),daughters(2,:);sons(3,:),daughters(3,:)];

        % Order offspring by descending fitness
        [temp, order] = sort(offspring(2,:),'descend');
        offspring = offspring(:,order);
   
        
        % If number of offspring is larger than max population size, select fittest P
        % regardless of sex
        if Pn>P   
            newgen=offspring(1:3,1:P);
            else
            newgen=offspring;
        end
      
        % Separate males and females
        is_female = newgen(3,:)== 1;
        newfem=newgen( :, is_female );
        is_male = newgen(3,:)== 0;
        newmale=newgen( :, is_male );

        % New generation
        n=n+1;
        
        % Sons become fathers
        cfs=newmale(1,:);
        S=length(newmale(1,:));
        
        % Daugthers become mothers
        cms=newfem(1,:)+0.14;
        F=length(newfem(1,:));
        mothers=[cms;arrayfun(@(x) invs(x),cms);zeros(1,F);zeros(1,F)];

    
        
    
end




% clean & save results and plot!
fi=max(results(1,:));
nexttile
plot(results(1,1:fi),results(6,1:fi),'b-',results(1,1:fi),results(7,1:fi),'k-',results(1,1:fi),results(8,1:fi),'LineWidth',1.5)
legend('Mean c_m','Variance c_m','Sex Ratio')
xlabel('Generations')
phi=phi+iter;
title(sprintf('results%d.mat', phi))
xlim([1 max(results(1,:))+0.01])
ylim([0 1])


end
iter=iter+1;

end
