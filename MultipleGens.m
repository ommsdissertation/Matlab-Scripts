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


%% MultipleGens.m script


% Initial population size
P=500;
N=50;

% choose type of fitness functions (1: linear, 2: male convex, 3: logistic)
Func=3;

% Initial number of females and males 
F=P/2;
S=P/2;


% Initial array of maternal condition 
cm=normrnd(0.5,0.1,1,F);


mothers=[cm; arrayfun(@(x) invs(x),cm);zeros(1,F);zeros(1,F);zeros(1,F)];
n=0;

results=zeros(8,gens);
cms=zeros(1,gens);
alphas=zeros(1,gens);
D=1;
S=1;


while  F>0 && S>0 && n<N
    n=n+1;
    u=1;
    results(1,n)=n;
    results(2,n)=F;
    results(6,n)=mean(cm);
    results(7,n)=var(cm);
    results(8,n)=F/P;
    sex_ratio=zeros(1,F);
    inv_ratio=zeros(1,F);
    
    
    figure(1)
    nexttile
    h1=histogram(mothers(1,:),'Normalization', 'pdf');
    hold on
    sons=zeros(3,100*P);
    daughters=zeros(3,100*P);

    for i=1:F
        [s,i2,fval]=optimal(mothers(1,i),Func);
        mothers(2,i)=s;
        mothers(3,i)=binornd(2,s); % number of sons
        mothers(4,i)=2-mothers(3,i); % number of daughters
        mothers(5,i)=i2(1); % investment to sons
        sex_ratio(i)=s;
        inv_ratio(i)=i2(1)/invs(mothers(1,i));
        I=invs(mothers(1,i));
    
        if mothers(3,i)==2 % if all-male
            sons(1,i)=i;
            sons(2,i)=co(mothers(1,i),I/2); % son 1 condition
            sons(3,i)=f_male(sons(2,i)); % son 1 fitness
            sons(1,100*P-i)=100*P-i;
            sons(2,100*P-i)=co(mothers(1,i),I/2); % son 2 condition
            sons(3,100*P-i)=f_male(sons(2,i)); % son 2 fitness
        
        elseif mothers(3,i)==1 % if mixed-brood
            sons(1,i)=i;
            sons(2,i)=co(mothers(1,i), mothers(5,i)); % son 1 condition
            sons(3,i)=f_male(sons(2,i)); % son 1 fitness
            daughters(1,i)=i;
            daughters(2,i)=co(mothers(1,i),I- mothers(5,i)); % daughter 1 condition
            daughters(3,i)=f_female(daughters(2,i)); % daughter 1 fitness
         
        elseif mothers(3,i)==0
            daughters(1,i)=i;
            daughters(2,i)=co(mothers(1,i), I/2); % daughter 1condition
            daughters(3,i)=f_female(daughters(2,i)); %% daughter 1 fitness
            daughters(1,100*P-i)=100*P-i;
            daughters(2,100*P-i)=co(mothers(1,i),I/2); % daughter 2 condition
            daughters(3,100*P-i)=f_female(daughters(2,i)); % daughter 2 fitness
         end
        
    end
    
    
    S=sum(mothers(3,:)); % number of sons born
    D=sum(mothers(4,:)); % number of sons born
    Pn=S+D; % Newborns population size

    results(3,n)=S; %save
    results(4,n)=D; 
    results(5,n)=Pn;

% remove zero entries
        daughters(:,all(daughters==0))=[];
        sons(:,all(sons==0))=[];
% find sex of all offspring
        sex=[zeros(1,sum(mothers(3,:))),ones(1,sum(mothers(4,:)))];
        
%put all offspring together: their condition, fitness and sex

        offspring=[sons(2,:),daughters(2,:);sons(3,:),daughters(3,:);sex];

 % order offspring by increasing fitness
        [temp, order] = sort(offspring(2,:),'descend');
        offspring = offspring(:,order);
   
        
% if offspring generation is larger than maximal P, select fittest P
% regardless of sex
    if Pn>P   
    newgen=offspring(1:3,1:P);
    else
        newgen=offspring;
    end
    
% split surviving offspring into females and males
is_female = newgen(3,:)== 1;
newfem=newgen( :, is_female );
is_male = newgen(3,:)== 0;
newmale=newgen( :, is_male );

% find number of surviving female and male offspring
F=length(newfem(1,:));
S=length(newmale(1,:));

% make new mother condition array using daughters condition
cm=newfem(1,:);

% daugters become mothers of next generation!
mothers=[cm;arrayfun(@(x) invs(x),cm);zeros(1,F);zeros(1,F)]; 
    
end

% plot evolution of the population
figure(2)
results(:,all(results==0))=[]; % remove zero entries
plot(results(1,:),results(6,:),results(1,:),results(7,:),results(1,:),1-results(8,:),'b-','LineWidth',1.5)
legend('mean','variance','alpha')
xlabel('Generations')
title('Evolution of a population optimizing sex allocation, basic case')
  
