
function [x,i2,fval]=optimal_uncertainty(cm,PM,ct,phi) 

% inputs: 
%       cm maternal condition
%       PM probability transition matrix
%       ct current state
%       phi mean value of change in condition
% 
% outputs: 
%       x optimal sex ratio
%       i2 array with investment to male and female in the mixed
%       fval optimized value of maternal fitness 

 
% find investment capability
I=invs(cm);

% define fitness functions (choose logistic)
f_m = @(x) f_male_g(x);
f_f = @(x) f_female_g(x);

% probability and expected adult condition of if next state being better
pb=PM(ct,1);
cb= @(i)co(cm,i) +phi;

% probability and condition of next state being normal
pn=PM(ct,2);
cn= @(i) co(cm,i);

%probability and condition of next state being worse
pw=PM(ct,3);
cw= @(i)co(cm,i)-phi; %max([0,co(cm,i)-phi]);

% for each possible state, find optimal investment in mixed brood

b2 = @(x) -(f_m(cb(x(1)))+f_f(cb(x(2))));
lb = [0,0];
ub = [];
A = [1,1];
b = I;
Aeq = [];
beq = [];
x0=[0.5,0.5];
i2b = fmincon(b2,x0,A,b,Aeq,beq,lb,ub);


b2 = @(x) -(f_m(cn(x(1)))+f_f(cn(x(2))));
lb = [0,0];
ub = [];
A = [1,1];
b = I;
Aeq = [];
beq = [];
x0=[0.5,0.5];
i2n = fmincon(b2,x0,A,b,Aeq,beq,lb,ub);



b2 = @(x) -(f_m(cw(x(1)))+f_f(cw(x(2))));
lb = [0,0];
ub = [];
A = [1,1];
b = I;
Aeq = [];
beq = [];
x0=[0.5,0.5];
i2w = fmincon(b2,x0,A,b,Aeq,beq,lb,ub);

% find optimal sex ratio 

repval = @(x) -( pb* ( x^2*2*f_m(cb(I/2)) + ((1-x)^2*2*f_f(cb(I/2)))+   2*(1-x)*x*(f_m(cb(i2b(1)))+f_f(cb(i2b(2)))) ) + pn* ( x^2*2*f_m(cn(I/2)) + ((1-x)^2*2*f_f(cn(I/2)))+   2*(1-x)*x*(f_m(cn(i2n(1)))+f_f(cn(i2n(2)))) ) + pw* ( x^2*2*f_m(cw(I/2)) + ((1-x)^2*2*f_f(cw(I/2)))+   2*(1-x)*x*(f_m(cw(i2w(1)))+f_f(cw(i2w(2)))) )   );

lb = 0;
ub = 1;
A = [];
b = [];
Aeq = [];
beq = [];
x0=[0.5];
[x,fval] = fmincon(repval,x0,A,b,Aeq,beq,lb,ub);
i2=pb*i2b + pw*i2w + pn*i2n;
i2=i2(1);
fval=-fval;
end
