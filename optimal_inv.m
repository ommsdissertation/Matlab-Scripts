
function [i1,i2,fval]=optimal_inv(cm,a) 

% choose male-convex fitness functions
f_m = @(x) f_male_x(x)*(1-a);
f_f = @(x) f_female_x(x)*a;

% find maximal investment and define offspring condition

I=invs(cm);
c= @(i) co(cm,i);

% find 
b2 = @(x) -(f_m(c(x(1)))+f_f(c(x(2))));
lb = [0,0];
ub = [];
A = [1,1];
b = I;
Aeq = [];
beq = [];
x0=[0.5,0.5];
[y,fval] = fmincon(b2,x0,A,b,Aeq,beq,lb,ub);
i1=y(1);
i2=y(2);
end