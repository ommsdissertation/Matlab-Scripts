
function [x,i1,i2,i3,fval]=optimal(cm,Func) 



if Func==1
    f_m = @(x) f_male_l(x);
    f_f = @(x) f_female_l(x);
elseif Func==2    
    f_m = @(x) f_male_x(x);
    f_f = @(x) f_female_x(x);
elseif Func==3
    f_m = @(x) f_male_g(x);
    f_f = @(x) f_female_g(x);
end

I=invs(cm);
c= @(i) co(cm,i);

b1 = @(x) -(f_m(c(x(1)))+f_m(c(x(2))));
lb = [0,0];
ub = [];
A = [1,1];
b = I;
Aeq = [];
beq = [];
x0=[0.5,0.5];
i1 = fmincon(b1,x0,A,b,Aeq,beq,lb,ub);

b2 = @(x) -(f_m(c(x(1)))+f_f(c(x(2))));
lb = [0,0];
ub = [];
A = [1,1];
b = I;
Aeq = [];
beq = [];
x0=[0.5,0.5];
i2 = fmincon(b2,x0,A,b,Aeq,beq,lb,ub);

b3 = @(x) -(f_f(c(x(1)))+f_f(c(x(2))));
lb = [0,0];
ub = [];
A = [1,1];
b = I;
Aeq = [];
beq = [];
x0=[0.5,0.5];
i3 = fmincon(b3,x0,A,b,Aeq,beq,lb,ub);


repval = @(x) -( x^2*(f_m(c(i1(1)))+f_m(c(i1(2)))) + (1-x)^2*(f_f(c(i3(1)))+f_f(c(i3(2)))) + 2*(1-x)*x*(f_m(c(i2(1)))+f_f(c(i2(2)))));

lb = 0;
ub = 1;
A = [];
b = [];
Aeq = [];
beq = [];
x0=[0.5];
[x,fval] = fmincon(repval,x0,A,b,Aeq,beq,lb,ub);

end
