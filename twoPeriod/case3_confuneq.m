function [c,ceq] = case3_confuneq(x)
%global alpha xbar beta p iV2 iW2 ci0 cj0
global k beta p xbar iui2 iuj2 ci0 cj0
% Nonlinear inequality constraints
c(1)=-x(1)+xbar;
c(2)=-x(2)+xbar;
c(3)=-x(3)+xbar;
c(4)=log(cj0)+k*log(xbar)+beta*((1-p)*iui2(ci0)+p*iuj2(cj0))-...
     (log(x(2))+k*log(x(3))+beta*((1-p)*iui2(x(1))+p*iuj2(x(2))));
 ceq=[];
% Nonlinear equality constraints
% ceq= alpha*log(cj0)+(1-alpha)*log(xbar)+beta*((1-p)*iV2(ci0)+p*iW2(cj0))-...
%     (alpha*log(x(2))+(1-alpha)*log(x(3))+beta*((1-p)*iV2(x(1))+p*iW2(x(2))));
%ceq=[];

