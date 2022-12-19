%Two period example

clear all
clear all
clc


sameX = 0.1

bet = 0.96;
beta=bet;Y=1.3;
xbar_c=sameX;
yP=xbar_c;
yR=Y-yP;

s_vec=xbar_c:0.001:Y-2*xbar_c;
ns=length(s_vec);

xbar_g=xbar_c;
xbar=xbar_c;
Ytil=Y;

delta=yR-yP;
tauD=(delta+xbar_c)/2;
tauDP=yR-xbar;
tauRP=tauD;

s_L=(Y-xbar_c)/2;
s_H=Y/2;
entD=(delta-xbar)/2;

tauvec=yR-s_vec;
evec=s_vec-yP;

%%%%%%%%%%%%%%%
% c-rule
%%%%%%%%%%%%%%

%%%%%%%%%%%%%%
%second period
%%%%%%%%%%%%%%

% in terms of tau

% in terms of cR
for i=1:ns
    tau=tauvec(i);
    
    if tau<delta/2
        
        g2(i)=Y/2;
        ci2(i)=(Y^2-4*xbar_c*(yR-tau))/(2*Y);
        cj2(i)=Y-ci2(i)-g2(i);
        tau2(i)=yR-2*xbar_c*(yR-tau)/Y;
        ent2(i)=delta/2-2*xbar_c*(yR-tau)/Y;
        
        %tau2b(i)=yR-cj2(i);
        %ent2b(i)=ci2(i)-yP;
    elseif tau>=delta/2 && tau<tauD
        
        g2(i)=yR-tau;
        ci2(i)=Y-g2(i)-xbar_c;
        cj2(i)=Y-ci2(i)-g2(i);
        tau2(i)=yR-xbar;
         ent2(i)=tau-xbar;
         %tau2b(i)=yR-cj2(i);
        %ent2b(i)=ci2(i)-yP;
    else
        g2(i)=(Y-xbar_c)/2;
        ci2(i)=(Y-xbar_c)/2;
        cj2(i)=Y-ci2(i)-g2(i);
       tau2(i)=tauDP;
        ent2(i)=entD;
        %tau2b(i)=yR-cj2(i);
        %ent2b(i)=ci2(i)-yP;
    end
     
 uPP2(i)=log(g2(i))+ log(ci2(i)); %P second period utilit when P is in power
 uRP2(i)=log(g2(i))+ log(cj2(i)); %R second period utilit when P is in power
   
end


% in terms of e, R in power in t=2

for i=1:ns
    ebar=evec(i);
    
    if ebar>=delta/2
        
        tau2(i)=delta/2+2*xbar_c*(yP-ebar)/Y;
        ent2(i)=2*xbar_c*(yP-ebar)/Y-yP;
        g2(i)=Y/2;
    
    elseif ebar<delta/2 && ebar>=entD
        
         tau2(i)=xbar-ebar;
         ent2(i)=entD;
        g2(i)=yP+entD;
        
    else
       tau2(i)=tauD;
        ent2(i)=entD;
        g2(i)=(Y-xbar_c)/2;
    end
        cR2(i)=yR-tau2(i);
        cP2(i)=yP+ent2(i);
     
uPR2(i)=log(g2(i))+ log(cP2(i)); %P second period utilit when R is in power
uRR2(i)=log(g2(i))+ log(cR2(i)); %R second period utilit when R is in power
   
end

%plot allocations

X1=tauvec;
YMatrix1=[g2',ci2',cj2'];
%create_oneshot_2020(X1, YMatrix1)

%plot policies
%figure(2)
%YMatrix1=[ent2',ent2',tau2'];
%createTaxEnt2(X1, YMatrix1)


%%%%%%%%%%%%%%%
%first period
%%%%%%%%%%%%%%
%status quo
ebar=0;
taubar=xbar;
cR_sq=yR-taubar;
g1_sq=taubar-ebar;
u1R_sq=log(cR_sq)+log(g1_sq);
u2R_sq=2*log(s_L);

Welf1=zeros(ns,ns);

lambda_L=2*xbar_c/Y;
lambda_H=(Y-2*xbar_c)/Y;
lambda=0.5;%lambda_H;
p=0;

for ie=1:length(evec)
    e1=evec(ie);
    for itau=1:length(tauvec)
        tau1=tauvec(itau);
        
        cR1(ie,itau)=yR-tau1;
        cP1(ie,itau)=yP-e1;
        g1(ie,itau)=tau1-e1;
        
        u1P(ie,itau)=log(cP1(ie,itau))+log(g1(ie,itau)); %P first period utilit when P is in power
        u1R(ie,itau)=log(cR1(ie,itau))+log(g1(ie,itau)); %R first period utilit when P is in power
        
        V1PP(ie,itau)=u1P(ie,itau)+beta*(p* uPP2(itau)+(1-p)*uPR2(ie)); %P lifetime welfare when P is in power
        V1RP(ie,itau)=u1R(ie,itau)+beta*((1-p)* uRR2(itau)+p*uRP2(ie)); %R lifetime welfare when P is in power
        
        %welfare of R under status quo
        V1RP_sq(ie,itau)=u1R_sq+beta * u2R_sq;
    end
end

return


for jp=1:length(pvec);
        p=pvec(jp);
   
    ci_crule(jp)=(1+beta*(1-p))*(Y-xbar_c)/(2+beta*(1-p));
    cj_crule(jp)=xbar_c;
     g_crule(jp)=Y-ci_crule(jp)-cj_crule(jp);
     
         
       ctilL_pol(jp)=xbar_c;
      ctilH_pol(jp)=(Ytil-xbar_c)/2;
      g_pol(jp)=(Ytil-xbar_c)/2;
      
      if lambda<=lambda_L
        
        gtil_PL(jp)=(Ytil-xbar_c)/(2-lambda);
        ctilH_PL(jp)=xbar_c;
        ctilL_PL(jp)=(1-lambda)/(2-lambda)*(Ytil-xbar_c);
        
    elseif lambda>=lambda_H
        
        gtil_PL(jp)=(Ytil-xbar_c)/(1+lambda);
        ctilH_PL(jp)=(Ytil-xbar_c)*lambda/(1+lambda);
        ctilL_PL(jp)=xbar_c;
        
    else
        gtil_PL(jp)=Ytil/2;
        ctilH_PL(jp)=lambda*Ytil/2;
        ctilL_PL(jp)=(1-lambda)*Ytil/2;
      end
    
     
end
%createfigure_grule(pvec, [ci_crule' ctilH_pol' ctilH_PL'], [g_crule' g_pol' gtil_PL'])
%createTaxEnt1(pvec, [ci_crule' ctilH_pol' ], [g_crule' g_pol' gtil_PL'])

URD=log(xbar)+(1+2*beta)*log((Y-xbar)/2);
URR=log(xbar)+log(g_crule(1))+beta*(log(ci2(1))+log(g2(1)));
Gain=URR-URD

return

% 1st period in terms of cR
for i=1:ns
    s=s_vec(i);
    
    if s<s_L
        
        g2(i)=(Y-xbar_c)/2;
        ci2(i)=(Y-xbar_c)/2;
        cj2(i)=xbar_c;
        
    elseif s>=s_L && s<s_H
        
        g2(i)=s;
        ci2(i)=Y-s-xbar_c;
        cj2(i)=xbar_c;
        
        
    else
        g2(i)=Y/2;
        ci2(i)=(Y^2-4*xbar_c*s)/(2*Y);
        cj2(i)=2*xbar_c*s/Y;
    end
     
    V(i)=  log(ci2(i))  + log(g2(i));
    W(i)=  log(cj2(i))  + log(g2(i));
    
   
end
%plot(s_vec,g2,s_vec,ci2,s_vec,cj2 )
%create_oneshot_2020(X1, YMatrix1)

%%%%%%%%%%%%%%%
% g-rule
%%%%%%%%%%%%%%
% pvec=0:0.01:1;
% 
% pstar=1-2*(1+beta)*xbar_c/beta/Y;
% 
% for j=1:length(pvec);
%     p=pvec(j);
%     
%     
%     if p>=pstar
%     g1algebra(j)=(1+beta)*(Y-xbar_c)/(2+beta*(1+p));
%     
%     else
%         A=4*xbar_c*(2+beta);
%         B=Y^2*(2+beta*(1-p))+4*xbar_c*(1+beta)*(Y-xbar_c);
%         C=(Y-xbar_c)*Y^2*(1+beta*(1-p));
%          g1algebra(j)=(B-sqrt(B^2-4*A*C))/(2*A);
%     end
%     ci_grule(j)=Y-xbar_c-g1algebra(j);
%     cj_grule(j)=xbar_c;
%     g_grule(j)=g1algebra(j);
%     
%     
% end
% 
%   ctilL_pol=xbar_c;
%       ctilH_pol=(Ytil-xbar_c)/2;
%       g_pol=(Ytil-xbar_c)/2;


%createfigure_grule(pvec, [ci_grule' ctilH_pol' ctilH_PL'], [g_grule' g_pol' gtil_PL'])

% figure(1)
% ax1=subplot(1,2,1)
% plot(ax1, pvec, ci_grule, 'g',pvec,ctilH_pol,'k',pvec,ctilL_PL,'b')
% title('Proposer Consumption')
% 
% ax2=subplot(1,2,2)
% plot(ax2, pvec, g_grule, 'g', pvec,g_pol,'k',pvec,gtil_PL,'b')
% axis([ax1 ax2],[0 1 0.35 0.85])
% title('Public Good')
% 
% figure(2)
% ax1=subplot(1,2,1)
% plot(ax1,pvec, ci_crule, 'r',pvec,ctilH_pol,'k',pvec,ctilL_PL,'b')
% title('Proposer Consumption')
% 
% ax2=subplot(1,2,2)
% plot(ax2,pvec, g_crule, 'r', pvec,g_pol,'k',pvec,gtil_PL,'b')
% axis([ax1 ax2],[0 1 0.35 0.85])
% title('Public Good')
% 
%  createfigure_grule(pvec, [ci_crule' ctilH_pol' ctilL_PL'], [g_crule' g_pol' gtil_PL'])
% axis([ax1 ax2],[0 1 0.35 0.85])
