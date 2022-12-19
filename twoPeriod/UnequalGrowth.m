%Two period example

clear all
clear all
clc


sameX = 0.1;

bet = 0.96;
beta=bet;
theta=1;

%initial case
growth=0;
Y0=1.3;

Y=Y0+growth;
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
s_L=(Y-xbar_c)/2;
s_H=Y/2;
entD=(delta-xbar)/2;

tauvec=yR-s_vec;

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
        
        tau2b(i)=yR-cj2(i);
        ent2b(i)=ci2(i)-yP;
    elseif tau>=delta/2 && tau<tauD
        
        g2(i)=yR-tau;
        ci2(i)=Y-g2(i)-xbar_c;
        cj2(i)=Y-ci2(i)-g2(i);
        tau2(i)=tauDP;
         ent2(i)=tau-xbar;
         tau2b(i)=yR-cj2(i);
        ent2b(i)=ci2(i)-yP;
    else
        g2(i)=(Y-xbar_c)/2;
        ci2(i)=(Y-xbar_c)/2;
        cj2(i)=Y-ci2(i)-g2(i);
       tau2(i)=tauDP;
        ent2(i)=entD;
        tau2b(i)=yR-cj2(i);
        ent2b(i)=ci2(i)-yP;
    end
     
 
   
end
%plot allocations

X1=tauvec;
YMatrix1_alloc=[g2',ci2',cj2'];
%create_oneshot_2020(X1, YMatrix1_alloc)

%plot policies

YMatrix1_policy=[ent2',ent2',tau2'];
%createTaxEnt2(X1, YMatrix1_policy)

%%%%%%%%%%%%%%%%%%%
%unequal growth case
%%%%%%%%%%%%%%%%%%%%5

%initial case
growth=0.2;

Y=Y0+growth;
xbar_c=sameX;
yP=xbar_c;
yR=Y-yP;

%s_vec=xbar_c:0.001:Y-2*xbar_c;
%ns=length(s_vec);

xbar_g=xbar_c;
xbar=xbar_c;
Ytil=Y;

delta=yR-yP;
tauD=(delta+xbar_c)/2;
tauDP=yR-xbar;
s_L=(Y-xbar_c)/2;
s_H=Y/2;
entD=(delta-xbar)/2;

tauvec=yR-s_vec;

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
        
        tau2b(i)=yR-cj2(i);
        ent2b(i)=ci2(i)-yP;
    elseif tau>=delta/2 && tau<tauD
        
        g2(i)=yR-tau;
        ci2(i)=Y-g2(i)-xbar_c;
        cj2(i)=Y-ci2(i)-g2(i);
        tau2(i)=tauDP;
         ent2(i)=tau-xbar;
         tau2b(i)=yR-cj2(i);
        ent2b(i)=ci2(i)-yP;
    else
        g2(i)=(Y-xbar_c)/2;
        ci2(i)=(Y-xbar_c)/2;
        cj2(i)=Y-ci2(i)-g2(i);
       tau2(i)=tauDP;
        ent2(i)=entD;
        tau2b(i)=yR-cj2(i);
        ent2b(i)=ci2(i)-yP;
    end
     
 
   
end
%plot allocations
figure(3)
X1=tauvec;
YMatrix1_alloc_g=[g2',ci2',cj2'];
%create_oneshot_2020(X1, YMatrix1_alloc_g)

%plot policies
figure(4)
YMatrix1_policy_g=[ent2',ent2',tau2'];
%createTaxEnt2(X1, YMatrix1_policy_g)

ent=YMatrix1_policy(:,1);
ent_g=YMatrix1_policy_g(:,1);

tax=YMatrix1_policy(:,3);
tax_g=YMatrix1_policy_g(:,3);


plot(X1,ent','r',X1,ent_g','b') 

plot(X1,tax','r',X1,tax_g','b') 