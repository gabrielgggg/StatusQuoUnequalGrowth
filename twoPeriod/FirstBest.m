  clear; clc;
clear all

sameX = 0.1;


bet = 0.96;
beta=bet;

xbar_c=0.1%;0.1;
xbar_g=xbar_c;

Y=1.3;

Ytil=1.3;
yP=xbar_c;
yR=Y-yP;

lambda_L=2*xbar_c/Ytil;
lambda_H=(1-lambda_L);
n=50;

%lambda_vec=lambda_L:(lambda_H-lambda_L)/n:lambda_H;
lambda_vec=0:0.05:1;
nlam=length(lambda_vec);
%Y+2*xbar_c+xbar_g;
pvec=lambda_vec;

for i=1:nlam
    lambda=lambda_vec(i);
    
 
        gtil_PL(i)=Ytil/2;
        ctilH_PL(i)=lambda*Ytil/2;
        ctilL_PL(i)=(1-lambda)*Ytil/2;
    
   %decentralize


        etil_PL1(i)=ctilL_PL(i)-yP;
        tautil_PL1(i)=yR-ctilH_PL(i);
        share1(i)=etil_PL1(i)/(etil_PL1(i)+ gtil_PL(i))
        
    check(i)=Ytil-gtil_PL(i)-ctilH_PL(i)-ctilL_PL(i);
    
    V_PL(i)=1/(1-beta)*(lambda*log(ctilH_PL(i))+(1-lambda)*log(ctilL_PL(i))+log(gtil_PL(i)));
    VH_PL(i)=1/(1-beta)*(log(ctilH_PL(i))+log(gtil_PL(i)));
    VL_PL(i)=1/(1-beta)*(log(ctilL_PL(i))+log(gtil_PL(i)));
end

        
for i=1:length(pvec)

    p=pvec(i);

    VH_pol(i)=1/(1-beta*p)*(2*log((Ytil-xbar_c)/2)+beta*(1-p)*(1-beta*p+2*(1-p)*beta)/(1-2*beta*p-beta^2*(1-2*p))*(log((Ytil-xbar_c)/2))...
        +beta*(1-p)*(1-beta*p)/(1-2*beta*p-beta^2*(1-2*p))*log(xbar_c));
    
    VH2_pol(i)=2*log((Ytil-xbar_c)/2)/(1-beta*p)+ beta*(1-p)*(log(xbar_c)+log((Ytil-xbar_c)/2)...
               +(1-p)*beta*2*log((Ytil-xbar_c)/2)...
               /(1-beta*p))*1/((1-beta*p)^2-(1-p)^2*beta^2);
               
    VL_pol(i)=(log(xbar_c)+log((Ytil-xbar_c)/2)+beta*(1-p)*VH_pol(i))/(1-beta*p);
    WH_pol(i)= VL_pol(i);
    
    %z(i)=2*log((Ytil-xbar_c)/2)+beta*(p*VH_pol(i)+(1-p)*WH_pol(i))-VH_pol(i);
    ctilL_pol(i)=xbar_c;
      ctilH_pol(i)=(Ytil-xbar_c)/2;
      g_pol(i)=(Ytil-xbar_c)/2;
end

        gCE=xbar_c;
        cRCE=yR-gCE;
        cPCE=yP;
VR_CE=(log(cRCE)+log(gCE));
    VP_CE=(log(cPCE)+log(gCE));
    
n1=1;
n2=nlam;
n2p=11;
n2p=nlam
n1p=n1;
%createfigureFB_2020_bis(VL_PL(n1p:n2p)*(1-beta), VH_PL(n1p:n2p)*(1-beta), VL_pol*(1-beta), VH_pol*(1-beta))

[VR_CE VP_CE]

VL_PL1=VL_PL;
VH_PL1=VH_PL;
VL_pol1=VL_pol;
VH_pol1=VH_pol;
VR_CE1=VR_CE ;
VP_CE1=VP_CE;


%growth

Y=1.3*1.6;

Ytil=Y;
yP=xbar_c;
yR=Y-yP;

lambda_L=2*xbar_c/Ytil;
lambda_H=(1-lambda_L);
n=50;

%lambda_vec=lambda_L:(lambda_H-lambda_L)/n:lambda_H;
lambda_vec=0:0.05:1;
nlam=length(lambda_vec);
%Y+2*xbar_c+xbar_g;
pvec=lambda_vec;

for i=1:nlam
    lambda=lambda_vec(i);
    
 
        gtil_PL(i)=Ytil/2;
        ctilH_PL(i)=lambda*Ytil/2;
       ctilL_PL(i) =(1-lambda)*Ytil/2;
    %decentralize


        etil_PL(i)=ctilL_PL(i)-yP;
        tautil_PL(i)=yR-ctilH_PL(i);

         share(i)=etil_PL(i)/(etil_PL(i)+ gtil_PL(i))
    check(i)=Ytil-gtil_PL(i)-ctilH_PL(i)-ctilL_PL(i);
    
    V_PL(i)=1/(1-beta)*(lambda*log(ctilH_PL(i))+(1-lambda)*log(ctilL_PL(i))+log(gtil_PL(i)));
    VH_PL(i)=1/(1-beta)*(log(ctilH_PL(i))+log(gtil_PL(i)));
    VL_PL(i)=1/(1-beta)*(log(ctilL_PL(i))+log(gtil_PL(i)));
end
        
for i=1:length(pvec)

    p=pvec(i);

    VH_pol(i)=1/(1-beta*p)*(2*log((Ytil-xbar_c)/2)+beta*(1-p)*(1-beta*p+2*(1-p)*beta)/(1-2*beta*p-beta^2*(1-2*p))*(log((Ytil-xbar_c)/2))...
        +beta*(1-p)*(1-beta*p)/(1-2*beta*p-beta^2*(1-2*p))*log(xbar_c));
    
    VH2_pol(i)=2*log((Ytil-xbar_c)/2)/(1-beta*p)+ beta*(1-p)*(log(xbar_c)+log((Ytil-xbar_c)/2)...
               +(1-p)*beta*2*log((Ytil-xbar_c)/2)...
               /(1-beta*p))*1/((1-beta*p)^2-(1-p)^2*beta^2);
               
    VL_pol(i)=(log(xbar_c)+log((Ytil-xbar_c)/2)+beta*(1-p)*VH_pol(i))/(1-beta*p);
    WH_pol(i)= VL_pol(i);
    
    %z(i)=2*log((Ytil-xbar_c)/2)+beta*(p*VH_pol(i)+(1-p)*WH_pol(i))-VH_pol(i);
    ctilL_pol(i)=xbar_c;
      ctilH_pol(i)=(Ytil-xbar_c)/2;
      g_pol(i)=(Ytil-xbar_c)/2;
end

        gCE=xbar_c;
        cRCE=yR-gCE;
        cPCE=yP;
VR_CE=(log(cRCE)+log(gCE));
    VP_CE=(log(cPCE)+log(gCE));
    
n1=1;
n2=nlam;
n2p=11;
n2p=nlam
n1p=n1;
%createfigureFB_2020_bis(VL_PL(n1p:n2p)*(1-beta), VH_PL(n1p:n2p)*(1-beta), VL_pol*(1-beta), VH_pol*(1-beta))

[VR_CE VP_CE]
plot(lambda_vec,share1,'r',lambda_vec,share,'b')

return




figure(3)
plot(pvec(n1p:n2p),(ctilH_pol(n1p:n2p)),'r', lambda_vec(n1:n2),(ctilH_PL(n1:n2)),'b',lambda_vec(n1:n2), gtil_PL(n1:n2),'k',pvec(n1p:n2p), g_pol(n1p:n2p),'k:')

return
%%%%%%%%%%%%%%%
% Static Model
%%%%%%%%%%%%%%
svec=lambda_vec;
for i=1:nlam
    s=lambda_vec(i);
    
    if s<=(Ytil-xbar_c)/2
        
        gtil_b(i)=(Ytil-xbar_c)/2;
        ctilH_b(i)=xbar_c;
        ctilL_b(i)=(Ytil-xbar_c)/2;
        
    elseif (s>(Ytil-xbar_c)/2 && s<=(Ytil)/2)
        
        gtil_b(i)=s;
        ctilL_b(i)=Ytil-s-xbar_c;
        ctilH_b(i)=xbar_c;
        
    else
         gtil_b(i)=(Ytil)/2;
         
        ctilL_b(i)=(Ytil^2-4*s*xbar_c)/(2*Ytil);
         ctilH_b(i)=(2*s*xbar_c)/(Ytil);
    end
    
 VH_b(i)=(log(ctilH_b(i))+log(gtil_b(i)));   
  VL_b(i)=(log(ctilL_b(i))+log(gtil_b(i)));     
  
end
figure(4)

plot(lambda_vec,gtil_b,'r', lambda_vec,ctilL_b,'b',lambda_vec,ctilH_b,'k')

VH_pl1=VH_PL*(1-beta);
VL_pl1=VL_PL*(1-beta);

plot(VL_PL, VH_PL, VL_pl1, VH_pl1)

