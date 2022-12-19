%Two period example

clear all
close all
clc

%Colors
color_erule=(1/255)*[0,153,255];
color_mrule=(1/255)*[0,120,0];
color_no=(1/255)*[255,0,0];
color_no_R=[0.635294117647059 0.07843137254902 0.184313725490196];
color_no_P=[0.243137254901961 0.368627450980392 0.074509803921569];
color_ci=(1/255)*[0,102,204];
color_cj=(1/255)*[255,102,102];
color_g=(1/255)*[0,153,76];


sameX = 0.1;

bet = 0.96;
beta=bet;Y=1.3;
xbar_c=0.1;
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
     
    ui2(i)=log(ci2(i))+log(g2(i));
    uj2(i)=log(cj2(i))+log(g2(i));
    
 
   
end
%plot allocations

X1=tauvec;
YMatrix1=[g2',ci2',cj2'];
%create_oneshot_2020(X1, YMatrix1)

%plot policies
YMatrix1=[ent2',ent2',tau2'];
%createTaxEnt2(X1, YMatrix1)


iui2= griddedInterpolant(s_vec, ui2, 'spline');
iuj2= griddedInterpolant(s_vec, uj2, 'spline');

%%%%%%%%%%%%%%%
%first period
%%%%%%%%%%%%%%
pvec=0:0.01:1;

Welf1=zeros(ns,ns);

lambda_L=2*xbar_c/Y;
lambda_H=(Y-2*xbar_c)/Y;
lambda=0.5;%lambda_H;

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
    
     
%Lifetime utility: P in power in the first period
%Note: P=i
uP(jp)=log(ci_crule(jp))+log(g_crule(jp))+beta*p*iui2(cj_crule(jp))+...
    beta*(1-p)*iuj2(ci_crule(jp));
uR(jp)=log(cj_crule(jp))+log(g_crule(jp))+beta*(1-p)*iui2(cj_crule(jp))+...
    beta*p*iuj2(ci_crule(jp));

end
%createfigure_grule(pvec, [ci_crule' ctilH_pol' ctilH_PL'], [g_crule' g_pol' gtil_PL'])
%createTaxEnt1(pvec, [ci_crule' ctilH_pol' ], [g_crule' g_pol' gtil_PL'])


%Pareto Frontier
lambdal=100;
lambdabar=linspace(0,1,lambdal);
pbarl=length(pvec);
pbar=pvec;

    alpha=0.5;
    for iL=1:lambdal
        lambda=lambdabar(iL);
        g_po(iL)=(1-alpha)*Y;
        cR_po(iL)=(1-lambda)*alpha*Y;
        cP_po(iL)=lambda*alpha*Y;
        
        uR_po(iL)=(1+beta)*(log(cR_po(iL))+log(g_po(iL)));
        uP_po(iL)=(1+beta)*(log(cP_po(iL))+log(g_po(iL)));
        
        if uR_po(iL)==-Inf
            uR_po(iL)=NaN;
        end
                
        if uP_po(iL)==-Inf
            uP_po(iL)=NaN;
        end
    end %end loop for lambda
    
    %Discretion
    cin_no=alpha*(Y-xbar);
    cout_no=xbar;
    g_no=(1-alpha)*(Y-xbar);
    
    for iP_plot=1:pbarl
        p=pbar(iP_plot);
        uR_R_no(iP_plot)=log(cin_no)+log(g_no)+...
            beta*p*(log(cin_no)+log(g_no))+...
            beta*(1-p)*(log(cout_no)+log(g_no));
        uP_R_no(iP_plot)=log(cout_no)+log(g_no)+...
            beta*p*(log(cout_no)+log(g_no))+...
            beta*(1-p)*(log(cin_no)+log(g_no));
        
       
        uR_P_no(iP_plot)=log(cout_no)+log(g_no)+...
            beta*p*(log(cout_no)+log(g_no))+...
            beta*(1-p)*(log(cin_no)+log(g_no));
        uP_P_no(iP_plot)=log(cin_no)+log(g_no)+...
            beta*p*(log(cin_no)+log(g_no))+...
            beta*(1-p)*(log(cout_no)+log(g_no));
        
        
        uR_veil(iP_plot)=(uR_R_no(iP_plot)+uR_P_no(iP_plot))/2;
        uP_veil(iP_plot)=(uP_R_no(iP_plot)+uP_P_no(iP_plot))/2;

    end %end p loop


%45degree
u45_R=uR_po;
    for iL=1:lambdal
if iL>=50
    u45_R(iL)=uR_po(iL);
else
 u45_R(iL)=NaN;
end
    end

uR_po_plot=uR_po(:);
uP_po_plot=uP_po(:);
uR_R_no_plot=uR_R_no(:);
uP_R_no_plot=uP_R_no(:);
uR_P_no_plot=uR_P_no(:);
uP_P_no_plot=uP_P_no(:);
uR_veil_plot=uR_veil(:);
uP_veil_plot=uP_veil(:);


Gain_P=uP_P_no_plot


%Plot PF
h=figure(301)
plot(uP_po_plot,uR_po_plot,'k','LineWidth', 3)
hold on
plot(u45_R,u45_R,'k','LineWidth', 3)
hold on
plot(uP_P_no_plot,uR_P_no_plot,'LineStyle',':','Color',color_no_R,'LineWidth', 3)
%hold on
%plot(uP_P_no_plot,uR_P_no_plot,'LineStyle',':','Color',color_no_P,'LineWidth', 3)
hold on
% plot(uP,uR,'LineStyle','-','Color','b','LineWidth', 3,...
% 'Marker','d','MarkerSize',11,...
%     'MarkerEdgeColor','b',...
%     'MarkerFaceColor','b')
plot(uP,uR,'LineStyle','-.','Color','b','LineWidth',3)
hold on
%title(['Pareto Frontier'])
xlabel('Welfare P')
set(gca,'ytick',[])
set(gca,'xtick',[])
ylabel('Welfare R')
xlim([min(uP_po_plot)+3,max(uP_po_plot)+0.5])
ylim([min(uR_po_plot)+3,max(uR_po_plot)+0.5])
annotation('textarrow',[0.588571428571429 0.613571428571429],...
    [0.517142857142857 0.545714285714286],'String','q=0','Interpreter','latex',...
    'FontSize',12);
annotation('arrow',[0.573571428571429 0.602142857142856],...
    [0.518095238095238 0.594285714285715]);
annotation('textbox',[0.15 0.83 0.3 0.1],'String','Pareto Frontier',...
    'Interpreter','latex',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation('textbox',[0.15 0.15 0.2 0.2],'String','Equitable Line',...
    'Interpreter','latex',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation('textbox',...
    [0.604999999999997 0.586666666666667 0.162857142857145 0.0485714285714332],...
    'Color',[0.635294117647059 0.07843137254902 0.184313725490196],...
    'String','Unconstrained',...
    'Interpreter','latex',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation('textbox',...
    [0.667857142857142 0.370476190476191 0.2 0.2],...
    'Color', 'b',...
    'String','Budget Rules',...
    'Interpreter','latex',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none');
box off
hold off
fname = 'pf_twoperiod';
savefig(sprintf('%s.fig', fname));
%saveas(gca,'test.pdf');
%system('pdfcrop test.pdf test.pdf');
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf, sprintf('%s.pdf', fname), '-dpdf');



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
