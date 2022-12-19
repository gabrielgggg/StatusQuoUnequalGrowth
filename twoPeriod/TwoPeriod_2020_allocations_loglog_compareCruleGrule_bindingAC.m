%Two period example

clear all
close all
clc

global Y xbar beta k cj0 ci0 iui2 iuj2 p

%Utility
k=1;
%NOTE: for k>((Y-2*xbar)/2), we end up with the non intersting case such that
%g2=Y-2*xbar. We should assume that k<((Y-2*xbar)/2), as BCE do in their paper
%(page 2945). Codewise, we can assume any value for k.
fu=@(c,g) log(c)+k*log(g);

%Colors
color_crule=(1/255)*[0,153,255];
color_allrule=(1/255)*[0,204,0];
color_grule=(1/255)*[0,120,0];
color_no=(1 /255)*[255,0,0];
color_no_R=[0.635294117647059 0.07843137254902 0.184313725490196];
color_no_P=[0.243137254901961 0.368627450980392 0.074509803921569];
color_ci=(1/255)*[0,102,204];
color_cj=(1/255)*[255,102,102];
color_g=(1/255)*[0,153,76];


beta=0.96;
Y=1.3;
xbar=0.1;
yP=0.1;
yR=Y-yP;

kmax=((Y-2*xbar)/2);


ns=20;
%s_vec=linspace(xbar,Y-2*xbar,ns);
s_vec=linspace(xbar,Y-xbar,ns);

ps=1000;
pvec=0.5;
%linspace(0,1,ps);

lambdal=100;
lambdabar=linspace(0,1,lambdal);
pbarl=length(pvec);

delta=yR-yP;
tauD=(delta+xbar)/2;
tauDP=yR-xbar;
s_L=(Y-xbar)/2;
s_H=Y/2;
entD=(delta-xbar)/2;

tauvec=yR-s_vec;
%evec=yP-s_vec;
evec=s_vec-yP;

%%%%%%%%%%%%%%%
% c-rule
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
%second period
%%%%%%%%%%%%%%
% in terms of allocations
sq2_low=((k^k)*((Y-xbar)^k))/(((1+k)^k)*(xbar^(k-1)));
sq2_high=((k*Y)^k)/((xbar^(k-1))*((1+k)^k));
sq2_high2=(((Y/(1+k))-xbar)*(((k*Y)^k))/((xbar^k)*((1+k)^k)));
fun1=@cj2case4loglog;
x0=0.6;
opts = optimset('Diagnostics','off', 'Display','off');
for i=1:ns
    sq2=s_vec(i);
    
    if sq2<sq2_low
        g2(i)=k*(Y-xbar)/(1+k);
        cj2(i)=xbar;
        ci2(i)=Y-g2(i)-cj2(i);
    elseif sq2>=sq2_low&& sq2<sq2_high
        g2(i)=(xbar*(sq2^(1/k)))/(xbar^(1/k));
        cj2(i)=xbar;
        ci2(i)=Y-g2(i)-cj2(i);
    elseif sq2>=sq2_high&& sq2<sq2_high2
        g2(i)=k*Y/(1+k);
        cj2(i)=((xbar^k)*((1+k)^k)*sq2) /((k*Y)^k);
        ci2(i)=Y-g2(i)-cj2(i);
    else
        ci2(i)=xbar;
        cj2(i)=fsolve(fun1,x0,opts);
        g2(i)=Y-ci2(i)-cj2(i);
    end
    
    ui2(i)=fu(ci2(i),g2(i));
    uj2(i)=fu(cj2(i),g2(i));
    
end
iui2= griddedInterpolant(s_vec, ui2, 'spline');
iuj2= griddedInterpolant(s_vec, uj2, 'spline');
ig2 = griddedInterpolant(s_vec, g2, 'spline');
ici2 = griddedInterpolant(s_vec, ci2, 'spline');

options=optimset('Display','off');

%%%%%%%%%%%%%%%
%first period
%%%%%%%%%%%%%%
for jp=1:length(pvec)
    p=pvec(jp);
    
    for jc=1:ns
        cj0=s_vec(jc);
        ci0=max(xbar,min(0.6,Y-cj0-k*(Y-xbar)));
        
        %     if ((xbar^(k-1)))/(k^k)<((Y-xbar)/(1+k))^(k-1)
        %         g_crule(jp)=(k*(Y-xbar))/(1+k);
        %         cj_crule(jp)=xbar;
        %         ci_crule(jp)=Y-g_crule(jp)-cj_crule(jp);
        %     else
        %         g_crule(jp)=(k*(Y-xbar))/(1+k+beta*(1-p));
        %         cj_crule(jp)=xbar;
        %         ci_crule(jp)=Y-g_crule(jp)-cj_crule(jp);
        %     end
        
        
        %Kj2(jc)=fu(cj0,xbar);
        
        Aeq=[1,1,1];
        beq=Y;
        if jc==1
            x03=[ci0,xbar,Y-ci0-xbar];
        else
            x03=[ci1_case3(jc-1),cj1_case3(jc-1),g1_case3(jc-1)];
        end
        [x,fval3]=fmincon(@Fobj,x03,[],[],Aeq,beq,[],[],...
            @case3_confuneq,options);
        ci1_case3(jc)=x(1);
        cj1_case3(jc)=x(2);
        g1_case3(jc)=x(3);
        
        ci1_val=x(1);
        cj1_val=x(2);
        g1_val=x(3);
        
        cj_crule(jc,jp)=cj1_val;
        g_crule(jc,jp)=g1_val;
        ci_crule(jc,jp)=ci1_val;
        
        
        %Lifetime utility: P in power in the first period
        %Note: P=i
        uP(jc,jp)=fu(ci_crule(jc,jp),g_crule(jc,jp))+...
            beta*p*iui2(cj_crule(jc,jp))+...
            beta*(1-p)*iuj2(ci_crule(jc,jp));
        uR(jc,jp)=fu(cj_crule(jc,jp),g_crule(jc,jp))+...
            beta*p*iuj2(cj_crule(jc,jp))+...
            beta*(1-p)*iui2(ci_crule(jc,jp));
        
        g2_crule_val(jc,jp)=p*ig2(cj_crule(jc,jp))+(1-p)*ig2(ci_crule(jc,jp));
        ci2_crule_val(jc,jp)=p*ici2(cj_crule(jc,jp))+(1-p)*ici2(ci_crule(jc,jp));
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Discretion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cin_no=(Y-xbar)/(1+k);
cout_no=xbar;
g_no=(k*(Y-xbar))/(1+k);

for jp_plot=1:pbarl
    p=pvec(jp_plot);
    uR_R_no(jp_plot)=fu(cin_no,g_no)+...
        beta*p*(fu(cin_no,g_no))+...
        beta*(1-p)*(fu(cout_no,g_no));
    uP_R_no(jp_plot)=fu(cout_no,g_no)+...
        beta*p*(fu(cout_no,g_no))+...
        beta*(1-p)*(fu(cin_no,g_no));
    
    
    uR_P_no(jp_plot)=fu(cout_no,g_no)+...
        beta*p*(fu(cout_no,g_no))+...
        beta*(1-p)*(fu(cin_no,g_no));
    uP_P_no(jp_plot)=fu(cin_no,g_no)+...
        beta*p*(fu(cin_no,g_no))+...
        beta*(1-p)*(fu(cout_no,g_no));
    
    
    uR_veil(jp_plot)=(uR_R_no(jp_plot)+uR_P_no(jp_plot))/2;
    uP_veil(jp_plot)=(uP_R_no(jp_plot)+uP_P_no(jp_plot))/2;
    
end %end p loop

g_unconsvec=0.4839*ones(1,length(g_crule));

g_novec=g_no*ones(1,length(g_crule));
plot(tauvec,g_crule,tauvec,g_novec, tauvec,g_unconsvec)


plot(tauvec,yR-cj_crule,tauvec,ci_crule-yP)


return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pareto Frontier
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iL=1:lambdal
    lambda=lambdabar(iL);
    g_po(iL)=k*Y/(1+k);
    cR_po(iL)=(1-lambda)*(Y/(1+k));
    cP_po(iL)=lambda*(Y/(1+k));
    
    uR_po(iL)=(1+beta)*fu(cR_po(iL),g_po(iL));
    uP_po(iL)=(1+beta)*fu(cP_po(iL),g_po(iL));
    
    if uR_po(iL)==-Inf
        uR_po(iL)=NaN;
    end
    
    if uP_po(iL)==-Inf
        uP_po(iL)=NaN;
    end
end %end loop for lambda




%45degree
u45_R=uR_po;
for i=1:lambdal
    if uP_po(i)>uR_po(i)
        u45_R(i)=NaN;
    else
        u45_R(i)=uP_po(i);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOT G1 wrt p
X1=pvec';
YMatrix1=[g_crule',g_crule'];
%create_oneshot_2020(X1, YMatrix1)

% Create figure
figure1 = figure('Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.142450142450142 0.140794223826715 0.762549857549857 0.784205776173286]);
hold(axes1,'on');

% Create multjple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'LineWidth',3,'Parent',axes1);
set(plot1(1),'LineStyle','--',...
    'Color',[0, 0.6000, 1.0000]);
set(plot1(2),'LineStyle',':',...
    'Color',[0.270588235294118 0.450980392156863 0.03921568627451]);

% Create ylabel
ylabel('$g_{1}$','FontSize',14,'Interpreter','latex');
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))

% Create xlabel
xlabel('$q$','FontWeight','bold','FontSize',14,...
    'Interpreter','latex');

% Set the remaining axes properties
set(axes1,'FontSize',12);

% Create textbox
annotation(figure1,'textbox',...
    [0.28135714285714 0.688095238095243 0.213285714285714 0.135714285714286],...
    'Color',[0.270588235294118 0.450980392156863 0.03921568627451],...
    'String','g-rule',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.288499999999999 0.138095238095242 0.213285714285715 0.135714285714286],...
    'Color',[0 0.6 1],...
    'String','c-rule',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.2885 0.423333333333334 0.0932857142857143 0.0523809523809529],...
    'Color',[0.635294117647059 0.0784313725490196 0.184313725490196],...
    'String','Unc.',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% % Create line
% annotation(figure1,'line',[0.626934523809524 0.626934523809524],...
%     [0.14151461207536 0.946276516837266],'Color',[0.8 0.8 0.8],'LineWidth',2,...
%     'LineStyle',':');
%
% % Create textarrow
% annotation(figure1,'textarrow',[0.683928571428571 0.648586309523808],...
%     [0.378571428571429 0.305525886367009],...
%     'Color',[0.650980392156863 0.650980392156863 0.650980392156863],...
%     'String',{'q^{*}=0.64'},...
%     'LineStyle',':',...
%     'FontSize',12);

% box off
% hold off
% fname = 'letter_ql_g1';
% savefig(sprintf('%s.fig', fname));
% saveas(gca,'test.pdf');
% system('pdfcrop test.pdf test.pdf');
% set(figure1,'Units','Inches');
% pos = get(h,'Position');
% set(figure1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(gcf, sprintf('%s.pdf', fname), '-dpdf');


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %PLOT expectd G2 wrt p (optimal t=1 allocation)
%
% X1=pvec';
% YMatrix1=[g2_crule_val',g2_grule_val'];
% %create_oneshot_2020(X1, YMatrix1)
%
% % Create figure
% figure1 = figure('Color',[1 1 1]);
%
% % Create axes
% axes1 = axes('Parent',figure1,...
%     'Position',[0.142450142450142 0.140794223826715 0.762549857549857 0.784205776173286]);
% hold(axes1,'on');
%
% % Create multjple lines using matrix input to plot
% plot1 = plot(X1,YMatrix1,'LineWidth',3,'Parent',axes1);
% set(plot1(1),'LineStyle','--',...
%     'Color',[0, 0.6000, 1.0000]);
% set(plot1(2),'LineStyle',':',...
%     'Color',[0.270588235294118 0.450980392156863 0.03921568627451]);
%
% % Create ylabel
% ylabel('$g_{2}$','FontSize',14,'Interpreter','latex');
%
% % Create xlabel
% xlabel('$q$','FontWeight','bold','FontSize',14,...
%     'Interpreter','latex');
%
% % Set the remaining axes properties
% set(axes1,'FontSize',12);
% % Create textbox
% annotation(figure1,'textbox',...
%     [0.281357142857142 0.659523809523815 0.213285714285715 0.135714285714286],...
%     'Color',[0 0.6 1],...
%     'String','c-rule',...
%     'Interpreter','latex',...
%     'FontSize',14,...
%     'FitBoxToText','off',...
%     'EdgeColor','none');
%
% % Create textbox
% annotation(figure1,'textbox',...
%     [0.290285714285711 0.335714285714294 0.213285714285714 0.135714285714286],...
%     'Color',[0.270588235294118 0.450980392156863 0.03921568627451],...
%     'String','g-rule',...
%     'Interpreter','latex',...
%     'FontSize',14,...
%     'FitBoxToText','off',...
%     'EdgeColor','none');
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %PLOT G2 wrt sq
%
% X1=s_vec';
% YMatrix1=[g2',g2_grule'];
% %create_oneshot_2020(X1, YMatrix1)
%
% % Create figure
% figure1 = figure('Color',[1 1 1]);
%
% % Create axes
% axes1 = axes('Parent',figure1,...
%     'Position',[0.142450142450142 0.140794223826715 0.762549857549857 0.784205776173286]);
% hold(axes1,'on');
%
% % Create multjple lines using matrix input to plot
% plot1 = plot(X1,YMatrix1,'LineWidth',3,'Parent',axes1);
% set(plot1(1),'LineStyle','--',...
%     'Color',[0, 0.6000, 1.0000]);
% set(plot1(2),'LineStyle',':',...
%     'Color',[0.270588235294118 0.450980392156863 0.03921568627451]);
%
% % Create ylabel
% ylabel('$g_{2}$','FontSize',14,'Interpreter','latex');
%
% % Create xlabel
% xlabel('$\overline{c}$','FontWeight','bold','FontSize',14,...
%     'Interpreter','latex');
%
% % Set the remaining axes properties
% set(axes1,'FontSize',12);
% % Create textbox
% annotation(figure1,'textbox',...
%     [0.733142857142854 0.376190476190482 0.213285714285714 0.135714285714286],...
%     'Color',[0.270588235294118 0.450980392156863 0.03921568627451],...
%     'String','g-rule',...
%     'Interpreter','latex',...
%     'FontSize',14,...
%     'FitBoxToText','off',...
%     'EdgeColor','none');
%
% % Create textbox
% annotation(figure1,'textbox',...
%     [0.758142857142855 0.769047619047625 0.213285714285715 0.135714285714286],...
%     'Color',[0 0.6 1],...
%     'String','c-rule',...
%     'Interpreter','latex',...
%     'FontSize',14,...
%     'FitBoxToText','off',...
%     'EdgeColor','none');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOT Ci1 wrt p

X1=pvec';
YMatrix1=[ci_crule',ci_crule'];
%create_oneshot_2020(X1, YMatrix1)

% Create figure
figure1 = figure('Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.142450142450142 0.140794223826715 0.762549857549857 0.784205776173286]);
hold(axes1,'on');

% Create multjple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'LineWidth',3,'Parent',axes1);
set(plot1(1),'LineStyle','--',...
    'Color',[0, 0.6000, 1.0000]);
set(plot1(2),'LineStyle',':',...
    'Color',[0.270588235294118 0.450980392156863 0.03921568627451]);

% Create ylabel
ylabel('$c_{i,1}$','FontSize',14,'Interpreter','latex');
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))

% Create xlabel
xlabel('q','FontWeight','bold','FontSize',14,...
    'Interpreter','latex');

% Set the remaining axes properties
set(axes1,'FontSize',12);

% Create textbox
annotation(figure1,'textbox',...
    [0.402785714285711 0.790476190476198 0.213285714285715 0.135714285714286],...
    'Color',[0 0.6 1],...
    'String','c-rule',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.42778571428571 0.600000000000008 0.213285714285715 0.135714285714287],...
    'Color',[0.270588235294118 0.450980392156863 0.03921568627451],...
    'String','g-rule',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% % Create line
% annotation(figure1,'line',[0.632142857142857 0.632142857142857],...
%     [0.932333333333334 0.140476190476191],...
%     'Color',[0.650980392156863 0.650980392156863 0.650980392156863],...
%     'LineStyle',':');
%
% % Create textarrow
% annotation(figure1,'textarrow',[0.685714285714286 0.642857142857143],...
%     [0.668047619047619 0.773809523809524],...
%     'Color',[0.650980392156863 0.650980392156863 0.650980392156863],...
%     'String',{'q^{*}=0.64'},...
%     'FontSize',12);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %PLOT expected Ci2 wrt p (optimal t=1 allocation)
%
% X1=pvec';
% YMatrix1=[ci2_crule_val',ci2_grule_val'];
% %create_oneshot_2020(X1, YMatrix1)
%
% % Create figure
% figure1 = figure('Color',[1 1 1]);
%
% % Create axes
% axes1 = axes('Parent',figure1,...
%     'Position',[0.142450142450142 0.140794223826715 0.762549857549857 0.784205776173286]);
% hold(axes1,'on');
%
% % Create multjple lines using matrix input to plot
% plot1 = plot(X1,YMatrix1,'LineWidth',3,'Parent',axes1);
% set(plot1(1),'LineStyle','--',...
%     'Color',[0, 0.6000, 1.0000]);
% set(plot1(2),'LineStyle',':',...
%     'Color',[0.270588235294118 0.450980392156863 0.03921568627451]);
%
% % Create ylabel
% ylabel('$c_{i,2}$','FontSize',14,'Interpreter','latex');
%
% % Create xlabel
% xlabel('$q$','FontWeight','bold','FontSize',14,...
%     'Interpreter','latex');
%
% % Set the remaining axes properties
% set(axes1,'FontSize',12);
% % Create textbox
% annotation(figure1,'textbox',...
%     [0.299214285714284 0.338095238095245 0.213285714285715 0.135714285714287],...
%     'Color',[0 0.6 1],...
%     'String','c-rule',...
%     'Interpreter','latex',...
%     'FontSize',14,...
%     'FitBoxToText','off',...
%     'EdgeColor','none');
%
% % Create textbox
% annotation(figure1,'textbox',...
%     [0.274214285714282 0.678571428571437 0.213285714285714 0.135714285714287],...
%     'Color',[0.270588235294118 0.450980392156863 0.03921568627451],...
%     'String','g-rule',...
%     'Interpreter','latex',...
%     'FontSize',14,...
%     'FitBoxToText','off',...
%     'EdgeColor','none');
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %PLOT Ci2 wrt sq
%
% X1=s_vec';
% YMatrix1=[ci2',ci2_grule'];
% %create_oneshot_2020(X1, YMatrix1)
%
% % Create figure
% figure1 = figure('Color',[1 1 1]);
%
% % Create axes
% axes1 = axes('Parent',figure1,...
%     'Position',[0.142450142450142 0.140794223826715 0.762549857549857 0.784205776173286]);
% hold(axes1,'on');
%
% % Create multjple lines using matrix input to plot
% plot1 = plot(X1,YMatrix1,'LineWidth',3,'Parent',axes1);
% set(plot1(1),'LineStyle','--',...
%     'Color',[0, 0.6000, 1.0000]);
% set(plot1(2),'LineStyle',':',...
%     'Color',[0.270588235294118 0.450980392156863 0.03921568627451]);
%
% % Create ylabel
% ylabel('$c_{i,2}$','FontSize',14,'Interpreter','latex');
%
% % Create xlabel
% xlabel('$\overline{c}$','FontWeight','bold','FontSize',14,...
%     'Interpreter','latex');
%
% % Set the remaining axes properties
% set(axes1,'FontSize',12);
% % Create textbox
% annotation(figure1,'textbox',...
%     [0.634928571428569 0.285714285714292 0.213285714285715 0.135714285714286],...
%     'Color',[0 0.6 1],...
%     'String','c-rule',...
%     'Interpreter','latex',...
%     'FontSize',14,...
%     'FitBoxToText','off',...
%     'EdgeColor','none');
%
% % Create textbox
% annotation(figure1,'textbox',...
%     [0.650999999999997 0.580952380952389 0.213285714285714 0.135714285714287],...
%     'Color',[0.270588235294118 0.450980392156863 0.03921568627451],...
%     'String','g-rule',...
%     'Interpreter','latex',...
%     'FontSize',14,...
%     'FitBoxToText','off',...
%     'EdgeColor','none');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOT total U for first proposer P for all p´s (note: Pareto Frontier has a hole)
X1=pvec';
YMatrix1=[uP',uP_crule', uP_P_no'];
%create_oneshot_2020(X1, YMatrix1)

% Create figure
figure1 = figure('Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.142450142450142 0.140794223826715 0.762549857549857 0.784205776173286]);
hold(axes1,'on');

% Create multjple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'LineWidth',3,'Parent',axes1);
set(plot1(1),'LineStyle','--',...
    'Color',[0, 0.6000, 1.0000]);
set(plot1(2),'LineStyle',':',...
    'Color',[0.270588235294118 0.450980392156863 0.03921568627451]);

set(plot1(3),'LineStyle','-.',...
    'Color',[0.635294117647059 0.07843137254902 0.184313725490196]);

% Create ylabel
ylabel('$U_{P}$','FontSize',14,'Interpreter','latex');

% Create xlabel
xlabel('q','FontWeight','bold','FontSize',14,...
    'Interpreter','latex');

% Set the remaining axes properties
set(axes1,'FontSize',12);
% Create multjple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'LineWidth',3,'Parent',axes1);
set(plot1(1),'LineStyle','--','Color',[0 0.6 1]);
set(plot1(2),'LineStyle',':',...
    'Color',[0.270588235294118 0.450980392156863 0.03921568627451]);
set(plot1(3),'LineStyle','-.',...
    'Color',[0.635294117647059 0.07843137254902 0.184313725490196]);

% Create ylabel
ylabel('$U_{P}$','Interpreter','latex');

% Create xlabel
xlabel('q','FontWeight','bold','Interpreter','latex');

hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'FontSize',12);
% Create textbox
annotation(figure1,'textbox',...
    [0.190285714285713 0.702380952380958 0.213285714285715 0.135714285714286],...
    'Color',[0 0.6 1],...
    'String','c-rule',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.218857142857138 0.130952380952391 0.213285714285714 0.135714285714287],...
    'Color',[0.635294117647059 0.0784313725490196 0.184313725490196],...
    'String','discretionary-rule',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.211714285714281 0.295238095238105 0.213285714285714 0.135714285714287],...
    'Color',[0.270588235294118 0.450980392156863 0.03921568627451],...
    'String','g-rule',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOT total U for first respondent R for all p´s (note: Pareto Frontier has a hole)
X1=pvec';
YMatrix1=[uR',uR_grule',uR_R_no'];
%create_oneshot_2020(X1, YMatrix1)

% Create figure
figure1 = figure('Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.142450142450142 0.140794223826715 0.762549857549857 0.784205776173286]);
hold(axes1,'on');

% Create multjple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'LineWidth',3,'Parent',axes1);
set(plot1(1),'LineStyle','--',...
    'Color',[0, 0.6000, 1.0000]);
set(plot1(2),'LineStyle',':',...
    'Color',[0.270588235294118 0.450980392156863 0.03921568627451]);

set(plot1(3),'LineStyle','-.',...
    'Color',[0.635294117647059 0.07843137254902 0.184313725490196]);

% Create ylabel
ylabel('$U_{R}$','FontSize',14,'Interpreter','latex');

% Create xlabel
xlabel('q','FontWeight','bold','FontSize',14,...
    'Interpreter','latex');

% Set the remaining axes properties
set(axes1,'FontSize',12);
% Create textbox
annotation(figure1,'textbox',...
    [0.334928571428567 0.495238095238105 0.213285714285714 0.135714285714286],...
    'Color',[0.270588235294118 0.450980392156863 0.03921568627451],...
    'String','g-rule',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.370642857142856 0.25476190476191 0.213285714285715 0.135714285714287],...
    'Color',[0 0.6 1],...
    'String','c-rule',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.338499999999998 0.702380952380959 0.213285714285715 0.135714285714287],...
    'Color',[0.635294117647059 0.0784313725490196 0.184313725490196],...
    'String','discretionary-rule',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOT total dif U wrt all discretionary for proposer for all p´s (note: Pareto Frontier has a hole)
X1=pvec';
YMatrix1=[(uP-uP_P_no)',(uP_grule-uP_P_no)', zeros(1,100)'];
%create_oneshot_2020(X1, YMatrix1)

% Create figure
figure1 = figure('Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.142450142450142 0.140794223826715 0.762549857549857 0.784205776173286]);
hold(axes1,'on');

% Create multjple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'LineWidth',3,'Parent',axes1);
set(plot1(1),'LineStyle','--',...
    'Color',[0, 0.6000, 1.0000]);
set(plot1(2),'LineStyle',':',...
    'Color',[0.270588235294118 0.450980392156863 0.03921568627451]);

set(plot1(3),'LineStyle',':',...
    'Color',[216/255,216/255,216/255],...
    'LineWidth',0.5);

% Create ylabel
ylabel('$U_{P}-U^{no}_{P}$','FontSize',14,'Interpreter','latex');

% Create xlabel
xlabel('q','FontWeight','bold','FontSize',14,...
    'Interpreter','latex');

% Set the remaining axes properties
set(axes1,'FontSize',12);
% Create textbox
% Create textbox
annotation(figure1,'textbox',...
    [0.242071428571425 0.166666666666673 0.213285714285714 0.135714285714286],...
    'Color',[0.270588235294118 0.450980392156863 0.03921568627451],...
    'String','g-rule',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.247428571428569 0.733333333333339 0.213285714285715 0.135714285714286],...
    'Color',[0, 0.6000, 1.0000],...
    'String','c-rule',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOT total dif U wrt all discretionary for respondent R for all p´s (note: Pareto Frontier has a hole)
X1=pvec';
YMatrix1=[(uR-uR_R_no)',(uR_grule-uR_R_no)', zeros(1,100)'];
%create_oneshot_2020(X1, YMatrix1)

% Create figure
figure1 = figure('Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.142450142450142 0.140794223826715 0.762549857549857 0.784205776173286]);
hold(axes1,'on');

% Create multjple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'LineWidth',3,'Parent',axes1);
set(plot1(1),'LineStyle','--',...
    'Color',[0, 0.6000, 1.0000]);
set(plot1(2),'LineStyle',':',...
    'Color',[0.270588235294118 0.450980392156863 0.03921568627451]);

set(plot1(3),'LineStyle',':',...
    'Color',[216/255,216/255,216/255],...
    'LineWidth',0.5);

% Create ylabel
ylabel('$U_{R}-U^{no}_{R}$','FontSize',14,'Interpreter','latex');

% Create xlabel
xlabel('q','FontWeight','bold','FontSize',14,...
    'Interpreter','latex');

% Set the remaining axes properties
set(axes1,'FontSize',12);
% Create textbox
annotation(figure1,'textbox',...
    [0.306357142857139 0.719047619047627 0.213285714285714 0.135714285714286],...
    'Color',[0.270588235294118 0.450980392156863 0.03921568627451],...
    'String','g-rule',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.320642857142854 0.438095238095246 0.213285714285715 0.135714285714286],...
    'Color',[0 0.6 1],...
    'String','c-rule',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOT total dif U wrt all discretionary for respondent R for all p´s (note: Pareto Frontier has a hole)
X1=pvec';
Uavg_crule=(((uR-uR_R_no)+(uP-uP_P_no))/2);
Uavg_grule=(((uR_grule-uR_R_no)+(uP_grule-uP_P_no))/2);
YMatrix1=[Uavg_crule',Uavg_grule', zeros(1,100)'];
%create_oneshot_2020(X1, YMatrix1)

% Create figure
figure1 = figure('Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.142450142450142 0.140794223826715 0.762549857549857 0.784205776173286]);
hold(axes1,'on');

% Create multjple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'LineWidth',3,'Parent',axes1);
set(plot1(1),'LineStyle','--',...
    'Color',[0, 0.6000, 1.0000]);
set(plot1(2),'LineStyle',':',...
    'Color',[0.270588235294118 0.450980392156863 0.03921568627451]);

set(plot1(3),'LineStyle',':',...
    'Color',[216/255,216/255,216/255],...
    'LineWidth',0.5);

% Create ylabel
ylabel('$U-U^{no}$','FontSize',14,'Interpreter','latex');

% Create xlabel
xlabel('$q$','FontWeight','bold','FontSize',14,...
    'Interpreter','latex');

% Set the remaining axes properties
set(axes1,'FontSize',12);
% Create textbox
% Create textbox
annotation(figure1,'textbox',...
    [0.286714285714283 0.445238095238104 0.213285714285715 0.135714285714286],...
    'Color',[0 0.6 1],...
    'String','c-rule',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.279571428571425 0.661904761904771 0.213285714285714 0.135714285714286],...
    'Color',[0.270588235294118 0.450980392156863 0.03921568627451],...
    'String','g-rule',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOT PF for all p´s
%Plot PF
%h=figure(301)
figure('Color',[1 1 1])
plot(uP_po,uR_po,'k','LineWidth', 3)
hold on
plot(uP_P_no,uR_P_no,'LineStyle',':','Color',color_no_R,'LineWidth', 3)
plot(uP,uR,'LineStyle','-.','Color',color_crule,'LineWidth',3)
plot(uP_grule,uR_grule,'LineStyle','-.','Color',color_grule,'LineWidth',3)
xlabel('Welfare P','FontWeight','bold','FontSize',14,...
    'Interpreter','latex');
ylabel('Welfare R','FontWeight','bold','FontSize',14,...
    'Interpreter','latex');
% Create textbox
annotation('textbox',...
    [0.820642857142857 0.111904763394879 0.0928571406219687 0.0738095223194077],...
    'Color',[0.650980392156863 0.650980392156863 0.650980392156863],...
    'String',{'q=1'},...
    'FontSize',12,...
    'EdgeColor',[1 1 1]);

% Create textarrow
annotation('textarrow',[0.325 0.732142857142857],...
    [0.520428571428571 0.223809523809524],...
    'Color',[0.650980392156863 0.650980392156863 0.650980392156863],...
    'String',{'q=0'},...
    'FontSize',12);

% Create arrow
annotation('arrow',[0.325 0.507142857142857],...
    [0.523809523809524 0.50952380952381],...
    'Color',[0.650980392156863 0.650980392156863 0.650980392156863]);

% Create arrow
annotation('arrow',[0.326785714285714 0.5375],...
    [0.522809523809524 0.540476190476191],...
    'Color',[0.650980392156863 0.650980392156863 0.650980392156863]);
set(gca,'ytick',[])
set(gca,'xtick',[])
xlim([min(uP_po),max(uP_po)])
ylim([min(uR_po),max(uR_po)])
box off
hold off
%fname = 'pf_twoperiod';
%savefig(sprintf('%s.fig', fname));
%saveas(gca,'test.pdf');
%system('pdfcrop test.pdf test.pdf');
%set(h,'Units','Inches');
%pos = get(h,'Position');
%set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(gcf, sprintf('%s.pdf', fname), '-dpdf');


save('loglog_var.mat','g_crule','g_grule', 'ci_crule', 'ci_grule', 's_vec', 'pvec', ...
    'Uavg_crule', 'Uavg_grule')

