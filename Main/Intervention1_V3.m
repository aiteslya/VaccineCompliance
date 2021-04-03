%This is the script that plots time series for Intervention 1 which aims to
%maintain compliance of non-vaccinated individuals buy keeping compliance
%duration independent of vaccination coverage, in the scenarios of 
%circulation of the original and a more transmissible, B.1.1.7-like strain 
%with and without the vaccination with high compliance rise
% during the lockdown for COVID-19-Compliance-Vaccine model
% R0 for the old virus strain is 2.5 and for the new 3.75
% such that Re(0)=1.1 for the old strain and Re(0)=1.65 for the new strain
%prepare settings
%This script also produces figures that detail proportions of vaccinated
%and non-vaccinated individuals in the attack rate
clc;
clear variables;
close all;
format long;

%set the parameters

%define colors
l=5;
red = [0.01, 0.4, 0.76];
pink = [0.45,0.76,0.98];
medblue=(red+pink)/2;
novacc=[0.9,0.13,0.13];
yearcol=[0.17,0.17,0.17];
colors_p = [linspace(pink(1),red(1),l)', linspace(pink(2),red(2),l)', linspace(pink(3),red(3),l)'];
col_p=colors_p([2,4],:);

%contact rate of non-compliant before the lockdown at the start of the
%epidemic
chat=14.9;
%1/gamma duration of infectious period
gamma=1/7;
alpha=1/4;
R0arr=[2.5, 2.5*1.5];
%calculate epsilon
epsilonarr=R0arr*gamma/chat;

%set up initial data
%seroprevalence
SP=0.08;
%total population
N=1.7e7;
%percentage of compliant people
PerCompl=0.65;
N0=N*(1-PerCompl);
Nc0=N*PerCompl;
TotalInf=37706;
Incid=TotalInf*gamma;
mu0=1/30;
delta=4e-5;

exposed=Incid/alpha;
TotalRec=SP*N;
TotalS=N-TotalInf-TotalRec-exposed;
%setting up of initial data
S0=(1-PerCompl)*TotalS;
E0=(1-PerCompl)*exposed;
I0=(1-PerCompl)*TotalInf;
R0=(1-PerCompl)*TotalRec;
Sc0=PerCompl*TotalS;
Ec0=PerCompl*exposed;
Ic0=PerCompl*TotalInf;
Rc0=PerCompl*TotalRec;
V0=0;
SV0=0;
EV0=0;
IV0=0;
RV0=0;
TV0=SV0+EV0+IV0+RV0;
%compliance

infect0=E0+I0+R0+Ec0+Ic0+Rc0;
popOut=1e5;
infect0pop=popOut*infect0/N;

init=[S0,E0,I0,R0,Sc0,Ec0,Ic0,Rc0,V0,SV0,EV0,IV0,RV0,TV0,0,0,0];
%compliance

%define array of contact rates
r1num=60;
r1arr=linspace(0,1,r1num);
%define the array of r1
num=2e4;
carr=linspace(0,15,num);

%for each c in carr calculate r1
%set up equation Re=1

Res=nan(r1num,2);
counter=1;
for r1=r1arr
    beta=carr*epsilonarr(1);
    Re1=beta*S0./(gamma*(N0+r1*Nc0))+beta*r1*Sc0.*(mu0*(alpha+gamma+mu0)+alpha*gamma*r1)./(gamma*(alpha+mu0)*(gamma+mu0).*(N0+r1*Nc0))-1.1;
    nulleqn=Re1(1:num-1).*Re1(2:num);
    ind=find(nulleqn<0);
    if numel(ind)==1
        Res(counter,1)=carr(ind);
        Res(counter,2)=r1;
        counter=counter+1;
    elseif numel(ind)>0
        error('RcContact: more than one root');
    end
end


ACR0=5;

ACR=Res(:,1)*(1-PerCompl)+Res(:,2).*Res(:,1)*PerCompl-ACR0;
eqn=ACR(1:r1num-1).*ACR(2:r1num);
ind=find(eqn<0,1);
c=Res(ind,1);
r1=Res(ind,2);

%integration options
Atol=1e-11;
opts = odeset('RelTol',1e-12,'AbsTol',Atol);
%integrating time
T=800;

numout=2;
upsilonarr=[0,5.9e-4,4.9e-3];
k1=1;
k2=1;
omega=0.6;
%format of the legend
formatSpec = '%.2e';

r2=13.47/c;
frac=1/3;
mu1=0;

%legend for the vaccinated plot
counter=1;
numpoints=30;
popOut=1e5;
XpopOut=N/popOut;
%the following the counter for the type of variant
i1=1;
%container for vaccination cumulative infected at mark of year 1 and 2
CumContainer=zeros(4,2);
for epsilon=epsilonarr
    beta=epsilon*c;
    i2=1;
    for upsilon0=upsilonarr
        pars=[beta,r1,r2,delta,mu0,mu1,upsilon0,alpha,gamma,k1,k2,omega];
        [t,y]=ode45(@(t,y)COVIDVaccineRHS2(t,y,pars),[0,T], init,opts);
                
        infectious=(y(:,2)+y(:,3)+y(:,6)+y(:,7)+y(:,11)+y(:,12)+y(:,15)+y(:,16))/XpopOut;
        cont=c*(y(:,1)+y(:,2)+y(:,3)+y(:,4))/N +c*r1*(y(:,5)+y(:,6)+y(:,7)+y(:,8))/N+c*r2*(y(:,9)+y(:,10)+y(:,11)+y(:,12)+y(:,13)+y(:,15)+y(:,16)+y(:,17))/N;
        cum=(y(:,2)+y(:,3)+y(:,4)+y(:,6)+y(:,7)+y(:,8)+y(:,11)+y(:,12)+y(:,13)+y(:,15)+y(:,16)+y(:,17))/XpopOut-infect0pop;
        %new
        cumvacc=(y(:,15)+y(:,16)+y(:,17))/XpopOut;
        if upsilon0==0 %plot 
           
            %%%
            if i1==1
                figure(6);
                subplot(1,2,1);
                plot(t,infectious,'-.','LineWidth',4,'color',novacc);hold on;
                figure(7);
                subplot(1,2,1);
                plot(t,cont,'-.','LineWidth',4,'color',novacc);hold on;
            else
                figure(6);
                subplot(1,2,2);
                plot(t,infectious,'-.','LineWidth',4,'color',novacc);hold on;
                figure(7);
                subplot(1,2,2);
                plot(t,cont,'-.','LineWidth',4,'color',novacc);hold on;
            end
            %collect baseline cumulative infections at mark of 1 year and 2
            %years
            %find index of the first year
            ind1year=find(t>=365,1);
            baseCum(i1,1)=cum(ind1year);
            %find index of the second year
            ind2year=find(t>=730,1);
            baseCum(i1,2)=cum(ind2year);
        elseif upsilon0==5.9e-4
            if i1==1
                figure(6);
                subplot(1,2,1);
                plot(t,infectious,'LineWidth',4,'color',col_p(1,:));hold on;
                figure(7);
                subplot(1,2,1);
                plot(t,cont,'LineWidth',4,'color',col_p(1,:));hold on;
                
                figure(8);% v1, slow
                subplot(1,2,1);
                h8(2)=area(t,100*cum/popOut,'FaceColor',novacc,'edgecolor','none');hold on;
                h8(3)=area(t,100*cumvacc/popOut,'FaceColor',medblue,'edgecolor','none');hold on;
                xlim([0,800]);
            else
                figure(6);
                subplot(1,2,2);
                plot(t,infectious,'LineWidth',4,'color',col_p(1,:));hold on;
                figure(7);
                subplot(1,2,2);
                plot(t,cont,'LineWidth',4,'color',col_p(1,:));hold on;
                
                figure(8);% v2, slow
                subplot(1,2,2);
                h9(2)=area(t,100*cum/popOut,'FaceColor',novacc,'edgecolor','none');hold on;
                h9(3)=area(t,100*cumvacc/popOut,'FaceColor',medblue,'edgecolor','none');hold on;
                xlim([0,800]);
            end
        else
            if i1==1
                figure(6);
                subplot(1,2,1);
                plot(t,infectious,'LineWidth',4,'color',col_p(2,:));hold on;
                figure(7);
                subplot(1,2,1);
                plot(t,cont,'LineWidth',4,'color',col_p(2,:));hold on;
                
                figure(9);% v1, fast
                subplot(1,2,1);
                h8(4)=area(t,100*cum/popOut,'FaceColor',novacc,'edgecolor','none');hold on;
                h8(5)=area(t,100*cumvacc/popOut,'FaceColor',medblue,'edgecolor','none');hold on;
                xlim([0,800]);
            else
                figure(6);
                subplot(1,2,2);
                plot(t,infectious,'LineWidth',4,'color',col_p(2,:));hold on;
                figure(7);
                subplot(1,2,2);
                plot(t,cont,'LineWidth',4,'color',col_p(2,:));hold on;
                
                figure(9);% v2, fast
                subplot(1,2,2);
                h9(4)=area(t,100*cum/popOut,'FaceColor',novacc,'edgecolor','none');hold on;
                h9(5)=area(t,100*cumvacc/popOut,'FaceColor',medblue,'edgecolor','none');hold on;
                xlim([0,800]);
            end
        end
        if upsilon0>0
            ind1year=find(t>=365,1);
            CumContainer(2*i1-2+i2,1)=cum(ind1year);
            ind2year=find(t>=730,1);
            CumContainer(2*i1-2+i2,2)=cum(ind2year);
            i2=i2+1;
        end
        
    end
    i1=i1+1;
end
%figure settings

figure(6);
subplot(1,2,1);
xlim([0, 800]);

xline(365,'-.','color',yearcol,'LineWidth',2);
xline(730,'-.','color',yearcol,'LineWidth',2);

xlabel('Time (days)','interpreter','latex');
ylabel({'Number of infected';'individuals (1/100,000)'},'interpreter','latex');
title('Original variant');
legend('Baseline, no vaccination','Slow vaccination','Fast vaccination','location','southoutside');
set(gca,'FontSize',25);

%figure settings
subplot(1,2,2);
xlim([0, 800]);
xline(365,'-.','color',yearcol,'LineWidth',2);
xline(730,'-.','color',yearcol,'LineWidth',2);

xlabel('Time (days)','interpreter','latex');
title('B.1.1.7-like variant');
set(gca,'FontSize',25);

figure(7);
subplot(1,2,1);
xlim([0, 800]);
xline(365,'-.','color',yearcol,'LineWidth',2);
xline(730,'-.','color',yearcol,'LineWidth',2);

xlabel('Time (days)','interpreter','latex');
ylabel({'Average contact rate,';'contacts per day'},'interpreter','latex');
title('Original variant');
legend('Baseline, no vaccination','Slow vaccination','Fast vaccination','location','southoutside');
set(gca,'FontSize',25);

%figure settings
figure(7);
subplot(1,2,2);
xlim([0, 800]);
xline(365,'-.','color',yearcol,'LineWidth',2);
xline(730,'-.','color',yearcol,'LineWidth',2);

xlabel('Time (days)','interpreter','latex');
%ylabel({'Average contact rate,';'contacts per day'},'interpreter','latex');
title('B.1.1.7-like variant');
%legend(legComp);
set(gca,'FontSize',25);

%create bar charts for old and new variants cumulative infections
CumContainer(1,1)=100*(CumContainer(1,1)-baseCum(1,1))/baseCum(1,1);
CumContainer(1,2)=100*(CumContainer(1,2)-baseCum(1,2))/baseCum(1,2);
CumContainer(2,1)=100*(CumContainer(2,1)-baseCum(1,1))/baseCum(1,1);
CumContainer(2,2)=100*(CumContainer(2,2)-baseCum(1,2))/baseCum(1,2);

temp=CumContainer(1,2);
CumContainer(1,2)=CumContainer(2,1);
CumContainer(2,1)=temp;

CumContainer(3,1)=100*(CumContainer(3,1)-baseCum(2,1))/baseCum(2,1);
CumContainer(3,2)=100*(CumContainer(3,2)-baseCum(2,2))/baseCum(2,2);
CumContainer(4,1)=100*(CumContainer(4,1)-baseCum(2,1))/baseCum(2,1);
CumContainer(4,2)=100*(CumContainer(4,2)-baseCum(2,2))/baseCum(2,2);

temp=CumContainer(3,2);
CumContainer(3,2)=CumContainer(4,1);
CumContainer(4,1)=temp;

figure(5);
subplot(1,2,1);
p1=bar(CumContainer(1:2,:));
set(p1(1),'FaceColor',col_p(1,:));%[116 120 128]/255);
set(p1(2),'FaceColor',col_p(2,:));%[24 74 69]/255);

set(gca, 'XTick', [1 2])
set(gca, 'XTickLabel', {'Year 1' 'Year 2'})
title('Original variant');
ylabel('Excess infections (\%)','interpreter','latex');
ylim([-80,60]);
legend('Slow vaccination','Fast vaccination','Location','southoutside','interpreter','latex');
set(gca,'FontSize',25);

subplot(1,2,2);
p2=bar(CumContainer(3:4,:));
set(p2(1),'FaceColor',col_p(1,:));
set(p2(2),'FaceColor',col_p(2,:));
set(gca, 'XTick', [1 2])
set(gca, 'XTickLabel', {'Year 1' 'Year 2'})
title('B.1.1.7-like variant');
ylim([-80,60]);
set(gca,'FontSize',25);

%output settings for figures 8 and 9
figure(8);
set(gca,'FontSize',25);
subplot(1,2,1);
xline(365,'-.','color',yearcol,'LineWidth',2);
xline(730,'-.','color',yearcol,'LineWidth',2);

xlabel('Time (days)','interpreter','latex');
ylabel('Attack rate ($$\%$$)','interpreter','latex');
title('Original variant');
legend(h8(4:5),'Non-vaccinated','Vaccinated','NumColumns',2,'location','southoutside');
set(gca,'FontSize',25);
subplot(1,2,2);
xline(365,'-.','color',yearcol,'LineWidth',2);
xline(730,'-.','color',yearcol,'LineWidth',2);

xlabel('Time (days)','interpreter','latex');
title('B.1.1.7-like variant');

figure(9);
subplot(1,2,1);
xline(365,'-.','color',yearcol,'LineWidth',2);
xline(730,'-.','color',yearcol,'LineWidth',2);

xlabel('Time (days)','interpreter','latex');
ylabel('Attack rate ($$\%$$)','interpreter','latex');
title('Original variant');
set(gca,'FontSize',25);
subplot(1,2,2);
xline(365,'-.','color',yearcol,'LineWidth',2);
xline(730,'-.','color',yearcol,'LineWidth',2);
xlabel('Time (days)','interpreter','latex');
title('B.1.1.7-like variant');
legend(h9(4:5),'Non-vaccinated','Vaccinated','NumColumns',2,'location','southoutside');
set(gca,'FontSize',25);

%labeling of the panels
figure(5);annotation('textbox', [0.03, 0.98, 0, 0], 'string', 'c','FontWeight','bold','FontSize',30)
figure(5);annotation('textbox', [0.5, 0.98, 0, 0], 'string', 'd','FontWeight','bold','FontSize',30)
figure(6);annotation('textbox', [0.03, 0.98, 0, 0], 'string', 'a','FontWeight','bold','FontSize',30)
figure(6);annotation('textbox', [0.5, 0.98, 0, 0], 'string', 'b','FontWeight','bold','FontSize',30)
figure(8);annotation('textbox', [0.03, 0.98, 0, 0], 'string', 'a','FontWeight','bold','FontSize',30)
figure(8);annotation('textbox', [0.5, 0.98, 0, 0], 'string', 'b','FontWeight','bold','FontSize',30)
figure(9);annotation('textbox', [0.03, 0.98, 0, 0], 'string', 'c','FontWeight','bold','FontSize',30)
figure(9);annotation('textbox', [0.5, 0.98, 0, 0], 'string', 'd','FontWeight','bold','FontSize',30)