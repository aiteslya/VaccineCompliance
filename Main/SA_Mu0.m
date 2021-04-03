%This is the script that plots time series for old and
%new strain without the vaccination given three different compliance decay
%rates
% during the lockdown for COVID-19-Compliance-Vaccine model
% R0 for the old virus strain is 2.5 and for the new 3.75
% such that Re(0)=1.1 for the old strain and Re(0)=1.65 for the new strain
% the cumulative output is adjusted for the initial number of infections
% these time series are done as a part of sensitivity analysis with respect
% to initial data: prevalence settings
% altogether we will have 6 scenarios: 3 mu0 settings, 2 virus
% settings, 4 figures
%prepare settings
clc;
clear variables;
close all;
format long;

%set the parameters

%contact rate of non-compliant before the lockdown at the start of the
%epidemic
chat=14.9;
%1/gamma duration of infectious period
gamma=1/7;
alpha=1/4;
R0arr=[2.5, 2.5*1.5];
%calculate epsilon
epsilonarr=R0arr*gamma/chat;
yearcol=[0.17,0.17,0.17];

%set up initial data
% factor of detected/total
X=1;
%prevalence
TotalInf=37706;
%seroprevalence
SP=0.08;
%total population
N=1.7e7;
%percentage of compliant people
PerCompl=0.65;

N0=N*(1-PerCompl);
Nc0=N*PerCompl;
TotalRec=SP*N;
TotalS=N-TotalInf-TotalRec;
%setting up of initial data
S0=(1-PerCompl)*TotalS;
E0=(1-PerCompl)*TotalInf/2;
I0=(1-PerCompl)*TotalInf/2;
R0=(1-PerCompl)*TotalRec;
Sc0=PerCompl*TotalS;
Ec0=PerCompl*TotalInf/2;
Ic0=PerCompl*TotalInf/2;
Rc0=PerCompl*TotalRec;

infect0=E0+I0+R0+Ec0+Ic0+Rc0;
popOut=1e5;
infect0pop=popOut*infect0/N;

V0=0;
SV0=0;
EV0=0;
IV0=0;
RV0=0;

TV0=V0+SV0+EV0+IV0+RV0;
init=[S0,E0,I0,R0,Sc0,Ec0,Ic0,Rc0,V0,SV0,EV0,IV0,RV0,TV0,0,0,0];
%compliance
mu0arr=[1/60,1/30,1/15];

%define array of contact rates, for in-loop calculations
r1num=60;
r1arr=linspace(0,1,r1num);
%define the array of r1
num=2e4;
carr=linspace(0,15,num);

%integration options
Atol=1e-11;
opts = odeset('RelTol',1e-12,'AbsTol',Atol);
%integrating time
T=800;

delta=4e-5;
upsilon0=0;
omega=0.6;
k1=1;
k2=1;

frac=1/3;
mu1=0;%3./(frac*5.1e8);
mu1hat=mu1*N;
c=8.752187609380469;
r1=0.338983050847458;
% contact rate of vaccinated 13.47
r2=13.47/c;
XpopOut=N/popOut;

%figure counter
counterPer=1;
for mu0=mu0arr
    %strain counter
    i1=1;
    for epsilon=epsilonarr
        beta=epsilon*c;
        i2=1;
        pars=[beta,r1,r2,delta,mu0,mu1,upsilon0,alpha,gamma,k1,k2,omega];
        [t,y]=ode45(@(t,y)COVIDVaccineRHS2(t,y,pars),[0,T], init,opts);
        infectious=(y(:,2)+y(:,3)+y(:,6)+y(:,7)+y(:,11)+y(:,12)+y(:,15)+y(:,16))/XpopOut;
        cum=(y(:,2)+y(:,3)+y(:,4)+y(:,6)+y(:,7)+y(:,8)+y(:,11)+y(:,12)+y(:,13)+y(:,15)+y(:,16)+y(:,17))/XpopOut-infect0pop;
        compl=(y(:,5)+y(:,6)+y(:,7)+y(:,8))/XpopOut;
        figure(1);
        subplot(2,2,i1);
        plot(t,infectious,'LineWidth',4);hold on;
        subplot(2,2,2+i1);
        h1(i1,counterPer)=plot(t,compl,'LineWidth',4);hold on;
        i1=i1+1;
    end
    counterPer=counterPer+1;
end
%figure settings
figure(1)
subplot(2,2,1);
xlim([0, 800]);
xline(365,'-.','color',yearcol,'LineWidth',2);
xline(730,'-.','color',yearcol,'LineWidth',2);

xlabel('Time (days)','interpreter','latex');
ylabel('Infected individuals (1/100,000)','interpreter','latex');
title('Original variant','interpreter','latex');
set(gca,'FontSize',25);

subplot(2,2,2);
xlim([0, 800]);
xline(365,'-.','color',yearcol,'LineWidth',2);
xline(730,'-.','color',yearcol,'LineWidth',2);

xlabel('Time (days)','interpreter','latex');
title('B.1.1.7 variant','interpreter','latex');
set(gca,'FontSize',25);

subplot(2,2,3);
xlim([0, 800]);
xline(365,'-.','color',yearcol,'LineWidth',2);
xline(730,'-.','color',yearcol,'LineWidth',2);

xlabel('Time (days)','interpreter','latex');
ylabel('Compliant individuals (1/100,000)','interpreter','latex');
set(gca,'FontSize',25);

subplot(2,2,4);
xlim([0, 800]);
xline(365,'-.','color',yearcol,'LineWidth',2);
xline(730,'-.','color',yearcol,'LineWidth',2);
xlabel('Time (days)','interpreter','latex');
set(gca,'FontSize',25);

legend(h1(2,:),'Slow compliance loss','Medium compliance loss','Fast compliance loss','numColumns',4,'location','southoutside','interpreter','latex');
xlabel('Time (days)','interpreter','latex');
set(gca,'FontSize',25);

figure(1);annotation('textbox', [0.03, 0.999, 0, 0], 'string', 'a','FontWeight','bold','FontSize',25)
figure(1);annotation('textbox', [0.49, 0.999, 0, 0], 'string', 'b','FontWeight','bold','FontSize',25)
figure(1);annotation('textbox', [0.03, 0.53, 0, 0], 'string', 'c','FontWeight','bold','FontSize',25)
figure(1);annotation('textbox', [0.49, 0.53, 0, 0], 'string', 'd','FontWeight','bold','FontSize',25)