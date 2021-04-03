%This is the script that plots how much vaccinations are realized within a certain period of time
% with vaccination in place as well as how decay of compliance depends on the
% vaccinated population: 1 figure in total
% during the lockdown for COVID-19-Compliance-Vaccine model
% This is the code for the old variant with R0=2.5 where
% Re(0)=1.1 for the new strain
%prepare settings
clc;
clear variables;
close all;
format long;

%set the parameters
%define colors
l=5;
red = [0.01, 0.4, 0.76];%[1 0 0];
pink = [0.45,0.76,0.98];
colors_p = [linspace(pink(1),red(1),l)', linspace(pink(2),red(2),l)', linspace(pink(3),red(3),l)'];
col_p=colors_p([2,4],:);
%contact rate of non-compliant before the lockdown at the start of the
%epidemic
chat=14.9;
%1/gamma duration of infectious period
gamma=1/7;
alpha=1/4;
R0=2.5;
%calculate epsilon
epsilon=R0*gamma/chat;

%set up initial data
%seroprevalence
SP=0.08;
%total population
N=1.7e7;
TotalInf=37706;
Incid=TotalInf*gamma;
mu0=1/30;
delta=4e-5;
%percentage of compliant people
PerCompl=0.65;%delta*Incid/(delta*Incid+mu0);
N0=N*(1-PerCompl);
Nc0=N*PerCompl;

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
init=[S0,E0,I0,R0,Sc0,Ec0,Ic0,Rc0,V0,SV0,EV0,IV0,RV0,TV0];
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
    beta=carr*epsilon;
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

upsilonarr=[5.9e-4,4.9e-3];

k1=1;
k2=1;
omega=0.6;
%format of the legend
formatSpec = '%.2e';

r2=13.47/c;
%frac=[2/3,1/3];
frac=1/3;
mu1arr=[0,3./(frac*5.1e8)];

fighandle=2;

%legend for the vaccinated plot
legVacc=['Slow vaccination';'Fast vaccination';'Netherlands data';'UK data         '];

legComp=['Const decay     ';'Slow vaccination';'Fast vaccination'];
counter=1;
numpoints=30;
tdiscr=linspace(0,T,numpoints);
vaccRateCount=1;
popOut=1e5;
XpopOut=N/popOut;

for upsilon0=upsilonarr
    for mu1=mu1arr
        beta=epsilon*c;
        pars=[beta,r1,r2,delta,mu0,mu1,upsilon0,alpha,gamma,k1,k2,omega];
        [t,y]=ode45(@(t,y)COVIDVaccineRHS(t,y,pars),[0,T], init,opts);
        vacc=y(:,14);
        vaccdiscr=interp1(t,vacc,tdiscr);
        compl=y(:,5)+y(:,6)+y(:,7)+y(:,8);
        if mu1==0 %plot vaccinated density
            figure(1);subplot(1,2,1);
            h1(vaccRateCount)=plot(tdiscr,100*vaccdiscr/N,'LineWidth',4,'color',col_p(vaccRateCount,:));hold on;
        end
        mu=mu0+mu1*vaccdiscr;
        switch mu1
            case 0
%                 if upsilon0==5.9e-4
%                     %figure(2);
%                     figure(1);subplot(1,2,2);
%                     plot(tdiscr,1./mu,'-.','LineWidth',4);hold on;
%                 end
            case mu1arr(2)
                %figure(2);plot(tdiscr,1./mu,'-s','MarkerSize',8,'color',get(h1,'color'),'MarkerEdgeColor',get(h1,'color'));hold on;
                %figure(2);
                init1=[(1-PerCompl)*N,PerCompl*N,0];
                pars1=[delta,mu0,mu1,upsilon0,Incid];
                [t1,y1]=ode45(@(t1,y1)VaccComplReduced(t1,y1,pars1),[0,T], init1,opts);
                figure(1);subplot(1,2,2);
                plot(t1,100*y1(:,2)/(1.7e7),'LineWidth',4,'color',get(h1(vaccRateCount),'color'),'MarkerEdgeColor',get(h1(vaccRateCount),'color'));hold on;
            otherwise
                init1=[(1-PerCompl)*N,PerCompl*N,0];
                pars1=[delta,mu0,mu1,upsilon0,Incid];
                [t1,y1]=ode45(@(t1,y1)VaccComplReduced(t1,y1,pars1),[0,T], init1,opts);
                figure(1);subplot(1,2,2);
                plot(t1,100*y1(:,2)/(1.7e7),'LineWidth',4,'color',get(h1(vaccRateCount),'color'),'MarkerEdgeColor',get(h1(vaccRateCount),'color'));hold on;
        end
    end
    vaccRateCount=vaccRateCount+1;
end
%figure settings
figure(1);
subplot(1,2,1);
%ylim([0,1]);
xline(365,'-.','color',[0.17,0.17,0.17],'LineWidth',2);
xline(730,'-.','color',[0.17,0.17,0.17],'LineWidth',2);
xlim([0,T])
%grid on;
xlabel('Time (days)','interpreter','latex');

set(gca,'FontSize',25);

%add the Netherlands and UK data to the vaccination plot, Figure 1
%read the CSV files and load them into array
TblNed = readtable('Data/VaccNetherlands.csv');
TVPHNed = (TblNed.total_vaccinations_per_hundred);%/100; 
DateNed = TblNed.Date;%+6;
TblUK = readtable('Data/VaccUK.csv');
TVPHUK = (TblUK.total_vaccinations_per_hundred);%/100; 
DateUK = TblUK.date;%+3;

figure(1);
h(1)=plot(DateNed,TVPHNed,'*','MarkerSize',10);hold on;
set(h(1),'color',get(h1(1),'color'));
h(2)=plot(DateUK,TVPHUK,'*','MarkerSize',10);hold on;
set(h(2),'color',get(h1(2),'color'));

xlabel('Time (days)','interpreter','latex');
ylabel({'Vaccination coverage, ($$\%$$)'},'interpreter','latex');
legend([h1,h],'Slow vaccination','Fast vaccination','Netherlands data','UK data         ','location','southoutside','interpreter','latex','NumColumns',4);

%calculate slopes of the three data sets
mNether=DateNed\TVPHNed;
mUKDoses=DateUK\TVPHUK;

figure(1); subplot(1,2,2);
ylim([0,100]);
xline(365,'-.','color',[0.17,0.17,0.17],'LineWidth',2);
xline(730,'-.','color',[0.17,0.17,0.17],'LineWidth',2);
xlim([0,T]);
%grid on;
xlabel('Time (days)','interpreter','latex');
ylabel({'Compliant population ($$\%$$)'},'interpreter','latex');
%legend('Constant','Slow vaccination','Fast vaccination','interpreter','latex');
set(gca,'FontSize',25);

figure(1);annotation('textbox', [0.03, 0.98, 0, 0], 'string', 'a','FontWeight','bold','FontSize',30)
figure(1);annotation('textbox', [0.5, 0.98, 0, 0], 'string', 'b','FontWeight','bold','FontSize',30)
