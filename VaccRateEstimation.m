%this file estimates vaccination rates approximating the rates of
%vaccination observed to date in the Netherlands and UK. These will be
%dubbed, slow and fast vaccination rates, respectively

%prepare the space
clc;
clear variables;
close all;
format long;

%total population
N=1.7e7;
%read the CSV files and load them into array
TblNed = readtable('Data/VaccNetherlands.csv');
TVPHNed = (TblNed.total_vaccinations_per_hundred)/100; 
DateNed = TblNed.Date;%+6;
TblUK = readtable('Data/VaccUK.csv');
TVPHUK = (TblUK.total_vaccinations_per_hundred)/100; 
DateUK = TblUK.date;%+3;
TPVPHUK = (TblUK.people_vaccinated_per_hundred)/100;

figure(1);
h(1)=plot(DateNed,TVPHNed,'*','MarkerSize',6);hold on;
h(2)=plot(DateUK,TVPHUK,'*','MarkerSize',6);hold on;
h(3)=plot(DateUK,TPVPHUK,'o','MarkerSize',6);hold on;
set(h(3),'MarkerEdgeColor',get(h(2),'color'));
xlabel('Time, days','interpreter','latex');
ylabel('Density of vaccinated','interpreter','latex');

leg=['Netherlands vaccine doses  ';'UK vaccine doses           ';'UK people vaccinated       ';'Lin approx to NL doses data';'Lin approx to UK doses data';'Lin approx to NL peopl data'];

%calculate slopes of the three data sets
mNether=DateNed\TVPHNed
mUKDoses=DateUK\TVPHUK
mUKPeople=DateUK\TPVPHUK;

numpoints=10;
tdiscrNed=linspace(min(DateNed),max(DateNed),numpoints);
tdiscrUK=linspace(min(DateUK),max(DateUK),numpoints);
figure(1)
h(4)=plot(tdiscrNed,tdiscrNed*mNether);
set(h(4),'Color',get(h(1),'color'));
h(5)=plot(tdiscrUK,tdiscrUK*mUKDoses);
set(h(5),'Color',get(h(2),'color'));
h(6)=plot(tdiscrUK,tdiscrUK*mUKPeople,'-.');
set(h(6),'Color',get(h(2),'color'));

legend(h,leg,'Location','eastoutside');
set(gca,'FontSize',25);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is the script that plots how much vaccinations are realized within a certain period of time
% with vaccination in place as well as how decay of compliance depends on the
% vaccinated population: 1 figure in total
% during the lockdown for COVID-19-Compliance-Vaccine model
% This is the code for the old variant with R0=2.5 where
% Re(0)=1.1 for the new strain


%set the parameters
%define colors
numout=20;
red = [0.6350, 0.0780, 0.1840];%[1 0 0];
pink = [255, 192, 203]/255;
colors_p = [linspace(pink(1),red(1),numout)', linspace(pink(2),red(2),numout)', linspace(pink(3),red(3),numout)'];
%col_p=colors_p([2,4],:);
%contact rate of non-compliant before the lockdown at the start of the
%epidemic
chat=14.9;
%1/gamma duration of infectious period
gamma=1/5;
alpha=1/4;
R0=2.5;
%calculate epsilon
epsilon=R0*gamma/chat;

%set up initial data
% factor of detected/total
X=1;
%seroprevalence
SP=0.08;
%percentage of compliant people
PerCompl=0.67;
N0=N*(1-PerCompl);
Nc0=N*PerCompl;
TotalInf=43522;
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
V0=0;
SV0=0;
EV0=0;
IV0=0;
RV0=0;
TV0=SV0+EV0+IV0+RV0;
init=[S0,E0,I0,R0,Sc0,Ec0,Ic0,Rc0,V0,SV0,EV0,IV0,RV0,TV0];
%compliance
mu0=1/30;

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

ind=find(Res(:,2)>0.355,1);

 % % %
 
c=Res(ind,1);
r1=Res(ind,2);

%integration options
Atol=1e-11;
opts = odeset('RelTol',1e-12,'AbsTol',Atol);
%integrating time
T=365;

delta=4e-5;
mu0=1/30;
Vtarg=linspace(0.5e6,17e6,numout);
upsilonarr=[5.9e-4,4.9e-3];%log(1.7e7./(-Vtarg+1.7e7))/365;
k1=1;
k2=1;
omega=0.6;
%format of the legend
formatSpec = '%.2e';

r2=1.5;
%frac=[2/3,1/3];
frac=1/3;
mu1arr=[0,3./(frac*5.1e8)];

fighandle=2;

%legend for the vaccinated plot
legVacc=['Slow vaccination';'Fast vaccination';'Netherlands vacc';'Israel vacc     ';'Netherlands data';'Israel data     '];
%legComp=['Const decay                  ';'Slow  decay, slow vaccination';'Fast  decay, slow vaccination';'Slow  decay, fast vaccination';'Fast  decay, fast vaccination'];
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
        if mu1==0 %plot vaccinated density
            figure(1);
            h1(vaccRateCount)=plot(tdiscr,vaccdiscr/N,'LineWidth',4,'color',colors_p(vaccRateCount,:));hold on;
        end
        mu=mu0+mu1*vaccdiscr;
        switch mu1
            case 0
                if upsilon0==log(1.7e7./(-1e6+1.7e7))/365
                    figure(3);plot(tdiscr,1./mu,'-.','LineWidth',4);hold on;
                end
            case mu1arr(2)
                %figure(2);plot(tdiscr,1./mu,'-s','MarkerSize',8,'color',get(h1,'color'),'MarkerEdgeColor',get(h1,'color'));hold on;
                figure(2);plot(tdiscr,1./mu,'LineWidth',4,'color',get(h1(vaccRateCount),'color'),'MarkerEdgeColor',get(h1(vaccRateCount),'color'));hold on;
            otherwise
                figure(2);plot(tdiscr,1./mu,'-^','MarkerSize',8,'color',get(h1(vaccRateCount),'color'),'MarkerEdgeColor',get(h1(vaccRateCount),'color'));hold on;
        end
    end
    vaccRateCount=vaccRateCount+1;
end
%figure settings
figure(1);
%ylim([0,1]);
grid on;
xlabel('Time, days','interpreter','latex');
ylabel('Density','interpreter','latex');

set(gca,'FontSize',25);

figure(2);
xlim([0,T]);
grid on;
xlabel('Time, days','interpreter','latex');
ylabel({'Average duration';'of compliance, days'},'interpreter','latex');
%legend(legComp,'Location','northeastoutside');
legend(legComp);
set(gca,'FontSize',25);