%This is the script that plots time series for old and
%new strain with the vaccination with high compliance rise and
%decay rates 
% during the lockdown for COVID-19-Compliance-Vaccine model
% R0 for the old virus strain is 2.5 and for the new 3.75
% such that Re(0)=1.1 for the old strain and Re(0)=1.65 for the new strain
% the cumulative output is adjusted for the initial number of infections
% these time series are done as a part of sensitivity analysis with respect
% to uptake of vaccination by infectious people, k1
% altogether we will have 6 scenarios: 3 k1 settings, 2 virus, 2
% vaccinations
% settings, 8 figures
%prepare settings
clc;
clear variables;
close all;
format long;

%set the parameters

%define colors
l=5;
dblue = [0.01, 0.4, 0.76];
lblue = [0.45,0.76,0.98];
medblue=(dblue+lblue)/2;
novacc=[0.9,0.13,0.13];
yearcol=[0.17,0.17,0.17];
colors_p = [linspace(lblue(1),dblue(1),l)', linspace(lblue(2),dblue(2),l)', linspace(lblue(3),dblue(3),l)'];
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
mu0=1/30;

%define array of contact rates
r1num=60;
r1arr=linspace(0,1,r1num);
%define the array of r1
num=2e4;
carr=linspace(0,15,num);

%for each c in carr calculate r1
%set up equation Re=1.1

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
%set up the array of initial average contacts:
%c*ProportionNonCompl+c*r1*ProportionCompl
temparr=Res(:,1)*(1-PerCompl)+Res(:,2).*Res(:,1)*PerCompl;
%find index where the average contact rate is closest to 5
temparr0=temparr-5;
eqn=temparr0(2:r1num).*temparr0(1:(r1num-1));
ind=find(eqn<0,1);

c=Res(ind,1);
r1=Res(ind,2);

%integration options
Atol=1e-11;
opts = odeset('RelTol',1e-12,'AbsTol',Atol);
%integrating time
T=800;

delta=4e-5;
upsilonarr=[5.9e-4,4.9e-3];
k1arr=[0,0.5,1];
k2=1;

%format of the legend
formatSpec = '%.2e';

% contact rate of vaccinated 13.47
r2=13.47/c;
frac=1/3;
mu1=3./(frac*5.1e8);

XpopOut=N/popOut;
omega=0.6;
%outer loop: circles through omega
%figure counter
counterPer=1;
for k1=k1arr
    %strain counter
    i1=1;
    for epsilon=epsilonarr
        beta=epsilon*c;
        i2=1;
        for upsilon0=upsilonarr
            pars=[beta,r1,r2,delta,mu0,mu1,upsilon0,alpha,gamma,k1,k2,omega];
            [t,y]=ode45(@(t,y)COVIDVaccineRHS2(t,y,pars),[0,T], init,opts);
            infectious=(y(:,2)+y(:,3)+y(:,6)+y(:,7)+y(:,11)+y(:,12)+y(:,15)+y(:,16))/XpopOut;
            compl=(y(:,5)+y(:,6)+y(:,7)+y(:,8))/XpopOut;

            if upsilon0==5.9e-4
                figure(i1);
                subplot(2,2,1);
                %plot(t,infectious,'LineWidth',2+(counterPer-1)*2,'color',col_p(1,:));hold on;
                plot(t,infectious,'LineWidth',4);hold on;
              
                subplot(2,2,3);
                %h1(i1,counterPer)=plot(t,compl,'LineWidth',2+(counterPer-1)*2,'color',col_p(1,:));hold on;
                %h1(i1,2)=plot(t,compl,'LineWidth',4);hold on;
                plot(t,compl,'LineWidth',4);hold on;
            else
                figure(i1);
                subplot(2,2,2);
                %plot(t,infectious,'LineWidth',2+(counterPer-1)*2,'color',col_p(2,:));hold on;
                plot(t,infectious,'LineWidth',4);hold on;
                subplot(2,2,4);
                %plot(t,compl,'LineWidth',2+(counterPer-1)*2,'color',col_p(2,:));hold on;
                %h1(i1,3)=plot(t,compl,'LineWidth',4);hold on;
                h1(i1,counterPer)=plot(t,compl,'LineWidth',4);hold on;
            end
           
        end
        i1=i1+1;
    end
    %figure settings
    for j=1:2
        figure(j);
        subplot(2,2,1);
        xlim([0, 800]);
        xline(365,'-.','color',yearcol,'LineWidth',2);
        xline(730,'-.','color',yearcol,'LineWidth',2);
        xlabel('Time (days)','interpreter','latex');
        ylabel('Infected individuals (1/100,000)','interpreter','latex');
        title('Slow vaccination','interpreter','latex');
        set(gca,'FontSize',25);

        subplot(2,2,2);
        xlim([0, 800]);
        xline(365,'-.','color',yearcol,'LineWidth',2);
        xline(730,'-.','color',yearcol,'LineWidth',2);
        xlabel('Time (days)','interpreter','latex');
        title('Fast vaccination','interpreter','latex');
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
        
        legend(h1(j,:),'No vaccine uptake','Low vaccine uptake','High vaccine uptake','numColumns',3,'location','southoutside','interpreter','latex');
        xlabel('Time (days)','interpreter','latex');
        set(gca,'FontSize',25);

    end 
    counterPer=counterPer+1;

end

figure(1);annotation('textbox', [0.03, 0.999, 0, 0], 'string', 'a','FontWeight','bold','FontSize',25)
figure(1);annotation('textbox', [0.49, 0.999, 0, 0], 'string', 'b','FontWeight','bold','FontSize',25)
figure(1);annotation('textbox', [0.03, 0.53, 0, 0], 'string', 'c','FontWeight','bold','FontSize',25)
figure(1);annotation('textbox', [0.49, 0.53, 0, 0], 'string', 'd','FontWeight','bold','FontSize',25)

figure(2);annotation('textbox', [0.03, 0.999, 0, 0], 'string', 'a','FontWeight','bold','FontSize',25)
figure(2);annotation('textbox', [0.49, 0.999, 0, 0], 'string', 'b','FontWeight','bold','FontSize',25)
figure(2);annotation('textbox', [0.03, 0.53, 0, 0], 'string', 'c','FontWeight','bold','FontSize',25)
figure(2);annotation('textbox', [0.49, 0.53, 0, 0], 'string', 'd','FontWeight','bold','FontSize',25)
