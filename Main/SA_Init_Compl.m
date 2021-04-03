%This is the script that plots time series for old and
%new strain with and without the vaccination with high compliance rise and
%decay rates 
% during the lockdown for COVID-19-Compliance-Vaccine model
% R0 for the old virus strain is 2.5 and for the new 3.75
% such that Re(0)=1.1 for the old strain and Re(0)=1.65 for the new strain
% the cumulative output is adjusted for the initial number of infections
% these time series are done as a part of sensitivity analysis with respect
% to initial data: compliance settings
% altogether we will have 6 scenarios: 3 compliance settings, 2 virus
% settings, 12 figures
%prepare settings
clc;
clear variables;
close all;
format long;

%set the parameters
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
%compliance
mu0=1/30;
c=8.752187609380469;
r1=0.338983050847458;
delta=4e-5;
upsilonarr=[0,5.9e-4,4.9e-3];
k1=1;
k2=1;
omega=0.6;
%format of the legend
formatSpec = '%.2e';

% contact rate of vaccinated 13.47
r2=13.47/c;
frac=1/3;
mu1=3./(frac*5.1e8);
    
%set up initial data
% factor of detected/total
X=1;
%seroprevalence
SP=0.08;
%total population
N=1.7e7;
popOut=1e5;
XpopOut=N/popOut;
%percentage of compliant people
PerComplArr=[0.2,0.67, 0.9];

%integration options
Atol=1e-11;
opts = odeset('RelTol',1e-12,'AbsTol',Atol);
%integrating time
T=800;
%outer loop: circles through initial conditions
%figure counter
counterPer=1;
for PerCompl=PerComplArr
    N0=N*(1-PerCompl);
    Nc0=N*PerCompl;
    TotalInf=37706;
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
    infect0pop=popOut*infect0/N;

    V0=0;
    SV0=0;
    EV0=0;
    IV0=0;
    RV0=0;

    TV0=SV0+EV0+IV0+RV0;
    init=[S0,E0,I0,R0,Sc0,Ec0,Ic0,Rc0,V0,SV0,EV0,IV0,RV0,TV0,0,0,0];
        
    %container for vaccination cumulative infected at mark of year 1 and 2
    CumContainer=zeros(4,2);
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

            if upsilon0==0 %plot 
                figure(i1);
                subplot(2,3,1);
                %plot(t,infectious,'LineWidth',2+(counterPer-1)*2,'color',[0,113,188]/255);hold on;
                plot(t,infectious,'LineWidth',4);hold on;
                subplot(2,3,4);
                %plot(t,compl,'LineWidth',2+(counterPer-1)*2,'color',[0,113,188]/255);hold on;
                plot(t,compl,'LineWidth',4);hold on;
            elseif upsilon0==5.9e-4
                figure(i1);
                subplot(2,3,2);
                %plot(t,infectious,'LineWidth',2+(counterPer-1)*2,'color',col_p(1,:));hold on;
                plot(t,infectious,'LineWidth',4);hold on;
                subplot(2,3,5);
                %h1(i1,counterPer)=plot(t,compl,'LineWidth',2+(counterPer-1)*2,'color',col_p(1,:));hold on;
                h1(i1,counterPer)=plot(t,compl,'LineWidth',4);hold on;
            else
                figure(i1);
                subplot(2,3,3);
                %plot(t,infectious,'LineWidth',2+(counterPer-1)*2,'color',col_p(2,:));hold on;
                plot(t,infectious,'LineWidth',4);hold on;
                subplot(2,3,6);
                %plot(t,compl,'LineWidth',2+(counterPer-1)*2,'color',col_p(2,:));hold on;
                plot(t,compl,'LineWidth',4);hold on;
            end
           
        end
        
        
        disp(['Given proportion of compliant ',num2str(PerCompl*100),'%, the contact rates are c=',num2str(c),' and cr1=',num2str(c*r1)]);
              
        i1=i1+1;
    end
    %figure settings
    for j=1:2
        figure(j);
        subplot(2,3,1);
        xlim([0, 800]);
        xline(365,'-.','color',yearcol,'LineWidth',2);
        xline(730,'-.','color',yearcol,'LineWidth',2);
        xlabel('Time (days)','interpreter','latex');
        ylabel('Infected individuals (1/100,000)','interpreter','latex');
        title('No vaccination','interpreter','latex');
        set(gca,'FontSize',25);

        subplot(2,3,2);
        xlim([0, 800]);
        xline(365,'-.','color',yearcol,'LineWidth',2);
        xline(730,'-.','color',yearcol,'LineWidth',2);
        xlabel('Time (days)','interpreter','latex');
        title('Slow vaccination','interpreter','latex');
        set(gca,'FontSize',25);

        subplot(2,3,3);
        xlim([0, 800]);
        xline(365,'-.','color',yearcol,'LineWidth',2);
        xline(730,'-.','color',yearcol,'LineWidth',2);
        xlabel('Time (days)','interpreter','latex');
        title('Fast vaccination','interpreter','latex');
        set(gca,'FontSize',25);
        
        subplot(2,3,4);
        xlim([0, 800]);
        xline(365,'-.','color',yearcol,'LineWidth',2);
        xline(730,'-.','color',yearcol,'LineWidth',2);
        xlabel('Time (days)','interpreter','latex');
        ylabel('Compliant individuals (1/100,000)','interpreter','latex');
        set(gca,'FontSize',25);

        subplot(2,3,5);
        xlim([0, 800]);
        xline(365,'-.','color',yearcol,'LineWidth',2);
        xline(730,'-.','color',yearcol,'LineWidth',2);
        legend(h1(j,:),'20\% compliant','67\% compliant','90\% compliant','numColumns',3,'location','southoutside','interpreter','latex');
        xlabel('Time (days)','interpreter','latex');
        set(gca,'FontSize',25);

        subplot(2,3,6);
        xlim([0, 800]);
        xline(365,'-.','color',yearcol,'LineWidth',2);
        xline(730,'-.','color',yearcol,'LineWidth',2);
        xlabel('Time (days)','interpreter','latex');
        set(gca,'FontSize',25);
    end 
    counterPer=counterPer+1;

end
figure(1);annotation('textbox', [0.03, 0.999, 0, 0], 'string', 'a','FontWeight','bold','FontSize',25)
figure(1);annotation('textbox', [0.35, 0.999, 0, 0], 'string', 'b','FontWeight','bold','FontSize',25)
figure(1);annotation('textbox', [0.65, 0.999, 0, 0], 'string', 'c','FontWeight','bold','FontSize',25)

figure(1);annotation('textbox', [0.03, 0.499, 0, 0], 'string', 'd','FontWeight','bold','FontSize',25)
figure(1);annotation('textbox', [0.35, 0.499, 0, 0], 'string', 'e','FontWeight','bold','FontSize',25)
figure(1);annotation('textbox', [0.65, 0.499, 0, 0], 'string', 'f','FontWeight','bold','FontSize',25)

figure(2);annotation('textbox', [0.03, 0.999, 0, 0], 'string', 'a','FontWeight','bold','FontSize',25)
figure(2);annotation('textbox', [0.35, 0.999, 0, 0], 'string', 'b','FontWeight','bold','FontSize',25)
figure(2);annotation('textbox', [0.65, 0.999, 0, 0], 'string', 'c','FontWeight','bold','FontSize',25)

figure(2);annotation('textbox', [0.03, 0.499, 0, 0], 'string', 'd','FontWeight','bold','FontSize',25)
figure(2);annotation('textbox', [0.35, 0.499, 0, 0], 'string', 'e','FontWeight','bold','FontSize',25)
figure(2);annotation('textbox', [0.65, 0.499, 0, 0], 'string', 'f','FontWeight','bold','FontSize',25)