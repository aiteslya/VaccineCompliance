% The script produces sensitivity analysis with respect to per capita rate
% of moving to compliant state, delta for the population which is partially
% compliant with physical distancing measures and where vaccination takes
% place, the original avriant circulates
% initialize space
clc;
format long;
clear variables;

deltaarr=linspace(1e-8,5e-5,4000);

%recovery rate in the infectious stage
Prev=37706;
gamma=1/7; %1/day
Incid=Prev*gamma;
mu0=1/30;

func=100*deltaarr*Incid./(deltaarr*Incid+mu0);
figure(1);plot(deltaarr,func,'LineWidth',4);hold on;
xlabel('Per capita rate of moving to compliant state, $$\delta$$','Interpreter','latex');
ylabel('Compliant population ($$\%$$)','interpreter','latex');
figure(1);linehandle=xline(4e-5,'-.','LineWidth',2,'color',[0.40,0.62,0.11]);
figure(1);linehandle2=xline(6.19781e-6,'-.','LineWidth',2,'color',[0.93,0.69,0.13]);
figure(1);linehandle2=xline(1.54758e-6,'-.','LineWidth',2,'color',[0.64,0.08,0.18]);
figure(1);linehandle2=yline(20,'-.','LineWidth',2,'color',[0.64,0.08,0.18]);
figure(1);linehandle2=yline(50,'-.','LineWidth',2,'color',[0.93,0.69,0.13]);
figure(1);linehandle2=yline(86.6056,'-.','LineWidth',2,'color',[0.40,0.62,0.11]);
xlim([0,5e-5])
ylim([0,100])
set(gca,'FontSize',25);

%array to be used in the subsequent production of figures
deltaArr=[1.55e-6,6.2e-6,4e-5];
upsilonarr=[0,5.9e-4,4.9e-3];

%set initial conditions

c=8.752187609380469;
r1=0.338983050847458;

%set up the rest of the parameters
chat=14.9;
R0=2.5;
%calculate epsilon
epsilon=R0*gamma/chat;
beta=epsilon*c;

r2=13.47/c;
frac=1/3; %vaccination coverage when the compliance lasts for a week
mu1=3./(frac*5.1e8);
alpha=1/4;
k1=1;
k2=1;
omega=0.6;

%set up untegration options
T=800;
init=1.0e7*[0.55532617 0.00075412 0.00131971 0.0476 1.012748601428572 0.001400508571429 0.00245089 0.0884 0 0 0 0 0 0 0 0 0];

Atol=1e-11;
opts = odeset('RelTol',1e-12,'AbsTol',Atol);

%set up ploting options
popOut=1e5;
N=1.7e7;
XpopOut=N/popOut;
%set up colors
cols=[0.9,0.13,0.13;
    0.45,0.76,0.98;
    0.01, 0.4, 0.76];
yearcol=[0.17,0.17,0.17];
deltacounter=1;
for delta=deltaArr
    %counter for vaccination rates
    upscounter=1;
    for upsilon0=upsilonarr
        pars=[beta,r1,r2,delta,mu0,mu1,upsilon0,alpha,gamma,k1,k2,omega];
        [t,y]=ode45(@(t,y)COVIDVaccineRHS2(t,y,pars),[0,T], init,opts);
        %collect outputs
        infected=(y(:,2)+y(:,3)+y(:,6)+y(:,7)+y(:,11)+y(:,12)+y(:,15)+y(:,16))/XpopOut;
        compl=(y(:,5)+y(:,6)+y(:,7)+y(:,8))/XpopOut;
        %plot infected
        figure(2);
        subplot(1,3,deltacounter);
        %collect line for legend
        plot(t,infected,'color',cols(upscounter,:),'LineWidth',4);hold on;
        xlabel('Time (days)','interpreter','latex');
        xlim([0,T]);
        ylim([0,1000]);
        set(gca,'FontSize',25);
        %plot compliant
        figure(3);
        subplot(1,3,deltacounter);
        if delta==deltaArr(1)
            h1(upscounter)=plot(t,compl,'color',cols(upscounter,:),'LineWidth',4);hold on;
        else
            plot(t,compl,'color',cols(upscounter,:),'LineWidth',4);hold on;
        end
        xlabel('Time (days)','interpreter','latex');
        xlim([0,T]);
        ylim([0,1e5]);
        set(gca,'FontSize',25);
        upscounter=upscounter+1;
    end
    deltacounter=deltacounter+1;
end
figure(2);
subplot(1,3,1);
xline(365,'-.','color',yearcol,'LineWidth',2);
xline(730,'-.','color',yearcol,'LineWidth',2);
ylabel('Infected individuals (1/100,000)','interpreter','latex');
title('Slow compliance acquisition','interpreter','latex');
subplot(1,3,2);
xline(365,'-.','color',yearcol,'LineWidth',2);
xline(730,'-.','color',yearcol,'LineWidth',2);
title('Medium compliance acquisition','interpreter','latex');
subplot(1,3,3);
xline(365,'-.','color',yearcol,'LineWidth',2);
xline(730,'-.','color',yearcol,'LineWidth',2);
title('Fast compliance acquisition','interpreter','latex');
figure(3);
subplot(1,3,1);
xline(365,'-.','color',yearcol,'LineWidth',2);
xline(730,'-.','color',yearcol,'LineWidth',2);
ylabel('Compliant individuals (1/100,000)','interpreter','latex');
subplot(1,3,2);
xline(365,'-.','color',yearcol,'LineWidth',2);
xline(730,'-.','color',yearcol,'LineWidth',2);
subplot(1,3,3);
xline(365,'-.','color',yearcol,'LineWidth',2);
xline(730,'-.','color',yearcol,'LineWidth',2);
legend(h1,'Baseline, no vaccination','Slow vaccination','Fast vaccination','location','southoutside','NumColumns',3);
figure(2);annotation('textbox', [0.03, 0.999, 0, 0], 'string', 'a','FontWeight','bold','FontSize',25)
figure(2);annotation('textbox', [0.35, 0.999, 0, 0], 'string', 'b','FontWeight','bold','FontSize',25)
figure(2);annotation('textbox', [0.65, 0.999, 0, 0], 'string', 'c','FontWeight','bold','FontSize',25)

figure(3);annotation('textbox', [0.03, 0.999, 0, 0], 'string', 'd','FontWeight','bold','FontSize',25)
figure(3);annotation('textbox', [0.35, 0.999, 0, 0], 'string', 'e','FontWeight','bold','FontSize',25)
figure(3);annotation('textbox', [0.65, 0.999, 0, 0], 'string', 'f','FontWeight','bold','FontSize',25)