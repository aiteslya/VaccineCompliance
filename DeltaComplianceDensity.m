%This is the script that outputs compliance curves for different
%delta's and find the maximum compliance achieved for each delta
% This is the code for the old Variant with R0=2.5, without lockdown

%prepare settings
clc;
clear variables;
close all;
format long;

%set the parameters
%contact rate of non-compliant before the lockdown at the start of the
%epidemic
c=14.9;
%1/gamma duration of infectious period
gamma=1/5;
alpha=1/4;
R0=2.5;
%calculate epsilon
epsilon=R0*gamma/c;

%set up initial data
%total population
N=1.7e7;
N0=N;
Nc0=0;
S0=N0-1;
E0=1;
I0=0;
R0=0;
Sc0=0;
Ec0=0;
Ic0=0;
Rc0=0;
V0=0;
SV0=0;
EV0=0;
IV0=0;
RV0=0;
TV0=SV0+EV0+IV0+RV0;
init=[S0,E0,I0,R0,Sc0,Ec0,Ic0,Rc0,V0,SV0,EV0,IV0,RV0,TV0];
%compliance
mu0=1/30;
mu1=0;

c=14.9;
% we are taking r1 from setting Re=1.1, c=8.164,r1=0.356
%this r1 is different from the one above since we still need to get the
%same contcat rate for compliant with and without lockdown
r1=3.2/c;
%r1=1;
%initial simulations
%integration options
Atol=1e-11;
opts = odeset('RelTol',1e-12,'AbsTol',Atol);
%integrating time
T=200;
%T=40;

%pure compliance process: looking for the compliance spread rate
%define parameter array
% modification in contact rate of the vaccinated individuals, for now set
% to 1.5
r2=1;
%deltanum=20;
%deltarr=10.^(linspace(-8,-4,deltanum));
deltarr=[4.28e-8,3.8e-6,4e-5];
deltanum=numel(deltarr);

k1=1;
k2=1;

%format of the legend
formatSpec = '%.2e';

iter=1;
fig_handle=2;
%number of points for discretization
n_points=100;
beta=c*epsilon;
omega=0;
upsilon0=0; 
%array containing peak compliant
PC=zeros(1,deltanum);
count=1;
popOut=100000;
for delta=deltarr
    pars=[beta,r1,r2,delta,mu0,mu1,upsilon0,alpha,gamma,k1,k2,omega];
    [t,y]=ode45(@(t,y)COVIDVaccineRHS(t,y,pars),[0,T], init,opts);
    %collect compliant
    compl=sum(y(:,5:8),2);    
    %collect infectious
    infectious=y(:,3)+y(:,7);
    %collect exposed
    expos=y(:,2)+y(:,6);
    %compare rates of infection and compliance spread
    figure(fig_handle);
   % subplot(1,3,1);
    %plot(t,alpha*delta*expos,'LineWidth',4);hold on;
   
    %subplot(1,3,2);
    %plot(t,beta*infectious/N,'LineWidth',4);hold on;
    %xlabel('Time, days','interpreter','latex');
    %ylabel('Infection rate');
    %subplot(1,3,3);
    h11(count)=plot(t,popOut*(expos)/N,'LineWidth',4);hold on;
    h12(count)=plot(t,popOut*compl/N,'-.','LineWidth',4);
    ylabel('Individuals (1/100,000)', 'interpreter','latex');
    set(h12(count),'Color',get(h11(count), 'color'));
    xlabel('Time, days','interpreter','latex');
    
    figure(fig_handle+1);
    h21(count)=semilogy(t,(expos)/N,'LineWidth',4);hold on;
    h22(count)=semilogy(t,compl/N,'-.','LineWidth',4);
    xlim([0,30]);
    ylabel('$$\log$$ Density', 'interpreter','latex');
    set(h22(count),'Color',get(h21(count), 'color'));
    xlabel('Time, days','interpreter','latex');
    
    ind30=find(t>=30,1);
    figure(fig_handle+2);
    h3(count)=plot(expos(1:ind30)/N,(compl(1:ind30))/N,'LineWidth',4);hold on;
    ylabel('Compliant', 'interpreter','latex');
    xlabel('Exposed','interpreter','latex');
    
    PC(count)=max(compl);
    
    count=count+1;
end
figure(fig_handle);
%legend(h11,['$$\delta=',num2str(deltarr(1),formatSpec),'$$'],['$$\delta=',num2str(deltarr(2),formatSpec),'$$'],['$$\delta=',num2str(deltarr(3),formatSpec),'$$'],'interpreter','latex');
legend(h11,['Slow  ';'Medium';'Fast  '],'interpreter','latex');
grid on;
set(gca,'FontSize',25);

figure(fig_handle+1);
legend(h21,['$$\delta=',num2str(deltarr(1),formatSpec),'$$'],['$$\delta=',num2str(deltarr(2),formatSpec),'$$'],['$$\delta=',num2str(deltarr(3),formatSpec),'$$'],'interpreter','latex');
grid on;
set(gca,'FontSize',25);

figure(fig_handle+2);
legend(h3,['$$\delta=',num2str(deltarr(1),formatSpec),'$$'],['$$\delta=',num2str(deltarr(2),formatSpec),'$$'],['$$\delta=',num2str(deltarr(3),formatSpec),'$$'],'interpreter','latex');
grid on;
set(gca,'FontSize',25);


fig_handle=fig_handle+3;
figure(fig_handle);
h1=semilogx(deltarr,PC/N,'o','MarkerSize',10);hold on;
set(h1,'MarkerFaceColor',get(h1,'color'));
deltax=linspace(1e-8,1e-5,20);
grid on
% semilogx(deltax,deltax*0+PC(1)/N,'--r');
% semilogx(deltax,deltax*0+PC(2)/N,'--r');
% semilogx(deltax,deltax*0+PC(3)/N,'--r');
ylim([0,1]);
xlabel('Logarithm of compliance rise rate, $$\log \delta$$','interpreter','latex');
ylabel({'Density of the peak'; 'of compliance population'},'interpreter','latex');
set(gca,'FontSize',25);