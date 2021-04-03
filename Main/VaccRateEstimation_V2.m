% this script loads the cleaned up vaccination uptake data (vaccination doses per 100 people) obtained from http://ourworldindata.org
% and calculated the vaccination rates approximating the rates of
% vaccination observed Between January 7-February, 2021 in the Netherlands and UK. These will be
% dubbed, slow and fast vaccination rates, respectively
% the estimation of rates is based on the assumption that at the start of
% the vaccination the growth rate of the vaccination population is
% V(t)=N(1-exp(-upsilon t)), which can be approximated by a line passing
% through the origin \hat{V}(t)=N upsilon t
% Since this was approximately the start of the vaccination campaign, then
% we assume that the number of vaccination doses approximates the number of
% people infected

%prepare the space
clc;
clear variables;
close all;
format long;

%total population, set to the size of the Netherlands
N=1.7e7;
%read the CSV files and load them into array
TblNed = readtable('Data/VaccNetherlands.csv');
%convert from percent to fraction and pad with (0,0) as the first point in
%the series, both the number of vaccine doses data and the time data
TVPHNed = [0;(TblNed.total_vaccinations_per_hundred)/100]; 
DateNed = [0;TblNed.Date];%
TblUK = readtable('Data/VaccUK.csv');
TVPHUK = [0;(TblUK.total_vaccinations_per_hundred)/100]; 
DateUK = [0;TblUK.date];%
TPVPHUK = [0;(TblUK.people_vaccinated_per_hundred)/100];

%Output the data points
figure(1);
h(1)=plot(DateNed,TVPHNed,'*','MarkerSize',6);hold on;
h(2)=plot(DateUK,TVPHUK,'*','MarkerSize',6);hold on;
h(3)=plot(DateUK,TPVPHUK,'o','MarkerSize',6);hold on;
set(h(3),'MarkerEdgeColor',get(h(2),'color'));
xlabel('Time, days','interpreter','latex');
ylabel('Density of vaccinated','interpreter','latex');

leg=['Netherlands vaccine doses';'UK vaccine doses         ';'UK people vaccinated     ';'Approx to NL doses data  ';'Approx to UK doses data  ';'Approx to NL peopl data  '];

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


figure(1);

xlabel('Time, days','interpreter','latex');
ylabel('Vaccinated individuals (density)','interpreter','latex');

set(gca,'FontSize',25);