% -------------------------------------------------------------------------
% The MIT License (MIT)
% 
% Copyright (c) [2015] [Matthew Crooks]
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% This code is a wrapper script that randomises the model parameters and
% runs the multiple mode equilibrium partitioning theory code to calculate
% the equilibrium condensed masses. Plots are produced at the end comparing
% the organic mass fraction and condensed concentrations of the two
% approximations against the direct numerical calculation of the solution
% using fsolve.
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Code dependencies:
% svp.m
% SOA_condensed_multiple_mode_fast_from_given_vapour
% -------------------------------------------------------------------------

% define a matrix with 3 columns to store the 3 values of the organic
% fraction calculated from the model. 
OF=zeros(0,3);
% specify the number of datat points on the graphs
num_runs=100;

% This plots the line y=x and is used in plotting Figure 3. 
plot([0 1], [0 1],'--')

% Repeats the calculations for different parameters num_runs times
while length(OF(:,1))<num_runs
    % an index to keep track of the numebr of runs
    i=length(OF(:,1))+1;
    % randomly chose an 
    Org_frac=rand(1,1);
    
    % randomly choose a number of modes from 2 to 6
    num_modes=ceil(1+rand(1,1)*5);
    % randomly choose a number concentration for a total from 50 to 500 pcc
    n_aerosols=(50e6+450e6*rand(1,num_modes))/num_modes;
    % randomly choose a diameter for the core aerosol from 50 to 500 nm
    sz=[50+450*rand(1,num_modes)]*1e-9;
    % set geometric standard deviation to 0.5
    sig_aer=0*n_aerosols+0.5;
    % set relative humidity to a random value from 0 to 1 (x100%)
    RH=rand(1,1);
     
     
    % scaling sv content by between 0 and 100
    Morg_in_bin=100*rand(1,1)*[0.005 0.01 0.02 0.03 0.06 0.08 0.16 0.3 0.42 0.8];
    
    % Material properties of core aerosol and organics. First num_modes
    % entries are core and then following 10 are the organucs
    % density:
    rhoa=[500+1500*rand(1,length(n_aerosols)),1000+1000*rand(1,length(Morg_in_bin))];
    % molecular weight
    Ma=[0.1+0.3*rand(1,length(n_aerosols)),0.1+0.3*rand(1,length(Morg_in_bin))];
    % van't Hoff factor
    nu=[1+2*rand(1,length(n_aerosols)),rand(1,length(Morg_in_bin))];
     

    % solve equilibrium partitioning theory
    [d_aer,sig_aer_out,Cij,org_frac,errors,Cond_mass_array]=...
    SOA_condensed_multiple_mode_fast_from_given_vapour(0,n_aerosols,[sz,sig_aer],...
    Morg_in_bin,rhoa,Ma,nu,RH);

    % collate calculated organic mass fractions of the first mode of each
    % simulation. First column is the numerical solution, the second is the 
    % quasi-leading order and the third column is the perturbation solution 
    OF=[OF;errors.of(1),errors.of_lead(1),errors.of_first(1)];

   


    
    
    
    % Collate condensed concentrations of organics on the first mode. num 
    % is the numerical solution, bar is the quasi-leading order and pert is
    % the perturbation solution
    num(i,:)=Cij.num(1,:)*1e9;
    bar(i,:)=Cij.Cij_bar(1,:);
    pert(i,:)=Cij.Cij_pert(1,:);
    




end


% -------------------------------------------------------------------------
% Plots organic fraction comparisons
% -------------------------------------------------------------------------
plot(OF(:,1),OF(:,2),'ko')
hold on
plot(OF(:,1),OF(:,3),'ro')
plot([0 1], [0 1],'--')
xlim([0 1])
ylim([0 1])
axis square

% -------------------------------------------------------------------------
% Plots condensed concentrations comparisons
% -------------------------------------------------------------------------
% loglog(num,bar,'ob')
figure
% loglog(num,pert,'.r')
xlim([0 100])
ylim([0 100])

loglog([min(min(num))/100 100],[min(min(num))/100 100] ,'--')
hold on

loglog1 = loglog(num(:,1),[bar(:,1),pert(:,1)],'LineStyle','none');
set(loglog1(1),'Marker','o','Color',[0.0431372560560703 0.517647087574005 0.780392169952393]);
set(loglog1(2),'MarkerFaceColor',[0.847058832645416 0.160784319043159 0],...
    'MarkerSize',20,...
    'Marker','.',...
    'Color',[0.847058832645416 0.160784319043159 0]);



% Create multiple lines using matrix input to loglog
loglog2 = loglog(num(:,2),[bar(:,2),pert(:,2)],'Marker','square','LineStyle','none');
set(loglog2(1),'MarkerSize',5,'Color',[0.0431372560560703 0.517647087574005 0.780392169952393]);
set(loglog2(2),'MarkerSize',5,'MarkerFaceColor',[0.847058832645416 0.160784319043159 0],...
    'Color',[0.847058832645416 0.160784319043159 0]);

% Create multiple lines using matrix input to loglog
loglog3 = loglog(num(:,3),[bar(:,3),pert(:,3)],'Marker','x','LineStyle','none');
set(loglog3(1),'MarkerSize',5,'Color',[0.0431372560560703 0.517647087574005 0.780392169952393]);
set(loglog3(2),'MarkerSize',5,'MarkerFaceColor',[0.847058832645416 0.160784319043159 0],...
    'Color',[0.847058832645416 0.160784319043159 0]);

% Create multiple lines using matrix input to loglog
loglog4 = loglog(num(:,4),[bar(:,4),pert(:,4)],'Marker','*','LineStyle','none');
set(loglog4(1),'MarkerSize',5,'Color',[0.0431372560560703 0.517647087574005 0.780392169952393]);
set(loglog4(2),'MarkerSize',5,'MarkerFaceColor',[0.847058832645416 0.160784319043159 0],...
    'Color',[0.847058832645416 0.160784319043159 0]);

% Create multiple lines using matrix input to loglog
loglog5 = loglog(num(:,5),[bar(:,5),pert(:,5)],'Marker','*','LineStyle','none');
set(loglog5(1),'MarkerSize',5,'Color',[0.0431372560560703 0.517647087574005 0.780392169952393]);
set(loglog5(2),'MarkerSize',5,'MarkerFaceColor',[0.847058832645416 0.160784319043159 0],...
    'Color',[0.847058832645416 0.160784319043159 0]);

% Create multiple lines using matrix input to loglog
loglog6 = loglog(num(:,6),[bar(:,6),pert(:,6)],'Marker','diamond','LineStyle','none');
set(loglog6(1),'MarkerSize',5,'Color',[0.0431372560560703 0.517647087574005 0.780392169952393]);
set(loglog6(2),'MarkerSize',5,'MarkerFaceColor',[0.847058832645416 0.160784319043159 0],...
    'Color',[0.847058832645416 0.160784319043159 0]);

% Create multiple lines using matrix input to loglog
loglog7 = loglog(num(:,7),[bar(:,7),pert(:,7)],'Marker','hexagram','LineStyle','none');
set(loglog7(1),'MarkerSize',5,'Color',[0.0431372560560703 0.517647087574005 0.780392169952393]);
set(loglog7(2),'MarkerSize',5,'MarkerFaceColor',[0.847058832645416 0.160784319043159 0],...
    'Color',[0.847058832645416 0.160784319043159 0]);

% Create multiple lines using matrix input to loglog
loglog8 = loglog(num(:,8),[bar(:,8),pert(:,8)],'Marker','<','LineStyle','none');
set(loglog8(1),'MarkerSize',5,'Color',[0.0431372560560703 0.517647087574005 0.780392169952393]);
set(loglog8(2),'MarkerSize',5,'MarkerFaceColor',[0.847058832645416 0.160784319043159 0],...
    'Color',[0.847058832645416 0.160784319043159 0]);

% Create multiple lines using matrix input to loglog
loglog9 = loglog(num(:,9),[bar(:,9),pert(:,9)],'Marker','^','LineStyle','none');
set(loglog9(1),'MarkerSize',5,'Color',[0.0431372560560703 0.517647087574005 0.780392169952393]);
set(loglog9(2),'MarkerSize',5,'MarkerFaceColor',[0.847058832645416 0.160784319043159 0],...
    'Color',[0.847058832645416 0.160784319043159 0]);

% Create multiple lines using matrix input to loglog
loglog10 = loglog(num(:,10),[bar(:,10),pert(:,10)],'Marker','v','LineStyle','none');
set(loglog10(1),'MarkerSize',5,'Color',[0.0431372560560703 0.517647087574005 0.780392169952393]);
set(loglog10(2),'MarkerSize',5,'MarkerFaceColor',[0.847058832645416 0.160784319043159 0],...
    'Color',[0.847058832645416 0.160784319043159 0]);


xlabel('x')
ylabel('y')
axis square

