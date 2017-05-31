function [Y_rescaled, dm_dry, sig_aer_cb] = dyn_partitioning(RHinit, ...
    n_aerosols, sz, sig_aer, Morg_in_bin, Ma1, rhoa1, nu1, t)
% Function version of the wrapper_multiple_cond_clean_trapz_logn_2016_param
% which just solves the dynamic partitioning solution at cloud base for use
% in the cloud activation parameterisation. This _clean version is a
% cleaned up faster version in which the notation matches the paper Crooks
% et al. 2016

global global_wind_speed

GRAV = 9.8;
cp = 1005;
Lv = 2.5e6;
Rv = 461.;    % specific gas constant for water vapour
Ra = 287;     % specific gas constant for dry air
sigma = 72e-3;% surface tension
R_gas = 8.314;% universal gas constant
Mw = 18./1000;% molecular weight of water
rhow = 1000.; % density of water
P = 95000;    % Initial pressure
T = 293.15;   % Initial Temperature
Mw = 18./1000;% molecular weight of water
rhow =1000; % denisty of water



% define temperature and pressure
PInit = 95000;
TInit = 283.15;

%     set the RH
RH = RHinit;





% scaling sv content to give 90% organics
% Morg_in_bin=1*[0.005 0.01 0.02 0.03 0.06 0.08 0.16 0.3 0.42 0.8];
num_org = length(Morg_in_bin);
      
     
log_10_cstar = [4-num_org:1:3];
delta_h_vap = 150.*ones(num_org);     
     
     

     num_modes = length(sz);
     
     

     
     
     
     

     









d.Morg_in_bin = Morg_in_bin;
d.RH = RH;
d.n_aerosols = n_aerosols;
d.sz = sz;
d.d_aer = sz;
d.sig_aer = sig_aer;
d.rhoa1 = rhoa1;
d.Ma1 = Ma1;
d.nu1 = nu1;
d.n_aer = n_aerosols;
d.d(5) = 95000;
d.d(6) = 293.15;
d.d(7) = RH;
d.num_org = num_org;
d.num_modes = num_modes;






        



% Calculate the moles of core aerosol
core = pi./6.*rhoa1(1:num_modes).*exp(3.*log(sz)+4.5.*sig_aer.^2).*n_aerosols;
moles_core = core./Ma1(1:num_modes).*nu1(1:num_modes); % moles / m^3

tic
[d_aer_IC,sig_aer_IC,Cij,org_frac,errors,Cond_mass_array]=...
    SOA_condensed_multiple_mode_fast_from_given_vapour(0,n_aerosols,[sz,sig_aer],...
    Morg_in_bin,rhoa1,Ma1,nu1,RH);


% if global_wind_speed == 0.01
%     vapour = Morg_in_bin - sum(Cond_mass_array*1e9);
%     save ICs d_aer_IC sig_aer_IC vapour
% end

% Calculate the moles of water
for i = 1:num_modes
    
   
    % Calculate the moles of water on a lognormal mode assuming it remains
    % a lognormal mode at equilibrium
    
    % calculate size of sz parles with water at eqm
    moles_core_med(i) = pi./6.*rhoa1(i).*sz(i).^3./Ma1(i).*nu1(i).*n_aerosols(i);
    Mwater(i) = fzero(@(x)water_eqm(x,Cond_mass_array(i,:)./1.1292,d,sz(i),...
        n_aerosols(i),0,moles_core_med(i),i),[0 1]);
    out = water_eqm(Mwater(i),Cond_mass_array(i,:)./1.1292,d,sz(i),n_aerosols(i),...
        0,moles_core_med(i),i,1);
    Kw(i) = out(1);
    dm(i) = out(2);
    

end
toc



Diam = dm;

% radius, obvs!
rad = Diam./2;

% Total number of moles in the aerosol at initial RH = RHinit%
moles_tot = Mwater./Mw+moles_core;
for i = 1:num_org
    moles_tot = moles_tot + Cond_mass_array(:,i)'./(P./Ra./T)./Ma1(num_modes+i).*nu1(num_modes+i);
end

% initial condensed mass
Cond_mass = sum(Cond_mass_array);
% mixing ratio
mix_ratio = (Morg_in_bin*1e-9-Cond_mass)./(P./Ra./T);
% mix_ratio0 is the mixing ratio when all the organics are in the condensed
% phase
% mix_ratio = Morg_in_bin*1e-9./(P./Ra./T);


% Define the Kelvin factor for the organics by extracting the exponent from
% the Kelvin term for water
if abs(sum(Ma1(num_modes+1:end))./num_org - Ma1(num_modes+1)) > 1e-4
    error('Molecular weights are not the same for each organic')
end
Kf = exp(Ma1(num_modes+1)./Mw.*log(Kw));


for j = 1:num_org;
    
    Molw_org = Ma1(j+num_modes);
    

    num_moles_org = mix_ratio*((P./Ra./T)/(1d-9));
    num_moles_org = num_moles_org*(1d-6)/(Molw_org*1d3);
    Porg = ((num_moles_org*8.314*T)); %in Pa
    
    %now need to have the saturation vapour pressure for each compound
    
    % k should be between 1 and 10 (NB: the input is i-1 and i runs from 2 to 11)
    c_star = 10.^(log_10_cstar(j));
    P_i_o_atm = (c_star.*298.15*8.2057d-5)/((Molw_org*1d3)*1d6);
    SVP_org = P_i_o_atm*(exp(delta_h_vap(j)/8.314d-3*(1/(298.15)-1/(T))))*1d5;
    RH_new(j) = Porg(j)/SVP_org;
    p0(j) = Porg(j);
    
    
    % Don't use dd as a function: Quicker...
    dd = 2.11e-5.*(T./273.15).^1.94.*(101325./P);
    
    alpha(j) = 2.*pi.*dd.*Ma1(j+num_modes)./(R_gas.*T);
    
    
end

k = p0./mix_ratio;
D0 = Diam;
C0 = moles_core;
Cw = Mwater;











     

    


                      


% humidity at cloud base
RH = 1;


% [d_aer_shifted,sig_aer_shifted,Cij2,org_frac2,errors2,Cond_mass_array2]=...
%     SOA_condensed_multiple_mode_fast_from_given_vapour(2,n_aerosols,[sz,sig_aer],...
%     Morg_in_bin,rhoa1,Ma1,nu1,RH);
% 
% 
% mass_inf = Cond_mass_array2;
% Cond_mass_array1 = Cond_mass_array2;
% 
% moles_core = moles_core_med./n_aerosols;
 
% -------------------------------------------------------------------------
% Calculate the bulk solution
% -------------------------------------------------------------------------





gamma = (dm*n_aerosols');

% Calculate the bulk solution at the initial time, ie F_i in the notes and
% at cloud base (which isn't necessary, its just nice to have it there to
% check)
Fi = sum(Cond_mass_array);
Fi(2,:) = Fi+p0./k.*(1-exp(-alpha.*gamma.*k.*t));


% -------------------------------------------------------------------------



Y=zeros(num_modes,num_org);

% -------------------------------------------------------------------------
% Now calculate the individual parle solutions
% -------------------------------------------------------------------------

for i=1:num_org
    for j=1:num_modes
        
            Y(j,i) = Cond_mass_array(j,i)+n_aerosols(j)*dm(j)*p0(i)./gamma./k(i)*(1-exp(-alpha(i)*gamma*k(i)*t));
        
    end
end

% -------------------------------------------------------------------------





% volume of core
dry_volume = pi./6.*n_aerosols'.*sz'.^3.*exp(4.5*sig_aer'.^2);

% Rescale individual condensed masses so that all SVOCs are in the 
% condensed phase at cb. ie we use Y to indicate the proportional split 
% across the two modes rather than taking it as the actual condensed mass

% [d_aer_IC1,sig_aer_IC1,Cij1,org_frac1,errors1,CB_Cond_mass]=...
%     SOA_condensed_multiple_mode_fast_from_given_vapour(2,n_aerosols,[sz,sig_aer],...
%     Morg_in_bin,rhoa1,Ma1,nu1,0.999999);
% CB_Cond_mass = sum(CB_Cond_mass);
CB_Cond_mass = Morg_in_bin*1e-9;

for i = 1:num_org
 Y_rescaled(:,i) = 0*Cond_mass_array(:,i)+(Y(:,i)-0*Cond_mass_array(:,i))...
     ./sum(Y(:,i)-0*Cond_mass_array(:,i))...
     .*(CB_Cond_mass(i));
%  Y_rescaled(:,i) = Cond_mass_array(:,i)+(Y(:,i)-Cond_mass_array(:,i))...
%      ./sum(Y(:,i)-Cond_mass_array(:,i))...
%      .*(CB_Cond_mass(i)-sum(Cond_mass_array(:,i)));
end


% Calculate volume at CB with the SVOCs but no water
for i = 1:num_org
    dry_volume = dry_volume + Y_rescaled(:,i)./rhoa1(num_modes+i);
end

% Initial arithme SD
SD = d_aer_IC.*sqrt(exp(sig_aer_IC.^2)).*sqrt(exp(sig_aer_IC.^2)-1);

% Calculate new geomer standard deviation


for j = 1:num_modes;
    sig_aer_cb(j)=fzero(@(x)shift_sig(x,dry_volume(j,end),SD(j),n_aerosols(j)),sig_aer_IC(j));
end

% New median diameter
dm_dry =(dry_volume(:,end)'*6/pi./n_aerosols.*exp(-4.5*sig_aer_cb.^2)).^(1/3);

Y_rescaled = reshape(Y_rescaled,num_modes,num_org);



core_mass = pi./6.*n_aerosols'.*sz'.^3.*exp(4.5*sig_aer'.^2).*rhoa1(1:num_modes)';
sum(Y_rescaled,2)./(sum(Y_rescaled,2)+core_mass);

Morg_in_bin;
sum(sum(Y_rescaled))./(sum(sum(Y_rescaled))+sum(core_mass));





end



function f=func(s,x)

% calculates gamma(a+1)*sum_k (b^k/gamma(a+2+k))

k = 0;
t = zeros(100,length(x));
t(1,:) = x.^(k)./prod(s:1:s+k);

for k=1:length(t(:,1))
    if abs(t(k))>1e-12
%         t(k+1,:) = x.^(k)./prod(s:1:s+k);
        bottom = [s+1:1:s+k]';
        top = x+0*bottom;
        t(k+1,:) = 1./s.*prod(top./bottom,1);
    end
end

t(isnan(t/inf))=0;
f=sum(t,1);

end

function F=water_eqm(Mwater,Cond_mass_array,d,d_aer,n_aer,sig_aer,moles_core,i,varargin)
% This function is used to iterates to find the number of moles of water on
% an aerosol parle with a given core and condensed mass of organics.
% This is very costly and is only used to check the approximate-iterative
% scheme implemented in func.

% Material parameters
Ra=287;
sigma=72e-3;
R_gas=8.314;
Mw=18./1000;
rhow=1000.;

% Extract variables from d
% d_aer = d.d_aer;
% n_aer = d.n_aer;
PInit = d.d(5);
TInit = d.d(6);
RHInit = d.d(7);
% moles_core = d.moles_core;
num_org = d.num_org;
num_modes = d.num_modes;
Ma1 = d.Ma1;
rhoa1 = d.rhoa1;
nu1 = d.nu1;
% sig_aer = d.sig_aer;


% Aerosol mass (in moles per m^3)
COA = moles_core + sum(Cond_mass_array./Ma1(num_modes+1:end).*nu1(num_modes+1:end));


% Mass of wet aerosol with condensed organics and water (in kg per m^3).
% Mass of water used is from previous iteration
mass_wet_aerosol = moles_core.*Ma1(i)./nu1(i)...
    +Mwater.*Mw + sum(Cond_mass_array);

% Volume of wet aerosol
Vol_wet_aerosol = moles_core.*Ma1(i)./nu1(i)./rhoa1(i)...
    +Mwater.*Mw./rhow +sum(Cond_mass_array./rhoa1(num_modes+1:end));

% Diameter of wet aerosol assuming sig_aer does not change
dm_wet = (Vol_wet_aerosol*6/pi./n_aer.*exp(-4.5*sig_aer.^2)).^(1/3);

% Kelvin factor for water
Kelvin_factor_water = exp(4.*Mw./(mass_wet_aerosol./Vol_wet_aerosol).*sigma./(R_gas.*TInit.*dm_wet));





% Calculate new mass of water
condensed_water_moles = RHInit.*COA./(Kelvin_factor_water-RHInit);


% Want this to be zero...
F = condensed_water_moles-Mwater;

% if we add another input variable (ie some sort f flag) then output the
% Kelvin factors instead
if nargin > 8
    F = [Kelvin_factor_water,dm_wet];
end

end

function f=shift_sig(sig_aer0,Vol_wet_aerosol,sd,n_aer)
% Shifts geometric standard deviation while keeping the arithme sd the
% same

% Calculate wet diameter from the given total volume
dm_wet=(Vol_wet_aerosol*6/pi./n_aer.*exp(-4.5*sig_aer0.^2)).^(1/3);

% iterate until the arithme standard deviation with sig_aer0 and dm_wet
% is the same as sd
f=sd-exp(log(dm_wet)+0.5.*sig_aer0.^2).*sqrt(exp(sig_aer0.^2)-1);


end

function f = trapz2(x,y)

f = [0;tril(ones(length(x)-1))*((y(2:end)+y(1:end-1))/2.*(x(2:end)-x(1:end-1)))']';
    
end



