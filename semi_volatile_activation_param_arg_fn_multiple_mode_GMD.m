function [act_frac,Smax] = semi_volatile_activation_param_arg_fn(varargin)
% attempt to parameterise activation of ccn with SVOCS using a lognormal 
% distribution and a value of cstar. Multiple aerosol particles can be used
% and the condensed mass of SVOCs are calculated using dynamic partitioning

% varargin
% 1: w, wind speed
% 2: n_aer, number concentration of aerosol (non-volatile), vector 
%       - component for each non-volatile mode, total number /kg 
% 3: flag, 0 for no co-condensation, 1 to include co-condensation
% 4: d_aer, median diameter of core aerosol
% 5: sig_aer, standard deviation of aerosol
% 6: org_content, organic mass at time t=0
% 7: Morg_in_bin, organic mass at cloud base
% 8: org_frac
% 9: rhoa
% 10:Ma
% 11:nu

global RHinit;
global Tinit;
global Pinit;
global EPS;
global kappa;
global Molw_org;
global Rorg;
global density_dummy;
global mass_dummy;
global n_aer;
global sig_aer;
global d_aer;
global sigma;
global Mair;
global Ra;
global Rv;
global cp;
global R_gas;
global GRAV;
global Lv;
global dd;
global ka;
global Mw;
global num_org;

density_dummy=[];
mass_dummy=[];
d_aer=[];



% For F-N scheme:
% Need Tcb, Pcb, ... A, B, sgi, w:
global Tcb;global Pcb;global A;global B;global sgi;global w;global rhow;
% Thermodynamic constants.
Ra=287;
Rv=461;
cp=1005;
R_gas=8.314;

EPS=Ra./Rv;
kappa=Ra./cp;
Mw=18e-3;
rhow=1000.;
sigma=72e-3;
Mair=29e-3;
Lv=2.5e6;

dd=inline('2.11e-5.*(T./273.15).^1.94.*(101325./P)','T','P'); % Diffusivity of water vapour in air
% ka=inline('(4.40+0.071.*T).*1e-3','T');
ka=inline('(5.69+0.017.*(T-273.15)).*1e-3.*4.187','T'); % thermal conductivity of water vapour.

% other constants
GRAV=9.8;

% dummy index
i=1;

%initial conditions
RHinit=0.95; %
Pinit=95000; % Pascals
Tinit=293.15; % K
% semi-volatiles: abundance and c* and van hoft factor.
log_c_star=[-6:1:3];
% initial organic content kg/kg

org_content=[0.005 0.01 0.02 0.03 0.06 0.08 0.16 0.3 0.42 0.8].*1e-9./(Pinit./Ra./Tinit);
% org_content=[];

if(nargin>5)
    org_content=varargin{6}./(Pinit./Ra./Tinit);% kg/kg - the amount of 
        % organic condensed at t=0
    Morg_in_bin=varargin{7}./(Pinit./Ra./Tinit); % kg/kg - the total mass
        % of SVOCs
end

num_org = length(org_content);

%--------------------------

if(nargin>0)
    w=varargin{1}; % updraft speed
else
    w=0.5;
end
semi_volatile_flag=1; % If 1 then consider semi volatiles, if 0 then do not.
if(nargin>2)
    semi_volatile_flag=varargin{3};
end
method1 = 2; % method flag for activation: 1=ARG; 2=F-N.

%initial aerosol size distribution parameters
n_aer=1000e6; % total number - #/kg
num_modes=1;
if(nargin>1)
    n_aer=varargin{2};
    num_modes=length(n_aer);
end

org_frac = 0.9*ones(1,num_modes);






% Deletes organic content data if this is provided but co-condensation is
% switched off
if semi_volatile_flag==0
    org_content=0*org_content;
end

nu = [3*ones(1,num_modes) ones(1,length(org_content))];
Molw_org = 200./1000.*ones(1,length(org_content));
Ma = [132./1000*ones(1,num_modes) Molw_org];
Rorg = R_gas./Molw_org; % specific gas constant for organic

if(nargin>8)
    rhoa = varargin{9};
end

if(nargin>9)
    Ma = varargin{10};
end
if(nargin>10)
    nu = varargin{11};
end
if(nargin>7)
    t = varargin{8};
end


sig_aer=0.5; % standard dev - care! defined for dN/dr in the paper.
d_aer=60e-9; % median diameter
if(nargin>3)
    d_aer = varargin{4};
end
if(nargin>4)
    sig_aer = varargin{5};
end


density_init = rhoa(1:num_modes);
mass_initial = pi./6.*n_aer.*d_aer.^3.*exp(4.5*sig_aer.^2)...
    .*rhoa(1:num_modes);

% steps are as follows
% (1) Calculate the temperature and pressure of cloud-base
% (2) Calculate the mass of organic condensed
% (3) Calculate the new sig_aer and d_aer parameters for the lognormal - so
% that number and mass are conserved.
% (4) Calculate the size vs critical supersaturation - hence number vs peak
% supersat.
% (5) Use Twomey-type approach to get number vs updraft.

% (1) Calculate the T&P of cloud base++++++++++++++++++++++++++++++++++++++
%         cwaitbar([4 j./length(w_expli)]);
%         save /tmp/co_cond.mat d1 d2 w_expli n_aerosols;
Pcb = fzero(@dry_potential,Pinit);
Tcb = Tinit.*(Pcb./Pinit).^kappa; % dry potential temp.
%--------------------------------------------------------------------------






% This is the mass of condensed SVOCs in the core. We assume there is none.
mass_org_init = 0*org_content;


% organic mass at t=0, ie at RHinit% humidity
mass_org_condensed = org_content;
%--------------------------------------------------------------------------





% do the first stage - at t=0 (ie RHinit% rel. hum.), do not shift sig_aer!
if(semi_volatile_flag == 1)
    % Therefore the aerosol mass at the point of activation:
    % Total mass of aerosol
%     mass_final = mass_initial+sum(mass_org_condensed,2)';
%     % And the density of the aerosol through mixing rule...
%     % Total volume
%     volume_final = mass_initial./density_core+sum(mass_org_condensed./density_org,2)';
%     % Mass/volume - with organics
%     density_final = mass_final./volume_final;
% 
%     % Now calculate what d_aer needs to be
%     density_dummy = density_final;
%     mass_dummy = mass_final;
% 
% 
%     % This should be the same as the value of sz in wrapper03 and
%     % Sweet_wrapper since we are adding back on the condensed organics that 
%     % SOA_condensed in wrapper scripts "removes":
%     d_aer_final = (mass_dummy.*6./pi./density_dummy./n_aer./exp(4.5*sig_aer.^2)).^(1/3);
%     d_aer_final_t_0 = d_aer_final;
    
    
    
    
    % now do the second stage - condensed at cloud base - shift sig_aer!
    % NB: SOA_condensed outputs mass per unit volume
    % If Morg_in_bin = org_content1 then this just recalculates the diameter
    % of aerosol with organics on and outputs mass_org_condensed

%     [red1,red2,red3,red4,red5,mass_org_condensed_cb] = ...
%         SOA_condensed_multiple_mode_fast_from_given_vapour(3,n_aer,...
%         [d_aer,sig_aer],Morg_in_bin.*(Pinit./Ra./Tinit),...
%         [density_core,density_org(1,:)],Ma,nu,0.99999);

% New dynamic partitioning model for multiple modes. Note that this also
% calculates the median diameter and standard deviation assuming constant
% arithmetic SD.
if sum(Morg_in_bin) > 0
    % find condensed mass at cloud base if there are SVOCs to condense

[mass_org_condensed_cb, d_aer_final, sig_aer] = dynamic_partitioning_param_GMD...
    (RHinit,n_aer, d_aer, sig_aer, Morg_in_bin.*(Pinit./Ra./Tinit),...
    Ma,rhoa,nu,t);


% sig_aer = sig_aer*1.5;
% d_aer_final = d_aer_final*0.7;
else
    % use initial values of core aerosol if there is no co-condensation
    mass_org_condensed_cb = zeros(num_modes,num_org);
    d_aer_final = d_aer;
end

    % mass per m^3?:
    mass_org_condensed_cb = mass_org_condensed_cb./(Pinit./Ra./Tinit);

    % If the condensed mass exceeds the total amount (due to slight
    % rounding errors) then rescale the condensed masses so that the total
    % is equal to the total abundance of the SVOC
    for i = 1:num_org
        if sum(mass_org_condensed_cb(:,i)) > Morg_in_bin(i)
            mass_org_condensed_cb(:,i) = mass_org_condensed_cb(:,i)...
                ./sum(mass_org_condensed_cb(:,i)).*Morg_in_bin(i);
        end
    end
    
    mass_final = mass_initial+sum(mass_org_condensed_cb,2)';
    for i = 1:num_modes
        volume_final(i) = mass_initial(i)./rhoa(i)+sum(mass_org_condensed_cb(i,:)...
            ./rhoa(num_modes+1:end),2);
    end
    density_final = mass_final./volume_final;
    
    
    % Now calculate what d_aer needs to be
    density_dummy = density_final; % density at cloud base
    mass_dummy = mass_final; % mass at cloud base
    cc = 1; % flag for changing standard deviation.
    d_aer = d_aer_final;
    
    % The arithmetic (= constant) standard deviation at the start.
%     sd=exp(log(d_aer)+0.5.*sig_aer.^2).*sqrt(exp(sig_aer.^2)-1);
    


    
%     d_aer_final
    

    % Shift the geometric standard deviation
%     for i = 1:num_modes
%         sig_aer(i) = fzero(@(x)shift_sig(x,volume_final(i),sd(i),n_aer(i)),sig_aer);
%     end
%     
%     % new median diameter
%     d_aer_final = (volume_final*6/pi./n_aer.*exp(-4.5*sig_aer.^2)).^(1/3);
    

else
   mass_org_condensed = zeros(num_modes,length(density_org));
   mass_org_init = zeros(num_modes,length(density_org));
   mass_final = mass_initial;
   density_final = density_core;
   d_aer_final = d_aer;
end
%--------------------------------------------------------------------------

if(method1==1)
    % Now apply Abdul-Razzak-Ghan (1998, JGR) - treat as internal mixture 
    % of sulphate and organics+++++++++++++++++++++++++++++++++++++++++++++
    % Solve dS/dt = 0. for maximum supersaturation.
    tau=sigma; % hang on, shouldn't tau be the surface tension?
    A=2.*tau.*Mw./(rhow.*R_gas.*Tcb); % eq. 
    B=zeros(size(n_aer));
    % Want a different B for each aerosol
%     B=([density_final]./(sum([Ma(1:num_modes)',repmat(Ma(num_modes+1:end),num_modes,1)]...
%         .*[mass_initial' mass_org_condensed]...
%         ./[nu(1:num_modes)',repmat(nu(num_modes+1:end),num_modes,1)],2)' ...
%         ./mass_final))./(rhow./Mw);
B=([density_final].*(sum([mass_initial' mass_org_condensed_cb]...
        .*[nu(1:num_modes)',repmat(nu(num_modes+1:end),num_modes,1)]...
        ./[Ma(1:num_modes)',repmat(Ma(num_modes+1:end),num_modes,1)],2)'...
        ./mass_final))./(rhow./Mw);
        
    % B=([density_final]./(2.*132./1000))./(rhow./Mw);  
    % eq. 6

    % eq. 8:
    Sm=2./sqrt(B).*(A./(3.*d_aer_final./2)).^1.5;

    alpha_sup=GRAV.*Mw.*Lv./(cp.*R_gas.*Tcb.^2)-GRAV.*Mair./(R_gas.*Tcb); % eq. 11
    sigma_sup=R_gas.*Tcb./(svp(Tcb,'buck2','liq').*Mw) +...
        Mw.*Lv.^2./(cp.*Pcb.*Mair.*Tcb); % eq. 12

    % Equation 16:
    G=rhow.*R_gas.*Tcb./(svp(Tcb,'buck2','liq').*dd(Tcb,Pcb).*Mw) + ...
        Lv.*rhow./(ka(Tcb).*Tcb).*(Lv.*Mw./(R_gas.*Tcb)-1);
    G=1./G;

    % Equation 22:
    eta=(alpha_sup.*w./G).^1.5./(2.*pi.*rhow.*sigma_sup.*n_aer);

    % Equation 23:
    chi=(2./3).*(alpha_sup.*w./G).^0.5.*A;

    % Equation 28:
%     f1=1.5.*exp(2.25.*sig_aer.^2);
    f1=0.5.*exp(2.5.*sig_aer.^2); % eq 7 of 2000 paper.
    % Equation 29:
    f2=1+0.25.*sig_aer;


%     % Fraction activated (eq. 30):
%     act_frac=...
%         0.5.*erfc(log(f1.*(chi./eta).^1.5+f2.*(Sm.^2./(eta+3.*chi)).^0.75)./...
%         (3.*sqrt(2).*sig_aer));
% 
    % maximum supersaturation (eq. 31):
    Smax=(f1.*(chi./eta).^1.5+f2.*(Sm.^2./eta).^0.75).^0.5;
    Smax=sum(Sm./Smax);
    
    % eq 6 of 2000 paper
    Smax=sum(1./Sm.^2.*(f1.*(chi./eta).^1.5+f2.*(Sm.^2./(eta+3.*chi)).^0.75)).^0.5;
    Smax=1./Smax;
    
    % eq 13 of 2000 paper
    act_frac=1./sum(n_aer).*sum(n_aer.*0.5.*(1.-erf(2.*log(Sm./Smax)./(3.*sqrt(2).*sig_aer) )));
    
    %----------------------------------------------------------------------
elseif(method1 == 2) % Fountoukis and Nenes (2005)
    tau = sigma; % hang on, shouldn't tau be the surface tension?
    A = 4.*tau.*Mw./(rhow.*R_gas.*Tcb); % there seem to be typos in F-N paper
    
%     B = ([density_final]./(sum([Ma(1:num_modes)',repmat(Ma(num_modes+1:end),num_modes,1)]...
%         .*[mass_initial' mass_org_condensed_cb]...
%         ./[nu(1:num_modes)',repmat(nu(num_modes+1:end),num_modes,1)],2)' ...
%         ./mass_final))./(rhow./Mw)
    
    B=([density_final].*(sum([mass_initial' mass_org_condensed_cb]...
        .*[nu(1:num_modes)',repmat(nu(num_modes+1:end),num_modes,1)]...
        ./[Ma(1:num_modes)',repmat(Ma(num_modes+1:end),num_modes,1)],2)'...
        ./mass_final))./(rhow./Mw);
    
    % B=([density_final]./(2.*132./1000))./(rhow./Mw);  

    sgi = sqrt(4.*A.^3./27./(B.*d_aer_final.^3)); % eq 17

    % calculate s_max by iteration
    Smax = max(fzero(@fountoukis_nenes,0.01),1e-20);
    
    % now we know Smax, calculate the activated fraction - eq 9 - note we
    % can have multiple lognormals
%     act_frac=0.5.*erfc(2.*log(sgi./Smax)./(3.*sqrt(2).*sig_aer)); % eq 8 and 9
    % eq 13 of 2000 paper
    act_frac = 1./sum(n_aer).*sum(n_aer.*0.5.*(1.-erf(2.*log(sgi./Smax)./(3.*sqrt(2).*sig_aer) )));
erf(2.*log(sgi./Smax)./(3.*sqrt(2).*sig_aer) );
end


function delta=fountoukis_nenes(Smax)
global EPS;
global kappa;
global Molw_org;
global Rorg;
global density_dummy;
global mass_dummy;
global n_aer;
global sig_aer;
global sigma;
global Mair;
global Ra;
global Rv;
global cp;
global R_gas;
global GRAV;
global Lv;
global dd;
global ka;
global Mw;

% Need Tcb, Pcb, ... A, B, sgi, w, rhow:
global Tcb;global Pcb;global A;global B;global sgi;global w;global rhow;

Smax=max(Smax,1e-20);

% This is just copied from ARG code...
a=GRAV.*Mw.*Lv./(cp.*R_gas.*Tcb.^2)-GRAV.*Mair./(R_gas.*Tcb); % eq. 11
gamma=R_gas.*Tcb./(svp(Tcb,'buck2','liq').*Mw) +...
    Mw.*Lv.^2./(cp.*Pcb.*Mair.*Tcb); % eq. 11

% Equation 12:
G=rhow.*R_gas.*Tcb./(svp(Tcb,'buck2','liq').*dd(Tcb,Pcb).*Mw) + ...
    Lv.*rhow./(ka(Tcb).*Tcb).*(Lv.*Mw./(R_gas.*Tcb)-1);
G=4./G;

% now calculate del - used for s_part - section 3.4
% del=Smax.^4-16.*A.^2.*a.*w./9./G;
del=1-16.*A.^2.*a.*w./(9.*G)./Smax.^4;
if(del>=0) % discriminant criterion
    s_part=Smax.*(0.5.*(1+(1-16.*A.^2.*a.*w./(9.*G.*Smax.^4)).^0.5)).^0.5;
else
    xic=(16.*A.^2.*a.*w./(9.*G)).^(0.25);
    s_part=Smax.*min((2e7.*A.*Smax.^-0.3824)./3,1.);
    s_part=Smax.*min((2e7.*A.*(Smax.^-0.3824-xic.^-0.3824))./3+1./sqrt(2),1.);
end
I=0.;

for i=1:length(sig_aer)
    % now calculate the integrals - eqs 18-20 - not sure how this can be done
    % for more than one mode...
    % equation 20:
    u_part=2.*log(sgi(i)./s_part)./(3.*sqrt(2).*sig_aer(i));
    UPART(i)=sgi(i)./s_part;
    u_max=2.*log(sgi(i)./Smax)./(3.*sqrt(2).*sig_aer(i));
    UMAX(i)=sgi(i)./Smax;
    
    % Equation 18:
    I1=0.5.*n_aer(i).*sqrt(G./a./w).*Smax.*...
        (erfc(u_part)-0.5.*(sgi(i)./Smax).^2.*exp(9.*sig_aer(i).^2./2).*...
        erfc(u_part+3.*sig_aer(i)./sqrt(2))) ;
    

    % Equation 19:
    I2=A.*n_aer(i)./3./sgi(i).*exp(9.*sig_aer(i).^2./8).*...
        (erf(u_part-3.*sig_aer(i)./(2.*sqrt(2)))-erf(u_max-3.*sig_aer(i)./(2.*sqrt(2))));

    % Note I don't believe the analytical integrations above
    % c=[(G./a./w).^0.5.*2.*n_aer./(3.*sqrt(2.*pi).*sig_aer) sgi 2.*sig_aer.^2 Smax];
    % I1=quad(@(x)integral1_fn(x,c),0.,s_part,1e-10);
    % c=[2.*A./3.*2.*n_aer./(3.*sqrt(2.*pi).*sig_aer) sgi 2.*sig_aer.^2];
    % I2=quad(@(x)integral2_fn(x,c),s_part,Smax,1e-10);
    % I=I1+I2; % equation 14
    
    %GCCN
    I3=1./sqrt(3).*A.*n_aer(i)./3./sgi(i).*exp(9.*sig_aer(i).^2./8).*...
        (1-erf(u_part-3.*sig_aer(i)./(2.*sqrt(2))));

    c=[2.*A./3 (G./a./w) 2.*n_aer(i)./(3.*sqrt(2.*pi).*sig_aer(i)) sgi(i) 2.*sig_aer(i).^2 Smax];
    I=I+I1+I2+I3;
    % I=I+quad(@(x)integral3_fn(x,c),0.,Smax,1e-3)
end
UPART;
UMAX;

% now cost function:
delta=2.*a.*w./(pi.*gamma.*rhow)-G.*Smax.*I; % equation to be minimised - eq 10

function I1=integral1_fn(s,c) % equation 15
I1=c(1).*(c(4).^2-s.^2).^0.5./s.*exp(-log((c(2)./s).^(2./3)).^2./(c(3)));

function I2=integral2_fn(s,c) % equation 16
I2=c(1).*(1./s).^2.*exp(-log((c(2)./s).^(2./3)).^2./(c(3)));

function I3=integral3_fn(s,c) % equation 14
I3=((c(1)./s).^2+c(2).*(c(6).^2-s.^2)).^0.5.*c(3)./s.*exp(-log((c(4)./s).^(2./3)).^2./(c(5)));

function x=dry_potential(P)
global RHinit;
global Tinit;
global Pinit;
global EPS;
global kappa;

total_water1=RHinit.*EPS.*svp(Tinit,'buck2','liq')./(Pinit-svp(Tinit,'buck2','liq'));

Tcalc=Tinit.*(P./Pinit).^kappa;

total_water2=EPS.*svp(Tcalc,'buck2','liq')./(P-svp(Tcalc,'buck2','liq'));

x=total_water2-total_water1;


function svp_org=vapour_pressure_org(T,log_c_star)
global Molw_org;
delta_h_vap=150;
c_star=10.^log_c_star;
P_i_o_atm=(c_star.*298.15.*8.2057d-5)./((Molw_org.*1d3).*1d6);
svp_org=P_i_o_atm.*(exp(delta_h_vap./8.314d-3.*(1./(298.15)-1./(T))))*1d5;

function y=lognorm1(x,c)
% This integrates mass.
y=pi.*x.^2./sqrt(2.*pi.*c(1).^2).*exp(-log(x./c(2)).^2./(2.*c(1).^2))./6;

function mass=integrate_mass(d_aer1,cc,i)
% i is number of the mode
global density_dummy;
global n_aer;
global sig_aer;
global d_aer;
global mass_dummy;

d_aer1=abs(d_aer1);
sig_aer1=sig_aer(i);
size_max=2e-6;
if(cc==1)
    % The standard deviation at the start.
    sd=exp(log(d_aer(i))+0.5.*sig_aer(i).^2).*sqrt(exp(sig_aer(i).^2)-1)
    % Root find
    c=[sd d_aer1];
    sig_aer1=fzero(@(x) find_sig_aer(x,c),sig_aer(i),optimset('TolX',1e-100,'TolFun',1e-10,'MaxFunEvals',100000));
    size_max=3e-6;
end
c=[(sig_aer1) abs(d_aer1)];
% mass=density_dummy.*n_aer.*quad(@(x)lognorm1(x,c(1:2)),0.,size_max,1e-30)-mass_dummy;
mass=density_dummy(i).*n_aer(i).*exp(3.*log(d_aer1)+4.5.*sig_aer1.^2).*pi./6-mass_dummy(i);




%broydn
function F=find_delta(delta,sd)

global density_dummy;
global n_aer;
global d_aer;
global mass_dummy;

delta=abs(delta);

% X_i=exp(sig_aer^2)
X_i=0.5+0.5*sqrt(1+4*sd.^2./(d_aer+delta).^2);

% geometric standard deviation
sig_i=sqrt(log(X_i));

% mass of aerosol per unit volume
mass_i=density_dummy.*n_aer.*exp(3*log(d_aer+delta)+4.5*(sig_i.^2))*pi/6;

F=sum(mass_i)-sum(mass_dummy);






function F=find_sig_aer(x,c)

sd=exp(log(c(2))+0.5.*x.^2).*sqrt(exp(x.^2)-1);


F=c(1)-sd;

function f = shift_sig(sig_aer0,Vol_wet_aerosol,sd,n_aer)
% Shifts geometric standard deviation while keeping the arithmetic sd the
% same

% Calculate wet diameter from the given total volume
dm_wet = (Vol_wet_aerosol*6/pi./n_aer.*exp(-4.5*sig_aer0.^2)).^(1/3);

% iterate until the arithmetic standard deviation with sig_aer0 and dm_wet
% is the same as sd
f = sd-exp(log(dm_wet)+0.5.*sig_aer0.^2).*sqrt(exp(sig_aer0.^2)-1);









function F=org_kohler(mass_vap)
global tot_org;
global mass_d;
global Mw;
global mol_wgt;
global smr_organic;
global n_aer;
mass_cond=abs(tot_org-mass_vap)./n_aer;
ns=mass_cond./mol_wgt;

RHeq=... exp(4.*Molw_org.*sigma./R_gas./T./rhoat./Diam).*...
    (ns)./(mass_d./Mw+ns);

F=RHeq-mass_vap./smr_organic;

function pi=svp(T,flag,flag2)
%usage: SVPIce(T,flag,flag2);
% flag='marti'
% flag='teten'
% flag2='ice'
% flag2='liq'


switch lower(flag2)
case 'ice'
    switch lower(flag)
    case 'buck2'
        pi = 100.*6.1115 .* exp((23.036 - (T-273.15)./ 333.7) .*(T-273.15) ./ (279.82 + (T-273.15)))  ;
    case 'buck'
        pi = 100.*6.1115 .* exp(22.452 .* (T-273.15) ./ (272.55+(T-273)));  

    case 'goff'
        pi =  100.*10.^(-9.09718.* (273.16./T - 1)       ...                                 
                   - 3.56654 .*log10(273.16./ T) ...
                   + 0.876793 .*(1 - T./ 273.16) ...
                   + log10(6.1071) );
        
    case 'marti'
        pi = 10.^((-2663.5 ./ T) + 12.537 );   
    case 'teten'
        pi = 100.*10.^(9.5 .*(T-273.15) ./ (T-273.15+265.5) + 0.7858  ) ;  
    case 'hyland'
        pi =  exp(-0.56745359e4 ./ T     ...                                            
              + 0.63925247e1 ...
              - 0.96778430e-2 .*T ...
             + 0.62215701e-6 .*T.^2 ...
             + 0.20747825e-8 .*T.^3 ...
             - 0.94840240e-12 .*T.^4 ...
             + 0.41635019e1 .*log(T) );
    case 'murphy'
        pi = exp(9.554605 - 5722.796./T + 3.5291623.*log(T) - 0.00727374.*T);
    case 'magnus'
        pi=610.7.*exp((22.44.*T-6.1186e3)./(T-0.75));
    case 'clausius'
        pi=611.73.*exp(2.501e6./461.5.*(1./273.16-1./T));
    otherwise
        error('invalid');
    end
case 'liq'
    switch lower(flag)
    case 'goff'
        pi =  100.*10.^(-7.90298 .*(373.16./T-1)    ...                       
                    + 5.02808 .*log10(373.16./T) ...
                    - 1.3816e-7 .*(10.^(11.344 .*(1-T./373.16))  -1) ...
                   + 8.1328e-3 .*(10.^(-3.49149 .*(373.16./T-1))  -1) ...
                   + log10(1013.246) );
        
    case 'bolton'
        pi = 100.*6.112 .*exp(17.67 .* (T-273.15) ./ (T-273.15+243.5));
    case 'roger'
        pi = 2.53e11 * exp(-5.42e3./(T));
    case 'buck2'
        pi = 100.*6.1121  .*exp((18.678 - (T-273.15)./ 234.5).* (T-273.15) ./ (257.14 + (T-273.15)));
    case 'buck1'
        pi = 100.*6.1121 .*exp(17.502 .*(T-273.15)./ (240.97 + T-273.15));

    case 'wmo'
        pi = 100.*10.^( 10.79574 .*(1-273.16./T)        ...                       
                    - 5.02800 .*log10(T./273.16) ...
                    + 1.50475e-4 .*(1 - 10.*(-8.2969.*(T./273.16-1))) ...
                    + 0.42873e-3 .*(10.*(+4.76955.*(1-273.16./T)) - 1) ...
                    + 0.78614 );
    case 'hyland'
        pi =  exp(-0.58002206e4 ./ T  ...                                     
              + 0.13914993e1 ...
              - 0.48640239e-1 .* T ...
              + 0.41764768e-4 .* T.^2 ...
              - 0.14452093e-7 .* T.^3 ...
              + 0.65459673e1 .* log(T)); 
            
    case 'sonntag'
        pi =  100.*exp(-6096.9385 ./ T  ...                        
                 + 16.635794 ...
                 - 2.711193e-2 .* T ...
                 + 1.673952e-5 .* T.^2 ... 
                 + 2.433502 .* log(T)); 
    case 'teten'
        pi = 100.*10.^(7.5 .*(T-273.15) ./ (T-273.15+237.3) + 0.7858  ) ;  
    case 'clausius'
        pi=611.73.*exp(2.834e6./461.5.*(1./273.16-1./T));
    case 'magnus'
        pi=610.7.*exp((17.38.*T-4.7473e3)./(T-34.15));
    otherwise
        error('invalid');
    end
end    

