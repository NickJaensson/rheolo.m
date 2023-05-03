close all; clear

global model flowtype rate lam alpha eps G alam tauy Kfac nexp

only_startup = 0;

model = 1; % 1:UCM, 2:Giesekus, 3:PTTlin, 4:PTTexp
flowtype = 1; % 1: shear, 2: planar extension, 3: uniaxial extension
lam  = 5.0; %
alpha = 0.1;
eps = 0.1;
G = 100.0;
alam = 3; % 0: no adapted alam  2: SRM1 model  3: SRM2 model

% if SRM1 or SRM2
tauy = 10.0; % yield stress

% if SRM2
Kfac = 100.0; % consistency factor of power law
nexp = 0.5;   % shear thinning index

rates = logspace(-3,2);
numsteps = 1000;
deltat = 100*max(lam)/numsteps;
rate = rates(1); % or use anither rate if only_startup == 1

% estimate first solution from transient
c0 = [1 0 0 1 0 1];
cn = c0;
visc = zeros(1,numsteps);
time = deltat*([1:numsteps]-1);

% explicit Euler scheme for first solution
for n=1:numsteps
  dcdt = rhs_viscoelastic(cn);
  cnp1 = cn + deltat*dcdt;
  taun = stress_viscoelastic_3D(cn);
  visc(n) = taun(2)/rate;
  cn = cnp1;
end

figure;
plot(time,visc);

if only_startup == 0

    % find steady state values
    options = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');
    c0 = cnp1; % initial guess from transient
    visc = zeros(1,length(rates));
    for i=1:length(rates)
        rate = rates(i);
        sol = fsolve(@rhs_viscoelastic,c0,options);
        visc(i) = G*sol(2)/rate;
        c0 = sol;
    end
    figure
    loglog(rates,visc);

%     % Giesekus solution for checking
%     if model == 2 && alam == 0
%         eta = G*lam;
%         chik = (((1+16*alpha*(1-alpha)*(lam*rate)^2)^(0.5) - 1) / ...
%                       (8*alpha*(1-alpha)*(lam*rate)^2))^0.5;
%         fk = (1-chik)/(1+(1-2*alpha)*chik);
%         visc_an = (eta*(1-fk)^2)/(1+(1-2*alpha)*fk)
%     end
end