close all; clear

only_startup = 0;

user.model = 1; % 1:UCM, 2:Giesekus, 3:PTTlin, 4:PTTexp
user.flowtype = 1; % 1: shear, 2: planar extension, 3: uniaxial extension
user.lam  = 5.0; %
user.alpha = 0.1;
user.eps = 0.1;
user.G = 100.0;
user.alam = 3; % 0: no adapted alam  2: SRM1 model  3: SRM2 model

% if SRM1 or SRM2
user.tauy = 10.0; % yield stress

% if SRM2
user.Kfac = 100.0; % consistency factor of power law
user.nexp = 0.5;   % shear thinning index

rates = logspace(-3,2);
numsteps = 100;
deltat = 100*max(user.lam)/numsteps; % do startup phase for 100*lambda
user.rate = rates(1); % or use another rate if only_startup == 1

% check if transient similation is at first rate for the steady simulations
if only_startup == 0 && user.rate ~=rates(1)
    error('Performing steady simulations but rate ~= to rates(1)')
end

% estimate first solution from transient
c0 = [1 0 0 1 0 1];
cn = c0;
visc = zeros(1,numsteps+1);
time = deltat*([1:numsteps+1]-1);

% store the viscosity
taun = stress_viscoelastic_3D(cn,user);
visc(1) = taun(2)/user.rate;

% time stepping with 2nd-order Runge-Kutta (Heun's method)
for n=1:numsteps

    % calculate k1 in Heun's method
    k1 = rhs_viscoelastic(cn,user);

    % calculate k2 in Heun's method
    k2 = rhs_viscoelastic(cn+deltat*k1,user);

    % do step
    cnp1 = cn + deltat*(k1+k2)/2;

    % store the viscosity
    taun = stress_viscoelastic_3D(cnp1,user);
    visc(n+1) = taun(2)/user.rate;
  
    % save old values
    cn = cnp1;

end

figure;
plot(time,visc);

if only_startup == 0

    % options for fsolve
    options = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');

    c0 = cnp1; % initial guess from transient
    
    visc = zeros(1,length(rates)); % initialize to store viscosity

    for i=1:length(rates)

        % update the current rate
        user.rate = rates(i);

        % anonymous function to pass extra parameters to rhs_viscoelastic
        % https://nl.mathworks.com/help/optim/ug/passing-extra-parameters.html)
        f = @(cvec)rhs_viscoelastic(cvec,user);

        % find solution for the current rate
        sol = fsolve(f,c0,options);
        
        % store the viscosity
        visc(i) = user.G*sol(2)/user.rate;
        
        c0 = sol; % store solution als initial guess for next rate

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