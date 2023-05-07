close all; clear

addpath('subs/')

only_startup = 0; % if 0: stop after performing the startup simulation
flowtype = 1; % 1: shear, 2: planar extension, 3: uniaxial extension

vemodel.model = 1; % 1:UCM, 2:Giesekus, 3:PTTlin, 4:PTTexp
vemodel.lam  = 5.0; %
vemodel.alpha = 0.1;
vemodel.eps = 0.1;
vemodel.G = 100.0;
vemodel.alam = 3; % 0:no adapted alam 1:elastic 2:SRM1 model  3:SRM2 model
vemodel.eta_s = 0.0; % solvent viscosity

% if SRM1 or SRM2
vemodel.tauy = 10.0; % yield stress

% if SRM2
vemodel.Kfac = 100.0; % consistency factor of power law
vemodel.nexp = 0.5;   % shear thinning index

rheodata.rates = logspace(-3,2);
rheodata.rate_for_startup = rheodata.rates(1); % or use another rate if only_startup == 1

% parameters for the startup problem
numsteps = 100;
deltat = 100*max(vemodel.lam)/numsteps; % do startup phase for 100*lambda

% check if transient similation is at first rate for the steady simulations
if only_startup == 0 && rheodata.rate_for_startup ~=rheodata.rates(1)
    error('Performing steady simulations but rate ~= to rates(1)')
end

% estimate first solution from transient
c0 = [1 0 0 1 0 1];
cn = c0;
rheodata.stress = zeros(6,numsteps+1);
rheodata.time = deltat*([1:numsteps+1]-1);

% store the stress
taun = stress_viscoelastic_3D(cn,vemodel);
solventstress = stress_solvent_3D(vemodel,rheodata.rate_for_startup,flowtype);
rheodata.stress(:,1) = taun+solventstress;

% time stepping with 2nd-order Runge-Kutta (Heun's method)
for n=1:numsteps

    % calculate k1 in Heun's method
    L = fill_L(vemodel,rheodata.rate_for_startup,flowtype);
    k1 = rhs_viscoelastic(cn,L,vemodel);

    % calculate k2 in Heun's method
    L = fill_L(vemodel,rheodata.rate_for_startup,flowtype);
    k2 = rhs_viscoelastic(cn+deltat*k1,L,vemodel);

    % do step
    cnp1 = cn + deltat*(k1+k2)/2;

    % store the stress
    taun = stress_viscoelastic_3D(cnp1,vemodel);
    solventstress = stress_solvent_3D(vemodel,rheodata.rate_for_startup,flowtype);
    rheodata.stress(:,n+1) = taun+solventstress;
  
    % save old values
    cn = cnp1;

end

rheoplot('startup',rheodata,vemodel);

rheodata.stress = zeros(6,length(rheodata.rates));

if only_startup == 0

    % options for fsolve
    options = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');

    c0 = cnp1; % initial guess from transient
    
    visc = zeros(1,length(rheodata.rates)); % initialize to store viscosity

    for i=1:length(rheodata.rates)

        % anonymous function to pass extra parameters to rhs_viscoelastic
        % https://nl.mathworks.com/help/optim/ug/passing-extra-parameters.html)
        L = fill_L(vemodel,rheodata.rates(i),flowtype);
        f = @(cvec)rhs_viscoelastic(cvec,L,vemodel);

        % find solution for the current rate
        cvec = fsolve(f,c0,options);
        
        % store the viscosity
        taun = stress_viscoelastic_3D(cvec,vemodel);
        solventstress = stress_solvent_3D(vemodel,rheodata.rates(i),flowtype);

        rheodata.stress(:,i) = taun+solventstress;
        
        c0 = cvec; % store solution als initial guess for next rate

    end

    rheoplot('steady',rheodata,vemodel);

    % Giesekus solution for checking
    if vemodel.model == 2 && vemodel.alam == 0
        eta = vemodel.G*vemodel.lam;
        chik = (((1+16*vemodel.alpha*(1-vemodel.alpha)*(vemodel.lam*rheodata.rates(end))^2)^(0.5) - 1) / ...
                      (8*vemodel.alpha*(1-vemodel.alpha)*(vemodel.lam*rheodata.rates(end))^2))^0.5;
        fk = (1-chik)/(1+(1-2*vemodel.alpha)*chik);
        error = abs ( (eta*(1-fk)^2)/(1+(1-2*vemodel.alpha)*fk)+vemodel.eta_s - ...
                           rheodata.stress(2,end)/rheodata.rates(end))
    end
end