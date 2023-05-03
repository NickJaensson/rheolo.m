close all; clear

numtimesteps1  = 40;    % number of time steps in zone 1
numtimesteps2  = 1000;  % number of time steps in zone 2
time1 = 4e-3;
deltat1 = time1/numtimesteps1;
deltat2 = 1e-3;

user.model = 1;     % 1:UCM, 2:Giesekus, 3:PTTlin, 4:PTTexp
user.flowtype = 1;  % 1: shear, 2: planar extension, 3: uniaxial extension
user.lam  = 5.0;    % relaxation time
user.alpha = 0.1;   % mobility in the Giesekus model 
user.eps = 0.1;     % epsilon in the PTT model
user.G = 10000.0;   % modulus
user.alam = 3;      % 0: no adapted alam  2: SRM1 model  3: SRM2 model
user.eta_s = 100.0; % solvent viscosity

% if SRM1 or SRM2
user.tauy = 2000.0; % yield stress

% if SRM2
user.Kfac = 100.0;  % consistency factor of power law
user.nexp = 0.5;    % shear thinning index

user.stress_imp = 2100; % imposed stress level

cvec = [1 0 0 1 0 1];
shearstrain = 0.0;

numsteps = numtimesteps1+numtimesteps2; % total number of steps

stress_all = zeros(1,numsteps);
strain_all = zeros(1,numsteps);
time_all = zeros(1,numsteps);

deltat = deltat1;
time = 0.0;

% time stepping with 2nd-order Runge-Kutta (Heun's method)
for n=1:numsteps

    % update the time
    time = time + deltat;

    % change the time step size
    if ( n == numtimesteps1+1 )
      deltat = deltat2;
    end

    % calculate k1 in Heun's method
    gdot1 = rate_for_stress(cvec,user);
    user.rate = gdot1;
    k1 = rhs_viscoelastic(cvec,user);

    % calculate k2 in Heun's method
    gdot2 = rate_for_stress(cvec+k1*deltat,user);
    user.rate = gdot2;
    k2 = rhs_viscoelastic(cvec+k1*deltat,user);

    % do step
    cvec = cvec + deltat * ( k1 + k2 ) / 2;
    shearstrain = shearstrain + deltat * ( gdot1 + gdot2 ) / 2;

    % get viscosity
    tau = stress_viscoelastic_3D(cvec,user);
    solventstress = user.eta_s * rate_for_stress(cvec,user);
    shearstress = solventstress + tau(2);

    % store the solutions
    stress_all(n) = shearstress;
    strain_all(n) = shearstrain;
    time_all(n) = time;

end

figure
plot(time_all,strain_all)

% function to calculate rate for a given stress
function [rate] = rate_for_stress(cvec,user)

    tauvec = stress_viscoelastic_3D(cvec,user);

    % calculate the rate of deformation
    if user.flowtype == 1
        rate = ( user.stress_imp - tauvec(2) ) / user.eta_s;
    elseif user.flowtype == 2
        rate = ( user.stress_imp  - (tauvec(1)-tauvec(4)) ) / ( 4*user.eta_s );
    elseif user.flowtype == 3
        rate = ( user.stress_imp  - (tauvec(1)-tauvec(4)) ) / ( 3*user.eta_s );
    end

end