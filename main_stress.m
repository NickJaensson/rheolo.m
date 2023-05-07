close all; clear

addpath('subs/')

numtimesteps1  = 40;    % number of time steps in zone 1
numtimesteps2  = 1000;  % number of time steps in zone 2
time1 = 4e-3;
deltat1 = time1/numtimesteps1;
deltat2 = 1e-3;

vemodel.model = 1;     % 1:UCM, 2:Giesekus, 3:PTTlin, 4:PTTexp
vemodel.flowtype = 1;  % 1: shear, 2: planar extension, 3: uniaxial extension
vemodel.lam  = 5.0;    % relaxation time
vemodel.alpha = 0.1;   % mobility in the Giesekus model 
vemodel.eps = 0.1;     % epsilon in the PTT model
vemodel.G = 10000.0;   % modulus
vemodel.alam = 3;      % 0: no adapted alam  2: SRM1 model  3: SRM2 model
vemodel.eta_s = 100.0; % solvent viscosity

% if SRM1 or SRM2
vemodel.tauy = 2000.0; % yield stress

% if SRM2
vemodel.Kfac = 100.0;  % consistency factor of power law
vemodel.nexp = 0.5;    % shear thinning index

vemodel.stress_imp = 2100; % imposed stress level

cvec = [1 0 0 1 0 1];
shearstrain = 0.0;

numsteps = numtimesteps1+numtimesteps2; % total number of steps

rheodata.stress_all = zeros(6,numsteps);
rheodata.strain_all = zeros(1,numsteps);
rheodata.time_all = zeros(1,numsteps);

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
    gdot1 = rate_for_stress(cvec,vemodel);
    vemodel.rate = gdot1;
    k1 = rhs_viscoelastic(cvec,vemodel);

    % calculate k2 in Heun's method
    gdot2 = rate_for_stress(cvec+k1*deltat,vemodel);
    vemodel.rate = gdot2;
    k2 = rhs_viscoelastic(cvec+k1*deltat,vemodel);

    % do step
    cvec = cvec + deltat * ( k1 + k2 ) / 2;
    shearstrain = shearstrain + deltat * ( gdot1 + gdot2 ) / 2;

    % get the stresses
    tau = stress_viscoelastic_3D(cvec,vemodel);
    vemodel.rate = rate_for_stress(cvec,vemodel);
    solventstress = stress_solvent_3D(vemodel);

    % store the solutions
    rheodata.stress_all(:,n) = tau+solventstress;
    rheodata.strain_all(n) = shearstrain;
    rheodata.time_all(n) = time;

end

rheoplot('transient_stress',rheodata,vemodel)

% function to calculate rate for a given stress
function [rate] = rate_for_stress(cvec,vemodel)

    tauvec = stress_viscoelastic_3D(cvec,vemodel);

    % calculate the rate of deformation
    if vemodel.flowtype == 1
        rate = ( vemodel.stress_imp - tauvec(2) ) / vemodel.eta_s;
    elseif vemodel.flowtype == 2
        rate = ( vemodel.stress_imp  - (tauvec(1)-tauvec(4)) ) / ( 4*vemodel.eta_s );
    elseif vemodel.flowtype == 3
        rate = ( vemodel.stress_imp  - (tauvec(1)-tauvec(4)) ) / ( 3*vemodel.eta_s );
    end

end