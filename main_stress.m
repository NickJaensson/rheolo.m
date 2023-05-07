close all; clear

addpath('subs/')

flowtype = 1;  % 1: shear, 2: planar extension, 3: uniaxial extension

numtimesteps1  = 40;    % number of time steps in zone 1
numtimesteps2  = 1000;  % number of time steps in zone 2
time1 = 4e-3;
deltat1 = time1/numtimesteps1;
deltat2 = 1e-3;

vemodel.model = 1;     % 1:UCM, 2:Giesekus, 3:PTTlin, 4:PTTexp
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

rheodata.stress_imp = 2100; % imposed stress level

cvec = [1 0 0 1 0 1];
shearstrain = 0.0;

numsteps = numtimesteps1+numtimesteps2; % total number of steps

rheodata.stress = zeros(6,numsteps);
rheodata.strain = zeros(1,numsteps);
rheodata.time = zeros(1,numsteps);

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
    gdot1 = rate_for_stress(cvec,vemodel,rheodata,flowtype);
    L = fill_L(vemodel,gdot1,flowtype);
    k1 = rhs_viscoelastic(cvec,L,vemodel);

    % calculate k2 in Heun's method
    gdot2 = rate_for_stress(cvec+k1*deltat,vemodel,rheodata,flowtype);
    L = fill_L(vemodel,gdot1,flowtype);
    k2 = rhs_viscoelastic(cvec+k1*deltat,L,vemodel);

    % do step
    cvec = cvec + deltat * ( k1 + k2 ) / 2;
    shearstrain = shearstrain + deltat * ( gdot1 + gdot2 ) / 2;

    % get the stresses
    tau = stress_viscoelastic_3D(cvec,vemodel);
    gdotnp1 = rate_for_stress(cvec,vemodel,rheodata,flowtype);
    solventstress = stress_solvent_3D(vemodel,gdotnp1,flowtype);

    % store the solutions
    rheodata.stress(:,n) = tau+solventstress;
    rheodata.strain(n) = shearstrain;
    rheodata.time(n) = time;

end

rheoplot('startup_stress',rheodata,vemodel,flowtype)

% function to calculate rate for a given stress
function [rate] = rate_for_stress(cvec,vemodel,rheodata,flowtype)

    tauvec = stress_viscoelastic_3D(cvec,vemodel);

    % calculate the rate of deformation
    if flowtype == 1
        rate = ( rheodata.stress_imp - tauvec(2) ) / vemodel.eta_s;
    elseif flowtype == 2
        rate = ( rheodata.stress_imp  - (tauvec(1)-tauvec(4)) ) / ( 4*vemodel.eta_s );
    elseif flowtype == 3
        rate = ( rheodata.stress_imp  - (tauvec(1)-tauvec(4)) ) / ( 3*vemodel.eta_s );
    end

end