close all; clear

addpath('subs/')

% general simulation parameters

flowtype = 1;  % 1: shear, 2: planar extension, 3: uniaxial extension
plottype = 'strain'; % 'strain': plot strain

numtimesteps1  = 40;   % number of time steps in zone 1
numtimesteps2  = 1000; % number of time steps in zone 2
time1 = 4e-3;          % duration of zone 1
time2 = 1.0;           % duration of zone 2

% viscoelastic model parameters

vemodel.model = 2;   % 1:UCM, 2:Giesekus, 3:PTTlin, 4:PTTexp
vemodel.alam = 2;    % 0:no adapted lambda 1:elastic 2:SRM1 model  3:SRM2 model
vemodel.lam  = 0.1;  % relaxation time
vemodel.G = 10.0;    % elastic modulus
vemodel.eta_s = 0.1; % solvent viscosity

% model specific material parameters

% if Giesekus (model == 2)
vemodel.alpha = 0.1; % alpha parameter for Giesekus

% if PTT (model == 3 or 4)
vemodel.eps = 0.1;   % epsilon parameter in PTT model

% if SRM1 or SRM2 (alam == 2 or 3)
vemodel.tauy = 10.0; % yield stress

% if SRM2 (alam == 3)
vemodel.Kfac = 10.0; % consistency factor of power law
vemodel.nexp = 0.5;  % shear thinning index

% imposed stress level
rheodata.stress_imp = 12; % stress_xy             for shear flow (flowtype 1)
                          % stress_xx - stress_yy for elongation flow (flowtype 2 & 3)

cvec = [1 0 0 1 0 1];
shearstrain = 0.0;

numsteps = numtimesteps1+numtimesteps2; % total number of steps

rheodata.stress = zeros(6,numsteps);
rheodata.strain = zeros(1,numsteps);
rheodata.rates = zeros(1,numsteps);
rheodata.time = zeros(1,numsteps);

deltat1 = time1/numtimesteps1;
deltat2 = time2/numtimesteps2;

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
    rate1 = rate_for_stress(cvec,vemodel,rheodata,flowtype);
    L = fill_L(vemodel,rate1,flowtype);
    k1 = rhs_viscoelastic(cvec,L,vemodel);

    % calculate k2 in Heun's method
    rate2 = rate_for_stress(cvec+k1*deltat,vemodel,rheodata,flowtype);
    L = fill_L(vemodel,rate1,flowtype);
    k2 = rhs_viscoelastic(cvec+k1*deltat,L,vemodel);

    % do step
    cvec = cvec + deltat * ( k1 + k2 ) / 2;
    shearstrain = shearstrain + deltat * ( rate1 + rate2 ) / 2;

    % get the stresses
    tau = stress_viscoelastic_3D(cvec,vemodel);
    ratenp1 = rate_for_stress(cvec,vemodel,rheodata,flowtype);
    solventstress = stress_solvent_3D(vemodel,ratenp1,flowtype);

    % store the solutions
    rheodata.stress(:,n) = tau+solventstress;
    rheodata.strain(n) = shearstrain;
    rheodata.rates(n) = ratenp1;
    rheodata.time(n) = time;

end

% plot the results
rheoplot('startup_stress',rheodata,vemodel,flowtype,plottype)

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