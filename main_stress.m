close all; clear

global model flowtype rate mode lam alpha eps G alam tauy Kfac nexp stress_imp eta_s

model = 1; % 1:UCM, 2:Giesekus, 3:PTTlin, 4:PTTexp
flowtype = 1; % 1: shear, 2: planar extension, 3: uniaxial extension
mode = 1; % current mode number
lam  = [5.0]; %
alpha = [0.1];
eps = [0.1];
G = 10000.0;
alam = 3; % 0: no adapted alam  2: SRM1 model  3: SRM2 model
eta_s = 100.0;

% if SRM1 or SRM2
tauy = [2000.0]; % yield stress

% if SRM2
Kfac = [100.0]; % consistency factor of power law
nexp = [0.5];   % shear thinning index

stress_imp = 2100; % imposed stress level

cvec = [1 0 0 1 0 1];
shearstrain = 0.0;

deltat = 1e-3;
numsteps = 1e4;

stress_all = zeros(1,numsteps);
strain_all = zeros(1,numsteps);
time = deltat*([1:numsteps]-1);

% explicit Euler scheme for first solution
for n=1:numsteps

    % initial step
    gdot1 = rate_for_stress(cvec);

    % calculate rhs
    rate = gdot1;
    rhs = rhs_viscoelastic(cvec);

    % do intermediate step
    k1 = deltat*rhs;
    gdot2 = rate_for_stress(cvec+k1);

    % calculate rhs
    rate = gdot2;
    rhs = rhs_viscoelastic(cvec+k1);

    % do step
    shearstrain = shearstrain + deltat * ( gdot1 + gdot2 ) / 2;
    cvec = cvec + ( k1 + deltat*rhs ) / 2;

    % get viscosity
    tau = stress_viscoelastic_3D(cvec);
    solventstress = eta_s * rate_for_stress(cvec);
    shearstress = solventstress + tau(2);

    stress_all(n) = shearstress;
    strain_all(n) = shearstrain;

end

figure
plot(time,strain_all)

% function to calculate rate for a given stress
function [rate] = rate_for_stress(cvec)

    global flowtype eta_s stress_imp

    tauvec = stress_viscoelastic_3D(cvec);

    % calculate the rate of deformation
    if flowtype == 1
        rate = ( stress_imp - tauvec(2) ) / eta_s;
    elseif flowtype == 2
        rate = ( stress_imp  - (tauvec(1)-tauvec(4)) ) / ( 4*eta_s );
    elseif flowtype == 3
        rate = ( stress_imp  - (tauvec(1)-tauvec(4)) ) / ( 3*eta_s );
    end

end