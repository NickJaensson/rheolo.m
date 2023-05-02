close all; clear

global model flowtype rate mode lam alpha eps G alam tauy Kfac nexp stress_imp eta_s

numtimesteps1  = 40;    % number of time steps in zone 1
numtimesteps2  = 1000;  % number of time steps in zone 2
time1 = 4e-3;
deltat1 = time1/numtimesteps1;
deltat2 = 1e-3;

model = 1;     % 1:UCM, 2:Giesekus, 3:PTTlin, 4:PTTexp
flowtype = 1;  % 1: shear, 2: planar extension, 3: uniaxial extension
mode = 1;      % current mode number
lam  = [5.0];  % relaxation time
alpha = [0.1]; % mobility in the Giesekus model 
eps = [0.1];   % epsilon in the PTT model
G = 10000.0;   % modulus
alam = 3;      % 0: no adapted alam  2: SRM1 model  3: SRM2 model
eta_s = 100.0; % solvent viscosity

% if SRM1 or SRM2
tauy = [2000.0]; % yield stress

% if SRM2
Kfac = [100.0]; % consistency factor of power law
nexp = [0.5];   % shear thinning index

stress_imp = 2100; % imposed stress level

cvec = [1 0 0 1 0 1];
shearstrain = 0.0;

numsteps = numtimesteps1+numtimesteps2; % total number of steps

stress_all = zeros(1,numsteps);
strain_all = zeros(1,numsteps);
time_all = zeros(1,numsteps);

deltat = deltat1;
time = 0.0;

% explicit Euler scheme for first solution
for n=1:numsteps

    % update the time
    time = time + deltat;

    % change the time step size
    if ( n == numtimesteps1+1 )
      deltat = deltat2;
    end

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

    % store the solutions
    stress_all(n) = shearstress;
    strain_all(n) = shearstrain;
    time_all(n) = time;

end

figure
plot(time_all,strain_all)

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