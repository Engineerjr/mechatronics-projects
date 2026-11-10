%% MAIN_CANTILEVER  Topology optimization of a 2D cantilever beam (SIMP + OC)
% Run this script. All other files are modular functions in the same folder.
clear; clc; close all;

%% ---- User parameters ----------------------------------------------------
nelx    = 160;     % elements in x
nely    = 100;     % elements in y
volfrac = 0.4;     % target volume fraction
penal   = 3;       % SIMP penalization power
rmin    = 1.5;     % filter radius (in elements)
ft      = 1;       % 1 = sensitivity filter, 2 = density filter
maxIter = 1000;    % max iterations
tolX    = 0.01;    % convergence tolerance on max design change

%% ---- Material -----------------------------------------------------------
mat.E0    = 1;
mat.Emin  = 1e-9;
mat.nu    = 0.3;
mat.penal = penal;

%% ---- Pre-processing -----------------------------------------------------
prob   = define_problem(nelx,nely);          % loads, supports, DOF sets
fe     = setup_fe(nelx,nely,mat.nu);         % KE, edofMat, sparse indices
[H,Hs] = build_filter(nelx,nely,rmin);       % filter matrix

%% ---- Initialization -----------------------------------------------------
x      = repmat(volfrac,nely,nelx);
xPhys  = x;
change = 1;
loop   = 0;
history = zeros(maxIter,1);

%% ---- Optimization loop --------------------------------------------------
while change > tolX && loop < maxIter
    loop = loop + 1;

    % 1) Finite element analysis
    U = fe_solve(xPhys,fe,prob,mat);

    % 2) Objective and sensitivities
    [c,dc,dv] = compute_sensitivities(xPhys,U,fe,mat);

    % 3) Filtering
    [dc,dv] = filter_sensitivities(dc,dv,x,H,Hs,ft);

    % 4) Design update (optimality criteria)
    [xnew,xPhys] = oc_update(x,dc,dv,H,Hs,volfrac,ft);

    change = max(abs(xnew(:)-x(:)));
    x = xnew;
    history(loop) = c;

    % 5) Output
    fprintf('It.:%4i  Obj.:%11.4f  Vol.:%7.3f  ch.:%7.3f\n', ...
        loop,c,mean(xPhys(:)),change);
    plot_design(xPhys);
end

history = history(1:loop);
figure; plot(history,'LineWidth',1.5); grid on;
xlabel('Iteration'); ylabel('Compliance'); title('Convergence history');