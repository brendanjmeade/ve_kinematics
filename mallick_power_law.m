function velocity_profile = mallick_power_law(ox, t_obs, Transition, etaval, powerval, Trecur)
    % Specify assumed parameters
    Nx = 0.2;
    M = 100;
    Vpl = 1e-9;
    
    % Call the spin physics and spin up function
    RUN_MAIN_ViscoEQ_imposedcycles(Nx, M, Transition, 0, powerval, etaval, Trecur, Vpl);
    
    % Surface observation points
    ox = ox(:);
    obs = [ox, zeros(length(ox), 1)];
    
    % Read in precomputed spun-up components
    odir = ['imposedviscocycles_' num2str(Vpl,'%1.e') '/Nx_' num2str(Nx) '_M_' num2str(M) '/power_' num2str(powerval) '_' num2str(etaval) '/Trec_' num2str(Trecur)];
    
    % Load results
    load([odir '/modelgeometry.mat'], 'ss', 'shz', 'evl', 'Teq')
    load([odir '/modeloutputs_f.mat'])
    load([odir '/modeloutputs_s.mat'])
    
    % Compute surface velocities    
    Gd = compute_displacementkernels(obs, ss, shz);
    deepx3 = max(shz.A(:, 2));
    vsurf_deep = (repmat(Vpl / pi .* atan2(ox, deepx3), 1, length(t)))';
    vsurf_v = (Gd.l12d * e12d' +  Gd.l13d * e13d')'; % no deep loading
    vsurf = vsurf_v + vsurf_deep;
  
    % Select velocity profile closed to observation time
    [~, index] = min(abs(t - t_obs));
    velocity_profile = vsurf(index, :) ./ Vpl;
end