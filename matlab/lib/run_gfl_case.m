function [t, dP] = run_gfl_case(A, B, C, D, u_mag, t_stop, rootDir)
%RUN_GFL_CASE  Simulate one GFL scenario: apply voltage sag, return ΔP(t).
%
%  Builds (or rebuilds) the Simulink state-space model se_gfl_physics_model,
%  injects a half-sine voltage-sag disturbance at t=0.10 s, runs the
%  simulation, and returns the active-power deviation matrix.
%
%  Inputs
%    A, B, C, D  – state-space matrices from build_gfl_state_space
%    u_mag       – sag magnitude (pu), e.g. 0.10
%    t_stop      – simulation end time (s), e.g. 1.0
%    rootDir     – path to the matlab/ directory (models stored here)
%
%  Outputs
%    t   – time vector  (T×1, seconds)
%    dP  – ΔP matrix    (T×N, per-unit active-power deviation)

dt    = 5e-4;
tv    = (0 : dt : t_stop)';
t0s   = 0.10;   % sag start (s)
t1s   = 0.12;   % sag end   (s)
u     = zeros(size(tv));
flt   = tv >= t0s & tv < t1s;
u(flt) = -u_mag * sin(pi*(tv(flt) - t0s)/(t1s - t0s)).^2;

modelName = 'se_gfl_physics_model';
modelPath = build_physics_simulink_model(rootDir, modelName, A, B, C, D);

assignin('base', 'ss_A',          A);
assignin('base', 'ss_B',          B);
assignin('base', 'ss_C',          C);
assignin('base', 'ss_D',          D);
assignin('base', 'fault_input',   timeseries(u, tv));
assignin('base', 'sim_stop_time', num2str(t_stop));

[modelDir, mname] = fileparts(modelPath);
oldDir = pwd;  cd(modelDir);
out = sim(mname, 'StopTime', num2str(t_stop));
cd(oldDir);

raw = out.get('simout_dP');
dP  = raw.signals.values;
if isvector(dP), dP = dP(:); end
t   = (0 : size(dP,1)-1)' * (t_stop / (size(dP,1)-1));
end
