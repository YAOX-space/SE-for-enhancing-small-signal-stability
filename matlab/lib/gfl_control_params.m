function ctrl = gfl_control_params()
%GFL_CONTROL_PARAMS  GFL converter control parameters (Table IX, Yuan et al. 2025).
%
% Returns a struct compatible with both:
%   build_gfl_state_space  (uses Kp_pll, Ki_pll, Kp_p, Ki_p)
%   gfl9_ode / run_gfl_full_nonlinear  (uses full 9-state ODE params)

ctrl.omega0  = 2*pi*50;              % base angular frequency (rad/s)

% PLL
ctrl.Kp_pll  = 26;
ctrl.Ki_pll  = 7800;

% Outer active-power loop
ctrl.Kp_p    = 0.5;
ctrl.Ki_p    = 5;

% Inner current loop
ctrl.Kp_i    = 1;
ctrl.Ki_i    = 10;

% LCL filter (per-unit reactance → physical: X_pu/omega0)
ctrl.Lf      = 0.05  / ctrl.omega0;
ctrl.Rf      = 0.002 / ctrl.omega0;
ctrl.Cf      = 0.05  / ctrl.omega0;  % used in 11-state ODE (gfl11_ode)

% GFF feedforward filter time constant (1/(1+tau*s), Table IX)
ctrl.tau_ff  = 0.01;
end
