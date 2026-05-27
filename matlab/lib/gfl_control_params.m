function ctrl = gfl_control_params()
%GFL_CONTROL_PARAMS  GFL converter control parameters (Table IX, Yuan et al. 2025).

ctrl.Kp_pll = 26;
ctrl.Ki_pll = 7800;
ctrl.Kp_p   = 0.5;
ctrl.Ki_p   = 5;
ctrl.Lf     = 0.05;
ctrl.Cf     = 0.05;
end
