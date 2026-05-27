% Code to plot simulation results from IEEE39BusSystem
%% Plot Description:
%
% This plot shows the rotor speed and terminal voltage of Generator-2 at 
% the Bus-31 swing bus and Generator-5 at the Bus-34 PV bus.

% Copyright 2024 The MathWorks, Inc.

% Generate new simulation results if they don't exist
if ~exist('simlog_IEEE39BusSystem', 'var')
    sim('IEEE39BusSystem')
end

% Reuse figure if it exists, else create new figure
if ~exist('h1_IEEE39BusSystem', 'var') || ...
        ~isgraphics(h1_IEEE39BusSystem, 'figure')
    h1_IEEE39BusSystem = figure('Name', 'IEEE39BusSystem');
end
figure(h1_IEEE39BusSystem)
clf(h1_IEEE39BusSystem)

plotData(simlog_IEEE39BusSystem)

% Plot rotor speeds and terminal voltages of generators
function plotData(simlog)

% Get simulation results
t = simlog.Generators.Gen1Bus39.Rotor_velocity.pu_output.series.time;
w_Gen2Bus31 = simlog.Generators.Gen2Bus31.Rotor_velocity.pu_output.series.values;
Vt_Gen2Bus31 = simlog.Generators.Gen2Bus31.Terminal_voltage.pu_output.series.values;
w_Gen5Bus34 = simlog.Generators.Gen5Bus34.Rotor_velocity.pu_output.series.values;
Vt_Gen5Bus34 = simlog.Generators.Gen5Bus34.Terminal_voltage.pu_output.series.values;

% Plot results
handles(1) = subplot(2, 1, 1);
plot(t, w_Gen2Bus31, 'LineWidth', 1)
hold on
plot(t, w_Gen5Bus34, 'LineWidth', 1)
hold off
grid on
title('Generator rotor speed and terminal voltage')
ylabel('Rotor speed (p.u.)')
legend({'Gen2@Bus31', 'Gen5@Bus34'},'Location','southwest');

handles(2) = subplot(2, 1, 2);
plot(t, Vt_Gen2Bus31, 'LineWidth', 1)
hold on
plot(t, Vt_Gen5Bus34, 'LineWidth', 1)
hold off
grid on
ylabel('Terminal voltage (p.u.)')
xlabel('Time (s)')
legend({'Gen2@Bus31', 'Gen5@Bus34'},'Location','southwest');
linkaxes(handles, 'x')

end