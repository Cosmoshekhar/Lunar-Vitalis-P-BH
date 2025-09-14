%% ballistic_hop_trajectory&heat
% Ballistic‐hop ΔV, propellant budget, preliminary trajectory,
% AND CBE/MEV heat power & energy (per hop & totals).

%% --- User Inputs --------------------------------------------------------
m0            = 500.0;                           % kg, initial wet mass
dry_mass      = 200.0;                           % kg, required dry mass at end
isp           = 317.0;                           % s, Isp (updated)
hop_heights   = [1500,  2500,  3500,  4500,  5500];        % m, apex height per hop
hop_distances = [2000,3000,4000,5000,6000];       % m, horizontal distance per hop
g             = 1.62;                            % m/s^2, lunar gravity
n_hops        = numel(hop_heights);

% Thermal / engine inputs for heat modeling
thrust_N      = 1780;    % N, TOTAL thrust (updated: 1.78 kN)
f_heat_CBE    = 0.05;    % fraction of exhaust/mech power entering vehicle (CBE)
f_heat_MEV    = 0.09;    % fraction for MEV 
g0            = 9.80665; % m/s^2, standard gravity (for Isp & ve)
ve            = isp * g0;                 % m/s, effective exhaust velocity
mdot          = thrust_N / ve;            % kg/s, mass flow at given thrust
Qdot_CBE_W    = 0.5 * thrust_N * ve * f_heat_CBE; % W, peak heat into vehicle
Qdot_MEV_W    = 0.5 * thrust_N * ve * f_heat_MEV; % W, peak heat into vehicle

if numel(hop_distances) ~= n_hops
    error('hop_heights and hop_distances must both have length = n_hops.');
end

%% --- Run calculation ----------------------------------------------------
res = ballistic_hop_calculator(m0, dry_mass, isp, hop_heights, hop_distances, g);

% Heat energy per hop using burn time = propellant / mdot
t_burn_per_hop_s = res.prop_per_hop_kg / mdot;           % s
E_CBE_J_per_hop  = Qdot_CBE_W * t_burn_per_hop_s;        % J
E_MEV_J_per_hop  = Qdot_MEV_W * t_burn_per_hop_s;        % J
E_CBE_J_total    = sum(E_CBE_J_per_hop);
E_MEV_J_total    = sum(E_MEV_J_per_hop);

% Per-hop percentage of total heat energy
pct_CBE          = 100 * (E_CBE_J_per_hop / max(E_CBE_J_total, eps));
pct_MEV          = 100 * (E_MEV_J_per_hop / max(E_MEV_J_total, eps));

%% --- Display summary ----------------------------------------------------
fprintf('\n=== Ballistic Hopping Summary ===\n');
for i = 1:n_hops
    fprintf('Hop %d → H=%4.0f m, D=%5.0f m, ΔV=%7.2f m/s, Prop=%7.2f kg, T=%5.1f s\n', ...
      i, hop_heights(i), hop_distances(i), ...
      res.dv_per_hop_mps(i), res.prop_per_hop_kg(i), res.time_per_hop_s(i));
end
fprintf('Total ΔV:           %.2f m/s\n', res.total_dv_mps);
fprintf('Total propellant:   %.2f kg\n', res.total_prop_kg);
fprintf('Mass at end:        %.2f kg\n', res.mass_end_kg);
fprintf('Dry mass margin:    %.2f kg\n', res.mass_margin_kg);

% Thrust-to-weight (Moon) sanity check
TW0   = thrust_N / (m0 * g);
TWend = thrust_N / (res.mass_end_kg * g);
fprintf('\n=== Thrust/Weight Check (Lunar g=%.2f m/s^2) ===\n', g);
fprintf('Initial T/W:  %.2f\n', TW0);
fprintf('End-of-mission T/W: %.2f\n', TWend);
if TW0 < 1
    warning('Initial T/W < 1: liftoff not possible under finite-thrust assumptions.');
elseif TW0 < 1.2
    warning('Initial T/W is marginal; gravity losses may be significant vs. impulsive model.');
end

% Heat results
fprintf('\n=== Peak Heat Power While Thrusting ===\n');
fprintf('CBE peak heat:  %.1f kW\n', Qdot_CBE_W/1e3);
fprintf('MEV peak heat:  %.1f kW\n', Qdot_MEV_W/1e3);

fprintf('\n=== Heat Energy Per Hop ===\n');
for i = 1:n_hops
    fprintf(['Hop %d:  CBE=%6.2f MJ (%5.1f%% of total CBE)  |  ' ...
             'MEV=%6.2f MJ (%5.1f%% of total MEV)  |  burn time=%5.1f s\n'], ...
        i, E_CBE_J_per_hop(i)/1e6, pct_CBE(i), E_MEV_J_per_hop(i)/1e6, pct_MEV(i), t_burn_per_hop_s(i));
end
fprintf('TOTAL:  CBE=%6.2f MJ  |  MEV=%6.2f MJ\n', E_CBE_J_total/1e6, E_MEV_J_total/1e6);

%% --- Plot trajectories --------------------------------------------------
figure('Name','Hop Trajectories','NumberTitle','off');
hold on;
colors = lines(n_hops);
for i = 1:n_hops
    traj = res.trajectories{i};
    plot(traj.x, traj.z, 'Color', colors(i,:), 'LineWidth', 1.5);
end
xlabel('Horizontal distance x (m)'); ylabel('Altitude z (m)');
title('Preliminary Ballistic Hop Trajectories');
legend(arrayfun(@(i) sprintf('Hop %d',i),1:n_hops,'UniformOutput',false),'Location','Best');
grid on; hold off;

%% --- Function Definitions -----------------------------------------------

function R = ballistic_hop_calculator(m0, dry_mass, isp, H, D, g)
% Returns ΔV, propellant, time‐of‐flight, and trajectory for each hop.

    n    = numel(H);
    dv   = zeros(n,1);
    mp   = zeros(n,1);
    tf   = zeros(n,1);
    traj = cell(n,1);
    mcur = m0;

    for i = 1:n
        % Vertical speed & flight time (apex at H(i))
        v_vert = sqrt(2 * g * H(i));
        tf(i)  = 2 * v_vert / g;

        % Horizontal speed to cover D(i) in time tf
        v_horiz = D(i) / tf(i);

        % Launch speed and ΔV per hop (launch + landing)
        v0    = hypot(v_vert, v_horiz);
        dv(i) = 2 * v0;

        % Propellant for this hop & mass update (equivalent single burn)
        mp(i) = propellant_mass(dv(i), mcur, isp);
        mcur  = mcur - mp(i);

        % Preliminary trajectory samples
        t_vec = linspace(0, tf(i), 100).';
        x_vec = v_horiz * t_vec;
        z_vec = v_vert  * t_vec - 0.5 * g * t_vec.^2;
        traj{i} = struct('t',t_vec, 'x',x_vec, 'z',z_vec);
    end

    R = struct( ...
      'dv_per_hop_mps',       dv, ...
      'total_dv_mps',         sum(dv), ...
      'prop_per_hop_kg',      mp, ...
      'total_prop_kg',        sum(mp), ...
      'time_per_hop_s',       tf, ...
      'total_flight_time_s',  sum(tf), ...
      'mass_end_kg',          mcur, ...
      'mass_margin_kg',       mcur - dry_mass, ...
      'trajectories',         {traj} ...
    );
end

function m_prop = propellant_mass(delta_v, m0, isp)
% Rocket eqn: m_prop = m0 - m0/exp(ΔV/(Isp*g0))
    g0      = 9.80665;   % m/s^2
    ratio   = exp(delta_v / (isp * g0));
    m_final = m0 / ratio;
    m_prop  = m0 - m_final;
end
