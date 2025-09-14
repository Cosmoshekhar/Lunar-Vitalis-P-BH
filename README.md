Lunar Hopper Mission Analysis: README

This MATLAB script, ballistic_hop_calculator.m, is a comprehensive tool for simulating and analysing a lunar hopping mission. It's designed to calculate crucial performance metrics for a vehicle undertaking a series of ballistic hops on the Moon's surface. The tool estimates the required delta-V (ΔV), propellant consumption, flight times, and thermal energy for a predefined mission profile.

The script is organised with a main section for user inputs and a core function (ballistic_hop_calculator) that handles the physics calculations. It's a great tool for preliminary mission design and trade studies.

Features

    ΔV & Propellant Budget: Calculates the ΔV required for each hop and the total mission, along with the corresponding propellant mass using the Tsiolkovsky Rocket Equation.

    Trajectory Simulation: Computes and plots the parabolic flight path for each hop, showing altitude versus horizontal distance.

    Thermal Analysis: Estimates the heat energy generated per hop and for the entire mission, using user-defined heat fractions for both a Current Best Estimate (CBE) and a Maximum Expected Value (MEV).

    Performance Metrics: Displays a summary of mission results, including:

        ΔV per hop and total.

        Propellant consumption per hop and total.

        Final vehicle mass and dry mass margin.

        Initial and end-of-mission thrust-to-weight (T/W) ratios.

How It Works

The script models each hop as a two-part impulsive burn: an ascent burn to impart the required initial velocity and a braking burn to nullify the velocity at landing. The trajectory during the coast phase is a simple parabola, assuming negligible atmospheric drag.

The core physics are based on:

    Vertical Velocity: vvert​=2gH​ (from energy conservation)

    Flight Time: tf​=g2vvert​​

    Horizontal Velocity: vhoriz​=tf​D​

    Launch Velocity: v0​=vvert2​+vhoriz2​​

    Per-Hop Delta-V: ΔVhop​=2v0​ (ascent + descent)

    Propellant Mass: mprop​=mcurrent​−eIsp​g0​ΔV​mcurrent​​ (Rocket Equation)

    Heat Power: Q˙​=0.5×Thrust×ve​×fheat​

Inputs

Users can configure the mission by modifying the following variables at the top of the script:

    m0: Initial wet mass (kg)

    dry_mass: Target dry mass (kg)

    isp: Specific impulse of the engine (s)

    hop_heights: An array of desired apex heights for each hop (m)

    hop_distances: An array of desired horizontal distances for each hop (m)

    g: Local gravitational acceleration (m/s²)

    thrust_N: Total engine thrust (N)

    f_heat_CBE: Heat fraction (Current Best Estimate)

    f_heat_MEV: Heat fraction (Maximum Expected Value)

Outputs

The script prints a detailed summary to the command window and generates a plot of the simulated trajectories.
