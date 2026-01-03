import control_and_sim_algorithms as csa
import matplotlib.pyplot as plt
import numpy as np
import time

#fig, ax = plt.subplots()
x0 = [1.2, 1.0]


start_time = time.perf_counter()  # Start the timer
result = csa.solve_exact_and_sim(30.0, 0.001, x0, csa.SystemType.UHO, csa.SimulationType.EC)
end_time = time.perf_counter()    # End the timer

print(f"Elapsed time: {end_time - start_time:.6f} seconds")

#print(result)
t = result[0]
sim = result[1]
sim_diff = result[2]
ex = result[3]

#ax.plot(t, sim, label='simulation')
#ax.plot(t, sim_diff, label='sim_diff')
#ax.plot(t, ex, label='exact solution')

#ax.set(xlabel='time (s)', ylabel='amplitude (m)',
       #title='Damped Harmonic Oscillator')
#ax.grid()
#ax.legend()
#plt.show()
def plot_sim_comparison(ax, t, sim, exact, method_name):
    """
    Calculates RMSE and plots simulation vs exact solution.
    """
    # Calculate RMSE
    sim_array = np.array(sim)
    exact_array = np.array(exact)
    rmse = np.sqrt(np.mean((sim_array - exact_array) ** 2))

    # Plotting
    ax.plot(t, exact, label='Exact', color='black', linestyle='--', alpha=0.7)
    ax.plot(t, sim, label=f'{method_name} (RMSE: {rmse:.5f})', color='tab:red', linewidth=1.2)

    ax.set_title(f'Method: {method_name}')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Amplitude (m)')
    ax.grid(True, which='both', linestyle='--', alpha=0.5)
    ax.legend(loc='upper right', fontsize='small')


# --- Simulation Configuration ---
x0 = [1.2, 1.0]
t_max = 30.0
dt = 0.008  # Increased slightly for better visual differentiation between solvers
sys_type = csa.SystemType.UHO

# Define the solvers to iterate through
solvers = [
    (csa.SimulationType.EC, "Euler-Cauchy"),
    (csa.SimulationType.MEC, "Modified Euler-Cauchy"),
    (csa.SimulationType.HN, "Heun"),
    #(csa.SimulationType.SN, "Simpson"),
    (csa.SimulationType.RK, "Runge-Kutta")
]

# Create a 4x2 grid
fig, axes = plt.subplots(4, 2, figsize=(12, 18), constrained_layout=True)
fig.suptitle(f'Damped Harmonic Oscillator: Solver Accuracy Analysis ($dt={dt}s$)', fontsize=16)

# --- Row 1: Ground Truth (Exact Solutions) ---
# We use the first simulation call to extract the analytical (exact) data
init_res = csa.solve_exact_and_sim(t_max, dt, x0, sys_type, solvers[0][0])
t = init_res[0]
ex_pos = init_res[3]
ex_vel = init_res[2]  # Assuming index 4 is exact derivative

axes[0, 0].plot(t, ex_pos, color='blue', label='Analytical $x(t)$')
axes[0, 0].set_title("Exact Solution: Position")
axes[0, 0].set_ylabel("Position (m)")
axes[0, 0].legend()
axes[0, 0].grid(True)

axes[0, 1].plot(t, ex_vel, color='green', label="Analytical $x'(t)$")
axes[0, 1].set_title("Exact Solution: Velocity")
axes[0, 1].set_ylabel("Velocity (m/s)")
axes[0, 1].legend()
axes[0, 1].grid(True)

# --- Rows 2-4: Comparison Plots ---
# Flatten the axes from the 2nd row onwards for easier iteration
comparison_axes = axes.flatten()[2:]

for i, (sim_type, name) in enumerate(solvers):
    start_time = time.perf_counter()
    res = csa.solve_exact_and_sim(t_max, dt, x0, sys_type, sim_type)
    end_time = time.perf_counter()

    print(f"Computed {name} in {end_time - start_time:.6f}s")

    # Passing data to helper function:
    # res[0]=time, res[1]=sim_pos, res[3]=exact_pos
    plot_sim_comparison(comparison_axes[i], res[0], res[1], res[3], name)

# Hide the 8th subplot (since we only have 5 solvers + 2 exact plots = 7 total)
if len(solvers) < len(comparison_axes):
    comparison_axes[-1].axis('off')

plt.show()