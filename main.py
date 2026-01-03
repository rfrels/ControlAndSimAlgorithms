import control_and_sim_algorithms as csa
import matplotlib.pyplot as plt
import time

fig, ax = plt.subplots()
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

ax.plot(t, sim, label='simulation')
ax.plot(t, sim_diff, label='sim_diff')
ax.plot(t, ex, label='exact solution')

ax.set(xlabel='time (s)', ylabel='amplitude (m)',
       title='Damped Harmonic Oscillator')
ax.grid()
ax.legend()
plt.show()