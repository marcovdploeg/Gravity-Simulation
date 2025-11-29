# A Python script to plot the output trajectories and energies from the C++ gravity simulation.
import matplotlib.pyplot as plt
import numpy as np

print("Plotting results...")

# Need these constants again
dt = 1e-4  # timestep
N_steps = 10000
times = np.arange(N_steps) * dt
N_objects = 4

filename = "output\\Gravity_simulation_output_cpp.csv"
data = np.loadtxt(filename, delimiter=",")

# Extract the proper arrays from the data
positions = np.zeros(shape=[N_steps, N_objects, 3])
kinetic_energy = np.zeros(N_steps)
potential_energy = np.zeros(N_steps)
energy = np.zeros(N_steps)
for i in range(N_steps):
    for j in range(N_objects):
        for k in range(3):
            positions[i][j][k] = data[i][3*j + k]
    kinetic_energy[i] = data[i][-3]
    potential_energy[i] = data[i][-2]
    energy[i] = data[i][-1]

# Plot of trajectories in xy-plane
plt.figure(figsize=(5,5))
for i in range(N_objects):
    obj_i_positions = positions[:,i,:]
    plt.plot(obj_i_positions[:,0], obj_i_positions[:,1], marker='.', label=f'Object {i+1}')
    # Object starts at green dot, ends at red dot
    plt.plot(obj_i_positions[0,0], obj_i_positions[0,1], marker='.', color='green', markersize=12)
    plt.plot(obj_i_positions[-1,0], obj_i_positions[-1,1], marker='.', color='red', markersize=12)

plt.xlabel('Position x')
plt.ylabel('Position y')
plt.title('Trajectories in the system')
plt.legend()
plt.xlim(-1.2, 1.2)
plt.ylim(-1.2, 1.2)
plt.savefig('output\\trajectories_xy_cpp.png', dpi=300)
plt.show()

# Plot the different energies over time
plt.figure()
plt.plot(times, energy, label="Total")
plt.plot(times, kinetic_energy, label="Kinetic")
plt.plot(times, potential_energy, label="Potential")
plt.title("Energies in the system")
plt.xlabel("Time")
plt.ylabel("Energy")
plt.legend()
plt.savefig('output\\energies_cpp.png', dpi=300)
plt.show()

print("Results plotted!")