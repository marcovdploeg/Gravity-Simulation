# A physical simulation of gravity between a chosen amount of objects.

import numpy as np
import matplotlib.pyplot as plt
import time

#####   START of set up of parameters and initial conditions   #####
# Constants of the simulation
G = 1
minimal_distance = 0.05  # to avoid infinities
dt = 1e-4  # timestep
N_steps = 10000
times = np.arange(N_steps) * dt
N_objects = 4

# Make an array with the mass of each object
masses = np.array([10, 10, 10, 10])

# Make an array 'positions' to mark the positions of the objects
# Here positions[i] are the positions of all objects at timestep i,
# positions[i][j] is the position of object j at timestep i,
# positions[i][j][k] is the k-coordinate of object j at timestep i.
positions = np.zeros(shape=[N_steps, N_objects, 3])  # 3 for x, y, z coordinates

# Choose initial positions
positions[0][0] = np.array([1,0,0])
positions[0][1] = np.array([0,1,0])
positions[0][2] = np.array([0,-1,0])
positions[0][3] = np.array([-1,0,0])

# Make an array 'velocities' to mark the velocities of the objects
# with the same structure as positions
velocities = np.zeros(shape=[N_steps, N_objects, 3])

# Choose initial velocities
velocities[0][0] = np.array([0,-1,0])
velocities[0][1] = np.array([1,0,0])
velocities[0][2] = np.array([-1,0,0])
velocities[0][3] = np.array([0,1,0])

#####   END of setup   #####

print("Starting gravity simulation...")
start_time = time.time()

# Make sure the user input is consistent for masses; if initial positions/velocities
# are not specified these just start at 0 so that won't cause issues
if len(masses) != N_objects:
    raise KeyboardInterrupt(f"The number of masses ({len(masses)}) does not match the number of objects ({N_objects})!")

def gravity_force(m_i, m_j, r_i, r_j):
    """
    Calculate the gravitational force exerted on object i by object j.
    Parameters:
        - m_i (float): Mass of object i
        - m_j (float): Mass of object j
        - r_i (np.array): Position vector of object i
        - r_j (np.array): Position vector of object j
    """
    r_ij = r_i - r_j
    distance = np.linalg.norm(r_ij)
    if distance < minimal_distance:
        distance = minimal_distance
    unit_r_ij = r_ij / distance
    return -G * m_i * m_j / distance**2 * unit_r_ij

def next_velocity(current_velocity, next_acceleration, current_acceleration):
    """
    Find the velocity at the next timestep using the Verlet integration method.
    Parameters:
        - current_velocity (np.array): Velocity at the current timestep
        - next_acceleration (np.array): Acceleration at the next timestep
        - current_acceleration (np.array): Acceleration at the current timestep
    """
    return current_velocity + 0.5 * dt * (next_acceleration + current_acceleration)

def next_position(current_position, previous_position, current_acceleration):
    """
    Find the position at the next timestep using the Verlet integration method.
    Parameters:
        - current_position (np.array): Position at the current timestep
        - previous_position (np.array): Position at the previous timestep
        - current_acceleration (np.array): Acceleration at the current timestep
    """
    return 2 * current_position - previous_position + current_acceleration * dt**2

# Also need to keep the previous acceleration for velocity
# Note we don't save these for all timesteps, only the previous one
previous_accelerations = np.zeros(shape=[N_objects, 3])

# Use Euler integration to find the positions at the second timestep
for j in range(N_objects):
    force = np.zeros(3)  # to add up all forces on object j
    for k in range(N_objects):
        # Consider all other objects k, so without k==j
        if k != j:
            force += gravity_force(masses[j], masses[k], positions[0][j], positions[0][k])
    acceleration = force / masses[j]
    previous_accelerations[j] = acceleration  # update these to contain step 0 accelerations

    # Now update positions with Euler
    positions[1][j] = positions[0][j] + velocities[0][j] * dt + 0.5 * acceleration * dt**2

# With the new positions of all the particles, calculate the acceleration at step 1
for j in range(N_objects):
    force = np.zeros(3)  # to add up all forces on object j
    for k in range(N_objects):
        # Consider all other objects k, so without k==j
        if k != j:
            force += gravity_force(masses[j], masses[k], positions[1][j], positions[1][k])
    acceleration = force / masses[j]

    # Now update velocities with Verlet
    velocities[1][j] = next_velocity(velocities[0][j], acceleration, previous_accelerations[j])

    # And update the previous acceleration for the next loop
    previous_accelerations[j] = acceleration

# Main simulation loop for all remaining steps
for i in range(2, N_steps):  # 0 set by initial conditions, 1 set by Euler
    for j in range(N_objects):
        # Update positions, with the already calculated previous acceleration
        positions[i][j] = next_position(positions[i-1][j], positions[i-2][j], previous_accelerations[j])
    
    # Once all positions are updated, calculate new accelerations and update velocities
    for j in range(N_objects):
        force = np.zeros(3)  # to add up all forces on object j
        for k in range(N_objects):
            # Consider all other objects k, so without k==j
            if k != j:
                force += gravity_force(masses[j], masses[k], positions[i][j], positions[i][k])
        acceleration = force / masses[j]

        # So update velocities
        velocities[i][j] = next_velocity(velocities[i-1][j], acceleration, previous_accelerations[j])
        # And update the previous acceleration for the next loop
        previous_accelerations[j] = acceleration

# Calculate energies
kinetic_per_particle = 0.5 * masses * np.sum(velocities**2, axis=2)
kinetic_total = np.sum(kinetic_per_particle, axis=1)

potential_total = np.zeros(shape=[N_steps])
for j in range(N_objects-1):  # don't need to consider the last object, is cancelled in next loop
    r_j = positions[:,j,:]
    for k in range(j+1, N_objects):  # start from the object after j
        r_k = positions[:,k,:]
        distance = np.linalg.norm(r_j - r_k, axis=1)  # distance between j and k at all timesteps
        
        # Filter distances to avoid infinities
        distance[ distance < minimal_distance ] = minimal_distance

        potential_jk = -G * masses[j] * masses[k] / distance
        potential_total += potential_jk

energy = kinetic_total + potential_total

print("Simulation finished!")
end_time = time.time()
print(f"Simulation took {end_time - start_time:.2f} seconds.")

print("Plotting results...")

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
plt.savefig('output\\trajectories_xy_python.png', dpi=300)
plt.show()

# Plot the different energies over time
plt.figure()
plt.plot(times, energy, label="Total")
plt.plot(times, kinetic_total, label="Kinetic")
plt.plot(times, potential_total, label="Potential")
plt.title("Energies in the system")
plt.xlabel("Time")
plt.ylabel("Energy")
plt.legend()
plt.savefig('output\\energies_python.png', dpi=300)
plt.show()

print("Results plotted!")