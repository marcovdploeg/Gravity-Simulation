# A physical simulation of gravity between a chosen amount of objects.
# Now also with collision physics and boundary conditions.

import numpy as np
import matplotlib.pyplot as plt
import time

#####   START of set up of parameters and initial conditions   #####
output_dir = 'example_output\\basic_test'
# Constants of the simulation
G = 1
dt = 1e-4  # timestep
N_steps = 10000
N_objects = 4
e = 1.0  # coefficient of restitution
box_length = None
periodic = False
e_wall = None
times = np.arange(N_steps) * dt

# Make an array with the mass and radius of each object
masses = np.array([10, 10, 10, 10])
radii = np.array([0.01, 0.01, 0.01, 0.01])

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

# Make sure the user input is consistent for masses/radii; 
# initial positions/velocities won't cause errors
if len(masses) != N_objects:
    raise KeyboardInterrupt(f"The number of masses ({len(masses)}) does not match the number of objects ({N_objects})!")
elif len(radii) != N_objects:
    raise KeyboardInterrupt(f"The number of radii ({len(radii)}) does not match the number of objects ({N_objects})!")

# Make sure boundary conditions are consistent
if box_length is None and periodic:
    print("You can't have periodic boundary conditions without a box length!")
    print("The simulation will now apply the minimal image convention but objects will not loop around the box.")
    raise KeyboardInterrupt("Please define box_length or set periodic to False.")
elif box_length is None and e_wall is not None:
    print("You can't have wall collisions without a box length!")
    print("The simulation will now ignore wall collisions, as if e_wall was set to None.")
    raise KeyboardInterrupt("Please define box_length or set e_wall to None.")

# Write the used parameters to a text file for later reference
with open(f'{output_dir}\\python_parameters.txt', 'w') as f:
    f.write(f'output_dir = {output_dir}\n')
    f.write(f'G = {G}\n')
    f.write(f'dt = {dt}\n')
    f.write(f'N_steps = {N_steps}\n')
    f.write(f'N_objects = {N_objects}\n')
    f.write(f'e = {e}\n')
    f.write(f'box_length = {box_length}\n')
    f.write(f'periodic = {periodic}\n')
    f.write(f'e_wall = {e_wall}\n')
    f.write(f'masses = {masses}\n')
    f.write(f'radii = {radii}\n')
    f.write(f'initial_positions = \n{positions[0]}\n')
    f.write(f'initial_velocities = \n{velocities[0]}\n')

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

def next_position(current_position, current_velocity, current_acceleration):
    """
    Find the position at the next timestep using the velocity Verlet integration method.
    Parameters:
        - current_position (np.array): Position at the current timestep
        - current_velocity (np.array): Velocity at the current timestep
        - current_acceleration (np.array): Acceleration at the current timestep
    """
    return current_position + current_velocity * dt + 0.5 * current_acceleration * dt**2

def galilean_transform(target_velocity, frame_velocity):
    """
    Perform a Galilean transformation.
    Parameters:
        - target_velocity (np.array): Velocity of the object that is to be transformed
        - frame_velocity (np.array): Velocity of the frame you want to transform to
    """
    return target_velocity - frame_velocity

def get_rotation_angles(position_i, position_j):
    """
    Calculate the rotation angles around the axes to align with 
    the line between the centres of the objects i and j.
    Parameters:
        - position_i (np.array): Position vector of object i
        - position_j (np.array): Position vector of object j
    Returns:
        - alpha (float): Angle around the x axis
        - beta (float): Angle around the y axis
        - gamma (float): Angle around the z axis
    """
    difference = position_i - position_j  # (delta_x, delta_y, delta_z)
    # Using arctan2 ensures the right quadrant sign
    alpha = np.arctan2(difference[2], difference[1])
    beta = np.arctan2(difference[0], difference[2])
    gamma = np.arctan2(difference[1], difference[0])
    return alpha, beta, gamma

def rotation(alpha, beta, gamma, target_vector):
    """
    Perform a rotation around the x, y and z axis.
    Parameters:
        - alpha (float): Angle around the x axis
        - beta (float): Angle around the y axis
        - gamma (float): Angle around the z axis
        - target_vector (np.array): Vector you want to rotate
    """
    rotation_matrix_x = np.array([
        [1, 0, 0],
        [0, np.cos(alpha), -np.sin(alpha)],
        [0, np.sin(alpha), np.cos(alpha)]
    ])
    rotation_matrix_y = np.array([
        [np.cos(beta), 0, np.sin(beta)],
        [0, 1, 0],
        [-np.sin(beta), 0, np.cos(beta)]
    ])
    rotation_matrix_z = np.array([
        [np.cos(gamma), -np.sin(gamma), 0],
        [np.sin(gamma), np.cos(gamma), 0],
        [0, 0, 1]
    ])
    combined_rotation = rotation_matrix_x @ rotation_matrix_y @ rotation_matrix_z  # Note xyz
    return combined_rotation @ target_vector

def inverse_rotation(alpha, beta, gamma, target_vector):
    """
    Perform an inverse rotation around the x, y and z axis,
    compared to the function above. Note no minus sign on the
    angles is needed relative to above!
    Parameters:
        - alpha (float): Angle around the x axis
        - beta (float): Angle around the y axis
        - gamma (float): Angle around the z axis
        - target_vector (np.array): Vector you want to rotate
    """
    rotation_matrix_x = np.array([
        [1, 0, 0],
        [0, np.cos(alpha), -np.sin(-alpha)],
        [0, np.sin(-alpha), np.cos(alpha)]
    ])
    rotation_matrix_y = np.array([
        [np.cos(beta), 0, np.sin(-beta)],
        [0, 1, 0],
        [-np.sin(-beta), 0, np.cos(beta)]
    ])
    rotation_matrix_z = np.array([
        [np.cos(gamma), -np.sin(-gamma), 0],
        [np.sin(-gamma), np.cos(gamma), 0],
        [0, 0, 1]
    ])
    combined_rotation = rotation_matrix_z @ rotation_matrix_y @ rotation_matrix_x  # Note zyx
    return combined_rotation @ target_vector

def velocity_transfer_collision(m_i, m_j, velocity_before_i):
    """
    Calculate the velocity of object i and j after an (in)elastic
    collision, where i moves and j is stationary.
    Parameters:
        - m_i (float): Mass of object i
        - m_j (float): Mass of object j
        - velocity_before_i (np.array): Velocity of object i before the collision
    Returns:
        - velocity_after_i (np.array): Velocity of object i after the collision
        - velocity_after_j (np.array): Velocity of object j after the collision
    """
    velocity_after_i = (m_i - e * m_j) / (m_i + m_j) * velocity_before_i
    velocity_after_j = (1 + e) * m_i / (m_i + m_j) * velocity_before_i
    return velocity_after_i, velocity_after_j

def objects_move_towards_each_other(position_i, position_j, velocity_i, velocity_j):
    """
    Determine if objects i and j are moving towards each other, so that a collision 
    can happen. This is the case if for any of the components of the position and 
    velocity vectors, the position component of i is smaller and the velocity 
    component of i is larger, or vice versa.
    Parameters:
        - position_i (np.array): Position vector of object i
        - position_j (np.array): Position vector of object j
        - velocity_i (np.array): Velocity vector of object i
        - velocity_j (np.array): Velocity vector of object j
    Returns:
        - bool: True if objects move towards each other
    """
    position_mask = position_i > position_j
    velocity_mask = velocity_i > velocity_j
    # We need position_i to be smaller and velocity_i larger or vice versa,
    # meaning these masks do not have the same value
    collisions = position_mask != velocity_mask
    # If any of these is True, there will be a collision
    return np.any(collisions)

def velocity_transfer_collision_wall(velocity_before_i):
    """
    Calculate the velocity of object i after an (in)elastic
    collision with the wall.
    Parameters:
        - velocity_before_i (np.array): Velocity of object i before the collision
    Returns:
        - velocity_after_i (np.array): Velocity of object i after the collision
    """
    velocity_after_i = -e_wall * velocity_before_i
    return velocity_after_i

def wall_collision_check_and_transfer(position_i, velocity_i):
    """
    Check if the object is colliding for each axis/component,
    and apply the velocity transfer if so.
    Parameters:
        - position_i (np.array): Position vector of object i
        - velocity_i (np.array): Velocity vector of object i
    Returns:
        - velocity_after_i (np.array): Velocity vector of object i after wall collision
    """
    walls = 0.5*np.array([box_length]*3)
    # Find where the object is outside the box
    where_cross = np.logical_or(position_i > walls, position_i < -walls)

    # Determine if a collision can happen, if the object moves towards the wall
    collisions = (position_i > 0) == (velocity_i > 0)

    # Apply velocity transfer where both conditions are met per component
    velocity_after_i = np.where(
        np.logical_and(where_cross, collisions),
        velocity_transfer_collision_wall(velocity_i),
        velocity_i
    )
    return velocity_after_i

def closest_copy_coordinates(position_i, positions):
    """
    Find the coordinates of the closest copy of each object relative to
    the object i according to the minimal image convention.
    Parameters:
        - position_i (np.array): Position vector of object i
        - positions (np.array): Array of position vectors of all objects
    Returns:
        - copy_positions (np.array): Array of position vectors of the closest copies of all objects
    """
    copy_positions = position_i -  (position_i - positions + 0.5*box_length) % box_length + 0.5*box_length
    return copy_positions

def apply_periodic_crossing(position_i):
    """
    Check if an object has crossed the periodic boundary, 
    and apply the correction in position if so.
    Parameters:
        - position_i (np.array): Position vector of object i
    Returns:
        - crossed_positions (np.array): Position vector of object i after applying periodic crossing correction
    """
    crossed_positions = np.copy(position_i)
    where_cross_positive = position_i > 0.5*box_length
    where_cross_negative = position_i < -0.5*box_length
    crossed_positions[where_cross_positive] -= box_length
    crossed_positions[where_cross_negative] += box_length
    return crossed_positions

# Also need to keep the previous acceleration for velocity
# Note we don't save these for all timesteps, only the previous one
previous_accelerations = np.zeros(shape=[N_objects, 3])

# Then also determine acceleration for step 0
for j in range(N_objects):
    force = np.zeros(3)  # to add up all forces on object j

    if periodic:  # get closest copy positions for object j
        copy_positions = closest_copy_coordinates(positions[0][j], positions[0])

    for k in range(N_objects):
        # Consider all other objects k, so without k==j
        if k != j:
            # Check for periodic box
            if not periodic:
                force += gravity_force(masses[j], masses[k], positions[0][j], positions[0][k])
            else:  # Apply minimal image convention
                force += gravity_force(masses[j], masses[k], positions[0][j], copy_positions[k])
    previous_accelerations[j] = force / masses[j]  # update these to contain step 0 accelerations

# Main simulation loop
for i in range(1, N_steps):  # 0 set by initial conditions
    for j in range(N_objects):
        # Update positions, with the already calculated previous acceleration
        positions[i][j] = next_position(positions[i-1][j], velocities[i-1][j], previous_accelerations[j])
    
    # Once all positions are updated, calculate new accelerations and update velocities
    for j in range(N_objects):
        force = np.zeros(3)  # to add up all forces on object j

        if periodic:  # get closest copy positions for object j
            copy_positions = closest_copy_coordinates(positions[i][j], positions[i])

        for k in range(N_objects):
            # Consider all other objects k, so without k==j
            if k != j:
                # Check for periodic box
                if not periodic:
                    force += gravity_force(masses[j], masses[k], positions[i][j], positions[i][k])
                else:  # Apply minimal image convention
                    force += gravity_force(masses[j], masses[k], positions[i][j], copy_positions[k])
        acceleration = force / masses[j]

        # So update velocities
        velocities[i][j] = next_velocity(velocities[i-1][j], acceleration, previous_accelerations[j])
        # And update the previous acceleration for the next loop
        previous_accelerations[j] = acceleration
    
    # Collisions
    for j in range(N_objects):
        position_j = positions[i][j]
        velocity_j = velocities[i][j]
        # Now only check for collisions between object pairs
        for k in range(j+1, N_objects):
            position_k = positions[i][k]
            velocity_k = velocities[i][k]
            # Check collision condition
            difference = position_j - position_k
            distance = np.linalg.norm(difference)
            if distance <= (radii[j] + radii[k]) and objects_move_towards_each_other(position_j, position_k, velocity_j, velocity_k):
                # 1 Galilean transform, with velocity_k as the frame velocity
                transformed_velocity = galilean_transform(velocity_j, velocity_k)

                # 2 Rotation
                alpha, beta, gamma = get_rotation_angles(position_j, position_k)
                rotated_transformed_velocity = rotation(alpha, beta, gamma, transformed_velocity)

                # 3 Velocity transfer
                velocity_after_j, velocity_after_k = velocity_transfer_collision(masses[j], masses[k], rotated_transformed_velocity)

                # 4 Inverse rotation
                velocity_after_j = inverse_rotation(alpha, beta, gamma, velocity_after_j)
                velocity_after_k = inverse_rotation(alpha, beta, gamma, velocity_after_k)

                # 5 Inverse Galilean transform
                velocity_after_j = galilean_transform(velocity_after_j, -velocity_k)
                velocity_after_k = galilean_transform(velocity_after_k, -velocity_k)

                # And finally assign
                velocities[i][j] = velocity_after_j
                velocities[i][k] = velocity_after_k
    
    # Boundary conditions
    if box_length is None:
        continue  # no boundary conditions, so skip to next timestep

    for j in range(N_objects):
        # Check if the object is outside the box
        if np.any(np.abs(positions[i][j]) > 0.5*box_length):
            # Check finite or periodic box
            if not periodic:
                velocities[i][j] = wall_collision_check_and_transfer(positions[i][j], velocities[i][j])
            else:
                positions[i][j] = apply_periodic_crossing(positions[i][j])

# Calculate energies
kinetic_per_particle = 0.5 * masses * np.sum(velocities**2, axis=2)
kinetic_total = np.sum(kinetic_per_particle, axis=1)

potential_total = np.zeros(shape=[N_steps])
for j in range(N_objects-1):  # don't need to consider the last object, is cancelled in next loop
    r_j = positions[:,j,:]
    for k in range(j+1, N_objects):  # start from the object after j
        if not periodic:
            r_k = positions[:,k,:]
            distance = np.linalg.norm(r_j - r_k, axis=1)  # distance between j and k at all timesteps
        else:  # Apply minimal image convention
            copy_r_k = closest_copy_coordinates(r_j, positions[:,k,:])
            distance = np.linalg.norm(r_j - copy_r_k, axis=1)
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
plt.grid()
#plt.xlim(-10.2, 10.2)
#plt.ylim(-10.2, 10.2)
plt.tight_layout()
plt.savefig(f'{output_dir}\\trajectories_xy_python.png', dpi=300)
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
plt.savefig(f'{output_dir}\\energies_python.png', dpi=300)
plt.show()

print("Results plotted!")