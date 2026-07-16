// A physical simulation of gravity between a chosen amount of objects.
// Now also with collision physics, boundary conditions and external forces.

#include "Gravity_constants_many_particles.hpp"
#include <iostream>
#include <vector>
#include <numeric>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <string>

/// @brief Calculate the distance and distance vector between two objects i and j.
/// @param r_i Position vector of object i
/// @param r_j Position vector of object j
/// @return The distance between object i and j
/// @return The distance vector from object i to j, so r_i[k] - r_j[k]
std::tuple<double, std::vector<double>> calculate_distance_vector(
    std::vector<double> r_i, 
    std::vector<double> r_j
) {
    double distance = 0.0;
    std::vector<double> r_ij(3);
    for (int k = 0; k < 3; k++) {
        r_ij[k] = r_i[k] - r_j[k];
        distance += r_ij[k] * r_ij[k];
    }
    return {std::sqrt(distance), r_ij};
}

/// @brief Calculate a gravitational force component exerted on object i by object j.
/// @param m_i Mass of object i
/// @param m_j Mass of object j
/// @param distance Distance between object i and j
/// @param r_ij Distance vector component from object i to j, so r_i[k] - r_j[k]
/// @return The gravitational force vector
double gravity_force(
    double m_i, 
    double m_j, 
    double distance, 
    double r_ij
) {
    // note r_ij / distance is the unit vector component
    return -G * m_i * m_j / (distance * distance) * r_ij / distance;
}

/// @brief Find the velocity at the next timestep using the Verlet integration method.
/// @param current_velocity Velocity vector component at the current timestep
/// @param next_acceleration Acceleration vector component at the next timestep
/// @param current_acceleration Acceleration vector component at the current timestep
/// @return The next velocity vector component
double next_velocity(
    double current_velocity,
    double next_acceleration,
    double current_acceleration
) {
    return current_velocity + 0.5 * dt * (next_acceleration + current_acceleration);
}

/// @brief Find the position at the next timestep using the velocity Verlet integration method.
/// @param current_position Position vector at the current timestep
/// @param current_velocity Velocity vector at the current timestep
/// @param current_acceleration Acceleration vector at the current timestep
/// @return The next position vector
std::vector<double> next_position(
    std::vector<double> current_position,
    std::vector<double> current_velocity,
    std::vector<double> current_acceleration
) {
    std::vector<double> next_position_vec(3);
    for (int k = 0; k < 3; k++) {
        next_position_vec[k] = current_position[k] + current_velocity[k] * dt + 0.5 * current_acceleration[k] * dt*dt;
    }
    return next_position_vec;
}

/// @brief Perform a Galilean transformation.
/// @param target_velocity Velocity vector of the object that is to be transformed
/// @param frame_velocity Velocity vector of the frame you want to transform to
/// @return The transformed velocity vector
std::vector<double> galilean_transform(
    std::vector<double> target_velocity,
    std::vector<double> frame_velocity
) {
    std::vector<double> transformed_velocity(3);
    for (int k = 0; k < 3; k++) {
        transformed_velocity[k] = target_velocity[k] - frame_velocity[k];
    }
    return transformed_velocity;
}

/// @brief Calculate the rotation angles around the axes to align with 
/// the line between the centres of the objects i and j.
/// @param position_i Position vector of object i
/// @param position_j Position vector of object j
/// @return alpha, beta, gamma; angles around the x, y and z axis respectively
std::vector<double> get_rotation_angles(
    std::vector<double> position_i,
    std::vector<double> position_j
) {
    std::vector<double> difference(3);
    for (int k = 0; k < 3; k++) {
        difference[k] = position_i[k] - position_j[k];
    }  // (delta_x, delta_y, delta_z)
    // Using atan2 ensures the right quadrant sign
    double alpha = std::atan2(difference[2], difference[1]);
    double beta = std::atan2(difference[0], difference[2]);
    double gamma = std::atan2(difference[1], difference[0]);
    return {alpha, beta, gamma};
}

/// @brief Multiply two 3x3 matrices
/// @param A The first 3x3 matrix to multiply
/// @param B The second 3x3 matrix to multiply
/// @return C; The resulting 3x3 matrix
std::vector<std::vector<double>> matrices_multiplication(
    std::vector<std::vector<double>> A,
    std::vector<std::vector<double>> B
) {
    std::vector<std::vector<double>> C(3, std::vector<double>(3));
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            C[i][j] = 0;
            for (int k = 0; k < 3; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

/// @brief Multiply a 3x3 matrix with a 3D vector
/// @param matrix The 3x3 matrix to multiply
/// @param vector The 3D vector to be multiplied
/// @return The resulting 3D vector
std::vector<double> matrix_vector_multiplication(
    std::vector<std::vector<double>> matrix,
    std::vector<double> vector
) {
    std::vector<double> result_vector(3);
    for (int i = 0; i < 3; i++) {
        result_vector[i] = 0;
        for (int j = 0; j < 3; j++) {
            result_vector[i] += matrix[i][j] * vector[j];
        }
    }
    return result_vector;
}

/// @brief Perform a rotation around the x, y and z axis.
/// @param alpha Angle around the x axis
/// @param beta Angle around the y axis
/// @param gamma Angle around the z axis
/// @param target_vector Vector you want to rotate
/// @return The rotated vector
std::vector<double> rotation(
    double alpha, double beta, double gamma,
    std::vector<double> target_vector
) {
    std::vector<std::vector<double>> rotation_matrix_x = {
        {1, 0, 0},
        {0, std::cos(alpha), -std::sin(alpha)},
        {0, std::sin(alpha), std::cos(alpha)}
    };
    std::vector<std::vector<double>> rotation_matrix_y = {
        {std::cos(beta), 0, std::sin(beta)},
        {0, 1, 0},
        {-std::sin(beta), 0, std::cos(beta)}
    };
    std::vector<std::vector<double>> rotation_matrix_z = {
        {std::cos(gamma), -std::sin(gamma), 0},
        {std::sin(gamma), std::cos(gamma), 0},
        {0, 0, 1}
    };
    std::vector<std::vector<double>> combined_rotation = matrices_multiplication(
        matrices_multiplication(rotation_matrix_x, rotation_matrix_y),
        rotation_matrix_z
    );  // note xyz
    return matrix_vector_multiplication(combined_rotation, target_vector);
}

/// @brief Perform an inverse rotation around the x, y and z axis,
/// compared to the function above. Note no minus sign on the
/// angles is needed relative to above!
/// @param alpha Angle around the x axis
/// @param beta Angle around the y axis
/// @param gamma Angle around the z axis
/// @param target_vector Vector you want to rotate
/// @return The inversely rotated vector
std::vector<double> inverse_rotation(
    double alpha, double beta, double gamma,
    std::vector<double> target_vector
) {
    std::vector<std::vector<double>> rotation_matrix_x = {
        {1, 0, 0},
        {0, std::cos(alpha), -std::sin(-alpha)},
        {0, std::sin(-alpha), std::cos(alpha)}
    };
    std::vector<std::vector<double>> rotation_matrix_y = {
        {std::cos(beta), 0, std::sin(-beta)},
        {0, 1, 0},
        {-std::sin(-beta), 0, std::cos(beta)}
    };
    std::vector<std::vector<double>> rotation_matrix_z = {
        {std::cos(gamma), -std::sin(-gamma), 0},
        {std::sin(-gamma), std::cos(gamma), 0},
        {0, 0, 1}
    };
    std::vector<std::vector<double>> combined_rotation = matrices_multiplication(
        matrices_multiplication(rotation_matrix_z, rotation_matrix_y),
        rotation_matrix_x
    );  // note zyx
    return matrix_vector_multiplication(combined_rotation, target_vector);
}

/// @brief Calculate the velocity of object i and j after an (in)elastic
/// collision, where i moves and j is stationary.
/// @param m_i Mass of object i
/// @param m_j Mass of object j
/// @param velocity_before_i Velocity of object i before the collision
/// @return velocity_after_i; Velocity of object i after the collision
/// @return velocity_after_j; Velocity of object j after the collision
std::tuple<std::vector<double>, std::vector<double>> velocity_transfer_collision(
    double m_i, double m_j,
    std::vector<double> velocity_before_i
) {
    std::vector<double> velocity_after_i(3);
    std::vector<double> velocity_after_j(3);
    for (int k = 0; k < 3; k++) {
        velocity_after_i[k] = (m_i - e * m_j) / (m_i + m_j) * velocity_before_i[k];
        velocity_after_j[k] = (1 + e) * m_i / (m_i + m_j) * velocity_before_i[k];
    }
    return {velocity_after_i, velocity_after_j};
}

/// @brief Determine if objects i and j are moving towards each other, so that a collision 
/// can happen. This is the case if for any of the components of the position and 
/// velocity vectors, the position component of i is smaller and the velocity 
/// component of i is larger, or vice versa.
/// @param position_i Position vector of object i
/// @param position_j Position vector of object j
/// @param velocity_i Velocity vector of object i
/// @param velocity_j Velocity vector of object j
/// @return any_true; true if objects move towards each other
bool objects_move_towards_each_other(
    std::vector<double> position_i,
    std::vector<double> position_j,
    std::vector<double> velocity_i,
    std::vector<double> velocity_j
) {
    bool any_true = false;
    for (int k = 0; k < 3; k++) {
        // Check to see if both conditions are true, make any_true also true if so
        any_true = any_true || (position_i[k] < position_j[k] && velocity_i[k] > velocity_j[k]);
        any_true = any_true || (position_i[k] > position_j[k] && velocity_i[k] < velocity_j[k]);
    }
    return any_true;
}

/// @brief Calculate the velocity of object i after an (in)elastic collision with the wall.
/// @param velocity_before_i Velocity component of object i before the collision
/// @return velocity_after_i; Velocity component of object i after the collision
double velocity_transfer_collision_wall(
    double velocity_before_i
) {
    double velocity_after_i;
    // dereference e_wall since it's an optional, but always has a value if this function is called
    velocity_after_i = -(*e_wall) * velocity_before_i;
    return velocity_after_i;
}

/// @brief Check if the object is colliding for each axis/component, 
/// and apply the velocity transfer if so.
/// @param position_i Position vector of object i
/// @param velocity_i Velocity vector of object i
/// @return velocity_after_i; Velocity vector of object i after wall collision
std::vector<double> wall_collision_check_and_transfer(
    std::vector<double> position_i,
    std::vector<double> velocity_i
) {
    std::vector<double> velocity_after_i(3);
    for (int k = 0; k < 3; k++) {
        // Find where the object is outside of the box
        bool cross;
        // dereference box_length since it's an optional, but always has a value if this function is called
        cross = (position_i[k] > 0.5*(*box_length)) || (position_i[k] < -0.5*(*box_length));
    
        // Determine if a collision can happen, if the object moves towards the wall
        bool collision;
        collision = (position_i[k] > 0) == (velocity_i[k] > 0);

        // Apply velocity transfer where both conditions are for component k
        if (cross && collision) {
            velocity_after_i[k] = velocity_transfer_collision_wall(velocity_i[k]);
        }
        else {
            velocity_after_i[k] = velocity_i[k];
        }
    }
    return velocity_after_i;
}

/// @brief A modulus function that works the same as in Python, like we need in 
/// the closest_copy_coordinates function, so that the result is always in [0, b),
/// even for negative a
/// @param a The first number to take the modulus of
/// @param b The second number to take the modulus of
/// @return The modulus of a and b, as in Python
double python_mod(double a, double b) {
    return std::fmod(std::fmod(a, b) + b, b);
}

/// @brief Find the coordinates of the closest copy of each object relative to
///        the object i according to the minimal image convention.
/// @param position_i Position vector of object i
/// @param positions Array of position vectors of all objects
/// @return copy_positions; Array of position vectors of the closest copies of all objects
std::vector<double> closest_copy_coordinates(
    std::vector<double> position_i,
    std::vector<double> positions
) {
    std::vector<double> copy_positions(N_objects * 3, 0.0);
    for (int j = 0; j < N_objects; j++) {
        for (int k = 0; k < 3; k++) {
            // dereference box_length since it's an optional, but always has a value if this function is called
            copy_positions[j*3 + k] = position_i[k] - python_mod(position_i[k] - positions[j*3 + k] + 0.5*(*box_length), *box_length)
                                    + 0.5*(*box_length);
        }
    }
    return copy_positions;
}

/// @brief Check if an object has crossed the periodic boundary, 
///        and apply the correction in position if so.
/// @param position_i Position vector of object i
/// @return crossed_positions; Position vector of object i after applying periodic crossing correction
std::vector<double> apply_periodic_crossing(
    std::vector<double> position_i
) {
    std::vector<double> crossed_positions(3);
    for (int k = 0; k < 3; k++) {
        // dereference box_length since it's an optional, but always has a value if this function is called
        if (position_i[k] > 0.5*(*box_length)) {
            crossed_positions[k] = position_i[k] - *box_length;
        }
        else if (position_i[k] < -0.5*(*box_length)) {
            crossed_positions[k] = position_i[k] + *box_length;
        }
        else {
            crossed_positions[k] = position_i[k];
        }
    }
    return crossed_positions;
}

int main() {
    std::cout << "Starting gravity simulation..." << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();

    // Make sure the user input is consistent
    if (masses.size() != N_objects) {
        throw std::runtime_error(
            "The number of masses (" + std::to_string(masses.size()) +
            ") does not match number of objects (" + std::to_string(N_objects) + ")!"
        );
    }
    else if (radii.size() != N_objects) {
        throw std::runtime_error(
            "The number of radii (" + std::to_string(radii.size()) +
            ") does not match number of objects (" + std::to_string(N_objects) + ")!"
        );
    }
    else if (initial_positions.size() != N_objects) {
        throw std::runtime_error(
            "The number of initial positions (" + std::to_string(initial_positions.size()) +
            ") does not match number of objects (" + std::to_string(N_objects) + ")!"
        );
    }
    else if (initial_velocities.size() != N_objects) {
        throw std::runtime_error(
            "The number of initial velocities (" + std::to_string(initial_velocities.size()) +
            ") does not match number of objects (" + std::to_string(N_objects) + ")!"
        );
    }

    // Make sure boundary conditions are consistent
    if (!box_length.has_value() && periodic) {
        std::cerr << "You can't have periodic boundary conditions without a box length!" << std::endl;
        std::cerr << "The simulation will now apply the minimal image convention but objects will not loop around the box."
        << std::endl;
        throw std::runtime_error("Please define box_length or set periodic to False.");
    }
    else if (!box_length.has_value() && e_wall.has_value()) {
        std::cerr << "You can't have wall collisions without a box length!" << std::endl;
        std::cerr << "The simulation will now ignore wall collisions, as if e_wall was set to None." << std::endl;
        throw std::runtime_error("Please define box_length or set e_wall to None.");
    }

    // Initialize positions and velocities
    // In C++ it is faster to use flat vectors than jagged ones,
    // meaning we only use a single index l that is related to the i,j,k
    // we use in Python with l = i * (N_objects * 3) + j * 3 + k.
    // It has the same length of N_steps*N_objects*3, and is initialized to 0.0
    std::vector<double> positions(N_steps * N_objects * 3, 0.0);
    std::vector<double> velocities(N_steps * N_objects * 3, 0.0);

    // Assign initial values (from jagged initial vectors)
    for (int j = 0; j < N_objects; j++) {
        for (int k = 0; k < 3; k++) {
            positions[j*3 + k] = initial_positions[j][k];
            velocities[j*3 + k] = initial_velocities[j][k];
        }
    }

    // Also need to keep the previous acceleration for velocity
    // Note we don't save these for all timesteps, only the previous one
    std::vector<double> previous_accelerations(N_objects * 3, 0.0);

    // To add up all forces on each object
    std::vector<double> forces(N_objects * 3, 0.0);

    // Loop for initial forces
    for (int j = 0; j < N_objects; j++) {
        // Get the slice of positions for object j (with i=0)
        std::vector<double> position_j(
            positions.begin() + j*3,
            positions.begin() + (j+1)*3
        );

        std::vector<double> copy_positions;
        if (periodic) {
            // Get the slice of all positions at step 0
            std::vector<double> positions_0(
                positions.begin(), 
                positions.begin() + N_objects*3
            );
            copy_positions = closest_copy_coordinates(position_j, positions_0);
        }

        for (int k = j+1; k < N_objects; k++) {
            // Consider remaining objects k, so  k>j
            // Check for periodic box
            if (!periodic) {
                // Get the slice of positions for object k (with i=0)
                std::vector<double> position_k(
                    positions.begin() + k*3,
                    positions.begin() + (k+1)*3
                );
                auto [distance, r_ij] = calculate_distance_vector(position_j, position_k);
                for (int l = 0; l < 3; l++) {
                    double force = gravity_force(masses[j], masses[k], distance, r_ij[l]);
                    forces[j*3 + l] += force;
                    forces[k*3 + l] -= force;
                }
            }
            else {  // Apply minimal image convention
                // Get the slice of copy positions for object k
                std::vector<double> copy_position_k(
                    copy_positions.begin() + k*3,
                    copy_positions.begin() + (k+1)*3
                );
                auto [distance, r_ij] = calculate_distance_vector(position_j, copy_position_k);
                for (int l = 0; l < 3; l++) {
                    double force = gravity_force(masses[j], masses[k], distance, r_ij[l]);
                    forces[j*3 + l] += force;
                    forces[k*3 + l] -= force;
                }
            }
        }

        // Calculate external force
        if (use_external_force) {
            std::vector<double> external_force_vector = external_force();  // masses[j], 0.0
            for (int l = 0; l < 3; l++) {
                forces[j*3 + l] += external_force_vector[l];
            }
        }
        
        for (int l = 0; l < 3; l++) {
            previous_accelerations[j*3 + l] = forces[j*3 + l] / masses[j];  // update these to contain step 0 accelerations
        }
    }
    
    // Main simulation loop
    for (size_t i = 1; i < N_steps; i++) {  // 0 set by initial conditions
        for (int j = 0; j < N_objects; j++) {
            // Get the slice of positions, velocities and accelerations for object j
            std::vector<double> previous_position_j(
                positions.begin() + (i-1)*N_objects*3 + j*3,
                positions.begin() + (i-1)*N_objects*3 + (j+1)*3
            );
            std::vector<double> previous_velocity_j(
                velocities.begin() + (i-1)*N_objects*3 + j*3,
                velocities.begin() + (i-1)*N_objects*3 + (j+1)*3
            );
            std::vector<double> previous_acceleration_j(
                previous_accelerations.begin() + j*3,
                previous_accelerations.begin() + (j+1)*3
            );
            // Update positions, with the already calculated previous acceleration
            std::vector<double> position_j = next_position(previous_position_j, previous_velocity_j, previous_acceleration_j);
            for (int k = 0; k < 3; k++) {
                positions[i*N_objects*3 + j*3 + k] = position_j[k];
            }
        }

        // Once all positions are updated, calculate new accelerations and update velocities
        std::vector<double> forces(N_objects * 3, 0.0);  // to add up all forces on each object
        for (int j = 0; j < N_objects; j++) {
            // Get the slice of positions for object j
            std::vector<double> position_j(
                positions.begin() + i*N_objects*3 + j*3,
                positions.begin() + i*N_objects*3 + (j+1)*3
            );

            std::vector<double> copy_positions;
            if (periodic) {
                // Get the slice of all positions at step i
                std::vector<double> positions_i(
                    positions.begin() + i*N_objects*3, 
                    positions.begin() + (i+1)*N_objects*3
                );
                copy_positions = closest_copy_coordinates(position_j, positions_i);
            }

            for (int k = j+1; k < N_objects; k++) {
                // Consider remaining objects k, so k>j
                // Check for periodic box
                if (!periodic) {
                    // Get the slice of positions for object k
                    std::vector<double> position_k(
                        positions.begin() + i*N_objects*3 + k*3,
                        positions.begin() + i*N_objects*3 + (k+1)*3
                    );
                    auto [distance, r_ij] = calculate_distance_vector(position_j, position_k);
                    for (int l = 0; l < 3; l++) {
                        double force = gravity_force(masses[j], masses[k], distance, r_ij[l]);
                        forces[j*3 + l] += force;
                        forces[k*3 + l] -= force;
                    }
                }
                else {  // Apply minimal image convention
                    // Get the slice of copy positions for object k 
                    // (note copies are defined at this step, so no i index)
                    std::vector<double> copy_position_k(
                        copy_positions.begin() + k*3,
                        copy_positions.begin() + (k+1)*3
                    );
                    auto [distance, r_ij] = calculate_distance_vector(position_j, copy_position_k);
                    for (int l = 0; l < 3; l++) {
                        double force = gravity_force(masses[j], masses[k], distance, r_ij[l]);
                        forces[j*3 + l] += force;
                        forces[k*3 + l] -= force;
                    }
                }
            }

            // Calculate external force
            if (use_external_force) {
                std::vector<double> external_force_vector = external_force();  // masses[j], i*dt
                for (int l = 0; l < 3; l++) {
                    forces[j*3 + l] += external_force_vector[l];
                }
            }

            for (int l = 0; l < 3; l++) {
                double acceleration = forces[j*3 + l] / masses[j];
                // So update velocities, per component
                velocities[i*N_objects*3 + j*3 + l] = next_velocity(velocities[(i-1)*N_objects*3 + j*3 + l], acceleration, previous_accelerations[j*3 + l]);
                
                // And update the previous acceleration for the next loop
                previous_accelerations[j*3 + l] = acceleration;
            }
        }

        // Collisions
        for (int j = 0; j < N_objects; j++) {
            std::vector<double> position_j(
                positions.begin() + i*N_objects*3 + j*3,
                positions.begin() + i*N_objects*3 + (j+1)*3
            );
            // Now only check for collisions between object pairs
            for (int k = j+1; k < N_objects; k++) {
                // (Fun sidenote: in Python, we can define velocity_j in the previous loop,
                // because when we change velocities[i][j] due to a collision later, this also changes
                // velocity_j. But in C++, that does not happen, so we need to define it here, so
                // that a potential change due to a collision with a previous object is not lost.
                // Before I fixed this, it altered the example run quite significantly.)
                std::vector<double> velocity_j(
                    velocities.begin() + i*N_objects*3 + j*3,
                    velocities.begin() + i*N_objects*3 + (j+1)*3
                );

                std::vector<double> position_k(
                    positions.begin() + i*N_objects*3 + k*3,
                    positions.begin() + i*N_objects*3 + (k+1)*3
                );
                std::vector<double> velocity_k(
                    velocities.begin() + i*N_objects*3 + k*3,
                    velocities.begin() + i*N_objects*3 + (k+1)*3
                );
                // Check collision condition
                double distance = 0.0;
                for (int l = 0; l < 3; l++) {
                    distance += (position_j[l] - position_k[l]) * (position_j[l] - position_k[l]);
                }
                distance = std::sqrt(distance);

                bool moving_towards_each_other = objects_move_towards_each_other(position_j, position_k, velocity_j, velocity_k);
                if ((distance < radii[j] + radii[k]) && moving_towards_each_other) {
                    // 1 Galilean transform, with velocity_k as the frame velocity
                    std::vector<double> transformed_velocity = galilean_transform(velocity_j, velocity_k);

                    // 2 Rotation
                    std::vector<double> angles = get_rotation_angles(position_j, position_k);
                    std::vector<double> rotated_transformed_velocity = rotation(angles[0], angles[1], angles[2], 
                                                                                transformed_velocity);

                    // 3 Velocity transfer
                    auto [velocity_after_j, velocity_after_k] = velocity_transfer_collision(masses[j], masses[k], 
                                                                                            rotated_transformed_velocity);
                    
                    // 4 Inverse rotation
                    velocity_after_j = inverse_rotation(angles[0], angles[1], angles[2], velocity_after_j);
                    velocity_after_k = inverse_rotation(angles[0], angles[1], angles[2], velocity_after_k);

                    // 5 Inverse Galilean transform
                    std::vector<double> negative_velocity_k = { -velocity_k[0], -velocity_k[1], -velocity_k[2] };
                    velocity_after_j = galilean_transform(velocity_after_j, negative_velocity_k);
                    velocity_after_k = galilean_transform(velocity_after_k, negative_velocity_k);

                    // And finally assign
                    for (int l = 0; l < 3; l++) {
                        velocities[i*N_objects*3 + j*3 + l] = velocity_after_j[l];
                        velocities[i*N_objects*3 + k*3 + l] = velocity_after_k[l];
                    }
                }
            }
        }

        // Boundary conditions
        if (!box_length.has_value()) {
            continue;  // no boundary conditions, so skip to next timestep
        }
        for (int j = 0; j < N_objects; j++) {
            // Check if the object is outside the box
            bool is_out = false;
            for (int k = 0; k < 3; k++) {
                if (std::abs(positions[i*N_objects*3 + j*3 + k]) > 0.5 * (*box_length)) {
                    is_out = true;
                    break;
                }
            }
            if (is_out) {
                std::vector<double> position_j(
                    positions.begin() + i*N_objects*3 + j*3,
                    positions.begin() + i*N_objects*3 + (j+1)*3
                );
                // Check finite or periodic box
                if (!periodic) {
                    std::vector<double> velocity_j(
                        velocities.begin() + i*N_objects*3 + j*3,
                        velocities.begin() + i*N_objects*3 + (j+1)*3
                    );
                    std::vector<double> new_velocity_j = wall_collision_check_and_transfer(position_j, velocity_j);
                    for (int k = 0; k < 3; k++) {
                        velocities[i*N_objects*3 + j*3 + k] = new_velocity_j[k];
                    }
                }
                else {
                    std::vector<double> new_position_j = apply_periodic_crossing(position_j);
                    for (int k = 0; k < 3; k++) {
                        positions[i*N_objects*3 + j*3 + k] = new_position_j[k];
                    }
                }
            }
        }
    }

    // Calculate energies
    std::vector<double> kinetic_total(N_steps, 0.0);
    std::vector<double> potential_total(N_steps, 0.0);
    std::vector<double> energy(N_steps, 0.0);
    for (size_t i = 0; i < N_steps; i++) {
        double kinetic_energy = 0.0;  // get the kinetic energy of all objects at this step
        for (int j = 0; j < N_objects; j++) {
            double sum_of_squares = 0.0;
            for (int k = 0; k < 3; ++k) {
                sum_of_squares += velocities[i*N_objects*3 + j*3 + k]*velocities[i*N_objects*3 + j*3 + k];
            }
            kinetic_energy += 0.5 * masses[j] * sum_of_squares;  // sum the kinetic energy for all particles
        }
        kinetic_total[i] = kinetic_energy;

        double potential_energy = 0.0;  // get the potential energy of all objects at this step
        for (int j = 0; j < N_objects-1; j++) {  // don't need to consider the last object, is cancelled in next loop
            std::vector<double> r_j(
                positions.begin() + i*N_objects*3 + j*3,
                positions.begin() + i*N_objects*3 + (j+1)*3
            );
            for (int k = j+1; k < N_objects; k++) {  // start from the object after j
                double distance = 0.0;
                if (!periodic) {
                    std::vector<double> r_k(
                        positions.begin() + i*N_objects*3 + k*3,
                        positions.begin() + i*N_objects*3 + (k+1)*3
                    );
                    for (int l = 0; l < 3; l++) {
                        distance += (r_j[l] - r_k[l]) * (r_j[l] - r_k[l]);
                    }
                    distance = std::sqrt(distance);
                }
                else {  // Apply minimal image convention
                    std::vector<double> positions_i(
                        positions.begin() + i*N_objects*3, 
                        positions.begin() + (i+1)*N_objects*3
                    );
                    std::vector<double> copy_positions = closest_copy_coordinates(r_j, positions_i);
                    std::vector<double> copy_r_k(
                        copy_positions.begin() + k*3,
                        copy_positions.begin() + (k+1)*3
                    );
                    for (int l = 0; l < 3; l++) {
                        distance += (r_j[l] - copy_r_k[l]) * (r_j[l] - copy_r_k[l]);
                    }
                    distance = std::sqrt(distance);
                }
                double potential_jk = -G * masses[j] * masses[k] / distance;
                potential_energy += potential_jk;  // sum the potential energy for all particles
            }
        }
        potential_total[i] = potential_energy;

        energy[i] = kinetic_energy + potential_energy;
    }
    std::cout << "Simulation finished!" << std::endl;
    auto end_time = std::chrono::high_resolution_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    // Print in seconds up to 2 decimal places
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Simulation took " << elapsed_time.count() << " seconds." << std::endl;

    // Save the positions and different energies to a file to plot in Python
    std::cout << "Writing output..." << std::endl;
    std::string filename = output_dir + "\\Gravity_simulation_output_cpp.csv";
    std::ofstream out_file(filename);
    for (size_t i = 0; i < N_steps; i++) {
        for (int j = 0; j < N_objects; j++) {
            for (int k = 0; k < 3; k++) {
                out_file << positions[i*N_objects*3 + j*3 + k] << ",";
            }
        }
        out_file << kinetic_total[i] << "," << potential_total[i] << "," << energy[i] << "\n";
    }
    out_file.close();
    std::cout << "Output written!" << std::endl;
}