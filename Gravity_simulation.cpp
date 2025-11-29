// A physical simulation of gravity between a chosen amount of objects.

#include "Gravity_constants.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <fstream>
#include <chrono>
#include <iomanip>

/// @brief Calculate the gravitational force exerted on object i by object j.
/// @param m_i Mass of object i
/// @param m_j Mass of object j
/// @param r_i Position vector of object i
/// @param r_j Position vector of object j
/// @return The gravitational force vector
std::vector<double> gravity_force(double m_i, double m_j, std::vector<double> r_i, std::vector<double> r_j) {
    std::vector<double> r_ij(3);
    for (int k = 0; k < 3; k++) {
        r_ij[k] = r_i[k] - r_j[k];
    }

    double distance = std::sqrt( std::inner_product(r_ij.begin(), r_ij.end(), r_ij.begin(), 0.0) );
    if (distance < minimal_distance) {
        distance = minimal_distance;
    }

    std::vector<double> force(3);
    for (int k = 0; k < 3; k++) {  // note r_ij / distance is the unit vector
        force[k] = -G * m_i * m_j / (distance * distance) * r_ij[k] / distance;
    }
    return force;
}

/// @brief Find the velocity at the next timestep using the Verlet integration method.
/// @param current_velocity Velocity vector at the current timestep
/// @param next_acceleration Acceleration vector at the next timestep
/// @param current_acceleration Acceleration vector at the current timestep
/// @return The next velocity vector
std::vector<double> next_velocity(
    std::vector<double> current_velocity,
    std::vector<double> next_acceleration,
    std::vector<double> current_acceleration
) {
    std::vector<double> next_velocity_vec(3);
    for (int k = 0; k < 3; k++) {
        next_velocity_vec[k] = current_velocity[k] + 0.5 * dt * (next_acceleration[k] + current_acceleration[k]);
    }
    return next_velocity_vec;
}

/// @brief Find the position at the next timestep using the Verlet integration method.
/// @param current_position Position vector at the current timestep
/// @param previous_position Position vector at the previous timestep
/// @param current_acceleration Acceleration vector at the current timestep
/// @return The next position vector
std::vector<double> next_position(
    std::vector<double> current_position,
    std::vector<double> previous_position,
    std::vector<double> current_acceleration
) {
    std::vector<double> next_position_vec(3);
    for (int k = 0; k < 3; k++) {
        next_position_vec[k] = 2 * current_position[k] - previous_position[k] + current_acceleration[k] * dt*dt;
    }
    return next_position_vec;
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

    // Initialize positions and velocities
    // Here positions[i] are the positions of all objects at timestep i,
    // positions[i][j] is the position of object j at timestep i,
    // positions[i][j][k] is the k-coordinate of object j at timestep i.
    std::vector<std::vector<std::vector<double>>> positions(
        N_steps,
        std::vector<std::vector<double>>(
            N_objects,
            std::vector<double>(3, 0.0)  // 3 for x, y, z coordinates, initialized to 0.0
        )
    );

    std::vector<std::vector<std::vector<double>>> velocities(
        N_steps,
        std::vector<std::vector<double>>(
            N_objects,
            std::vector<double>(3, 0.0)  // 3 for x, y, z coordinates, initialized to 0.0
        )
    );

    // Assign initial values
    for (int j = 0; j < N_objects; j++) {
        positions[0][j] = initial_positions[j];
        velocities[0][j] = initial_velocities[j];
    }

    // Also need to keep the previous acceleration for velocity
    // Note we don't save these for all timesteps, only the previous one
    std::vector<std::vector<double>> previous_accelerations(
        N_objects,
        std::vector<double>(3, 0.0)
    );

    // Use Euler integration to find the positions at the second timestep
    for (int j = 0; j < N_objects; j++) {
        std::vector<double> force(3, 0.0);  // to add up all forces on object j
        for (int k = 0; k < N_objects; k++) {
            // Consider all other objects k, so without k==j
            if (k != j) {
                // Loop for: force += gravity_force(masses[j], masses[k], positions[0][j], positions[0][k])
                std::vector<double> f_ij = gravity_force(masses[j], masses[k], positions[0][j], positions[0][k]);
                for (int l = 0; l < 3; l++) {
                    force[l] += f_ij[l];
                }
            }
        }
        
        for (int l = 0; l < 3; l++) {
            double acceleration = force[l] / masses[j];
            previous_accelerations[j][l] = acceleration;  // update these to contain step 0 accelerations

            // Now update positions with Euler
            positions[1][j][l] = positions[0][j][l] + velocities[0][j][l] * dt + 0.5 * acceleration * dt*dt;
        }
    }

    // With the new positions of all the particles, calculate the acceleration at step 1
    for (int j = 0; j < N_objects; j++) {
        std::vector<double> force(3, 0.0);  // to add up all forces on object j
        for (int k = 0; k < N_objects; k++) {
            // Consider all other objects k, so without k==j
            if (k != j) {
                // Loop for: force += gravity_force(masses[j], masses[k], positions[0][j], positions[0][k])
                std::vector<double> f_ij = gravity_force(masses[j], masses[k], positions[0][j], positions[0][k]);
                for (int l = 0; l < 3; l++) {
                    force[l] += f_ij[l];
                }
            }
        }

        std::vector<double> acceleration(3, 0.0);
        for (int l = 0; l < 3; l++) {
            acceleration[l] = force[l] / masses[j];
        }

        // Now update velocities with Verlet
        velocities[1][j] = next_velocity(velocities[0][j], acceleration, previous_accelerations[j]);

        // And update the previous acceleration for the next loop
        previous_accelerations[j] = acceleration;
    }
    
    // Main simulation loop for all remaining steps
    for (size_t i = 2; i < N_steps; i++) {  // 0 set by initial conditions, 1 set by Euler
        for (int j = 0; j < N_objects; j++) {
            // Update positions, with the already calculated previous acceleration
            positions[i][j] = next_position(positions[i-1][j], positions[i-2][j], previous_accelerations[j]);
        }

        // Once all positions are updated, calculate new accelerations and update velocities
        for (int j = 0; j < N_objects; j++) {
            std::vector<double> force(3, 0.0);  // to add up all forces on object j
            for (int k = 0; k < N_objects; k++) {
                // Consider all other objects k, so without k==j
                if (k != j) {
                    // Loop for: force += gravity_force(masses[j], masses[k], positions[i][j], positions[i][k])
                    std::vector<double> f_ij = gravity_force(masses[j], masses[k], positions[i][j], positions[i][k]);
                    for (int l = 0; l < 3; l++) {
                        force[l] += f_ij[l];
                    }
                }
            }

            std::vector<double> acceleration(3, 0.0);
            for (int l = 0; l < 3; l++) {
                acceleration[l] = force[l] / masses[j];
            }

            // So update velocities
            velocities[i][j] = next_velocity(velocities[i-1][j], acceleration, previous_accelerations[j]);
            // And update the previous acceleration for the next loop
            previous_accelerations[j] = acceleration;
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
                sum_of_squares += velocities[i][j][k]*velocities[i][j][k];
            }
            kinetic_energy += 0.5 * masses[j] * sum_of_squares;  // sum the kinetic energy for all particles
        }
        kinetic_total[i] = kinetic_energy;

        double potential_energy = 0.0;  // get the potential energy of all objects at this step
        for (int j = 0; j < N_objects-1; j++) {  // don't need to consider the last object, is cancelled in next loop
            std::vector<double> r_j = positions[i][j];
            for (int k = j+1; k < N_objects; k++) {  // start from the object after j
                std::vector<double> r_k = positions[i][k];
                std::vector<double> r_ij(3);
                for (int l = 0; l < 3; l++) {
                    r_ij[l] = r_j[l] - r_k[l];
                }

                double distance = std::sqrt( std::inner_product(r_ij.begin(), r_ij.end(), r_ij.begin(), 0.0) );
                if (distance < minimal_distance) {
                    distance = minimal_distance;
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
    std::string filename = "output\\Gravity_simulation_output_cpp.csv";
    std::ofstream out_file(filename);
    for (size_t i = 0; i < N_steps; i++) {
        for (int j = 0; j < N_objects; j++) {
            for (int k = 0; k < 3; k++) {
                out_file << positions[i][j][k] << ",";
            }
        }
        out_file << kinetic_total[i] << "," << potential_total[i] << "," << energy[i] << "\n";
    }
    out_file.close();
    std::cout << "Output written!" << std::endl;
}