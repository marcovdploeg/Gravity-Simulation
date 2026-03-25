// Constants and initial conditions for the gravity simulation.
#pragma once
#include <vector>
#include <optional>
#include <string>
#include <cmath>

const std::string output_dir = "example_output\\external_force_test";  // Directory to save the output files in

constexpr double G = 1.0;                                   // Gravitational constant, 1 in simulation units
constexpr double dt = 1e-3;                                 // Timestep
constexpr int N_steps = 5000;                               // Number of timesteps
constexpr int N_objects = 3;                                // Number of objects
constexpr double e = 1.0;                                   // Coefficient of restitution
constexpr std::optional<double> box_length = std::nullopt;  // Length of the box for boundary conditions, if any
constexpr bool periodic = false;                            // Whether to use periodic boundary conditions
constexpr std::optional<double> e_wall = std::nullopt;      // Coefficient of restitution for wall collisions, if any

inline std::vector<double> masses = {0.01, 0.01, 0.01};  // Masses of the objects
inline std::vector<double> radii = {0.01, 0.01, 0.01};   // Radii of the objects
inline std::vector<std::vector<double>> initial_positions = {
    {1.0, 0.0, 0.0},  // For each object, (x,y,z)
    {-1.0, 0.0, 0.0},
    {0.0, 0.0, 0.0}
};
inline std::vector<std::vector<double>> initial_velocities = {
    {1.0, 0.0, 0.0},  // For each object, (vx,vy,vz)
    {-1.0, 0.0, 0.0},
    {0.0, 0.0, 0.0}
};

// External force setup
constexpr bool use_external_force = true;
inline std::vector<double> external_force(
    double mass,
    double time,
    double initial_time=0.0,
    double force_constant=100.0
) {
    double used_time = time + initial_time;
    std::vector<double> direction_vector = {0.0, 1.0, 0.0};  // force in the y-direction
    std::vector<double> force_vector(3);
    for (int i = 0; i < 3; i++) {
        force_vector[i] = force_constant * mass * std::sin(2.0 * M_PI * used_time) * direction_vector[i];
    }
    return force_vector;
}