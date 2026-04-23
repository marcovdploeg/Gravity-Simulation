// Constants and initial conditions for the gravity simulation.
#pragma once
#include <vector>
#include <optional>
#include <string>
#include <cmath>
#include <random>

const std::string output_dir = "example_output\\many_particles_test";  // Directory to save the output files in

constexpr double G = 1.0;                                   // Gravitational constant, 1 in simulation units
constexpr double dt = 5e-3;                                 // Timestep
constexpr int N_steps = 1000;                               // Number of timesteps
constexpr int N_objects = 100;                              // Number of objects
constexpr double e = 1.0;                                   // Coefficient of restitution
constexpr std::optional<double> box_length = 10.0;          // Length of the box for boundary conditions, if any
constexpr bool periodic = true;                             // Whether to use periodic boundary conditions
constexpr std::optional<double> e_wall = std::nullopt;      // Coefficient of restitution for wall collisions, if any

inline std::vector<double> masses(N_objects, 1.0);    // Masses of the objects
inline std::vector<double> radii(N_objects, 0.1);     // Radii of the objects

// Generate random uniformly distributed initial positions
inline const std::vector<std::vector<double>> initial_positions = [] {
    std::mt19937 engine(42);  // Seed = 42
    std::uniform_real_distribution<double> dist(-0.5*box_length.value(), 0.5*box_length.value());

    std::vector<std::vector<double>> positions(N_objects, std::vector<double>(3));
    for (int i = 0; i < N_objects; i++) {
        for (int j = 0; j < 3; j++) {
            positions[i][j] = dist(engine);
        }
    }
    return positions;
}();

// Generate random normally distributed initial velocities
inline const std::vector<std::vector<double>> initial_velocities = [] {
    std::mt19937 engine(42);  // Seed = 42
    constexpr double temp = 10.0;
    constexpr double k_B = 1.0;  // Boltzmann constant
    double sigma = std::sqrt(k_B * temp / masses[0]);
    std::normal_distribution<double> dist(0.0, sigma);  // Mean=0, stddev=sigma

    std::vector<std::vector<double>> velocities(N_objects, std::vector<double>(3));
    for (int i = 0; i < N_objects; i++) {
        for (int j = 0; j < 3; j++) {
            velocities[i][j] = dist(engine);
        }
    }
    return velocities;
}();

// External force setup
constexpr bool use_external_force = false;
inline std::vector<double> external_force() {
    // placeholder
    return {0.0, 0.0, 0.0};
}