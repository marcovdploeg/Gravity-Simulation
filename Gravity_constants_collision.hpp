// Constants and initial conditions for the gravity simulation with collisions.
#pragma once
#include <vector>
#include <optional>
#include <string>

const std::string output_dir = "example_output\\collision_test";  // Directory to save the output files in

constexpr double G = 1.0;                                   // Gravitational constant, 1 in simulation units
constexpr double dt = 1e-3;                                 // Timestep
constexpr int N_steps = 20000;                              // Number of timesteps
constexpr int N_objects = 4;                                // Number of objects
constexpr double e = 0.8;                                   // Coefficient of restitution
constexpr std::optional<double> box_length = std::nullopt;  // Length of the box for boundary conditions, if any
constexpr bool periodic = false;                            // Whether to use periodic boundary conditions
constexpr std::optional<double> e_wall = std::nullopt;      // Coefficient of restitution for wall collisions, if any

inline std::vector<double> masses = {100.0, 10.0, 10.0, 1.0};  // Masses of the objects
inline std::vector<double> radii = {1.0, 0.5, 0.5, 0.1};       // Radii of the objects
inline std::vector<std::vector<double>> initial_positions = {
    {0.0, 0.0, 0.0},  // For each object, (x,y,z)
    {0.0, 10.0, 0.0},
    {0.0, -10.0, 0.0},
    {5.0, 5.0, 0.0}
};
inline std::vector<std::vector<double>> initial_velocities = {
    {0.0, 0.0, 0.0},  // For each object, (vx,vy,vz)
    {-4.0, -2.0, 0.0},
    {-4.0, 2.0, 0.0},
    {-1.0, -1.0, 0.0}
};