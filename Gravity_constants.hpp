// Constants and initial conditions for the gravity simulation.
#pragma once
#include <vector>
#include <optional>

constexpr double G = 1.0;                                   // Gravitational constant, 1 in simulation units
constexpr double dt = 1e-4;                                 // Timestep
constexpr int N_steps = 10000;                              // Number of timesteps
constexpr int N_objects = 4;                                // Number of objects
constexpr double e = 1.0;                                   // Coefficient of restitution
constexpr std::optional<double> box_length = std::nullopt;  // Length of the box for boundary conditions, if any
constexpr bool periodic = false;                            // Whether to use periodic boundary conditions
constexpr std::optional<double> e_wall = std::nullopt;      // Coefficient of restitution for wall collisions, if any

inline std::vector<double> masses = {10.0, 10.0, 10.0, 10.0};  // Masses of the objects
inline std::vector<double> radii = {0.01, 0.01, 0.01, 0.01};   // Radii of the objects
inline std::vector<std::vector<double>> initial_positions = {
    {1.0, 0.0, 0.0},  // For each object, (x,y,z)
    {0.0, 1.0, 0.0},
    {0.0, -1.0, 0.0},
    {-1.0, 0.0, 0.0}
};
inline std::vector<std::vector<double>> initial_velocities = {
    {0.0, -1.0, 0.0},  // For each object, (vx,vy,vz)
    {1.0, 0.0, 0.0},
    {-1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0}
};