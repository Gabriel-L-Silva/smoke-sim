#ifndef _GRID_FLUIDS_H
#define _GRID_FLUIDS_H

#include <Godot.hpp>
#include <Node2D.hpp>
#include <vector>


using namespace std;
using namespace godot;

class Vect{
public:
    double pressure;
    Vector2 pos;
    Vector2 vel;
};


class GridFluids: public Node2D 
{
    GODOT_CLASS(GridFluids, Node2D);

    // Exposed props
    double maxSpeed = 1000.0;
    int maxIterPoisson = 20;
    int subSteps = 10;
    double rho = 1.0;
    Vector2 tile_size = Vector2(1, 1);
    Vector2 grid_size = Vector2(800, 800);
    Vector2 vector_size = Vector2(0, 0);

public:
    static void _register_methods();
    void _init();

    // Functions 
    void update_field(double delta, Array grid, Vector2 externalForces);
    void update_grid(vector<vector<Vect>> &vectors, Array grid);
    void project(vector<vector<Vect>> &vectors);
    Vector2 gradient_at_point(int x, int y, vector<vector<double>> &grid);
    vector<vector<Vector2>> gradient(vector<vector<double>> &grid);
    double divergent_at_point(int x, int y, vector<vector<Vect>> &vectors);
    vector<vector<double>> divergent(vector<vector<Vect>> &vectors);
    void poisson_solver(vector<vector<double>> &div, vector<vector<double>> &x0, double tol);
    void advect(vector<vector<Vect>> &vectors, double timestep);
    void add_force(vector<vector<Vect>> &vectors, double delta, Vector2 force);
    void update_boundary(vector<vector<Vect>> &vectors);
    Vector2 bilinear_interpolation(vector<vector<Vect>> &vectors, Vector2 pos, bool pressure);
    Vector2 bilinear_interpolation_grid(Array grid, Vector2 pos, bool pressure);
    Vector2 get_minmax_velocity(Array grid);
    Vector2 get_minmax_pressure(Array grid);
    void update_particles(Array grid, Array particles, double delta);
};

#endif