#include <cmath>

#include "GridFluids.h"



void GridFluids::_register_methods()
{
    register_method("update_field", &GridFluids::update_field);
    
    register_property<GridFluids, double>("max_speed", &GridFluids::maxSpeed, 1000.0);
    register_property<GridFluids, int>("max_iter_poisson", &GridFluids::maxIterPoisson, 20);
    register_property<GridFluids, int>("num_sub_steps", &GridFluids::subSteps, 10);
    register_property<GridFluids, double>("rho_const", &GridFluids::rho, 1.0);
    register_property<GridFluids, Vector2>("tile_size", &GridFluids::tile_size, Vector2(1, 1));
    register_property<GridFluids, Vector2>("grid_size", &GridFluids::grid_size, Vector2(800, 800));
    register_property<GridFluids, Vector2>("vector_size", &GridFluids::vector_size, Vector2(0, 0));
}

void GridFluids::_init()
{
    
}

vector<vector<Vect>> copy_grid(Array grid)
{
    vector<vector<Vect>> vectors;


    for (int x = 0; x < grid.size(); x++){
        Array inner = grid[x];
        vector<Vect> aux;
        for (int y = 0; y < inner.size(); y++){
            Object *v = inner[y];
            Vect vec;
            vec.pos = v->get("pos");
            vec.pressure = v->get("pressure");
            vec.vel = v->get("velocity");

            aux.push_back(vec);
        }
        vectors.push_back(aux);
    }

    return vectors;
}

Array GridFluids::update_field(double delta, Array grid, Vector2 externalForces)
{
    vector<vector<Vect>> vectors = copy_grid(grid);
	add_force(vectors, delta, externalForces);
	update_boundary(vectors);
	advect(vectors, delta);
	update_boundary(vectors);
	// diffuse(vectors);
	project(vectors);
	update_boundary(vectors);
    return get_update_grid(vectors);
}

Array GridFluids::get_update_grid(vector<vector<Vect>> &vectors){
    Array grid_updated;
    for (int x=0; x < vector_size.x; x++){
        Array aux;
		for (int y=0; y < vector_size.y; y++){
            Array aux2;
            aux2.append(vectors[x][y].pressure);
            aux2.append(vectors[x][y].vel);
            aux.append(aux2);
        }
        grid_updated.append(aux);
    }
    return grid_updated;
}

void GridFluids::project(vector<vector<Vect>> &vectors)
{
    vector<vector<double>> x0(vector_size.x, vector<double> (vector_size.y, 0));
	
	poisson_solver(divergent(vectors), x0, 10e-5);
	vector<vector<Vector2>> grad_q = gradient(x0);
	for (int x=1; x < vector_size.x-1; x++){
		for (int y=1; y < vector_size.y-1; y++){
			vectors[x][y].vel -= grad_q[x][y];

            if(vectors[x][y].vel.x > maxSpeed) 
                vectors[x][y].vel.x = maxSpeed; 
            if (vectors[x][y].vel.x < -maxSpeed)
                vectors[x][y].vel.x = -maxSpeed;
            if (vectors[x][y].vel.y > maxSpeed)
                vectors[x][y].vel.y = maxSpeed; 
            if (vectors[x][y].vel.y < -maxSpeed)
                vectors[x][y].vel.y = -maxSpeed; 

			vectors[x][y].pressure = rho/(0.5*(tile_size.x + tile_size.y)) * x0[x][y];
        }
    }
}

Vector2 GridFluids::gradient_at_point(int x, int y, vector<vector<double>> &grid)
{
	double right = grid[x][y+1];
	double left = grid[x][y-1];
	double up = grid[x+1][y];
	double down = grid[x-1][y];
	
	return 0.5 * Vector2((right-left)/tile_size.x, (up-down)/tile_size.y);
}

vector<vector<Vector2>> GridFluids::gradient(vector<vector<double>> &grid)
{
    vector<vector<Vector2>> grad(vector_size.x, vector<Vector2>(vector_size.y, Vector2(0, 0)));
	for (int x=1; x < vector_size.x-1; x++){
		for (int y=1; y < vector_size.y-1; y++){
			grad[x][y] = gradient_at_point(x, y, grid);
        }
    }
	return grad;
}

double GridFluids::divergent_at_point(int x, int y, vector<vector<Vect>> &vectors)
{
	double right = vectors[x][y+1].vel.x; 
	double left = vectors[x][y-1].vel.x; 
	double up = vectors[x+1][y].vel.y; 
	double down = vectors[x-1][y].vel.y;
	
	return -0.5 * tile_size.x * (right - left) + (up - down);
}

vector<vector<double>> GridFluids::divergent(vector<vector<Vect>> &vectors)
{
    vector<vector<double>> div(vector_size.x, vector<double>(vector_size.y, 0));
	for (int x=1; x < vector_size.x-1; x++){
		for (int y=1; y < vector_size.y-1; y++){
			// iterate over columns which is equivalent to iterate over x axis
			div[x][y] = divergent_at_point(x, y, vectors);
        }
    }
	return div;
}

void GridFluids::poisson_solver(vector<vector<double>> &div, vector<vector<double>> &x0, double tol)
{
	for (int i=0; i < maxIterPoisson; i++){
		vector<vector<double>> old = x0;
		double accum = 0.0;
		for (int x=1; x < vector_size.x-1; x++){
			for (int y=1; y < vector_size.y-1; y++){
				x0[x][y] = 0.25*(x0[x+1][y] + x0[x-1][y] + x0[x][y+1] + x0[x][y-1] + div[x][y]);
				accum += abs(old[x][y] - x0[x][y]);
            }
        }
		if (accum < tol)
			break;
    }
}

void GridFluids::advect(vector<vector<Vect>> &vectors, double timestep)
{
    double s = timestep/subSteps;

    for (int x=1; x < vector_size.x-1; x++){
		for (int y=1; y < vector_size.y-1; y++){    
			Vector2 pos = vectors[x][y].pos;
			Vector2 vel = vectors[x][y].vel;
			for (int _s=0; _s < subSteps; _s++){
				pos -= s * vel;
				if (pos.x < 0)
                    pos.x = 0;
				if (pos.y < 0)
                    pos.y = 0;
				if (pos.x > grid_size.x)
                    pos.x = grid_size.x;
				if (pos.y > grid_size.y)
                    pos.y = grid_size.y;
				vel = bilinear_interpolation(vectors, pos.x, pos.y, false);
            }
            Vector2 p = bilinear_interpolation(vectors, pos.x, pos.y, true);
			vectors[x][y].pressure = p.x;        }
    }
}

void GridFluids::add_force(vector<vector<Vect>> &vectors, double delta, Vector2 force)
{   
    for(int i=1; i < vector_size.x-1; i++){
        for (int j=1; j < vector_size.y-1; j++){
            vectors[i][j].vel += delta * force;
        }
    }
}

void GridFluids::update_boundary(vector<vector<Vect>> &vectors)
{
    // vertical
    for (int x = 0; x < vector_size.x; x++) {
        vectors[x][0].vel = vectors[x][1].vel;
        vectors[x][0].pressure = vectors[x][1].pressure;

        vectors[x][vector_size.y-1].vel = vectors[x][vector_size.y-2].vel;
        vectors[x][vector_size.y-1].pressure = vectors[x][vector_size.y-2].pressure;
    }
    
    // horizontal
    for (int y = 0; y < vector_size.y; y++) {
        vectors[0][y].vel = vectors[1][y].vel;
        vectors[0][y].pressure = vectors[1][y].pressure;

        vectors[vector_size.x-1][y].vel = vectors[vector_size.x-2][y].vel;
        vectors[vector_size.x-1][y].pressure = vectors[vector_size.x-2][y].pressure;
    }
}

Vector2 GridFluids::bilinear_interpolation(vector<vector<Vect>> &vectors, int x, int y, bool pressure)
{
    int i = floor( (y - tile_size.y / 2.0) / tile_size.y) + 1;
	int j = floor( (x - tile_size.x / 2.0) / tile_size.x) + 1;
	
	Vect q11 = vectors[i][j];
	Vect q12 = vectors[i][j+1];
	Vect q21 = vectors[i+1][j];
	Vect q22 = vectors[i+1][j+1];
	
	double x1 = q11.pos.x;
	double y1 = q11.pos.y;
	
	double x2 = x1 + tile_size.x;
	double y2 = y1 + tile_size.y;

    if (pressure){
        double result = (q11.pressure * (x2 - x) * (y2 - y) +
            q21.pressure * (x - x1) * (y2 - y) +
            q12.pressure * (x2 - x) * (y - y1) +
            q22.pressure * (x - x1) * (y - y1)
            ) / ((x2 - x1) * (y2 - y1));
        
        return Vector2(result, 0);
    }else{
        return (q11.vel * (x2 - x) * (y2 - y) +
            q21.vel * (x - x1) * (y2 - y) +
            q12.vel * (x2 - x) * (y - y1) +
            q22.vel * (x - x1) * (y - y1)
            ) / ((x2 - x1) * (y2 - y1));
    }	
}
