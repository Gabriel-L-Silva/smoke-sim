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
    register_property<GridFluids, bool>("bilinear_interp", &GridFluids::bilinear_interp, true);
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


void GridFluids::update_field(double delta, Array grid, Array particles, Vector2 externalForces)
{

    if (bilinear_interp){
        initial = 1;
        end = 1;
    }else{
        initial = 2;
        end = 2;
    }

    vector<vector<Vect>> vectors = copy_grid(grid);
	add_force(vectors, delta, externalForces);
	update_boundary(vectors);
	advect(vectors, delta);
	update_boundary(vectors);
	// diffuse(vectors);
	project(vectors);
	update_boundary(vectors);
    update_grid(vectors, grid);
    update_particles(vectors, particles, delta);
}

void GridFluids::update_grid(vector<vector<Vect>> &vectors, Array grid){
    for (int x=0; x < vector_size.x; x++){
        Array aux = grid[x];
		for (int y=0; y < vector_size.y; y++){
            Node* vec = aux[y];
            vec->set("pressure", vectors[x][y].pressure); 
            vec->set("velocity", vectors[x][y].vel);
            vec->call("update");
        }
    }
}

void GridFluids::project(vector<vector<Vect>> &vectors)
{
    vector<vector<double>> x0 (vector_size.x, vector<double> (vector_size.y, 0.0));

	poisson_solver(divergent(vectors), x0, 10e-5);
	vector<vector<Vector2>> grad_q = gradient(x0);

	for (int x=initial; x < vector_size.x-end; x++){
		for (int y=initial; y < vector_size.y-end; y++){
			vectors[x][y].vel -= grad_q[x][y];

            if(vectors[x][y].vel.x > maxSpeed) 
                vectors[x][y].vel.x = maxSpeed; 
            if (vectors[x][y].vel.x < -maxSpeed)
                vectors[x][y].vel.x = -maxSpeed;
            if (vectors[x][y].vel.y > maxSpeed)
                vectors[x][y].vel.y = maxSpeed; 
            if (vectors[x][y].vel.y < -maxSpeed)
                vectors[x][y].vel.y = -maxSpeed;

			vectors[x][y].pressure = rho/(tile_size.x + tile_size.y) * x0[x][y];
        }
    }
}

Vector2 GridFluids::gradient_at_point(int x, int y, vector<vector<double>> &grid)
{
	double right = grid[x][y+1];
	double left = grid[x][y-1];
	double up = grid[x+1][y];
	double down = grid[x-1][y];
	
	return Vector2((right-left)/(tile_size.x+tile_size.y), (up-down)/(tile_size.x+tile_size.y));
}

vector<vector<Vector2>> GridFluids::gradient(vector<vector<double>> &grid)
{
    vector<vector<Vector2>> grad(vector_size.x, vector<Vector2>(vector_size.y, Vector2(0, 0)));
	for (int x=initial; x < vector_size.x-end; x++){
		for (int y=initial; y < vector_size.y-end; y++){
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
	
	return -1/(tile_size.x+tile_size.y) * ((right - left) + (up - down));
}

vector<vector<double>> GridFluids::divergent(vector<vector<Vect>> &vectors)
{
    vector<vector<double>> div(vector_size.x, vector<double>(vector_size.y, 0));
	for (int x=initial; x < vector_size.x-end; x++){
		for (int y=initial; y < vector_size.y-end; y++){
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
		for (int x=initial; x < vector_size.x-end; x++){
			for (int y=initial; y < vector_size.y-end; y++){
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

    vector<vector<Vect>> old_vector = vectors;

    for (int x=initial; x < vector_size.x-end; x++){
		for (int y=initial; y < vector_size.y-end; y++){    
			Vector2 pos = vectors[x][y].pos - tile_size;
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

                if (bilinear_interp)
				    vel = bilinear_interpolation(old_vector, pos, false);
                else
                    vel = catmull_rom_interp(old_vector, pos, false);

            }
            Vector2 p;
            if (bilinear_interp)
                p = bilinear_interpolation(old_vector, pos, true);
            else
                p = catmull_rom_interp(old_vector, pos, true);

			vectors[x][y].pressure = p.x;

            if (bilinear_interp) 
                vectors[x][y].vel = bilinear_interpolation(old_vector, pos, false);
            else
                vectors[x][y].vel = catmull_rom_interp(old_vector, pos, false);
        }
    }
}

void GridFluids::add_force(vector<vector<Vect>> &vectors, double delta, Vector2 force)
{   
    for(int i=initial; i < vector_size.x-end; i++){
        for (int j=initial; j < vector_size.y-end; j++){
            vectors[i][j].vel += delta * force;
        }
    }
}

void GridFluids::update_boundary(vector<vector<Vect>> &vectors)
{
    // vertical
    if (bilinear_interp){
        for (int x = 0; x < vector_size.x; x++) {
            vectors[x][0].vel = -vectors[x][1].vel;
            vectors[x][0].pressure = vectors[x][1].pressure;

            vectors[x][vector_size.y-1].vel = -vectors[x][vector_size.y-2].vel;
            vectors[x][vector_size.y-1].pressure = vectors[x][vector_size.y-2].pressure;
        }
    }else{
        for (int x = 0; x < vector_size.x; x++) {
            vectors[x][0].vel = -vectors[x][2].vel;
            vectors[x][0].pressure = vectors[x][2].pressure;
            vectors[x][1].vel = -vectors[x][2].vel;
            vectors[x][1].pressure = vectors[x][2].pressure;

            vectors[x][vector_size.y-1].vel = -vectors[x][vector_size.y-3].vel;
            vectors[x][vector_size.y-1].pressure = vectors[x][vector_size.y-3].pressure;
            vectors[x][vector_size.y-2].vel = -vectors[x][vector_size.y-3].vel;
            vectors[x][vector_size.y-2].pressure = vectors[x][vector_size.y-3].pressure;
        }
    }
    
    
    // horizontal
    if (bilinear_interp){
        for (int y = 0; y < vector_size.y; y++) {
            int fact = -1;
            
            vectors[0][y].pressure = vectors[1][y].pressure;
            vectors[vector_size.x-1][y].pressure = vectors[vector_size.x-2][y].pressure;

            if (y == 0)
                fact = 1;

            vectors[0][y].vel = vectors[1][y].vel * fact;
            vectors[vector_size.x-1][y].vel = vectors[vector_size.x-2][y].vel * fact;
        }
    }else{
        for (int y = 0; y < vector_size.y; y++) {
            int fact = -1;
            
            vectors[0][y].pressure = vectors[2][y].pressure;
            vectors[vector_size.x-1][y].pressure = vectors[vector_size.x-3][y].pressure;
            vectors[1][y].pressure = vectors[2][y].pressure;
            vectors[vector_size.x-2][y].pressure = vectors[vector_size.x-3][y].pressure;

            if (y == 0)
                fact = 1;

            vectors[0][y].vel = vectors[2][y].vel * fact;
            vectors[vector_size.x-1][y].vel = vectors[vector_size.x-3][y].vel * fact;
            vectors[1][y].vel = vectors[2][y].vel * fact;
            vectors[vector_size.x-2][y].vel = vectors[vector_size.x-3][y].vel * fact;
        }
    }
}

Vector2 GridFluids::bilinear_interpolation(vector<vector<Vect>> &vectors, Vector2 pos, bool pressure)
{
    // top left
    int i = floor( (pos.y - tile_size.y / 2.0) / tile_size.y) + 1;
	int j = floor( (pos.x - tile_size.x / 2.0) / tile_size.x) + 1;

    pos += tile_size;

    Vect q11 = vectors[i+1][j];
	Vect q12 = vectors[i][j];
	Vect q21 = vectors[i+1][j+1];
	Vect q22 = vectors[i][j+1];
	
	double x1 = q11.pos.x;
	double y1 = q11.pos.y;
	
	double x2 = q22.pos.x;
	double y2 = q22.pos.y;

    if (pressure){
        double result = (q11.pressure * (x2 - pos.x) * (y2 - pos.y) +
            q12.pressure * (x2 - pos.x) * (pos.y - y1) +
            q21.pressure * (pos.x - x1) * (y2 - pos.y) +
            q22.pressure * (pos.x - x1) * (pos.y - y1)
            ) / ((x2 - x1) * (y2 - y1));
        
        return Vector2(result, 0);
    }else{
        return (q11.vel * (x2 - pos.x) * (y2 - pos.y) +
            q12.vel * (x2 - pos.x) * (pos.y - y1) +
            q21.vel * (pos.x - x1) * (y2 - pos.y) +
            q22.vel * (pos.x - x1) * (pos.y - y1)
            ) / ((x2 - x1) * (y2 - y1));
    }	
}

void GridFluids::update_particles(vector<vector<Vect>> &vectors, Array particles, double delta)
{
    for (int x = 0; x < particles.size(); x++){
        Node* p = particles[x];

        Vector2 position = p->get("position");
        Vector2 vel = p->get("vel");
        double speed = p->get("speed");
        Vector2 p_size = p->get("p_size");

        Node* p_particles = p->get_node("Particles2D");

        // Vector2 vel = bilinear_interpolation(vectors, position, false);
        Vector2 pos = position;
        double s = delta/subSteps;
        for (int _s=0; _s < subSteps; _s++){
            pos -= s * vel * speed;
            if (pos.x < 0)
                pos.x = 0;
            if (pos.y < 0)
                pos.y = 0;
            if (pos.x > grid_size.x)
                pos.x = grid_size.x;
            if (pos.y > grid_size.y)
                pos.y = grid_size.y;
            
            if (bilinear_interp)
                vel = bilinear_interpolation(vectors, pos, false);
            else
                vel = catmull_rom_interp(vectors, pos, false);
        }
        
        if (bilinear_interp)
            vel = bilinear_interpolation(vectors, pos, false);
        else
            vel = catmull_rom_interp(vectors, pos, false);

        p->set("vel", vel);

        Vector3 s_speed = Vector3(-vel.x, -vel.y, 0) * speed;
        p_particles->set_indexed("process_material:gravity", s_speed);
        p_particles->set_indexed("process_material:angular_velocity", (vel.x + vel.y)/10.0);

        position += delta * vel * speed;
        if (position.x < 0 + p_size.x)
            position.x = p_size.x;
        if (position.x > grid_size.x - p_size.x)
            position.x = grid_size.x - p_size.x ;
        if (position.y < 0 + p_size.y)
            position.y = p_size.y;
        if (position.y > grid_size.y - p_size.y)
            position.y = grid_size.y - p_size.y;

        p->set("position", position);
        p->call("update");
    }
}


Vector2 catmull_rom_aux(float f, Vector2 f0, Vector2 f1, Vector2 f2, Vector2 f3){
    Vector2 d1 = (f2 - f0) / 2.0;
    Vector2 d2 = (f3 - f1) / 2.0;

    Vector2 D = f2 - f1;

    Vector2 a3 = d1 + d2 - 2*D;
    Vector2 a2 = 3 * D - 2 * d1 - d2;
    Vector2 a1 = d1;
    Vector2 a0 = f1;

    return a3 * pow(f, 3) + a2 * pow(f, 2) + a1*f + a0;
}

Vector2 GridFluids::catmull_rom_interp(vector<vector<Vect>> &vectors, Vector2 pos, bool pressure){
    // top left
    int i = floor( (pos.y - tile_size.y / 2.0) / tile_size.y) + 2;
	int j = floor( (pos.x - tile_size.x / 2.0) / tile_size.x) + 2;

    Vector2 res[4];

    for (int k = 0; k<4; k++){
        Vect q11 = vectors[i+k-1][j-1];
        Vect q12 = vectors[i+k-1][j];
        Vect q21 = vectors[i+k-1][j+1];
        Vect q22 = vectors[i+k-1][j+2];    

        float f = pos.x/grid_size.x;
        
        if (pressure)
            res[k] = catmull_rom_aux(f, Vector2(q11.pressure, 0), Vector2(q12.pressure, 0), Vector2(q21.pressure, 0), Vector2(q22.pressure, 0));
        else
            res[k] = catmull_rom_aux(f, q11.vel, q12.vel, q21.vel, q22.vel);
        
    }

    float f = pos.y/grid_size.y;

    return catmull_rom_aux(f, res[0], res[1], res[2], res[3]);
}
