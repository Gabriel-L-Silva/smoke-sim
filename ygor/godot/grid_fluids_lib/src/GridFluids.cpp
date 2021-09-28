#include "GridFluids.h"


void GridFluids::_register_methods()
{
    register_method("update_field", &GridFluids::update_field);
    register_method("bilinear_interpolation_grid", &GridFluids::bilinear_interpolation_grid);
    register_method("get_minmax_pressure", &GridFluids::get_minmax_pressure);
    register_method("get_minmax_velocity", &GridFluids::get_minmax_velocity);
    register_method("update_particles", &GridFluids::update_particles);
    register_method("get_density_primitive_vertex", &GridFluids::get_density_primitive_vertex);
    register_method("get_density_primitive_colors", &GridFluids::get_density_primitive_colors);
    register_method("get_density_primitive", &GridFluids::get_density_primitive);
    
    register_property<GridFluids, double>("max_speed", &GridFluids::maxSpeed, 1000.0);
    register_property<GridFluids, int>("max_iter_poisson", &GridFluids::maxIterPoisson, 20);
    register_property<GridFluids, int>("num_sub_steps", &GridFluids::subSteps, 10);
    register_property<GridFluids, double>("rho_const", &GridFluids::rho_const, 1.0);
    register_property<GridFluids, double>("diff_const", &GridFluids::diff_const, 0);
    register_property<GridFluids, double>("force_const", &GridFluids::force_const, 2.0);
    register_property<GridFluids, double>("source_const", &GridFluids::source_const, 2.0);
    register_property<GridFluids, Vector2>("tile_size", &GridFluids::tile_size, Vector2(1, 1));
    register_property<GridFluids, Vector2>("grid_size", &GridFluids::grid_size, Vector2(800, 800));
    register_property<GridFluids, Vector2>("vector_size", &GridFluids::vector_size, Vector2(0, 0));
    register_property<GridFluids, Vector2>("mouse_pos", &GridFluids::mouse_pos, Vector2(1, 1));
    register_property<GridFluids, Vector2>("prev_mouse_pos", &GridFluids::prev_mouse_pos, Vector2(-1, -1));
    register_property<GridFluids, Vector2>("source_pos", &GridFluids::source_pos, Vector2(-1,-1));
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
            vec.density = v->get("density");
            aux.push_back(vec);
        }
        vectors.push_back(aux);
    }

    return vectors;
}

//accesed at 08/09/21 11:51am https://stackoverflow.com/a/17299623
template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T>>& v) {
    std::size_t total_size = 0;
    for (const auto& sub : v)
        total_size += sub.size(); // I wish there was a transform_accumulate
    std::vector<T> result;
    result.reserve(total_size);
    for (const auto& sub : v)
        result.insert(result.end(), sub.begin(), sub.end());
    return result;
}

Vector2 GridFluids::get_minmax_velocity(Array grid)
{
    vector<vector<Vect>> vectors = copy_grid(grid);
    auto flat = flatten(vectors);
    const auto result = minmax_element(flat.begin(), flat.end(),
        [](const Vect& a, const Vect& b)
        {
            return a.vel.length() < b.vel.length();
        });
    Vect min = *(result.first);
    Vect max = *(result.second);

    return Vector2(min.vel.length(), max.vel.length());
}

Vector2 GridFluids::get_minmax_pressure(Array grid)
{
    vector<vector<Vect>> vectors = copy_grid(grid);
    auto flat = flatten(vectors);
    const auto result = minmax_element(flat.begin(), flat.end(),
        [](const Vect& a, const Vect& b)
        {
            return a.pressure < b.pressure;
        });
    Vect min = *(result.first);
    Vect max = *(result.second);
    return Vector2(min.pressure, max.pressure);
}

double GridFluids::check_divfree(vector<vector<Vect>>& vectors)
{
    auto accum = 0.0;
    for (int x = 1; x < vector_size.x - 1; x++) {
        for (int y = 1; y < vector_size.y - 1; y++) {
            accum += GridFluids::divergent_at_point(x, y, vectors);
        }
    }
    return accum;
}

Array GridFluids::get_density_primitive(Array grid)
{
    vector<vector<Vect>> vectors = copy_grid(grid);
    Array prim = Array();
    for (int x=0; x < vector_size.x-1; x++){
		for (int y=0; y < vector_size.y-1; y++)
        {
            Array vert = Array();
            Array colr = Array();
            Vector2 position = vectors[x][y].pos;
			Vector2 lb = position - tile_size;
			Vector2 lt = lb + Vector2(0,tile_size.y);
			Vector2 rb = lb + Vector2(tile_size.x, 0);
			Vector2 rt = rb + Vector2(0,tile_size.y);
            vert.append(lb);
            vert.append(rb);
            vert.append(rt);
            vert.append(lt);
			prim.append(PoolVector2Array(vert));
            colr.append(Color(1,1,1, vectors[x][y].density));
            colr.append(Color(1,1,1, vectors[x][y+1].density));
            colr.append(Color(1,1,1, vectors[x+1][y+1].density));
            colr.append(Color(1,1,1, vectors[x+1][y].density));
            prim.append(PoolColorArray(colr));
        }
    }
    return prim;
}
Array GridFluids::get_density_primitive_vertex(Array grid)
{
    vector<vector<Vect>> vectors = copy_grid(grid);
    Array vertex = Array();
    for (int x=0; x < vector_size.x-1; x++){
		for (int y=0; y < vector_size.y-1; y++)
        {
            Array arr = Array();
            Vector2 position = vectors[x][y].pos;
			Vector2 lb = position - tile_size;
			Vector2 lt = lb + Vector2(0,tile_size.y);
			Vector2 rb = lb + Vector2(tile_size.x, 0);
			Vector2 rt = rb + Vector2(0,tile_size.y);
            arr.append(lb);
            arr.append(rb);
            arr.append(rt);
            arr.append(lt);
			vertex.append(PoolVector2Array(arr));
        }
    }
    return vertex;
}

Array GridFluids::get_density_primitive_colors(Array grid)
{
    vector<vector<Vect>> vectors = copy_grid(grid);
    Array colors = Array();
    for (int x=0; x < vector_size.x-1; x++){
		for (int y=0; y < vector_size.y-1; y++)
        {
            Array arr = Array();
            arr.append(Color(1,1,1, vectors[x][y].density));
            arr.append(Color(1,1,1, vectors[x][y+1].density));
            arr.append(Color(1,1,1, vectors[x+1][y+1].density));
            arr.append(Color(1,1,1, vectors[x+1][y].density));
			colors.append(PoolColorArray(arr));
        }
    }
    return colors;
}

double GridFluids::update_field(double delta, Array grid, Vector2 externalForces)
{
    vector<vector<Vect>> vectors = copy_grid(grid);
	add_force(vectors, delta, externalForces);
	update_boundary(vectors);
	diffuse(vectors, delta, diff_const);
	update_boundary(vectors);
	advect(vectors, delta);
	update_boundary(vectors);
	project(vectors);
	update_boundary(vectors);
    update_grid(vectors, grid);
    // return check_divfree(vectors);
    return 0;
}

void GridFluids::update_grid(vector<vector<Vect>> &vectors, Array grid){
    for (int x=0; x < vector_size.x; x++){
        Array aux = grid[x];
		for (int y=0; y < vector_size.y; y++){
            Node* vec = aux[y];
            vec->set("pressure", vectors[x][y].pressure); 
            vec->set("velocity", vectors[x][y].vel);
            vec->set("density", vectors[x][y].density);
            vec->call("update");
        }
    }
}

void GridFluids::get_prev_grid(vector<vector<Vect>> &vectors, vector<vector<tuple<double,double,double,double>>> &x0)
{
    /*
        0 velx
        1 vely
        2 press
        3 dens
    */
    for (int x=0; x < vector_size.x; x++){
		for (int y=0; y < vector_size.y; y++){
            std::get<0>(x0[x][y]) = vectors[x][y].vel.x;
            std::get<1>(x0[x][y]) = vectors[x][y].vel.y;
            std::get<2>(x0[x][y]) = vectors[x][y].pressure;
			std::get<3>(x0[x][y]) = vectors[x][y].density;
        }
    }
}

void GridFluids::diffuse(vector<vector<Vect>> &vectors, double delta, double diff)
{
    vector<vector<tuple<double,double,double,double>>> x0(vector_size.x, vector<tuple<double,double,double,double>> (vector_size.y, make_tuple(0,0,0,0)));

    get_prev_grid(vectors, x0);

    double a = (double)(delta * diff * int(vector_size.x-2.0) * int(vector_size.y-2.0));
    for (int i=0; i < maxIterPoisson; i++){
		auto old = x0;
		for (int x=1; x < vector_size.x-1; x++){
			for (int y=1; y < vector_size.y-1; y++){
                auto vel = Vector2(std::get<0>(old[x][y]), std::get<1>(old[x][y]));
                auto velip1 = Vector2(std::get<0>(old[x+1][y]), std::get<1>(old[x+1][y]));
                auto velim1 = Vector2(std::get<0>(old[x-1][y]), std::get<1>(old[x-1][y]));
                auto veljp1 = Vector2(std::get<0>(old[x][y+1]), std::get<1>(old[x][y+1]));
                auto veljm1 = Vector2(std::get<0>(old[x][y-1]), std::get<1>(old[x][y-1]));

				vel = (vel + a*(velip1 + velim1 + veljp1 + veljm1))/(1+4*a);

                std::get<0>(x0[x][y]) = vel.x;
                std::get<1>(x0[x][y]) = vel.y;

                auto dens = std::get<3>(old[x][y]);
                auto densip1 = std::get<3>(old[x+1][y]);
                auto densim1 = std::get<3>(old[x-1][y]);
                auto densjp1 = std::get<3>(old[x][y+1]);
                auto densjm1 = std::get<3>(old[x][y-1]);

                dens = (dens + a*(densip1 + densim1 + densjp1 + densjm1))/(1+4*a);
                std::get<3>(x0[x][y]) = dens;
            }
        }
    }

    for (int x=1; x < vector_size.x-1; x++){
		for (int y=1; y < vector_size.y-1; y++){
            vectors[x][y].vel.x = std::get<0>(x0[x][y]);
            vectors[x][y].vel.y = std::get<1>(x0[x][y]);
            vectors[x][y].pressure = std::get<2>(x0[x][y]);
            vectors[x][y].density = std::get<3>(x0[x][y]);
        }
    }
}

void GridFluids::project(vector<vector<Vect>> &vectors)
{
    vector<vector<tuple<double,double,double,double>>> x0(vector_size.x, vector<tuple<double,double,double,double>> (vector_size.y, make_tuple(0,0,0,0)));

    get_prev_grid(vectors, x0);

	poisson_solver(divergent(vectors), x0, 10e-5);
	vector<vector<Vector2>> grad_q = gradient(x0);
	for (int x=1; x < vector_size.x-1; x++){
		for (int y=1; y < vector_size.y-1; y++){
			vectors[x][y].vel -= grad_q[x][y];

            /*if(vectors[x][y].vel.x > maxSpeed) 
                vectors[x][y].vel.x = maxSpeed; 
            if (vectors[x][y].vel.x < -maxSpeed)
                vectors[x][y].vel.x = -maxSpeed;
            if (vectors[x][y].vel.y > maxSpeed)
                vectors[x][y].vel.y = maxSpeed; 
            if (vectors[x][y].vel.y < -maxSpeed)
                vectors[x][y].vel.y = -maxSpeed; */

			vectors[x][y].pressure = rho_const/(tile_size.x + tile_size.y) * std::get<2>(x0[x][y]);
        }
    }
}

Vector2 GridFluids::gradient_at_point(int x, int y, vector<vector<tuple<double,double,double,double>>> &grid)
{
	double right = std::get<2>(grid[x][y+1]);
	double left = std::get<2>(grid[x][y-1]);
	double up = std::get<2>(grid[x+1][y]);
	double down = std::get<2>(grid[x-1][y]);
	
	return Vector2((right-left)/(tile_size.x+tile_size.y), (up-down)/(tile_size.x+tile_size.y));
}

vector<vector<Vector2>> GridFluids::gradient(vector<vector<tuple<double,double,double,double>>> &grid)
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
	
	return -1/(tile_size.x+tile_size.y) * ((right - left) + (up - down));
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

void GridFluids::poisson_solver(vector<vector<double>> &div, vector<vector<tuple<double,double,double,double>>> &x0, double tol)
{
	for (int i=0; i < maxIterPoisson; i++){
		auto old = x0;
		for (int x=1; x < vector_size.x-1; x++){
			for (int y=1; y < vector_size.y-1; y++){
				std::get<2>(x0[x][y]) = 0.25*(std::get<2>(old[x+1][y]) + std::get<2>(old[x-1][y]) + std::get<2>(old[x][y+1]) + std::get<2>(old[x][y-1]) + div[x][y]);
            }
        }
    }
}


void GridFluids::advect(vector<vector<Vect>> &vectors, double timestep)
{
    double s = timestep/subSteps;

    vector<vector<tuple<double,double,double,double>>> x0(vector_size.x, vector<tuple<double,double,double,double>> (vector_size.y, make_tuple(0,0,0,0)));

    get_prev_grid(vectors, x0);
    
    auto old = x0;
    for (int x=1; x < vector_size.x-1; x++){
		for (int y=1; y < vector_size.y-1; y++){    
			Vector2 pos = vectors[x][y].pos-tile_size;
			Vector2 vel = Vector2(std::get<0>(x0[x][y]) ,std::get<1>(x0[x][y]));
			for (int _s=0; _s < subSteps; _s++){
				pos -= s * vel;

                bool brk = false;
				if (pos.x < 0)
                {
                    pos.x = 0;
                    brk = true;
                }
				if (pos.y < 0)
                {
                    pos.y = 0;
                    brk = true;
                }
				if (pos.x > grid_size.x)
                {
                    pos.x = grid_size.x;
                    brk = true;
                }
				if (pos.y > grid_size.y)
                {
                    pos.y = grid_size.y;
                    brk = true;
                }
                if (brk)
                    break;
				vel = bilinear_interpolation(vectors, pos, false);
            }

            Vector2 p = bilinear_interpolation(vectors, pos, true);
            vel = bilinear_interpolation(vectors, pos, false);
            std::get<0>(x0[x][y]) = vel.x;
            std::get<1>(x0[x][y]) = vel.y;
            std::get<2>(x0[x][y]) = p.x;
            std::get<3>(x0[x][y]) = p.y;
        }
    }

    for (int x=1; x < vector_size.x-1; x++){
		for (int y=1; y < vector_size.y-1; y++){
            vectors[x][y].vel = Vector2(std::get<0>(x0[x][y]),std::get<1>(x0[x][y]));
            vectors[x][y].pressure = std::get<2>(x0[x][y]);
            vectors[x][y].density = std::get<3>(x0[x][y]);
        }
    }
}

Vector2 GridFluids::buoyancy(int i, int j)
{
    if (j*tile_size.y >= grid_size.y/4 - 2*tile_size.y &&
        j*tile_size.y <= grid_size.y/4 + 2*tile_size.y )
        return Vector2(0, -20);
    if (j*tile_size.y >= 3*grid_size.y/4 - 2*tile_size.y &&
        j*tile_size.y <= 3*grid_size.y/4 + 2*tile_size.y )
        return Vector2(0, -20);
    return Vector2(0,0);
}

Vector2 GridFluids::mouse_repellent(int j, int i, Vector2 pos, Vector2 vel)
{
    //TODO: força aplicada em todas as células no caminho
    
    int mouse_i = floor(mouse_pos.x / tile_size.x) + 1;
    int mouse_j = floor(mouse_pos.y / tile_size.y) + 1;
    // if (mouse_pos.x >= (i-2)*tile_size.x && mouse_pos.y >= (j-2)*tile_size.y && mouse_pos.x <= (i+2)*tile_size.x && mouse_pos.y <= (j+2)*tile_size.y)
    // if (mouse_i >= i-1 && mouse_j >= j-1 && mouse_i <= i+1 && mouse_j <= j+1)
    // {
    // cout << mouse_i << ", " << mouse_j <<endl;
    if(i>=mouse_i-1 && i<=mouse_i+1 && j>=mouse_j-1 && j<=mouse_j+1)
    {
        auto force =force_const * (mouse_pos - prev_mouse_pos);
        return force.length() != 0? force : vel;
    }
        // force = 100*f/f.length();
        // force += vel;
    // }
    // return pos-mouse_pos;
    return vel;
}

double GridFluids::get_source(int i, int j, double delta, Vector2 pos, double density)
{
    int source_i = floor(source_pos.y / tile_size.y) + 1;
    int source_j = floor(source_pos.x / tile_size.x) + 1;
    
    //cout << source_i << ", " << source_j <<endl;
    double dens = density;
    dens += delta * (i==source_i & j==source_j ? source_const : 0);
    
    return dens;
}
void GridFluids::add_force(vector<vector<Vect>> &vectors, double delta, Vector2 force)
{   
    for(int i=1; i < vector_size.x-1; i++){
        for (int j=1; j < vector_size.y-1; j++){
            // cout << vectors[i][j].pos.x << ", " << vectors[i][j].pos.y << endl;
            vectors[i][j].vel = (mouse_pos.x >= 0 && mouse_pos.y >= 0) && (prev_mouse_pos.x >= 0 && prev_mouse_pos.y >= 0)  ? mouse_repellent(i, j, vectors[i][j].pos, vectors[i][j].vel) : vectors[i][j].vel;
            // vectors[i][j].vel = mouse_repellent(i, j, vectors[i][j].pos, vectors[i][j].vel);//(force + buoyancy(i, j) + mouse_repellent(i, j, vectors[i][j].pos));
            // vectors[i][j].density += i*tile_size.x >= grid_size.x/2-3*tile_size.x & j*tile_size.y >= grid_size.y/2-3*tile_size.y & i*tile_size.x <= grid_size.x/2+3*tile_size.x & j*tile_size.y <= grid_size.y/2+3*tile_size.y ? delta : 0.0;
            vectors[i][j].density = (source_pos.x >= 0 && source_pos.y >= 0) ? get_source(i, j, delta, source_pos, vectors[i][j].density) : vectors[i][j].density;
        }
    }
}

void GridFluids::update_boundary(vector<vector<Vect>> &vectors)
{
    // vertical
    for (int x = 0; x < vector_size.x; x++) {

        vectors[x][1].vel.x = -vectors[x][2].vel.x ;
        vectors[x][vector_size.y-2].vel.x = -vectors[x][vector_size.y-3].vel.x;
        vectors[x][1].density = vectors[x][2].density;
        vectors[x][vector_size.y-2].density = vectors[x][vector_size.y-3].density;

        vectors[1][1].vel.x  = 0.5*(vectors[2][1].vel.x + vectors[1][2].vel.x);
        vectors[1][vector_size.y-2].vel.x = 0.5*(vectors[2][vector_size.y-2].vel.x + vectors[1][vector_size.y-3].vel.x); 
        vectors[vector_size.x-2][1].vel.x  = 0.5*(vectors[vector_size.x-3][1].vel.x + vectors[vector_size.x-2][2].vel.x); 
        vectors[vector_size.x-2][vector_size.y-2].vel.x  = 0.5*(vectors[vector_size.x-3][vector_size.y-2].vel.x + vectors[vector_size.x-2][vector_size.y-3].vel.x);
        // vectors[x][0].pressure = vectors[x][1].pressure;
        // vectors[x][vector_size.y-1].pressure = vectors[x][vector_size.y-2].pressure;
    }
    
    // horizontal
    for (int y = 0; y < vector_size.y; y++) {        
        vectors[1][y].vel.y = -vectors[2][y].vel.y;
        vectors[vector_size.x-2][y].vel.y = -vectors[vector_size.x-3][y].vel.y;
        vectors[1][y].density = vectors[2][y].density;
        vectors[vector_size.x-2][y].density = vectors[vector_size.x-3][y].density;

        vectors[1][1].vel.y  = 0.5*(vectors[2][1].vel.y + vectors[1][2].vel.y);
        vectors[1][vector_size.y-2].vel.y = 0.5*(vectors[2][vector_size.y-2].vel.y + vectors[1][vector_size.y-3].vel.y); 
        vectors[vector_size.x-2][1].vel.y  = 0.5*(vectors[vector_size.x-3][1].vel.y + vectors[vector_size.x-2][2].vel.y); 
        vectors[vector_size.x-2][vector_size.y-2].vel.y  = 0.5*(vectors[vector_size.x-3][vector_size.y-2].vel.y + vectors[vector_size.x-2][vector_size.y-3].vel.y);
        // vectors[0][y].pressure = vectors[1][y].pressure;
        // vectors[vector_size.x-1][y].pressure = vectors[vector_size.x-2][y].pressure;
    }
    // density    
    vectors[1][1].density  = 0.5*(vectors[2][1].density + vectors[1][2].density);
    vectors[1][vector_size.y-2].density = 0.5*(vectors[2][vector_size.y-2].density + vectors[1][vector_size.y-3].density); 
    vectors[vector_size.x-1][1].density  = 0.5*(vectors[vector_size.x-3][1].density + vectors[vector_size.x-2][2].density); 
    vectors[vector_size.x-2][vector_size.y-2].density  = 0.5*(vectors[vector_size.x-3][vector_size.y-2].density + vectors[vector_size.x-2][vector_size.y-3].density);
}

Vector2 GridFluids::bilinear_interpolation(vector<vector<Vect>> &vectors, Vector2 pos, bool pressure)
{
    // top left
    int i = floor( (pos.y - tile_size.y / 2.0) / tile_size.y) + 1;
	int j = floor( (pos.x - tile_size.x / 2.0) / tile_size.x) + 1;
	
    pos += tile_size;
    
	Vect q11 = vectors[i][j];
	Vect q12 = vectors[i+1][j];
	Vect q21 = vectors[i][j+1];
	Vect q22 = vectors[i+1][j+1];

	double x1 = q11.pos.x;
	double y1 = q11.pos.y;
	
	double x2 = q22.pos.x;
	double y2 = q22.pos.y;

    if (pressure){
        double press = (q11.pressure * (x2 - pos.x) * (y2 - pos.y) +
            q21.pressure * (pos.x - x1) * (y2 - pos.y) +
            q12.pressure * (x2 - pos.x) * (pos.y - y1) +
            q22.pressure * (pos.x - x1) * (pos.y - y1)
            ) / ((x2 - x1) * (y2 - y1));
        double dens = (q11.density * (x2 - pos.x) * (y2 - pos.y) +
            q21.density * (pos.x - x1) * (y2 - pos.y) +
            q12.density * (x2 - pos.x) * (pos.y - y1) +
            q22.density * (pos.x - x1) * (pos.y - y1)
            ) / ((x2 - x1) * (y2 - y1));

        return Vector2(press, dens);
    }else{
        return (q11.vel * (x2 - pos.x) * (y2 - pos.y) +
            q21.vel * (pos.x - x1) * (y2 - pos.y) +
            q12.vel * (x2 - pos.x) * (pos.y - y1) +
            q22.vel * (pos.x - x1) * (pos.y - y1)
            ) / ((x2 - x1) * (y2 - y1));
    }	
}

Vector2 GridFluids::bilinear_interpolation_grid(Array grid, Vector2 pos, bool pressure)
{
    // top left
    int i = (pos.y - tile_size.y/2.0)/tile_size.y + 1;
	int j = (pos.x - tile_size.x/2.0)/tile_size.x + 1;

    pos += tile_size;

    Array t;

    Array aux;
    Array aux2;

    aux = grid[i];
    aux2 = grid[i+1];

    Object* q11 = aux[j];
	Object* q12 = aux2[j];

	Object* q21 = aux[j+1];
	Object* q22 = aux2[j+1];

    Vector2 q11_pos = q11->get("pos");
    Vector2 q22_pos = q22->get("pos");

	double x1 = q11_pos.x;
	double y1 = q11_pos.y;
	
	double x2 = q22_pos.x;
	double y2 = q22_pos.y;

    if (pressure){
        double p11 = q11->get("pressure");
        double p12 = q12->get("pressure");
        double p21 = q21->get("pressure");
        double p22 = q22->get("pressure");
        
        double p = (p11 * (x2 - pos.x) * (y2 - pos.y) +
            p12 * (x2 - pos.x) * (pos.y - y1) +
            p21 * (pos.x - x1) * (y2 - pos.y) +
            p22 * (pos.x - x1) * (pos.y - y1)
            ) / ((x2 - x1) * (y2 - y1));
        
        p11 = q11->get("density");
        p12 = q12->get("density");
        p21 = q21->get("density");
        p22 = q22->get("density");
        double d = (p11 * (x2 - pos.x) * (y2 - pos.y) +
            p12 * (x2 - pos.x) * (pos.y - y1) +
            p21 * (pos.x - x1) * (y2 - pos.y) +
            p22 * (pos.x - x1) * (pos.y - y1)
            ) / ((x2 - x1) * (y2 - y1));
        
        return Vector2(p, d);
    }else{
        Vector2 v11 = q11->get("velocity");
        Vector2 v12 = q12->get("velocity");
        Vector2 v21 = q21->get("velocity");
        Vector2 v22 = q22->get("velocity");

        return (v11 * (x2 - pos.x) * (y2 - pos.y) +
            v12 * (x2 - pos.x) * (pos.y - y1) +
            v21 * (pos.x - x1) * (y2 - pos.y) +
            v22 * (pos.x - x1) * (pos.y - y1)
            ) / ((x2 - x1) * (y2 - y1));
    }	    
}

void GridFluids::update_particles(Array grid, Array particles, double delta)
{
    for(int x = 0; x < particles.size(); x++){
        Node* p = particles[x];
        Vector2 pos = p->get("position");
        double speed = p->get("speed");
        Vector2 p_size = p->get("p_size");
        
        Node* p_particles = p->get_node("Particles2D");

        Vector2 vel = bilinear_interpolation_grid(grid, pos, false);
        p->set("vel", vel);
        
        auto s_speed = Vector3(-vel.x, -vel.y, 0) * speed;
        p_particles->set_indexed("process_material:gravity", s_speed);
        
        pos += delta * vel * speed;

        if (pos.x - p_size.x < 0) 
            pos.x = p_size.x;
        if (pos.x + p_size.x > grid_size.x) 
            pos.x = grid_size.x - p_size.x;
        if (pos.y - p_size.y < 0) 
            pos.y = p_size.y;
        if (pos.y + p_size.y > grid_size.y) 
            pos.y = grid_size.y - p_size.y;

        p->set("position", pos);

        p->call("update");
    }
}
