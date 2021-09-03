extends Node2D

var grid_size 		= OS.get_window_size()
var squares_qtd 	= Vector2(16, 16)
var tile_size 		= Vector2(grid_size.x/squares_qtd.x, grid_size.y/squares_qtd.y)
var show_vectors 	= false
var show_grid 		= false
var timer = 0
var vel_min = pow(2, 63)
var vel_max = -vel_min

var rho = 1.0
var gravity = Vector2(0, 9.81)
var sub_steps = squares_qtd.x #random value, maybe be lowered for performance improvement

var Particle = preload("res://scenes/particle.tscn")
var Vector = preload("res://scenes/vector.tscn")

var grid_vectors = []
var particles = []

func get_velocity(pos):
#	return Vector2(pow(pos.x,2)*abs(cos(pos.x)), -pow(pos.y,2)*abs(sin(pos.x)))/100
	return Vector2(0,sin(pos.y))

func get_pressure(pos):
	return pos.x+pos.y

func _ready():
	for x in range(squares_qtd.x):
		grid_vectors.append([])
		for y in range(squares_qtd.y):
			# middle point
			var pos = Vector2(x*tile_size.x+tile_size.x/2, y*tile_size.y+tile_size.y/2)
			var vector = Vector.instance()
			vector.position = pos
			vector.velocity = get_velocity(pos)
			
			if vector.velocity.length() > vel_max:
				vel_max = vector.velocity.length()
			
			if vector.velocity.length() < vel_min:
				vel_min = vector.velocity.length() 
			
			vector.pressure = get_pressure(pos) 
			grid_vectors[x].append(vector)
			$vector_visualizer.add_child(vector)

func bilinear_interpolation(pos) -> Vector2:
	var dx = tile_size.x
	var dy = tile_size.y
	
	var x0 = 0.5 * dx
	var y0 = 0.5 * dy
	
	if pos.x - x0 < 0:
		pos.x = x0
	
	if pos.y - y0 < 0:
		pos.y = y0
	
	var i = int((pos.x-x0)/dx)
	var j = int((pos.y-y0)/dy)
	
	if i+1 >= squares_qtd.x:
		pos.x = grid_size.x - x0
		i = squares_qtd.x-2
	
	if j+1 >= squares_qtd.y:
		pos.y = grid_size.y - y0
		j = squares_qtd.y-2

	var xi = grid_vectors[i][j].position.x
	var yj = grid_vectors[i][j].position.y
	
	var x1 = xi
	var y1 = yj
	var q11 = grid_vectors[i][j]
	
	var x2 = xi + dx
	var _y1 = yj
	var q21 = grid_vectors[i+1][j]
	
	var _x1 = xi
	var y2 = yj + dy
	var q12 = grid_vectors[i][j+1]
	
	var _x2 = xi + dx
	var _y2 = yj + dy
	var q22 = grid_vectors[i+1][j+1]

	return (q11.velocity * (x2 - pos.x) * (y2 - pos.y) +
		q21.velocity * (pos.x - x1) * (y2 - pos.y) +
		q12.velocity * (x2 - pos.x) * (pos.y - y1) +
		q22.velocity * (pos.x - x1) * (pos.y - y1)
		) / ((x2 - x1) * (y2 - y1))

func add_particle():
	var pos = get_global_mouse_position()
	var new_particle = Particle.instance()
	new_particle.position = pos
	particles.append(new_particle)
	add_child(new_particle)
	print('particle added')

func external_forces():
	return gravity #+ bouancy (+ mouse_reppelant_force)

func add_force(w0: Array, delta):
	var w1 = w0.duplicate(true)
	for x in range(squares_qtd.x):
		for y in range(squares_qtd.y):
			w1[x][y].velocity = w1[x][y].velocity + delta * external_forces()
	return w1

func advect(w1: Array, timestep):
	var w2 = w1.duplicate(true)
	var s = timestep/sub_steps
	for x in range(squares_qtd.x):
		for y in range(squares_qtd.y):
			var pos = grid_vectors[x][y].position
			for _s in range(sub_steps):
				if pos.x > grid_size.x and pos.y > grid_size.y:
					break
				pos = pos - s*bilinear_interpolation(pos)
			w2[x][y].velocity = bilinear_interpolation(pos)
	return w2

func diffuse(w2):
	return w2 #Passthrough function because smoke has no viscosity

func poisson_solver(grid, x0, tol):
	var max_it = 10000
	var m = squares_qtd.x
	var n = squares_qtd.y
	var dx = tile_size.x
	var dy = tile_size.y
#	for idx_x in range(squares_qtd.x):
#		for idx_y in range(squares_qtd.y):
#			var x = vec_field.pos.x
#			var y = vec_field.pos.y
#			var p = vec_field.pressure
#			if idx_y == 0:
#				x0[idx_y, idx_x] = np.sin(np.pi * x)

	var A = grid.duplicate(true)
	var it = 0
	for _i in range(max_it):
		var old = x0.duplicate(true)
		for i in range(m):
			for j in range(n):
				#boundary condition
				var xpp = 0
				var ypp = 0
				var center = A[i][j].pressure
				var dx_sqr_inv = 1/pow(dx,2)
				var dy_sqr_inv = 1/pow(dy,2)
				
				if i == m-1 or i == 0:
					if i == 0:
						#esquerda
						var right = x0[i+1][j]
						
						xpp = (right - center)*dx_sqr_inv
					if i== m-1:
						#direita
						var left = x0[i-1][j]
						
						xpp = (left - center)*dx_sqr_inv
				else:
					var right = x0[i+1][j]
					var left = x0[i-1][j]
					xpp = (right + left - 2*center)*dx_sqr_inv
				
				if j == n-1 or j == 0:
					if j == 0:
						#baixo
						var up = x0[i][j+1]
						
						ypp = (up - center)*dy_sqr_inv
					if j == n-1:
						#cima
						var down = x0[i][j-1]
						
						ypp = (down - center)*dy_sqr_inv
				else:
					var up = x0[i][j+1]
					var down = x0[i][j-1]
					ypp = (up + down - 2*center)*dy_sqr_inv
					
				x0[i][j] = xpp + ypp
		
		var accum = 0.0
		print(old[0][0]," , ", x0[0][0]," , ", A[0][0].pressure )
		for x in range(m):
			for y in range(n):
				accum += abs(old[x][y]-x0[x][y])
		if accum < tol:
			it = _i
			break
		old = x0.duplicate(true)
	print(it)
	
	#Fazendo isso pq o godot nos odeia
	for i in range(squares_qtd.x):
		for j in range(squares_qtd.y):
			A[i][j].pressure = x0[i][j]
	return A

func divergent(grid):
	var m = squares_qtd.x
	var n = squares_qtd.y
	
	var div = grid.duplicate(true)
	
	for i in range(m):
		for j in range(n):
			#iterate over columns which is equivalent to iterate over x axis
			div[i][j].pressure = divergent_at_point(grid, i,j)
	
	return div

func divergent_at_point(grid, x, y):
	var m = squares_qtd.x
	var n = squares_qtd.y
	var dx = tile_size.x * (2 if x!= m-1 or x != 0 else 1)
	var dy = tile_size.y * (2 if y != n-1 or y != 0 else 1)
	
	var center = grid[x][y].velocity

	var right = grid[x+1][y].velocity.x if x != m-1 else -center.x
	var left = grid[x-1][y].velocity.x if x != 0 else -center.x
	var up = grid[x][y+1].velocity.y if y != n-1 else -center.y
	var down = grid[x][y-1].velocity.y if y != 0 else -center.y
	
	return (right - left) / dx + (up - down) / dy

func gradient(grid):
	var m = squares_qtd.x
	var n = squares_qtd.y
	
	var grad = grid.duplicate(true)
	for x in range(m):
		for y in range(n):
			#iterate over columns which is equivalent to iterate over x axis
			grad[x][y].velocity = gradient_at_point(grid, x, y)
	
	return grad

func gradient_at_point(grid, x, y):	
	var m = squares_qtd.x
	var n = squares_qtd.y
	
	var dx = tile_size.x * (2 if x != m-1 or x != 0 else 1)
	var dy = tile_size.y * (2 if y != n-1 or y != 0 else 1)
	
	var center = grid[x][y].pressure

	var right = grid[x+1][y].pressure if x != m-1 else center
	var left = grid[x-1][y].pressure if x != 0 else center
	var up = grid[x][y+1].pressure if y != n-1 else center
	var down = grid[x][y-1].pressure if y != 0 else center
	
	return Vector2((right-left)/dx, (up-down)/dy)

func project(w3):
	#Refactor make func
	var x0 = []
	for x in range(squares_qtd.x):
		x0.append([])
		for y in range(squares_qtd.y):
			x0[x].append(1)
	
	var q = poisson_solver(divergent(w3), x0, 1e-6)
	var acc_w3 = 0.0
	var acc_q = 0.0
	for x in range(squares_qtd.x):
		for y in range(squares_qtd.y):
			acc_w3 += w3[x][y].pressure
			acc_q += q[x][y].pressure
	if acc_w3 == acc_q:
		print("deu ruim")
	var grad_q = gradient(q)
	var w4 = w3.duplicate(true)
	for x in range(squares_qtd.x):
		for y in range(squares_qtd.y):
			w4[x][y].velocity = w3[x][y].velocity - grad_q[x][y].velocity
			w4[x][y].pressure = rho/(0.5*(tile_size.x+tile_size.y)) * q[x][y].pressure
			print(w4[x][y].velocity, w4[x][y].pressure)
	return w4

func update_field(delta):
	var w0 = grid_vectors.duplicate(true)
	var w1 = add_force(w0, delta)
	var w2 = advect(w1, delta)
	var w3 = diffuse(w2)
	var w4 = project(w3)
	return w4

func _process(delta):
	if get_parent().interface_visible:
		return
	
	timer += delta
	
	grid_vectors = update_field(delta/5.0) #division so animation is slower
	
	if Input.is_action_pressed("left_click") and not get_parent().interface_visible:
		if timer >= 0.05:
			add_particle()
			timer = 0
	
	for p in particles:
		if p.must_remove:
			particles.erase(p)
			remove_child(p)
		else:
			p.vel = bilinear_interpolation(p.position)

func _on_interface_show_grid_signal():
	show_grid = not show_grid
	$grid_visualizer.update()

func _on_interface_show_vectors_signal():
	show_vectors = not show_vectors
	if show_vectors:
		$vector_visualizer.show()
	else:
		$vector_visualizer.hide()
