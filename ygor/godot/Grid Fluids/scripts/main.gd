extends Node2D

var interface_visible = true
	
func _input(event):
	if event.is_action_pressed("show_interface"):
		if interface_visible:
			$interface.hide()
		else:
			$interface.show()
		interface_visible = not interface_visible

func check_inter_col(pos):
	return $interface/InterRect.get_rect().has_point(pos)

func _ready():
	var t1 = "Grid Dim: " + str($Grid.grid_size.x) + "X" + str($Grid.grid_size.y)
	$interface.set_grid_dim_label(t1)

func _process(delta):
	var t = "Num of Particles: " + str($Grid.particles.size())
	$interface.set_n_particles_label(t)
