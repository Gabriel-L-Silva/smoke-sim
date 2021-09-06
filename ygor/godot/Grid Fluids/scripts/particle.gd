extends Node2D


var vel = Vector2(0, -1)
var pressure = 0
var speed = 0.02

onready var grid = get_parent()

func update_velocity():
	$Particles2D.process_material.gravity = Vector3(-vel.x*speed, -vel.y*speed, 0)

func _process(delta):
	if (position.x > grid.grid_size.x or position.x <= 0 or 
		position.y > grid.grid_size.y or position.y <= 0):
			grid.particles.erase(self)
			grid.remove_child(self)
			return
			
	vel = grid.bilinear_interpolation_vel(position)
	update_velocity()
	position += vel * delta * speed
