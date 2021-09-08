extends Node2D


var vel = Vector2(0, -1)
var pressure = 0
var speed = 0.1

onready var grid = get_parent()
onready var p_size = ($Particles2D.texture.get_size() * $Particles2D.scale)/2

func update_velocity():
	$Particles2D.process_material.gravity = Vector3(-vel.x*speed, -vel.y*speed, 0)

func _process(delta):
	vel = grid.bilinear_interpolation_vel(position)
	update_velocity()
	position += vel * delta * speed
	
	if position.x - p_size.x < 0: 
		position.x = p_size.x
	if position.x + p_size.x > grid.grid_size.x: 
		position.x = grid.grid_size.x - p_size.x
	if position.y - p_size.y < 0: 
		position.y = p_size.y
	if position.y + p_size.y > grid.grid_size.y: 
		position.y = grid.grid_size.y - p_size.y
