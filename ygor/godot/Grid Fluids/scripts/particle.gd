extends Node2D


var vel = Vector2(0, -1)
var pressure = 0
var speed = 0.1
var must_remove = false

onready var grid = get_parent()

func update_velocity():
	$Particles2D.process_material.gravity = Vector3(-vel.x*speed, -vel.y*speed, 0)

func _process(delta):
	if (position.x > grid.grid_size.x or position.x < 0 or 
		position.y > grid.grid_size.y or position.y < 0):
			must_remove = true
			return
	
	update_velocity()
	self.position += vel * delta * speed
