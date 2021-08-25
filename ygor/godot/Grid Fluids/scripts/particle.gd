extends Node2D


var vel = Vector2(0, -1)
var pressure = 0
var speed = 0.1

func update_velocity():
	$Particles2D.process_material.gravity = Vector3(-vel.x*speed, -vel.y*speed, 0)

func _process(delta):
	update_velocity()
	
	self.position += vel * delta * speed
