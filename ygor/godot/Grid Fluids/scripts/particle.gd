extends Node2D


var vel = Vector2(0, -0.5)
var pressure = 0

func _ready():
	$Particles2D.process_material.gravity = Vector3(0, 100, 0)

func _process(delta):
	self.position += vel
