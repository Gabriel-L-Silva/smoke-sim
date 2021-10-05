extends Node2D


var speed = 8.0
var vel = Vector2(0, 0)
onready var p_size = ($Particles2D.texture.get_size() * $Particles2D.scale)/2
