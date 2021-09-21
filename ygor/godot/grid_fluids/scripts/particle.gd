extends Node2D


var vel = Vector2(0, -1)
var speed = 0.1

onready var p_size = ($Particles2D.texture.get_size() * $Particles2D.scale)/2
