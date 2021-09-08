extends Node2D


var vel = Vector2(0, -1)
var pressure = 0
var speed = 0.1

onready var grid = get_parent()
onready var p_size = ($Particles2D.texture.get_size() * $Particles2D.scale)/2
