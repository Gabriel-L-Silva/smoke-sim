extends Node2D


var velocity: Vector2
var pressure: float

onready var grid = get_parent().get_parent()

func _draw():
	var l = (velocity.length() - grid.vel_min) / (grid.vel_max - grid.vel_min)
	var c = Color.from_hsv(170/360.0, 1, l)
	$Sprite.modulate = c
	$Sprite.rotation = velocity.angle() + PI/2
