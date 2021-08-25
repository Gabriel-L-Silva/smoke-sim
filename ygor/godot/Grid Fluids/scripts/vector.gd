extends Node2D


var velocity
var pressure

func _draw():
	$Sprite.rotate(velocity.angle())
