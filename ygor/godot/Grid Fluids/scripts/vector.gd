extends Node2D


var velocity
var pressure

func _draw():
	$Sprite.rotation = velocity.angle() + PI/2
	$Sprite.scale.x += (abs(velocity.x) + abs(velocity.y)) / 10000
