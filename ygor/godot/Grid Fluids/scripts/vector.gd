extends Node2D


var velocity: Vector2
var pressure: float
var pos: Vector2

onready var grid = get_parent().get_parent()

func _process(delta):
	update()

func _ready():
	position.x = pos.x - grid.tile_size.x
	position.y = pos.y - grid.tile_size.y

func _draw():
	var l = velocity.length() / grid.MAX_VELOCITY
	var c = Color.from_hsv(170/360.0, 1, l)
	
	$Sprite.modulate = c
	$Sprite.rotation = velocity.angle() + PI/2
	$Sprite.scale = grid.tile_size/$Sprite.texture.get_size().length()
