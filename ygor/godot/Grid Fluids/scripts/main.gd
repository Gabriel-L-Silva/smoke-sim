extends Node2D


func _input(event):
	if event.is_action_pressed("show_interface"):
		if $interface.visible:
			$interface.hide()
		else:
			$interface.show()

