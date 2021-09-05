extends Node2D

var interface_visible = true
	
func _input(event):
	if event.is_action_pressed("show_interface"):
		if interface_visible:
			$interface.hide()
		else:
			$interface.show()
		interface_visible = not interface_visible
