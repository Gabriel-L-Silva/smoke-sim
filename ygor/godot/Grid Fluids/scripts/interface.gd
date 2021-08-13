extends Control


signal show_grid_signal
signal show_vectors_signal

func _on_show_grid_pressed():
	emit_signal("show_grid_signal")

func _on_show_vectors_pressed():
	emit_signal("show_vectors_signal")

