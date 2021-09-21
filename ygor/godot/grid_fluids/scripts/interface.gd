extends Control


signal show_grid_signal
signal show_vectors_signal

func set_grid_dim_label(t):
	$InterRect/VBoxContainer/Grid_dim.text = t
	
func set_n_particles_label(t):
	$InterRect/VBoxContainer/N_Particles.text = t

func _on_show_grid_pressed():
	emit_signal("show_grid_signal")

func _on_show_vectors_pressed():
	emit_signal("show_vectors_signal")
