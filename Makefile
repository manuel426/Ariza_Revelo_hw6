all: trayectoria.pdf

trayectoria.pdf: trayectoria_E_alpha.dat
	python plot3d.py
	python plot.py

trayectoria_E_alpha.dat: particle_in_field
	./particle_in_field 3 30

particle_in_field: particle_in_field.c
	gcc particle_in_field.c -o particle_in_field -lm

clean:
	rm trayectoria.pdf
	rm trayectoria_E_alpha.dat
