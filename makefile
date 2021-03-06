all: runge_kutta.x
 poblaciones_30_20.dat
 poblaciones_29_20.dat
 poblaciones_28_20.dat
 poblaciones_28_20.dat
 poblaciones_27_20.dat
 poblaciones_26_20.dat
 poblaciones_25_20.dat
 poblaciones_24_20.dat
 poblaciones_23_20.dat
 poblaciones_22_20.dat
 poblaciones_21_20.dat
 poblaciones_20_20.dat
 poblaciones_19_20.dat
 poblaciones_18_20.dat
 poblaciones_17_20.dat
 poblaciones_16_20.dat
 poblaciones_15_20.dat
 poblaciones_14_20.dat
 poblaciones_13_20.dat
 poblaciones_12_20.dat
 poblaciones_11_20.dat
 poblaciones_10_20.dat
 poblaciones_9_20.dat
 poblaciones_8_20.dat
 poblaciones_7_20.dat
 poblaciones_6_20.dat
 poblaciones_5_20.dat
 poblaciones_4_20.dat
 poblaciones_3_20.dat
 poblaciones_2_20.dat
 poblaciones_1_20.dat
 poblaciones_30_20.dat.pdf
 poblaciones_29_20.dat.pdf
 poblaciones_28_20.dat.pdf
 poblaciones_28_20.dat.pdf
 poblaciones_27_20.dat.pdf
 poblaciones_26_20.dat.pdf
 poblaciones_25_20.dat.pdf
 poblaciones_24_20.dat.pdf
 poblaciones_23_20.dat.pdf
 poblaciones_22_20.dat.pdf
 poblaciones_21_20.dat.pdf
 poblaciones_20_20.dat.pdf
 poblaciones_19_20.dat.pdf
 poblaciones_18_20.dat.pdf
 poblaciones_17_20.dat.pdf
 poblaciones_16_20.dat.pdf
 poblaciones_15_20.dat.pdf
 poblaciones_14_20.dat.pdf
 poblaciones_13_20.dat.pdf
 poblaciones_12_20.dat.pdf
 poblaciones_11_20.dat.pdf
 poblaciones_10_20.dat.pdf
 poblaciones_9_20.dat.pdf
 poblaciones_8_20.dat.pdf
 poblaciones_7_20.dat.pdf
 poblaciones_6_20.dat.pdf
 poblaciones_5_20.dat.pdf
 poblaciones_4_20.dat.pdf
 poblaciones_3_20.dat.pdf
 poblaciones_2_20.dat.pdf
 poblaciones_1_20.dat.pdf
runge_kutta.x: runge_kutta.c
cc runge_kutta.c -lm -o runge_kutta.x
poblaciones_30_20.dat: runge_kutta.x
./runge_kutta.x 30 20
poblaciones_29_20.dat: runge_kutta.x
./runge_kutta.x 29 20
poblaciones_28_20.dat: runge_kutta.x
./runge_kutta.x 28 20
poblaciones_27_20.dat: runge_kutta.x
./runge_kutta.x 27 20
poblaciones_26_20.dat: runge_kutta.x
./runge_kutta.x 26 20
poblaciones_25_20.dat: runge_kutta.x
./runge_kutta.x 25 20
poblaciones_24_20.dat: runge_kutta.x
./runge_kutta.x 24 20
poblaciones_23_20.dat: runge_kutta.x
./runge_kutta.x 23 20
poblaciones_22_20.dat: runge_kutta.x
./runge_kutta.x 22 20
poblaciones_21_20.dat: runge_kutta.x
./runge_kutta.x 21 20
poblaciones_20_20.dat: runge_kutta.x
./runge_kutta.x 20 20
poblaciones_19_20.dat: runge_kutta.x
./runge_kutta.x 19 20
poblaciones_18_20.dat: runge_kutta.x
./runge_kutta.x 18 20
poblaciones_17_20.dat: runge_kutta.x
./runge_kutta.x 17 20
poblaciones_16_20.dat: runge_kutta.x
./runge_kutta.x 16 20
poblaciones_15_20.dat: runge_kutta.x
./runge_kutta.x 15 20
poblaciones_14_20.dat: runge_kutta.x
./runge_kutta.x 14 20
poblaciones_13_20.dat: runge_kutta.x
./runge_kutta.x 13 20
poblaciones_12_20.dat: runge_kutta.x
./runge_kutta.x 12 20
poblaciones_11_20.dat: runge_kutta.x
./runge_kutta.x 11 20
poblaciones_10_20.dat: runge_kutta.x
./runge_kutta.x 10 20
poblaciones_9_20.dat: runge_kutta.x
./runge_kutta.x 9 20
poblaciones_8_20.dat: runge_kutta.x
./runge_kutta.x 8 20
poblaciones_7_20.dat: runge_kutta.x
./runge_kutta.x 7 20
poblaciones_6_20.dat: runge_kutta.x
./runge_kutta.x 6 20
poblaciones_5_20.dat: runge_kutta.x
./runge_kutta.x 5 20
poblaciones_4_20.dat: runge_kutta.x
./runge_kutta.x 4 20
poblaciones_3_20.dat: runge_kutta.x
./runge_kutta.x 3 20
poblaciones_2_20.dat: runge_kutta.x
./runge_kutta.x 2 20
poblaciones_1_20.dat: runge_kutta.x
./runge_kutta.x 1 20
poblaciones_30_20.dat.pdf:  poblaciones_30_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_30_20.dat
				open poblaciones_30_20.pdf
poblaciones_29_20.dat.pdf:  poblaciones_29_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_29_20.dat
				open poblaciones_29_20.pdf
poblaciones_28_20.dat.pdf:  poblaciones_28_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_28_20.dat
				open poblaciones_28_20.pdf
poblaciones_27_20.dat.pdf:  poblaciones_27_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_27_20.dat
				open poblaciones_27_20.pdf
poblaciones_26_20.dat.pdf:  poblaciones_26_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_26_20.dat
				open poblaciones_26_20.pdf
poblaciones_25_20.dat.pdf:  poblaciones_25_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_25_20.dat
				open poblaciones_25_20.pdf
poblaciones_24_20.dat.pdf:  poblaciones_24_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_24_20.dat
				open poblaciones_24_20.pdf
poblaciones_23_20.dat.pdf:  poblaciones_23_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_23_20.dat
				open poblaciones_23_20.pdf
poblaciones_22_20.dat.pdf:  poblaciones_22_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_22_20.dat
				open poblaciones_22_20.pdf
poblaciones_21_20.dat.pdf:  poblaciones_21_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_21_20.dat
				open poblaciones_21_20.pdf
poblaciones_20_20.dat.pdf:  poblaciones_20_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_20_20.dat
				open poblaciones_20_20.pdf
poblaciones_19_20.dat.pdf:  poblaciones_19_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_19_20.dat
				open poblaciones_19_20.pdf
poblaciones_18_20.dat.pdf:  poblaciones_18_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_18_20.dat
				open poblaciones_18_20.pdf
poblaciones_17_20.dat.pdf:  poblaciones_17_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_17_20.dat
				open poblaciones_17_20.pdf
poblaciones_16_20.dat.pdf:  poblaciones_16_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_16_20.dat
				open poblaciones_16_20.pdf
poblaciones_15_20.dat.pdf:  poblaciones_15_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_15_20.dat
				open poblaciones_15_20.pdf
poblaciones_14_20.dat.pdf:  poblaciones_14_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_14_20.dat
				open poblaciones_14_20.pdf
poblaciones_13_20.dat.pdf:  poblaciones_13_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_13_20.dat
				open poblaciones_13_20.pdf
poblaciones_12_20.dat.pdf:  poblaciones_12_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_12_20.dat
				open poblaciones_12_20.pdf
poblaciones_11_20.dat.pdf:  poblaciones_11_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_11_20.dat
				open poblaciones_11_20.pdf
poblaciones_10_20.dat.pdf:  poblaciones_10_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_10_20.dat
				open poblaciones_10_20.pdf
poblaciones_9_20.dat.pdf:   poblaciones_9_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_9_20.dat
				open poblaciones_9_20.pdf
poblaciones_8_20.dat.pdf:   poblaciones_8_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_8_20.dat
				open poblaciones_8_20.pdf
poblaciones_7_20.dat.pdf: poblaciones_7_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_7_20.dat
				open poblaciones_7_20.pdf
poblaciones_6_20.dat.pdf: poblaciones_6_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_6_20.dat
				open poblaciones_6_20.pdf
poblaciones_5_20.dat.pdf:poblaciones_5_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_5_20.dat
				open poblaciones_5_20.pdf
poblaciones_4_20.dat.pdf:   poblaciones_4_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_4_20.dat
				open poblaciones_4_20.pdf
poblaciones_3_20.dat.pdf:   poblaciones_3_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_3_20.dat
				open poblaciones_3_20.pdf
poblaciones_2_20.dat.pdf:   poblaciones_2_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_2_20.dat
				open poblaciones_2_20.pdf
poblaciones_1_20.dat.pdf:   poblaciones_1_20.dat plot_poblaciones.py
				python plot_poblaciones.py poblaciones_1_20.dat
				open poblaciones_1_20.pdf
