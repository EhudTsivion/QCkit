import molecule
from thermalDesorption import simulatedTDP

tiocatMg = molecule.from_xyz('./geometries/hydrogen_dimer.xyz')

threads_per_process = 1

simulatedTDP.TPD(molecule=tiocatMg,
                 basis="6-31g*",
                 exchange="b97-d3",
                 threads=str(threads_per_process),
                 high_T=155,
                 low_T=150,
                 temp_advance=1,  # advance temperature by 1 K
                 thermostat="nose_hoover",
                 # for simplicity, only N_H is supported now
                 nh_length=3,
                 nh_timescale=45,
                 time_step=10,
                 aimd_steps=20)
