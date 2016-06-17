from  QCkit import molecule
from QCkit.thermalDesorption import simulatedTDP

tiocatMg = molecule.from_xyz('./geometries/hydrogen_dimer.xyz')

threads_per_process = 1

simulatedTDP.TPD(molecule=tiocatMg,
                 basis="6-31g*",
                 exchange="b97-d3",
                 threads=str(threads_per_process),
                 low_T=10,
                 high_T=12,
                 temp_advance=1,    # advance temperature by 1 K
                 thermostat="langevin",
                 thermostat_timescale=120,  # fs
                 time_step=60,      # atomic units
                 aimd_steps=2)     # per step
