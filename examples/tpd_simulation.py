from  QCkit import molecule
from QCkit.thermalDesorption import simulatedTDP

tiocatMg = molecule.from_xyz('./geometries/hydrogen_dimer.xyz')

threads_per_process = 1

simulatedTDP.TPD(molecule=tiocatMg,
                 basis="6-31g*",
                 exchange="b97-d3",
                 threads=str(threads_per_process),
                 low_temperature=10,
                 high_temperature=50,
                 tpd_job_name='test1',
                 temp_advance=1,  # advance temperature by 1 K
                 thermostat="langevin",
                 thermostat_timescale=120,  # fs
                 aimd_step_duration=60,  # atomic units
                 aimd_num_steps=3)  # per step
