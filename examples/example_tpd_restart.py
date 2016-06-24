from QCkit.thermalDesorption.restartJob import TPDrestart

TPDrestart(pos_vel_file='./example.last_pos_vel',
           basis="6-31g*",
           exchange="b97-d3",
           threads='4',
           high_temperature=124,
           temp_advance=2,  # advance temperature by 2 K
           thermostat="langevin",
           thermostat_timescale=120,  # fs
           aimd_step_duration=60,  # atomic units
           aimd_num_steps=10)