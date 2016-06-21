from QCkit.thermalDesorption import trjProcessor

trjp = trjProcessor.TrjProcessor('./example_trj.trj')

trjp.plot_h2_metal_distance()
trjp.dump_data()