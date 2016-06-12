import molecule
import qcjob

# create the molecule from xyz file
h2 = molecule.from_xyz('geometries/hydrogen.xyz')

# generate Q-Chem job
job = qcjob.QCjob(h2)

job.run()
# h2.to_jmol()