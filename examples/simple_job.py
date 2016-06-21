import molecule
import qcjob

# create the molecule from xyz file
h2 = molecule.from_xyz('geometries/hydrogen.xyz')

# generate Q-Chem job
job = qcjob.QCjob(h2, rems={"gen_scfman": "true"})

job.run()

# if you wish to test the JMol viewer, uncomment to following line:
# h2.to_jmol()
