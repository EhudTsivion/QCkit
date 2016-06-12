import os
import datetime
import random
import string
from molecule import Molecule
from input_template import template_text
import subprocess


class QCjob:
    def __init__(self,
                 molecule,
                 comment="",
                 exchange='hf',
                 basis='6-31g*',
                 job_type='sp',
                 molden_format='true',
                 threads=None,
                 job_name=None,
                 rems={}):

        if isinstance(molecule, Molecule):

            self.molecule = molecule

        else:
            raise TypeError('please enter an instance of a Molecule object')

        self.exchange = exchange

        self.basis = basis

        self.job_type = job_type

        self.molden_format = molden_format

        self.charge = molecule.charge

        self.multiplicity = molecule.multiplicity

        self.extra_rems = rems

        self.comment = comment

        # set number of openmp threads
        if not threads:

            try:
                self.threads = os.environ['OMP_NUM_THREADS']

            except KeyError:
                self.threads = 1

        # set automatic job name
        if not job_name:
            self.job_name = 'qc-' + datetime.datetime.now().strftime('%Y-%m-%d-') \
                            + str(random.randint(1000, 9999))

        self.input_file = self.job_name + '.in'

        self.output_file = self.job_name + '.qchem'

    def write_file(self):
        """
        write the input file to the filesystem
        so that it'll be available for the QChem executable

        :return: None
        """

        content = string.Template(template_text)

        # generate text for additional rem parameters
        rem_text = str()

        if self.extra_rems:
            for key, value in self.extra_rems.iteritems():
                rem_text += '{} {}'.format(key, value)

        # fill the template
        content = content.substitute(charge=self.charge,
                                     multiplicity=self.multiplicity,
                                     geometry=self.molecule.xyz(),
                                     job_type=self.job_type,
                                     exchange=self.exchange,
                                     molden_format=self.molden_format,
                                     other_rems=rem_text,
                                     basis=self.basis,
                                     comment=self.comment)

        with open(self.job_name + ".in", 'w') as f:

            f.write(content)

    def run(self):

        self.write_file()

        subprocess.call(['qchem', '-nt', self.threads, self.input_file])
        