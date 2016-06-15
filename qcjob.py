import os
import datetime
import random
import string
from molecule import Molecule
from input_template import template_text
import subprocess
import logging


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
                 rems={},
                 other_things={}):

        self.molecule = molecule

        self.exchange = exchange

        self.basis = basis

        self.job_type = job_type

        self.molden_format = molden_format

        self.charge = molecule.charge

        self.multiplicity = molecule.multiplicity

        self.rems = rems

        self.comment = comment

        self.job_done = False

        self.other_things = other_things

        # set number of openmp threads
        if not threads:

            try:
                self.threads = os.environ['OMP_NUM_THREADS']

            except KeyError:
                self.threads = str(1)

        else:
            self.threads = threads

        # set automatic job name
        if not job_name:
            self.job_name = 'qc-' + datetime.datetime.now().strftime('%Y-%m-%d-') \
                            + str(random.randint(10000, 99999))

        else:
            self.job_name = job_name

        self.input_file = self.job_name + '.in'

        self.output_file = self.job_name + '.qchem'

        # logging
        logging.getLogger('qcrun')

        logging.info('Initialized new Q-Chem job')

        # write the file
        self.write_file()

    def write_file(self):
        """
        write the input file to the filesystem
        so that it'll be available for the QChem executable

        :return: None
        """

        content = string.Template(template_text)

        # generate text for additional rem parameters
        rem_text = str()

        if self.rems:
            for key, value in self.rems.items():
                rem_text += '{} {}\n'.format(key, value)

        things_text = " "

        if self.other_things:
            for key, value in self.other_things.items():
                things_text += '${}\n{}\n$end\n'.format(key, value)


        # fill the template
        content = content.substitute(charge=self.charge,
                                     multiplicity=self.multiplicity,
                                     geometry=self.molecule.xyz(),
                                     job_type=self.job_type,
                                     exchange=self.exchange,
                                     molden_format=self.molden_format,
                                     other_rems=rem_text,
                                     basis=self.basis,
                                     comment=self.comment,
                                     other_things=things_text)

        with open(self.input_file, 'w') as f:

            f.write(content)

    def add_rem(self, key, value):
        """
        remove rem value from extera rems

        :param value:
        :return:
        """
        self.rems[key] = value

        self.write_file()

    def rm_rem(self, key):
        """
        add a rem value to input file

        :return:
        """
        self.rems.pop(key, None)

    def run(self, save=None):

        self.write_file()

        if save == 0:

            subprocess.call(['qchem', '-nt', self.threads, self.input_file, self.output_file])

        # we assume that for aimd runs the default is to save
        elif save == 1 or self.job_type.lower() == "aimd":

            subprocess.call(['qchem', '-save', '-nt',
                             self.threads, self.input_file, self.output_file, self.job_name])

        self.job_done = True
