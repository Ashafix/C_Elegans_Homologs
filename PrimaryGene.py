import sys
import os
import json
import subprocess
from bs4 import BeautifulSoup


class PrimaryGene():
    def __init__(self,
                 uniprot_id='',
                 wormbase_id='',
                 protein_name='',
                 wormbase_name='',
                 homolog_finder=None,
                 overwrite=False):
        self.uniprot = uniprot_id
        self.wormbase_id = wormbase_id
        self.wormbase_name = wormbase_name
        self.protein_name = protein_name
        self.uniprot_id = uniprot_id
        self.homolog_finder = homolog_finder  # the outer homolog_finder object which stores all the settings
        self.homologs = list()  # stores a list of homolog objects
        self.compara = list()  # stores the info from compara database
        self.overwrite = False  # whether to overwrite existing data or not
        self.jackhmmer_filename = ''  # the filename where the jackhmmer output is located
        self.wormbase_homologs = ''
        self.overwrite = overwrite
        self.uniprot_ids_from_hmmer = list()
        if homolog_finder and uniprot_id:
            self.gene_id = homolog_finder.uniprot_to_genename(uniprot_id)
        else:
            self.gene_id = ''

    def find_homologs(self):
        self.get_homologs_from_compara()
        self.run_jackhmmer()
        self.get_uniprot_ids_from_hmmer()
        self.get_info_from_wormbase()

    def run_jackhmmer(self):
        """

        :return:
        """
        filename = '{0}.fasta'.format(self.uniprot_id)
        jackhmmer_output = 'jackhmmer_{0}.aln'.format(self.uniprot_id)

        if not os.path.isfile(filename):
            r = requests.get('{0}/{1}.fasta'.format(self.homolog_finder.urls['uniprot'], self.uniprot_id))
            with open(filename, 'w') as f:
                f.write(r.text)
        if not os.path.isfile(jackhmmer_output) or self.overwrite:
            args = [self.homolog_finder.jackhmmer,
                    '-A',
                    jackhmmer_output,
                    filename,
                    self.homolog_finder.sprot_database]
            subprocess.check_output(args)
        self.jackhmmer_filename = jackhmmer_output

    def get_uniprot_ids_from_hmmer(self):
        """
        gets all UniProt IDs from a hmmer alignment file
        does not return the original input ID
        :return: a list of IDs in order of appearance
        """

        ids = list()
        print(self.jackhmmer_filename)
        if not os.path.isfile(self.jackhmmer_filename):
            self.uniprot_ids_from_hmmer = list()
            return False

        with open(self.jackhmmer_filename, 'r') as f:
            lines = f.read().splitlines()
        print(self.jackhmmer_filename)
        print(lines)
        if not lines[0].startswith('# STOCKHOLM'):
            self.uniprot_ids_from_hmmer = list()
            return False

        for line in lines:
            if line.startswith('#=GR'):
                id = line.split('|')[1]
                if id not in ids:
                    ids.append(id)
            else:
                break
        self.uniprot_ids_from_hmmer = ids
        return True

    def get_homologs_from_compara(self,
                                  species='homo_sapiens'):
        filename = 'compara/{0}_{1}.txt'.format(species, self.wormbase_id)
        if not os.path.isfile(filename):
            r = requests.get(
                '{0}{1}?content-type=application/json&target_species={2}'.format(self.homolog_finder.urls['compara'],
                                                                                 self.wormbase_id,
                                                                                 species))
            data = r.text
            with open(filename, 'w') as f:
                f.write(data)
        else:
            # print(filename)
            with open(filename, 'r') as f:
                data = f.read()

        j = json.loads(data)
        if not 'error' in j.keys() and len(j['data']) > 0:
            for i in range(len(j['data'][0]['homologies'])):
                self.compara.append(j['data'][0]['homologies'][i]['target']['id'])

    def get_info_from_wormbase(self):
        filename = 'wormbase_protein/{0}.html'.format(self.wormbase_id)

        if not os.path.isfile(filename):
            r = requests.get('{0}/cds/{1}/overview'.format(self.homolog_finder.urls['wormbase'], self.wormbase_id))
            if r.status_code != 200:
                r = requests.get('http://www.wormbase.org/rest/widget/cds/{0}a/overview'.format(self.wormbase_id))
            data = r.text.encode('ascii', errors='ignore')
            with open(filename, 'wb') as f:
                f.write(data)
        else:
            with open(filename, 'rb') as f:
                data = f.read()

        filename = 'wormbase_protein/{0}_protein.html'.format(self.wormbase_id)
        if not os.path.isfile(filename):
            soup = BeautifulSoup(data)
            links = soup.find_all('a', {'class': 'protein-link'})
            if links:
                r = requests.get(
                    '{0}/protein/{1}/overview'.format(self.homolog_finder.urls['wormbase'],
                                                      links[0].get('href').split('/')[-1]))
                data = r.text
            else:
                data = ''
            with open(filename, 'w') as f:
                f.write(data)
        else:
            with open(filename, 'r') as f:
                data = f.read()
        soup = BeautifulSoup(data)
        links = soup.find_all('a', {'class': 'protein-link'})
        if links:
            self.wormbase_homologs = [links[0].get('href').split(':')[-1], links[0].text]
            return True
        else:
            self.wormbase_homologs = ['', '']
            return False