import sys
import os
import xlrd
import re
import requests
import json
import pprint
import bs4
import time
import subprocess
from Bio import AlignIO
import jinja2
import urllib.request, urllib.error, urllib.parse

class homolog_finder():
    def __init__(self, xls_filename="",
                 urlUniprotId='http://www.uniprot.org/uniprot/?query=id',
                 urlCompara='https://rest.ensembl.org/homology/symbol/caenorhabditis_elegans/',
                 urlBrain='http://api.brain-map.org/api/v2/data/query.json?criteria=',
                 urlWormbase='http://www.wormbase.org/rest/widget',
                 urlOverlap='http://rest.ensembl.org/overlap',
                 urlUniprot='http://www.uniprot.org/uniprot',
                 urlDisgenet='http://www.disgenet.org/oql',
                 main_folder='',
                 wormbase_name_column=4,
                 wormbase_id_column=2,
                 uniprot_id_column=3,
                 sprot_database='/media/ashafix/WesternDigital/data/uniprot_sprot.fasta',
                 jackhmmer='/media/ashafix/WesternDigital/tools/hmmer-3.1b2/src/jackhmmer'):
        self.xls_filename = xls_filename
        self.urls = dict()
        self.urls['xrefs'] = 'http://rest.ensembl.org/xrefs/id'
        self.urls['id'] = urlUniprotId
        self.urls['compara'] = urlCompara
        self.urls['brain'] = urlBrain
        self.urls['wormbase'] = urlWormbase
        self.urls['overlap'] = urlOverlap
        self.urls['uniprot'] = urlUniprot
        self.urls['disgenet'] = urlDisgenet
        self.required_folders = ['brain', 'compara', 'wormbase_protein', 'aliases', 'disgenet', 'uniprot']
        self.main_folder = main_folder
        self.disgene = dict() #stores all the gene to disease assocations
        self.disgene_raw = dict() #stores all raw data from disgenet
        self.disgene = dict() #stores parsed data from disgenet
        self.gene_ids = list()
        self.uniprot_ids = list()
        self.initialize_folders()
        self.primary_genes = list()
        self.homolog_gene = list()
        self.homologs = dict() #key: Uniprot ID or gene name, value: a list of homolog object
        self.sprot_database = sprot_database #location of the SwissProt database for HMMER
        self.jackhmmer = jackhmmer #location of the jackhmmer executable
        self.uniprot_genename = dict() #key: UniProt ID, value: Gene name
        self.expressed_in_brain = dict() #key: gene name, value: boolean (expressed in brain or not)
        self.wormbase_names = list() #a list of wormbase names
        if xls_filename:
            #if a xls filename is specified, read all the data
            #initalize primary genes for all entries and find possible homologs
            self.workbook = xlrd.open_workbook(self.xls_filename)
            #if not self.validate_inputfile():
            #    sys.exit(1)
            self.get_genes_from_xls_file(wormbase_id_column)
            self.gene_ids = self.get_ids_from_xls_file(wormbase_id_column)
            self.uniprot_ids = self.get_ids_from_xls_file(uniprot_id_column)
            self.uniprot_ids = self.uniprot_ids
            self.gene_ids = self.gene_ids
            self.wormbase_names = self.get_ids_from_xls_file(wormbase_name_column)
            if len(self.gene_ids) != len(self.uniprot_ids):
                print('Gene length and Uniprot length do not match: {} vs {}'.format(len(self.gene_ids), len(self.uniprot_ids)))
                sys.exit(1)
            self.find_homologs()
        else:
            self.workbook = None

    def find_homologs(self):
        for ii, gene in enumerate(self.gene_ids):
            self.primary_genes.append(primary_gene(uniprot_id=self.uniprot_ids[ii],
                                                   wormbase_id=gene,
                                                   homolog_finder=self,
                                                   wormbase_name=self.wormbase_names[ii]))
            self.primary_genes[-1].find_homologs()
        for primary in self.primary_genes:
            for uniprot in primary.uniprot_ids_from_hmmer:
                self.create_homolog(primary,
                                    uniprot=uniprot,
                                    source='hmmer')
            for compara in primary.compara:
                self.create_homolog(primary,
                                    gene_name=compara,
                                    source='compara')
            self.create_homolog_html(primary)

    def get_homolog(self, uniprot_id):
        return self.homologs.get(uniprot_id)

    def create_homolog(self,
                       primary,
                       gene_name='',
                       uniprot='',
                       source='',
                       double_check=True):
        """
        
        :param primary: 
        :param double_check: 
        :return: 
        """
        if double_check and self.get_homolog(primary.wormbase_name):
            if uniprot and uniprot in self.get_homolog(primary.wormbase_name):
                return None
            if gene_name and gene_name in self.get_homolog(primary.wormbase_name):
                return None

        new_homolog = homolog_gene(primary,
                                   uniprot_id=uniprot,
                                   gene_name=gene_name,
                                   source=source,
                                   homolog_finder=self
                                   )
        if not self.homologs.get(primary.wormbase_name):
            self.homologs[primary.wormbase_name] = list()
        self.homologs[primary.wormbase_name].append(new_homolog)

    def initialize_folders(self):
        """
        creates all input folders needed to run the script
        :return: None
        """

        for folder in self.required_folders:
            joinedFolder = os.path.join(folder, self.main_folder)
            if not os.path.isdir(joinedFolder):
                os.makedirs(joinedFolder)

    def validate_inputfile(self):
        """
        validates the input file
        :return:
        """
        # let's check if the gene matches the ID, UniProt and name
        xl_sheet = self.workbook.sheet_by_index(0)
        wb_genes = xl_sheet.col(0)
        print('validating genes')
        for i, wb_gene in enumerate(wb_genes):
            error = (None, '')
            print('Validating {0}'.format(wb_gene.value))
            r = requests.get(
                '{0}/{1}?content-type=application/json'.format(self.urls['xrefs'], wb_gene.value.strip(), timeout=10))
            if r.status_code == 200:
                data = json.loads(r.text)
            else:
                print('error retrieving data for {0}'.format(wb_gene))
                data = list()
            # UniProt check
            for d in data:
                if 'UniProt' in d.get('db_display_name'):
                    if not xl_sheet.cell_value(rowx=i, colx=3) == d.get('primary_id'):
                        if error[0] == None:
                            error = (True, d.get('primary_id'))
                    else:
                        error = (False, '')
            if error[0] == True:
                print('Error in row {0}: Given Uniprot ID {1} does not match retrieved Uniprot ID {2}'.format(i + 1,
                                                                                                              xl_sheet.cell_value(
                                                                                                                  rowx=i,
                                                                                                                  colx=3),
                                                                                                              error[1]))

            if error[0] == None:
                print('Warning in row {0}: No Uniprot ID found'.format(i + 1))
            # ID check
            error = (None, '')
            for d in data:
                if 'wormbase_gseqname' in d.get('dbname'):
                    if not xl_sheet.cell_value(rowx=i, colx=2) == d.get('display_id'):
                        if error[0] == None:
                            error = (True, d.get('display_id'))
                    else:
                        error = (False, '')
            if error[0] == True:
                print('Error in row {0}: Given display ID {1} does not match retrieved display ID {2}'.format(i + 1,
                                                                                                              xl_sheet.cell_value(
                                                                                                                  rowx=i,
                                                                                                                  colx=2),
                                                                                                              error[1]))

            if error[0] == None:
                print('Warning in row {0}: No Uniprot ID found'.format(i + 1))


            # let's check if the gene matches the ID, UniProt and name

    def get_genes_from_xls_file(self, col):
        """
        
        :param filename: 
        :param col: 
        :return: 
        """

        xl_sheet = self.workbook.sheet_by_index(0)
        cols = xl_sheet.col(col)

        p = re.compile('^.*\(([-.a-zA-Z0-9]+)\)')

        for c in cols:
            m = p.match(c.value)
            if m:
                self.gene_ids.append(m.group(1))
            else:
                self.gene_ids.append('')
        return self.gene_ids

    def get_ids_from_xls_file(self, col):
        """
        
        :param filename: 
        :param col: 
        :return: 
        """
        xl_sheet = self.workbook.sheet_by_index(0)
        cols = xl_sheet.col(col)

        ids = [c.value for c in cols]
        return ids

    def render(self, template_path, context):
        """
    
        :param template_path: location of the template
        :param context: a dictionary which is used for the template
        :return: the parsed template as a string
        """
        path, filename = os.path.split(template_path)
        return jinja2.Environment(
            loader=jinja2.FileSystemLoader(path or './')
        ).get_template(filename).render(context)

    def create_homolog_html(self, primary):
        """
        
        :param primary: 
        :return: 
        """
        input_dict = dict()
        input_dict['primary'] = primary.wormbase_name
        input_dict['homologs'] = self.homologs[primary.wormbase_name]
        table = list()
        table.append(list())
        table[0].append('Uniprot ID')
        table[0].append('Source')
        table[0].append('Expressed in brain')
        table[0].append('Alignment')
        filename = 'jackhmmer_{}.aln'.format(primary.uniprot_id)
        for homolog in self.homologs[primary.wormbase_name]:
            table.append(list())
            table[-1].append('<a href="http://www.uniprot.org/uniprot/{0}">{0}</a>'.format(homolog.uniprot_id))
            table[-1].append(homolog.source)
            if homolog.expressed_in_brain:
                table[-1].append('yes')
            elif homolog.expressed_in_brain == False:
                table[-1].append('no')
            else:
                table[-1].append('unknown')
            table[-1].append('<a href="{}.html" target="_blank">Open alignment</a>'.format(filename))
        result = self.render('table.template.html', {'table': table})
        print(primary.wormbase_name)
        with open('{}.html'.format(primary.wormbase_name.replace('/', '_')), 'w') as f:
            f.write(result)
        self.alignment(filename)


    def alignment(self, filename='align.txt'):
        """
        Creates an alignment file from an alignment, only human sequences are used
        :param filename: the file where the alignment is aved
        :return: None
        """
        alignment = AlignIO.read(filename, "stockholm")
        items = list()
        for i, a in enumerate(alignment):
            if 'HUMAN' in a.id or i == 0:
                items.append({'id': a.id, 'sequence': str(a.seq)})

        result = self.render('alignment.template', {'items': items})

        with open('{}.html'.format(filename), 'w') as f:
            f.write(result)


    def hmmer_to_table(self, filename='h.txt'):

        #the location of the individual rows
        ranges = [(0, 22), (22, 33), (33, 56), (56, 66), (67, 76), (77, 83), (84, 89), (90, 99), (100, 106), (107, 112), (114, 118), (118, 122), (122, 126), (126, 130), (130, 134), (134, 138), (138, 142), (142, 146), (146, 10**3)]
        lines = list()
        header = 0
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    header += 1
                if header == 2 or not line.startswith('#'):
                    lines.append(list())
                    for r in ranges:
                        lines[-1].append(line[r[0]:r[1]].strip())
                if header > 3:
                    break
        lines[0][1] = 'Expressed in the brain'
        for i, line in enumerate(lines[1:]):
            gene = self.uniprot_to_genename(line[0].split('|')[1])
            lines[i + 1][1] = expressed_in_brain(gene)

        result = self.render('table.template.html', {'table': lines})
        with open('table.html', 'w') as f:
            f.write(result)


    def uniprot_to_genename(self, uniprot):
        """
        gets the gene name for a UniProt identifier
        :param uniprot: a UniProt identifier
        :return: a gene name, in case of multiple gene names only the first is returned
        """

        if not self.uniprot_genename.get(uniprot):
            self.uniprot_genename[uniprot] = None
            filename = 'uniprot/id_{}.txt'.format(uniprot)
            if not os.path.isfile(filename):
                r = requests.get('{0}:{1}&columns=id,entry%20name,genes&format=json'.format(self.urls['id'],
                                                                                        uniprot))


                if r.status_code == 200:
                    with open(filename, 'w') as f:
                        f.write(r.text)
                    self.uniprot_genename[uniprot] = json.loads(r.text)[0]['genes'].split(' ')[0]
            else:
                with open(filename, 'r') as f:
                    self.uniprot_genename[uniprot] =json.load(f)[0]['genes'].split(' ')[0]

        return self.uniprot_genename[uniprot].split(';')[0]


class primary_gene():
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
        self.homolog_finder = homolog_finder #the outer homolog_finder object which stores all the settings
        self.homologs = list() #stores a list of homolog objects
        self.compara = list() #stores the info from compara database
        self.overwrite = False #whether to overwrite existing data or not
        self.jackhmmer_filename = '' #the filename where the jackhmmer output is located
        self.wormbase_homologs = ''
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
        :param hmmer_file: the filename which is searched
        :return: a list of IDs in order of appearance
        """

        ids = list()
        if not os.path.isfile(self.jackhmmer_filename):
            self.uniprot_ids_from_hmmer = list()
            return False

        with open(self.jackhmmer_filename, 'r') as f:
            lines = f.read().split('\n')
        if not lines[0].startswith('# STOCKHOLM'):
            self.uniprot_ids_from_hmmer = list()
            return False

        orig_id = lines[1].split('|')[1]
        for line in lines[3:]:
            if line.startswith('#'):
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
            print(filename)
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
            soup = bs4.BeautifulSoup(data)
            links = soup.find_all('a', {'class': 'protein-link'})
            if links:
                r = requests.get(
                    '{0}/protein/{1}/overview'.format(self.homolog_finder.urls['wormbase'], links[0].get('href').split('/')[-1]))
                data = r.text
            else:
                data = ''
            with open(filename, 'w') as f:
                f.write(data)
        else:
            with open(filename, 'r') as f:
                data = f.read()
        soup = bs4.BeautifulSoup(data)
        links = soup.find_all('a', {'class': 'protein-link'})
        if links:
            self.wormbase_homologs = [links[0].get('href').split(':')[-1], links[0].text]
            return True
        else:
            self.wormbase_homologs = ['', '']
            return False


class homolog_gene():
    """
    stores all the information about homologs
    """
    def __init__(self,
                 primary, #the orgin of the homology, a primary gene object
                 uniprot_id='', #the uniprot ID, used as primary identifier
                 gene_name='', #the gene name, used as primary identifier
                 jackhmmer_files='', #the filenames with the jackhmmer alignment
                 jackhmmer_ids=list(), #a list of uniprot IDs identified by jackhmmer
                 source='', #from where the homology was deduced, e.g. hmmer or compara
                 homolog_finder=None
                 ):
        self.primary = primary
        self.uniprot_id = uniprot_id

        self.source = source
        if type(jackhmmer_files) == list:
            self.jackhmmer_files = jackhmmer_files
        else:
            self.jackhmmer_files = list()
            self.jackhmmer_files.append(jackhmmer_files)

        self.jackhmmer_ids = jackhmmer_ids
        if not gene_name and homolog_finder and uniprot_id:
            self.gene_name = homolog_finder.uniprot_to_genename(uniprot_id)
        else:
            self.gene_name = gene_name.upper()
        if self.gene_name:
            self.expressed_in_brain = self.expressed_in_brain(self.gene_name)
        else:
            self.expressed_in_brain = None

    def add_info(self,
                 jackhmmer_file=None,
                 jackhmmer_ids=None,
                 compara=None,
                 primary=None):

        if primary:
            jackhmmer_ids = primary.uniprot_ids_from_hmmer
            jackhmmer_file = primary.jackhmmer_file
            compara = primary.compara
        if jackhmmer_file:
            if jackhmmer_file not in self.jackhmmer_files:
                self.jackhmmer_files.append(jackhmmer_file)
        if jackhmmer_ids:
            if type(jackhmmer_ids) != list:
                jackhmmer_ids = [jackhmmer_ids]
            for jackhmmer_id in jackhmmer_ids:
                if jackhmmer_id not in self.jackhmmer_ids:
                    self.jackhmmer_ids.append(jackhmmer_id)
        if compara:
            if type(compara) != list:
                compara = list(compara)
            for comp in compara:
                if comp not in self.compara:
                    self.compara.append(comp)

    def get_diseases_from_disgenet(self, gene, overwrite=False):
        """

        inspired from here:
        http://www.disgenet.org/ds/DisGeNET/scripts/disgenet_python3.py
        :param gene:
        :return:
        """

        filename = 'disgenet/{0}.txt'.format(gene)
        if not overwrite:
            if os.path.isfile(filename):
                with open(filename, 'r') as f:
                    return (f.read())

        query = """
                DEFINE
                    c0='/data/gene_disease_summary',
            c1='/data/diseases',
            c2='/data/genes',
            c3='/data/gene_roles',
            c4='/data/sources'
                ON
                   'http://www.disgenet.org/web/DisGeNET'
                SELECT
                    c2 (geneId, name, uniprotId, description, pathName, pantherName),
                    c1 (cui, name, diseaseClassName, STY, MESH, omimInt),
                   c3 (PI, PL),
                   c0 (score, pmids,  snps)

                FROM
                    c0
                WHERE
                    (
                        c2.name = '{0}'
                    AND
                        c4 = 'ALL'
                    )
                ORDER BY
                    c0.score DESC""".format(gene)
        query = query.encode('utf-8')
        req = urllib.request.Request(self.urls['disgenet'])
        res = urllib.request.urlopen(req, query)
        data = res.read().decode("utf-8")
        with open(filename, 'w') as f:
            f.write(data)
        self.parse_disgenet(data)
        return (data)

    def parse_disgenet(self, data):
        """

        :param gene:
        :return:
        """
        if self.disgenet.get(gene):
            return self.disgenet['gene']

        filename = 'disgenet/{}.txt'.format(gene)
        if not os.path.isfile(filename):
            data = self.get_diseases_from_disgenet(gene)
        else:
            with open(filename, 'r') as f:
                data = f.read()

        return data

    def expressed_in_brain(self,
                           gene,
                           donor='15496',
                           thresh_avg=3,
                           thresh_max=5):
        """Checks if a gene is expressed in the brain via the Allen Brain Atlas
        :param gene: a gene name according to HGNC
        :param donor: the id of the brain donor
        :param thresh_avg: threshold for average in all measurements
        :param thresh_max: threshold which needs to met at least once
        :return: True, if above thresholds, False, if below thresholds
        """
        gene = gene.upper()
        exp_max = 0
        exp_avg = 0
        probes = list()

        filename = 'brain/{0}_probes.txt'.format(gene.replace('/', '_'))
        if os.path.isfile(filename):
            with open(filename, 'r') as f:
                msg = json.load(f)
        else:
            r = requests.get(
                "http://api.brain-map.org/api/v2/data/query.json?criteria=model::Probe,rma::criteria,[probe_type$eq'DNA'],products[abbreviation$eq'HumanMA'],gene[acronym$eq'{0}'],rma::options[only$eq'probes.id']".format(
                    gene),
                timeout=10)

            if r.status_code != 200:
                return None
            msg = json.loads(r.text)['msg']
            with open(filename, 'w') as f:
                json.dump(msg, f)
        for id in msg:
            probes.append(str(id['id']))
        filename = 'brain/{0}_expression_donor_{1}.txt'.format(gene.replace('/', '_'), donor)
        if os.path.isfile(filename):
            with open(filename, 'r') as f:
                try:
                    probes = json.load(f)
                except:
                    return None
        else:
            r = requests.get(
                "http://api.brain-map.org/api/v2/data/query.json?criteria=service::human_microarray_expression[probes$eq{0}][donors$eq{1}]".format(
                    ','.join(probes), donor), timeout=10)
            if r.status_code == 200:
                print(r.text)
                if json.loads(r.text).get('success'):
                    probes = json.loads(r.text)['msg']['probes']
                    with open(filename, 'w') as f:
                        json.dump(probes, f)
                else:
                    with open(filename, 'w') as f:
                        f.write('')
                    return None
            else:
                return None
        for probe in probes:
            exp_level = [float(p) for p in probe['expression_level']]
            exp_max += max(exp_level)
            exp_avg += sum(exp_level) / len(exp_level)

        exp_max /= len(probes)
        exp_avg /= len(probes)

        if exp_avg > thresh_avg and exp_max > thresh_max:
            return True
        else:
            return False




new_homolog_finder = homolog_finder(xls_filename='57_top_candidates.xlsx')
for k in new_homolog_finder.homologs.keys():
    pass



