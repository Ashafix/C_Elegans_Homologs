import sys
import os
import xlrd
import re
import requests
import json
import pprint
import time
import jinja2
import urllib.request, urllib.error, urllib.parse
import xml.etree.ElementTree as ET
import pronto
import argparse
from HomologGene import HomologGene
from PrimaryGene import PrimaryGene


class HomologFinder():
    def __init__(self, xls_filename='',
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
                 sprot_database='uniprot_sprot.fasta',
                 jackhmmer='jackhmmer',
                 disease_ontology_filename='HumanDO.obo'):
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
        self.required_folders = ['brain', 'compara', 'wormbase_protein', 'aliases', 'disgenet',
                                 'uniprot', 'jackhmmer', 'fasta']
        self.main_folder = main_folder
        self.disgenet = dict() #stores all the gene to disease assocations
        self.disgenet_raw = dict() #stores all raw data from disgenet
        self.disgenet = dict() #stores parsed data from disgenet
        self.gene_ids = list() #stores a list of gene IDs for which homologs are searched
        self.uniprot_ids = list() #tores a list of UniProt IDs for which homologs are searched
        self.initialize_folders()
        self.primary_genes = list()
        self.homolog_gene = list()
        self.homologs = dict() #key: Uniprot ID or gene name, value: a list of homolog object
        self.sprot_database = sprot_database #location of the SwissProt database for HMMER
        self.jackhmmer = jackhmmer #location of the jackhmmer executable
        self.uniprot_genename = dict() #key: UniProt ID, value: Gene name
        self.expressed_in_brain = dict() #key: gene name, value: boolean (expressed in brain or not)
        self.wormbase_names = list() #a list of wormbase names
        self.neuro_diseases = list() #a list of all diseases associated with brain or nerves
        self.disease_synonyms = dict() #a dict with synonyms for neuro_diseases, key: a disease ontology entry, value: the corresponding MESH value
        self.disease_ontology_filename = disease_ontology_filename #the filename where the human disease ontology is stored
        self.uniprot_info = dict() #stores info about Uniprot entries which were retrieved, key: Uniprot ID, value: a dict with information
        self.disgenet_map = None #a placeholder for a dict, key: UniProtID, value: GeneID
        self.disgenet_map_filename ='mapa_geneid_4_uniprot_crossref.tsv'
        self.diseases = None #a placeholder for a list with all diseases
        self.disease_stats = None #a placeholder for a dict with all disease statistics
        self.alignments = dict() #a dict storing all the alignment objects, needed for performance

        
        if xls_filename:
            #if a xls filename is specified, read all the data
            #initalize primary genes for all entries and find possible homologs
            self.workbook = xlrd.open_workbook(self.xls_filename)
            #if not self.validate_inputfile():
            #    sys.exit(1)
            self.get_genes_from_xls_file(wormbase_id_column)
            self.gene_ids = self.get_ids_from_xls_file(wormbase_id_column)
            self.uniprot_ids = self.get_ids_from_xls_file(uniprot_id_column)
            self.gene_ids = self.gene_ids
            self.wormbase_names = self.get_ids_from_xls_file(wormbase_name_column)
            if len(self.gene_ids) != len(self.uniprot_ids):
                print('Gene length and Uniprot length do not match: {0} vs {1}'.format(len(self.gene_ids), len(self.uniprot_ids)))
                sys.exit(1)
            self.find_homologs()
        else:
            self.workbook = None


            #parse and store the disease ontology
        if self.disease_ontology_filename:
            self.onto = pronto.Ontology(self.disease_ontology_filename)
            # self.disease_parents =
            # parents = [self.onto['DOID:150'], self.onto['DOID:863']] #mental and nerve diseases
            self.disease_parents = [self.onto['DOID:9120'], self.onto['DOID:1289']]  # amyloidosis and neurodegenerative diseases
            #self.disease_parents = [self.onto['DOID:9120']]  # amyloidosis
            # self.disease_parents = [self.onto['DOID:1289']]  #neurodegenerative diseases

            for parent in self.disease_parents:
                for c, child in enumerate(parent.children):
                    self.neuro_diseases.append(child)
                    while child.children:
                        child = child.children
                        self.neuro_diseases.extend(child)

            for neuro_disease in self.neuro_diseases:
                if neuro_disease.other.get('xref'):
                    for term in neuro_disease.other['xref']:
                        if term.startswith('MESH:'):
                            self.disease_synonyms[term] = neuro_disease
                        elif term.startswith('OMIM:'):
                            self.disease_synonyms[term] = neuro_disease

    def calculate_scores(self):
        """
        
        :return: 
        """


        for primary_gene in self.homologs:
            for homolog in self.homologs[primary_gene]:
                homolog.calculate_score(homolog_finder=self)

        return None

    def get_all_diseases(self):
        """

        :return:
        """
        if self.diseases:
            return self.diseases
        self.calculate_scores()
        self.get_human_homologs(verbose=False)
        diseases = list()
        for primary_gene in self.human_homologs:
            for homolog in self.human_homologs[primary_gene]:
                if homolog.disgenet.get('diseases'):
                    for disease in homolog.disgenet['diseases']:
                        diseases.append(disease['disease'])
        self.diseases = list(set(diseases))
        return self.diseases

    def get_disease_statistics(self):
        """

        :return:
        """
        if self.disease_stats:
            return self.disease_stats
        self.get_human_homologs(verbose=False)
        disease_stats = dict()
        for primary_gene in self.human_homologs:
            disease_stats[primary_gene] = 0
            for homolog in self.human_homologs[primary_gene]:
                if homolog.disgenet and homolog.disgenet.get('diseases'):
                    disease_stats[primary_gene] = len(homolog.disgenet['diseases'])
        self.disease_stats = disease_stats

    def get_homologs(self):
        """
        
        :return: 
        """
        return [homolog for primary_gene in self.homologs.values() for homolog in primary_gene]

    def get_homologs_expressed_in_brain(self):
        """

        :return:
        """
        expressed = list()
        for homolog in self.get_homologs():
            homolog.calculate_score()
            self.expressed_in_brain[homolog.name] = homolog.expressed_in_brain
            if homolog.expressed_in_brain:
                expressed.append(homolog)

        return expressed

    def print_homologs(self):
        """
        prints information about homologs
        :return: None
        """

        for primary_gene in self.homologs:
            print('Primary gene: {}'.format(primary_gene))
            for homolog in self.homologs[primary_gene]:
                print(homolog)
                print('=' * 20)

    def get_human_homologs(self, verbose=True):
        """
        returns all human homologs
        :return: 
        """

        self.human_homologs = dict()
        for primary_gene in self.homologs:
            if verbose:
                print('#' * 40)
                print('Primary gene: {}'.format(primary_gene))
            self.human_homologs[primary_gene] = list()
            for homolog in self.homologs[primary_gene]:
                if 'homo' in homolog.organism.lower() or 'human' in homolog.organism.lower():
                    self.human_homologs[primary_gene].append(homolog)
                    if verbose:
                        print(homolog)
                        print('=' * 20)

        counter = dict()
        counter['no homolog'] = 0
        counter['potential homolog'] = 0
        counter['total'] = 0
        for primary_gene in self.human_homologs:
            if verbose:
                print('Primary gene {} has {} potential homologs'.format(primary_gene,
                                                                         len(self.human_homologs[primary_gene])
                                                                         )
                      )
            if len(self.human_homologs[primary_gene]) > 0:
                counter['potential homolog'] += 1
            else:
                counter['no homolog'] += 1
            counter['total'] += len(self.human_homologs[primary_gene])
            counter[primary_gene] = len(self.human_homologs[primary_gene])
        return counter

    def calculate_brain_expression(self):
        """
        calculates how many of the human homologs are expressed in the brain
        :return: a dict with three keys, expressed, not_expressed, unknown
        """
        self.calculate_scores()
        self.get_human_homologs(verbose=False)

        expression = dict(expressed=0, not_expressed=0, unknown=0)
        for human_homolog in self.human_homologs:
            for homolog in self.human_homologs[human_homolog]:
                if homolog.expressed_in_brain:
                    expression['expressed'] += 1
                elif not homolog.expressed_in_brain:
                    expression['not_expressed'] += 1
                else:
                    expression['unknown'] += 1

        return expression


    def get_info_from_uniprot(self, uniprot_id):
        """
        retrieves information from Uniprot and stores it in uniprot_info attribute
        :param uniprot_id: a Uniprot ID, e.g. Q9TU53
        :return: a dictionary with attributes
        """

        attributes = dict()
        attributes['organism'] = ''
        attributes['fullname'] = ''
        attributes['sequence'] = ''
        attributes['gene_id'] = ''
        filename = 'uniprot/{}.xml'.format(uniprot_id)
        if not os.path.isfile(filename):
            r = requests.get('{0}/{1}.xml'.format(self.urls['uniprot'], uniprot_id))
            with open(filename, 'w') as f:
                if r.status_code == 200:
                    f.write(r.text)
                else:
                    #create the file but leave it empty
                    print('Error retrieving Uniprot {}'.format(uniprot_id))
                    return attributes

        #with open(filename, 'r') as f:
        #    uniprot_xml = f.read()
        try:
            root = ET.parse(filename).getroot()
        except:
            print('could not get root for UniProt ID: {}'.format(uniprot_id))
            return attributes
        entries = [node for node in root.getchildren() if 'entry' in node.tag]
        if len(entries) != 1:
            print('malformed XML for UniProt ID: {}\nexpected 1 entry'.format(uniprot_id))


        organism_root = [node for node in entries[0].getchildren() if node.tag.endswith('organism')]
        if len(organism_root) != 1:
            print('malformed XML for UniProt ID: {}\nexpected 1 organism'.format(uniprot_id))
        else:
            organism = [node for node in organism_root[0].getchildren() if
                    (node.tag.endswith('name') and node.attrib.get('type') == 'common')]
            if len(organism) == 1:
                attributes['organism'] = organism[0].text
            else:
                organism = [node for node in organism_root[0].getchildren() if
                        (node.tag.endswith('name') and node.attrib.get('type') == 'scientific')]
                if len(organism) == 1:
                    attributes['organism'] = organism[0].text
                else:
                    print('malformed XML for UniProt ID: {}\nexpected 1 common organism name'.format(uniprot_id))

        seq = [node for node in root[0].getchildren() if node.tag.endswith('sequence')]
        if len(seq) != 1:
            print('malformed XML for UniProt ID: {}\nexptected 1 sequence'.format(uniprot_id))
        else:
            attributes['sequence'] = seq[0].text.strip().replace('\n', '')


        protein = [node for node in root[0].getchildren() if node.tag.endswith('protein')]
        if len(protein) != 1:
            print('malformed XML for UniProt ID: {}\nexptected 1 protein'.format(uniprot_id))
            attributes['fullname'] = ''
        else:
            submittedName = [node for node in protein[0].getchildren() if node.tag.endswith('submittedName')]
            if len(submittedName) != 1:
                submittedName = [node for node in protein[0].getchildren() if node.tag.endswith('recommendedName')]
            if len(submittedName) != 1:
                print('malformed XML for UniProt ID: {}\nexpected 1 submittedName or recommendedName'.format(uniprot_id))
            else:
                attributes['fullname'] = submittedName[0].getchildren()[0].text

        gene_id_root = [node for node in entries[0].getchildren() if node.tag.endswith('gene')]
        gene_found = False
        if len(gene_id_root) != 1:
            #print('malformed XML for UniProt ID: {}\nexpected 1 gene\nTrying via API'.format(uniprot_id))
            pass
        else:
            #print(gene_id_root)
            gene_id = [node for node in gene_id_root[0].getchildren() if
                       (node.tag.endswith('name') and node.attrib.get('type') == 'primary')]
            if len(gene_id) == 1:
                gene_found = True
                attributes['gene_id'] = gene_id[0].text
        if not gene_found:
            #print('Could not find GeneName in XML for UniProt ID: {}\nTrying via API'.format(uniprot_id))
            attributes['gene_id'] = self.uniprot_to_genename(uniprot_id)

        return attributes
    def find_homologs(self):
        for ii, gene in enumerate(self.gene_ids):
            self.primary_genes.append(PrimaryGene(uniprot_id=self.uniprot_ids[ii],
                                                  wormbase_id=gene,
                                                  homolog_finder=self,
                                                  wormbase_name=self.wormbase_names[ii]))
            self.primary_genes[-1].find_homologs()
        #self.calculate_scores()
        #self.get_diseases_for_genes()
        #self.calculate_scores()
        #self.get_disease_statistics()

        for primary in self.primary_genes:
            for uniprot in primary.uniprot_ids_from_hmmer:
                self.create_homolog(primary,
                                    uniprot=uniprot,
                                    source='hmmer')
            for compara in primary.compara:
                self.create_homolog(primary,
                                    gene_name=compara,
                                    source='compara')
            #self.create_homolog_html(primary)

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

        new_homolog = HomologGene(primary,
                                  uniprot_id=uniprot,
                                  gene_name=gene_name,
                                  source=source,
                                  homolog_finder=self,
                                  jackhmmer_files=primary.jackhmmer_filename)
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
                        if error[0] is None:
                            error = (True, d.get('primary_id'))
                    else:
                        error = (False, '')
            if error[0]:
                print('Error in row {0}: Given Uniprot ID {1} does not match retrieved Uniprot ID {2}'.format(i + 1,
                                                                                                              xl_sheet.cell_value(
                                                                                                                  rowx=i,
                                                                                                                  colx=3),
                                                                                                              error[1]))

            if error[0] is None:
                print('Warning in row {0}: No Uniprot ID found'.format(i + 1))
            # ID check
            error = (None, '')
            for d in data:
                if 'wormbase_gseqname' in d.get('dbname'):
                    if not xl_sheet.cell_value(rowx=i, colx=2) == d.get('display_id'):
                        if error[0] is None:
                            error = (True, d.get('display_id'))
                    else:
                        error = (False, '')
            if error[0]:
                print('Error in row {0}: Given display ID {1} does not match retrieved display ID {2}'.format(i + 1,
                                                                                                              xl_sheet.cell_value(
                                                                                                                  rowx=i,
                                                                                                                  colx=2),
                                                                                                              error[1]))

            if error[0] is None:
                print('Warning in row {0}: No Uniprot ID found'.format(i + 1))


            # let's check if the gene matches the ID, UniProt and name

    def get_genes_from_xls_file(self, col):
        """
        Reads the gene IDs from the Excel file
        updates the gene_ids attribute
        :param col: 
        :return: a list of gene IDs
        """

        xl_sheet = self.workbook.sheet_by_index(0)
        cols = xl_sheet.col(col)

        p = re.compile('^.*\(([-.a-zA-Z0-9]+)\)')

        for c in cols:
            m = p.match(c.value)
            if m is not None:
                self.gene_ids.append(m.group(1))
            else:
                self.gene_ids.append('')
        return self.gene_ids

    def get_ids_from_xls_file(self, col):
        """ 
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

    def create_homolog_html(self, primary, only_human=True):
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
        table[0].append('Name')
        table[0].append('Expressed in brain')
        table[0].append('Alignment')
        table[0].append('Diseases')
        table[0].append('Alignment length [%]')
        table[0].append('Total Score')
        filename = 'jackhmmer_{}.aln'.format(primary.uniprot_id)
        for homolog in self.homologs[primary.wormbase_name]:
            if only_human and ('human' not in homolog.organism.lower() and 'homo' not in homolog.organism.lower()):
                continue
            table.append(list())
            table[-1].append('<a href="http://www.uniprot.org/uniprot/{0}">{0}</a>'.format(homolog.uniprot_id))
            table[-1].append(homolog.name)
            if homolog.expressed_in_brain:
                table[-1].append('yes')
            elif not homolog.expressed_in_brain:
                table[-1].append('no')
            else:
                table[-1].append('unknown')
            table[-1].append('<a href="{}.html?uniprot_id={}" target="_blank">Open alignment</a>'.format(filename, homolog.uniprot_id))
            if homolog.disgenet and homolog.disgenet.get('diseases'):
                table[-1].append('<br>'.join(list(set([dd['disease'] for dd in homolog.disgenet['diseases']]))))
            else:
                table[-1].append('no associated diseases')
            table[-1].append(round(homolog.alignment_score * 100, 2))
            table[-1].append(homolog.score)
        result = self.render('table.template.html', {'table': table})
        print(primary.wormbase_name)
        with open('{}.html'.format(primary.wormbase_name.replace('/', '_')), 'w') as f:
            f.write(result)
        self.alignment(filename)


    def alignment(self, filename='align.txt'):
        """
        Creates an alignment file from an alignment, only human sequences are used
        :param filename: the file where the alignment is saved
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


    def uniprot_to_genename_json(self, uniprot):
        """
        gets the gene name for a UniProt identifier
        :param uniprot: a UniProt identifier
        :return: a gene name, in case of multiple gene names only the first is returned
        """

        if not self.uniprot_genename.get(uniprot):
            self.uniprot_genename[uniprot] = None
            filename = 'uniprot/id_{}.txt'.format(uniprot)
            if not os.path.isfile(filename):
                r = requests.get('{0}:{1}&columns=id,entry name,genes&format=json'.format(self.urls['id'],
                                                                                          uniprot))


                if r.status_code == 200:
                    with open(filename, 'w') as f:
                        f.write(r.text)
                    self.uniprot_genename[uniprot] = json.loads(r.text)[0]['genes'].split(' ')[0]
            else:
                with open(filename, 'r') as f:
                    self.uniprot_genename[uniprot] = json.load(f)[0]['genes'].split(' ')[0]

        return self.uniprot_genename[uniprot].split(';')[0]

    def uniprot_to_genename(self, uniprot):
        """
        gets the gene name for a UniProt identifier
        :param uniprot: a UniProt identifier
        :return: a gene name, in case of multiple gene names only the first is returned
        """

        if not self.uniprot_genename.get(uniprot):
            self.uniprot_genename[uniprot] = None
            filename = 'uniprot/id_{}.tsv'.format(uniprot)
            if not os.path.isfile(filename):
                r = requests.get('{0}:{1}&columns=id,entry name,genes&format=tab'.format(self.urls['id'],
                                                                                         uniprot))


                if r.status_code == 200:
                    with open(filename, 'w') as f:
                        f.write(r.text)
                else:
                    raise RuntimeError('could not load UniProt entry: {}'.format(uniprot))
        with open(filename, 'r') as f:
            lines = f.read().splitlines()
        assert len(lines) == 2, 'Uniprot {} has {} lines instead of 2'.format(uniprot, len(lines))
        uniprot_id = lines[1].split('\t')[0].lower()
        assert uniprot_id == uniprot.lower(), 'Uniprot {} has wrong id: '.format(uniprot, uniprot_id)
        assert len(lines[1].split('\t')) == 3, 'Uniprot {} has wrong number of columns: '.format(uniprot)
        self.uniprot_genename[uniprot] = lines[1].split('\t')[-1].split(' ')[0]

        return self.uniprot_genename[uniprot].split(';')[0]

    def get_diseases_for_genes(self):
        """
        
        :return: 
        """
        for gene in self.get_homologs():
            self.get_disease_for_gene(gene.gene_name)


    def get_disease_for_gene(self, gene):
        """
        
        :param gene: a string with GeneID, e.g. SEZ6 
        :return: 
        """
        pass


    def read_disgenet_map(self):
        """
        
        :return: 
        """

        with open(self.disgenet_map_filename, 'r') as f:
            lines = f.read().splitlines()

        self.disgenet_map = dict()
        for line in lines:
            cells = line.split('\t')
            self.disgenet_map[cells[0]] = cells[1]

        return None


    def get_diseases_from_disgenet(self, uniprot_id, filename='', overwrite=False):
        """

        inspired by:
        http://www.disgenet.org/ds/DisGeNET/scripts/disgenet_python3.py
        :param gene:
        :return:
        """

        if not self.disgenet_map:
            self.read_disgenet_map()
        gene_id = self.disgenet_map.get(uniprot_id)
        if not gene_id:
            print('no gene ID found for {}'.format(uniprot_id))
            return None
        if filename == '':
            filename = 'disgenet/{}'.format(uniprot_id)
        print(filename)
        if not overwrite and os.path.isfile(filename):
            with open(filename, 'r') as f:
                return f.read()

        query = """
        DEFINE
          	c0='/data/gene_disease_summary',
	c1='/data/diseases',
	c2='/data/genes',
	c4='/data/sources'
        ON
           'http://www.disgenet.org/web/DisGeNET'
        SELECT
         	c1 (diseaseId, name, diseaseClassName, STY, MESH, OMIM, type ),
	c2 (geneId, symbol,   uniprotId, description, pantherName ),
	c0 (score, EI, Npmids, Nsnps)
           
        FROM
            c0
        WHERE
            (
                c2.geneId = '{0}'
            AND
                c4 = 'ALL'
            )
        ORDER BY
            c0.score DESC""".format(gene_id)
        query = query.encode('utf-8')
        req = urllib.request.Request(self.urls['disgenet'])
        res = urllib.request.urlopen(req, query)
        data = res.read().decode("utf-8").strip()
        if overwrite or not os.path.isfile(filename):
            with open(filename, 'w') as f:
                f.write(data)
        return data

    def get_disgenet(self, uniprot_id, overwrite=False):
        """

        :param uniprot_id: The uniprot ID for which associated diseases should be found
        :return:
        """
        if self.disgenet.get(uniprot_id) is None or overwrite:
            filename = 'disgenet/{}.txt'.format(uniprot_id)
            data = self.get_diseases_from_disgenet(uniprot_id, filename)
            self.disgenet[uniprot_id] = self.parse_disgenet(self.filter_disgenet(data))

        return self.disgenet[uniprot_id]

    def filter_disgenet(self, data):
        """

        :param self: 
        :param data: a result string from a DisGenet query 
        :return: all lines which are associated with nerve disease, taken from self.neuro_diseases
        """
        filtered_diseases = list()
        if not data:
            return filtered_diseases

        lines = data.split('\n')
        for line in lines:
            cells = line.split('\t')
            if len(cells) < 6:
                continue
            if self.disease_synonyms.get('MESH:{}'.format(cells[4])):
                filtered_diseases.append(line)
            if self.disease_synonyms.get('OMIM:{}'.format(cells[5])):
                filtered_diseases.append(line)
        #print(filtered_diseases)
        return filtered_diseases

    def parse_disgenet(self, data):
        """
        
        :param data: 
        :return: 
        """

        parsed_values = dict()
        parsed_values['total score'] = 0
        parsed_values['diseases'] = list()
        for line in data:
            cells = line.split('\t')
            parsed_values['total score'] += float(cells[14])
            parsed_values['diseases'].append(dict(disease=cells[1], score=float(cells[14])))
        return parsed_values


def report(homolog_finder, disease_threshold=[0, 10**6]):
    print('PrimaryGene\tHomolog\tAssociatedDisease')
    for homo in homolog_finder.homologs:
        for x in homolog_finder.homologs[homo]:
            if isinstance(disease_threshold, list):
                if x.disgenet.get('diseases')in disease_threshold:
                    print('{}\t{}\t{}\{}'.format(x.primary.wormbase_name, x.name, '\n\t\t'.join(list(set([dd['disease'] for dd in x.disgenet.get('diseases')]))), homologs.x.expressed_in_brain))
            elif disease_threshold is None:
                if not x.disgenet.get('diseases'):
                    print('{}\\t{}\{}'.format(x.primary.wormbase_name, x.name, homologs.x.expressed_in_brain))


def parse_args(args):
    parser = argparse.ArgumentParser(description='Finds homologs for C. Elegans genes in humans')
    parser.add_argument('input_filename', help='an Excel file which contains the C. Elegans genes')
    parser.add_argument('--jackhmmer',
                        help='location of the jackhmmer executable',
                        default='jackhmmer')
    parser.add_argument('--uniprot',
                        help='location of the UniProt Swiss-Pro FASTA file',
                        default='uniprot_sprot.fasta')
    return parser.parse_args(args)


def validate_args(args):
    if not os.path.isfile(args.input_filename):
        return (False, 'Input Excel file is missing')
    if not os.path.isfile(args.jackhmmer) or not os.access(args.jackhmmer, os.X_OK):
        return (False, 'jackHMMER executable is not found not or cannot be executed')
    if not os.path.isfile(args.uniprot):
        return (False, 'UniProt Swiss-Pro FASTA file is missing')

    return True, ''


if __name__ == '__main__':
    
    args = parse_args(sys.argv[1:])
    if not validate_args(args)[0]:
        print(validate_args(args)[1])
        sys.exit(1)
    homolog_finder = HomologFinder(xls_filename=args.input_filename,
                                       jackhmmer=args.jackhmmer)

    for i in range(len(homolog_finder.homologs.keys())):
        k = list(homolog_finder.homologs.keys())[i]
        print(homolog_finder.homologs[k][0])

    homolog_finder.calculate_scores()
    homolog_finder.get_disease_statistics()

    for homo in homolog_finder.homologs:
        for x in homolog_finder.homologs[homo]:
            if x.disgenet and x.disgenet.get('diseases'):
                print(x)
