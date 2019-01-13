import sys
import os
import json
import requests
from Bio import AlignIO


class HomologGene():
    """
    stores all the information about homologs
    """

    def __init__(self,
                 primary,  # the origin of the homology, a primary gene object
                 uniprot_id='',  # the uniprot ID, used as primary identifier
                 gene_name='',  # the gene name, used as primary identifier
                 jackhmmer_files='',  # the filenames with the jackhmmer alignment
                 jackhmmer_ids=list(),  # a list of uniprot IDs identified by jackhmmer
                 source='',  # from where the homology was deduced, e.g. hmmer or compara
                 homolog_finder=None,
                 score=None,
                 # the best (for multiple homologs) weighted score taking alignment, expression, etc. into account
                 description='',  # description, taken from UniProt
                 name='',  # human readable name of the homolog
                 organism='',
                 sequence=''):
        self.primary = primary
        self.uniprot_id = uniprot_id
        self.score = score
        self.alignment_length = 0
        self.source = source
        self.description = description
        self.name = name
        self.organism = organism
        self.sequence = sequence
        self.disgenet = None  # a placeholder for info from disgenet
        if isinstance(jackhmmer_files, list):
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

        if uniprot_id and homolog_finder:
            uniprot_info = homolog_finder.get_info_from_uniprot(uniprot_id)
            # print(uniprot_info)
            if not self.name:
                self.name = uniprot_info['fullname']
            # if not self.description:
            #    self.description = uniprot_info['description']
            if not self.sequence:
                self.sequence = uniprot_info['sequence']
            if not self.organism:
                self.organism = uniprot_info['organism']

    def __repr__(self):
        return 'Primary: {0}\nName: {1}\nUniprot: {2}\nOrganism: {3}\nScore: {4}\nExpressed in Brain: {5}\n'.format(
            self.primary.wormbase_name, self.name, self.uniprot_id, self.organism, self.score, self.expressed_in_brain)

    def add_info(self,
                 jackhmmer_file=None,
                 jackhmmer_ids=None,
                 compara=None,
                 primary=None):

        if primary:
            jackhmmer_ids = primary.uniprot_ids_from_hmmer
            jackhmmer_file = primary.jackhmmer_filename
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
                # print(r.text)
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

    def calculate_score(self,
                        weights=dict(brain=dict(values=[-1, 0, 1],
                                                weight=1),
                                     hmmer=dict(weight=1)),
                        homolog_finder=None):
        values = (False, None, True)
        if self.score is None:
            # brain
            if self.expressed_in_brain in values:
                self.score = weights['brain']['values'][values.index(self.expressed_in_brain)]
                self.score *= weights['brain']['weight']
            # diseases
            if self.disgenet:
                if self.disgenet.get('total score'):
                    self.score += self.disgenet.get('total score')
            else:
                self.disgenet = self.primary.homolog_finder.get_disgenet(self.uniprot_id)
            # alignment length
            self.calculate_alignment_score()
            self.score += self.alignment_score

        return self.score

    def calculate_alignment_score(self):

        score = 0
        for filename in self.jackhmmer_files:
            if filename in self.primary.homolog_finder.alignments:
                alignment = self.primary.homolog_finder.alignments[filename]
            else:
                alignment = AlignIO.read(filename, "stockholm")
                self.primary.homolog_finder.alignments[filename] = alignment
            for r, record in enumerate(alignment):
                if r == 0:
                    reference = str(record.seq)
                    if '_' in reference:
                        print('argh, damm it')
                        assert True == False, 'breaking here'
                elif self.uniprot_id in record.id:
                    # print(record.id, record.seq)
                    for pos, aa in enumerate(str(record.seq)):
                        if aa != '-':
                            reference = reference[0:pos] + '_' + reference[pos + 1:]
            score += reference.count('_') / len(reference)
        self.alignment_score = score
        return score