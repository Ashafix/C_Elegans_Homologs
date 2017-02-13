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

results = dict()
results['wormbase'] = dict()
results['compara'] = dict()

#to do
#all urls in variables

urls = dict()
urls['xrefs'] = 'http://rest.ensembl.org/xrefs/id'

def initialize():
    """
    creates all input folders needed to run the script
    :return: None
    """
    folders = ['brain', 'compara', 'wormbase_protein']
    for folder in folders:
        if not os.path.isdir(folder):
            os.makedirs(folder)

def validate(filename='57_top_candidates.xlsx'):
    """
    validates the input file
    :return:
    """
    #let's check if the gene matches the ID, UniProt and name
    xl_workbook = xlrd.open_workbook(filename)
    xl_sheet = xl_workbook.sheet_by_index(0)
    wb_genes = xl_sheet.col(0)
    print('validating genes')
    for i, wb_gene in enumerate(wb_genes):
        error = (None, '')
        print('Validating {0}'.format(wb_gene.value))
        r = requests.get('{0}/{1}?content-type=application/json'.format(urls['xrefs'], wb_gene.value.strip(), timeout=10))
        if r.status_code == 200:
            data = json.loads(r.text)
        else:
            print('error retrieving data for {0}'.format(wb_gene))
            data = list()
        #UniProt check
        for d in data:
            if 'UniProt' in d.get('db_display_name'):
                if not xl_sheet.cell_value(rowx=i, colx=3) == d.get('primary_id'):
                    if error[0] ==  None:
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
        #ID check
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

def get_genes_from_xls_file(filename, col):

    xl_workbook = xlrd.open_workbook(filename)

    xl_sheet = xl_workbook.sheet_by_index(0)
    cols = xl_sheet.col(col)

    p = re.compile('^.*\(([-.a-zA-Z0-9]+)\)')
    genes = list()

    for c in cols:
        m = p.match(c.value)
        if m:
            genes.append(m.group(1))
    return genes

def get_ids_from_xls_file(filename, col):
    xl_workbook = xlrd.open_workbook(filename)

    xl_sheet = xl_workbook.sheet_by_index(0)
    cols = xl_sheet.col(col)
    return [c.value for c in cols]

def expressed_in_brain(gene, donor='15496', thresh_avg=3, thresh_max=5):
    """Checks if a gene is expressed in the brain via the Allen Brain Atlas

    :param gene: a gene name according to HGNC
    :param donor: the id of the brain donor
    :param thresh_avg: threshold for average in all measurements
    :param thresh_max: threshold which needs to met at least once
    :return: True, if above thresholds, False, if below thresholds
    """
    exp_max = 0
    exp_avg = 0
    probes = list()
    filename = 'brain/{0}_probes.txt'.format(gene)
    if os.path.isfile(filename):
        with open(filename, 'r') as f:
            msg = json.load(f)
    else:
        r = requests.get("http://api.brain-map.org/api/v2/data/query.json?criteria=model::Probe,rma::criteria,[probe_type$eq'DNA'],products[abbreviation$eq'HumanMA'],gene[acronym$eq'{0}'],rma::options[only$eq'probes.id']".format(gene),
                     timeout=10)

        if r.status_code != 200:
            return None
        msg = json.loads(r.text)['msg']
        with open(filename, 'w') as f:
            json.dump(msg, f)
    for id in msg:
        probes.append(str(id['id']))
    filename = 'brain/{0}_expression_donor_{1}.txt'.format(gene, donor)
    if os.path.isfile(filename):
        with open(filename, 'r') as f:
            probes = json.load(f)
    else:
        r = requests.get("http://api.brain-map.org/api/v2/data/query.json?criteria=service::human_microarray_expression[probes$eq{0}][donors$eq{1}]".format(','.join(probes), donor), timeout=10)
        if r.status_code == 200:
            probes = json.loads(r.text)['msg']['probes']
            with open(filename, 'w') as f:
                json.dump(probes, f)
        else:
            return None
    for probe in probes:
        exp_level = [float(p) for p in probe['expression_level']]
        exp_max += max(exp_level)
        exp_avg += sum(exp_level)/len(exp_level)

    exp_max /= len(probes)
    exp_avg /= len(probes)

    if exp_avg > thresh_avg and exp_max > thresh_max:
        return True
    else:
        return False

def render(tpl_path, context):
    """

    :param tpl_path: location of the template
    :param context: a dictionary which is used for the template
    :return: the parsed template as a string
    """
    path, filename = os.path.split(tpl_path)
    return jinja2.Environment(
        loader=jinja2.FileSystemLoader(path or './')
    ).get_template(filename).render(context)

def alignment(filename='align.txt'):
    """
    Creates an alignment file from an alignment, only human sequences are used
    :param filename: the file where the alignment is aved
    :return: None
    """
    alignment = AlignIO.read(filename, "stockholm")
    items = list()
    for a in alignment:
        if 'HUMAN' in a.id:
            items.append({'id': a.id, 'sequence': str(a.seq)})

    result = render('alignment.template', {'items': items})

    with open('output.html', 'w') as f:
        f.write(result)


def hmmer_to_table(filename='h.txt'):

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
        gene = uniprot_to_genename(line[0].split('|')[1])
        lines[i + 1][1] = expressed_in_brain(gene)

    result = render('table.template.html', {'table': lines})
    with open('table.html', 'w') as f:
        f.write(result)


def uniprot_to_genename(uniprot):
    r = requests.get('http://www.uniprot.org/uniprot/?query=id:{0}&columns=id,entry%20name,genes&format=json'.format(uniprot))
    if r.status_code == 200:
        return json.loads(r.text)[0]['genes'].split(' ')[0]
    else:
        return None

def main():
    genes = get_genes_from_xls_file('57_top_candidates.xlsx', 4)
    ids = get_ids_from_xls_file('57_top_candidates.xlsx', 2)

    ortho_info = dict()

    s = requests.Session()
    s.mount('http://rest.ensembl.org', requests.adapters.HTTPAdapter(max_retries=10))

    for gene in genes:
        filename = 'compara/{0}.txt'.format(gene)
        if not os.path.isfile(filename):
            r = requests.get('https://rest.ensembl.org/homology/symbol/caenorhabditis_elegans/{0}?content-type=application/json&target_species=homo_sapiens'.format(gene))
            ortho_info[gene] = r.text
            with open(filename, 'w') as f:
                f.write(r.text)
        else:
            with open(filename.format(gene), 'r') as f:
                ortho_info[gene] = f.read()



    for id in ids:
        filename = 'wormbase_protein/{0}.html'.format(id)

        if not os.path.isfile(filename):
            r = requests.get('http://www.wormbase.org/rest/widget/cds/{0}/overview'.format(id))
            if r.status_code != 200:
                r = requests.get('http://www.wormbase.org/rest/widget/cds/{0}a/overview'.format(id))
            data = r.text.encode('ascii', errors='ignore')
            with open(filename, 'wb') as f:
                f.write(data)
        else:
            with open(filename, 'rb') as f:
                data = f.read()

        filename = 'wormbase_protein/{0}_protein.html'.format(id)
        if not os.path.isfile(filename):
            soup = bs4.BeautifulSoup(data)
            links = soup.find_all('a', {'class': 'protein-link'})
            if links:
                r = requests.get('http://www.wormbase.org/rest/widget/protein/{0}/overview'.format(links[0].get('href').split('/')[-1]))
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
            results['wormbase'][id] = [links[0].get('href').split(':')[-1], links[0].text]
        else:
            results['wormbase'][id] = ['', '']

    ortho_human = dict()
    for ortho in ortho_info:
        j = json.loads(ortho_info[ortho])
        results['compara'][ortho] = dict(ids=list())
        if not 'error' in j.keys():
            for i in range(len(j['data'][0]['homologies'])):
                results['compara'][ortho]['ids'].append(j['data'][0]['homologies'][i]['target']['id'])

    aliases = dict()
    for id in results['wormbase']:
        aliases[id] = dict()

        filename = 'aliases/{0}.txt'.format(id)
        if os.path.isfile(filename):
            with open(filename, 'r') as f:
                aliases.update(json.load(f))
        else:
            r = requests.get('http://rest.ensembl.org/overlap/translation/{0}?content-type=application/json'.format(results['wormbase'][id][0]), timeout=10)

            if r.status_code == 200:
                r = requests.get('http://rest.ensembl.org/overlap/id/{0}?feature=gene;content-type=application/json'.format(json.loads(r.text)[0]['Parent']), timeout=10)
            else:
                r = requests.get('http://rest.ensembl.org/overlap/translation/{0}?content-type=application/json'.format(id))

            if r.status_code == 200 and 'gene_id' in json.loads(r.text)[0].keys():
                aliases[id]['ensg'] = json.loads(r.text)[0]['gene_id']
                r = requests.get('{0}/{1}?content-type=application/json&external_db=hgnc'.format(urls['xrefs'], aliases[id]['ensg']), timeout=10)
                if r.status_code == 200 and len(json.loads(r.text)) > 0:
                    aliases[id]['hgnc'] = json.loads(r.text)[0]['display_id']

            with open(filename, 'w') as f:
                json.dump({id: aliases[id]}, f)

    for gene in results['compara']:
        if len(results['compara'][gene]['ids']):
            for id in results['compara'][gene]['ids']:
                filename = 'aliases/{0}.txt'.format(id)
                if os.path.isfile(filename):
                    with open(filename, 'r') as f:
                        aliases.update(json.load(f))
                else:
                    r = requests.get('{0}/{1}?content-type=application/json&external_db=hgnc'.format(urls['xrefs'], id), timeout=10)
                    aliases[id] = dict()
                    if r.status_code == 200:
                        aliases[id]['hgnc'] = json.loads(r.text)[0]['display_id']
                    with open(filename, 'w') as f:
                        json.dump({id: aliases[id]}, f)
    hgnc = set()
    for alias in aliases:
        if 'hgnc' in aliases[alias].keys():
            hgnc.add(aliases[alias]['hgnc'])

def phmmer(uniprot):
    filename = '{0}.fasta'.format(uniprot)
    phmmer_filename = '/media/ashafix/WesternDigital/tools/hmmer-3.1b2/src/phmmer'
    database = '/media/ashafix/WesternDigital/data/uniprot_sprot.fasta'
    if not os.path.isfile(filename):
        r = requests.get('http://www.uniprot.org/uniprot/P90732.fasta')
        with open(filename, 'w') as f:
            f.write(r.text)
    args = [phmmer_filename, filename, database]
    p = subprocess.check_output(args)

