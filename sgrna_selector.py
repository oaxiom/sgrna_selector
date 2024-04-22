
import sys, os
from glbase3 import *

if len(sys.argv) != 4:
    print('Usage: sgrna_selector.py genenames_file.txt libraryname num_sgrnas_per_gene')
    sys.exit()

all_genes = genelist(sys.argv[1], format={'name': 0})

max_sgrnas = int(sys.argv[3])
library_name = sys.argv[2]

gf = {'force_tsv': True,
    'sgrna': {'code': 'column[0][:-3]'}, # 0,
    'name': 3}

df = {'force_tsv': True,
    'sgrna': 7, # {'code': '"{}{}".format(column[7], column[9])'},
    'name': 2}

tko3 = {'force_tsv': True,
    'sgrna': 1,
    'name': 0}

Source_sgRNA_table = {'force_tsv': True,
    'sgrna': {'code': 'column[3][:-3]'}, # 0,
    'original_library': 1,
    'name': 2}

crisprab = {'force_tsv': True,
    'name': 0,
    'sgrna': 1}

# Some gene names are the old names in the KO libs;
aliases = {
    'ABRAXAS1': 'FAM175A',
    'BABAM2': 'BRE',
    'ELP1': 'IKBKAP',
    'EMSY': 'C11orf30',
    'KAT14': 'CSRP2BP',
    'KMT5A' : 'SETD8',
    'KMT5B': 'SUV420H1',
    'KMT5C': 'SUV420H2',
    'NSD2': 'WHSC1',
    'NSD3': 'WHSC1L1',
    'OGA': 'MGEA5',
    'PPP4R3A': 'SMEK1',
    'PPP4R3B': 'SMEK2',
    #'RIOX1': '',
    'RIOX2': 'MINA',
    'RTL5': 'RGAG4',
    'SGF29': 'CCDC101',
    'SHTN1': 'KIAA1598',
    'SLF1': 'ANKRD32',
    'TOMM70': 'TOMM70A',

    #'RYBP': '', # Seems impossible to generate sgRNAs?
    'PWWP3A': 'MUM1',
    'AC006064.6': '', # NOD2 CHD4 fusion transcript
    'AC012184.2': '', # DDX19A and DDX19B fusion transcript;
    'AC106886.5': 'SRCAP',
    'AC118549.1': '', #?
    'AL022318.4': '', # APOBECC D fusion transcript;
    'AL031681.2': '', # SRSF6/L3MBTL1 fusion
    'AL133500.1': '', # CITED1/HDAC8
    'AL360181.3': '', # PAOX/MTG1 fusion
    'AL451062.3': '', #?
    'FP565260.4': '', #?
    'KDM4F': '', # single exon transcript, looks kind of odd;
    'PPP4R3C': '', # predicted noncoding;
    'TSPY9P': '', # predicted noncoding;
    'WASHC1': '', #?
    }

# The list of genes that are lethal, from : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5555476/
core_essential_genes = genelist('datasets/CEG2.txt.gz', format={'force_tsv': True, 'name': 0}, gzip=True)

newl = []
for gn in all_genes:
    if '-' in gn['name']:
        continue
    newl.append(gn)
newl.append({'name': 'CTRL LacZ', 'source': 'Control'}) # will grab two ctrls from TKO3;
newl.append({'name': 'CTRL luciferase', 'source': 'Control'})
newl.append({'name': 'CTRL NONTARG', 'source': 'Control'})
all_genes.load_list(newl)

minlibcas9 = genelist('datasets/sgRNA_libraries/MinLibCas9.txt.gz', format=gf, gzip=True)
minlibcas9 = minlibcas9.addEmptyKey('original_library', 'MinLibCas9')
brunello = genelist('datasets/sgRNA_libraries/Brunello.txt.gz', format=df, gzip=True)
brunello = brunello.addEmptyKey('original_library', 'Brunello')
tko3 = genelist('datasets/sgRNA_libraries/tkov3_guide_sequence.txt.gz', format=tko3, gzip=True)
tko3 = tko3.addEmptyKey('original_library', 'TKOv3')
super_table = genelist('datasets/sgRNA_libraries/Source_sgRNA_table.txt.gz', format=Source_sgRNA_table, gzip=True) # Table S2 from the Minlabcas9.
super_table = super_table.addEmptyKey('original_library', 'super_table')
crisprab = genelist('datasets/sgRNA_libraries/CRISPRab.packed.tsv.gz', format=crisprab, gzip=True)
crisprab = crisprab.addEmptyKey('original_library', 'CRISPRA/B')

# rank by confidence:
minlibcas9 = minlibcas9.addEmptyKey('rank', 1)
brunello = brunello.addEmptyKey('rank', 2)
tko3 = tko3.addEmptyKey('rank', 3)
super_table = super_table.addEmptyKey('rank', 3)
crisprab = crisprab.addEmptyKey('rank', 3)

all_sgrnas = minlibcas9 + brunello + tko3 + super_table + crisprab

sgrna_hits = all_genes.map(genelist=all_sgrnas, key='name', greedy=True)

# Do a rescue for the old name genes:
for gene_name in aliases:
    if not aliases[gene_name]:
        continue
    r = all_sgrnas.get(key='name', value=aliases[gene_name])
    if r:
        print('Rescued {}/{}'.format(gene_name, aliases[gene_name]))
        for g in r:
            g['name'] = gene_name
        sgrna_hits += r

sgrna_hits = sgrna_hits.removeDuplicates('sgrna')
sgrna_hits.sort('name')
sgrna_hits.saveTSV('sgrna_{}-all.tsv'.format(library_name), key_order=['name', 'source', 'original_library', 'sgrna'])
sgrna_hits_gene_names = sgrna_hits.removeDuplicates('name').getColumns(['name'])
sgrna_hits_gene_names.saveTSV('sgrna_{}-genes.tsv'.format(library_name), key_order=['name'])

num_genes_with_1_sgrna = len(sgrna_hits.removeDuplicates('name'))

ceg_hits = core_essential_genes.map(genelist=sgrna_hits, key='name').removeDuplicates('name')

print()
for g in ceg_hits:
    print('CEG2 match: {}'.format(g['name']))

num_hits = {}
for sgrna in sgrna_hits:
    if sgrna['name'] not in num_hits:
        num_hits[sgrna['name']] = 0
    num_hits[sgrna['name']] += 1

__less_than_2sgRNAs = 0
__more_than_5sgRNAs = 0
print()
for g in sorted(num_hits):
    if num_hits[g] <= 2:
        print('{} = {} sgRNAs'.format(g, num_hits[g]))
        __less_than_2sgRNAs += 1
    elif num_hits[g] > 5:
        __more_than_5sgRNAs += 1

print()
for g in all_genes:
    if g['name'] not in num_hits:
        print('{} = 0 sgRNAs!'.format(g['name']))

print()
print('Full set:')
print('  All genes                        : {}'.format(len(all_genes)))
print('  Genes with at least 1 sgRNA      : {}'.format(num_genes_with_1_sgrna))
print('  Number of sgRNAs                 : {}'.format(len(sgrna_hits)))
print('  Average sgRNAs/gene              : {:.1f}'.format(len(sgrna_hits) / num_genes_with_1_sgrna))
print('  Number of genes with <= 2 sgRNAs : {}'.format(__less_than_2sgRNAs))
print('  Number of genes with > 5 sgRNAs  : {}'.format(__more_than_5sgRNAs))
print('  CEG2 hits                        : {}'.format(len(ceg_hits)))
print()

# Filter to max_sgrnas based on the original libraries
max_sgrnas = int(sys.argv[3])

sgrna_hits.sort('rank')

by_gene = {}
for hit in sgrna_hits:
    if hit['name'] not in by_gene:
        by_gene[hit['name']] = []
    by_gene[hit['name']].append(hit)

final_set = []
for gene in by_gene:
    if 'CTRL ' in gene:
        final_set += by_gene[gene][0:59]
        continue # add all CTRLs;
    final_set += by_gene[gene][0:max_sgrnas]

# get final sgrna counts per genes;

gl = genelist()
gl.load_list(final_set)
gl.sort('name')
gl.saveTSV('sgrna_{}-final.tsv'.format(library_name), key_order=['name', 'source', 'original_library', 'sgrna'])
gl = gl.renameKey('sgrna', 'seq')
gl = gl.joinKey('name_sgrna', '{}-{}', 'name', 'seq', keep_originals=True)
gl.saveFASTA('sgrna_{}-final.fasta'.format(library_name), name='name_sgrna')
gl_names = len(set(gl['name']))
num_ctrls = [g['name'] for g in gl if 'CTRL ' in g['name']]

print()
print('Final set:')
print('  Number of sgRNAs                 : {}'.format(len(gl)))
print('  Average sgRNAs/gene              : {:.1f}'.format(len(gl) / gl_names))
print('  Total number of CTRL sgRNAs      : {}'.format(len(num_ctrls)))
print()
