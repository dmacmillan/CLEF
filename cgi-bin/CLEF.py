#!/home3/dmacmill/software/anaconda2/bin/python
import cgi
import sys
import math
from scipy import stats
from Codon import *
from Epitope import *
import logging
from itertools import *

logging.basicConfig(filename='log', level='DEBUG', filemode='w')

## Global variables
output_columns = (
    'PATIENT_ID',
    'HIV_PROTEIN',
    'HLA_ALLELE',
    'CTL_EPITOPE',
    'PATIENT_CTL_SEQUENCE',
    'HLA_RESTRICTION',
    'EPITOPE_COORDINATES',
    'EPITOPE_SOURCE',
    'EXPANDED_HLA_DEFINITION'
)

null_char = 'NAN'

epitopes_file = '../epitopes_v1.0.1.txt'
sys.stderr = open("../error-cgi.log", "w")

form = cgi.FieldStorage()
patients = form.getvalue("patients_input")
protein = form.getvalue("protein_selection")

def printHtmlHeaders(content_type='text/html'):
    print "Content-Type: {}".format(content_type)
    print
    print """<!DOCTYPE html><html><head>
    <script src="//ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js"></script>
    <script src="../cgi-bin/script.js"></script>
    <link rel="stylesheet" href="../css/style.css"></head><body>
    <!-- Bootstrap core CSS -->
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0-alpha.6/css/bootstrap.min.css" integrity="sha384-rwoIResjU2yc3z8GV/NPeZWAv56rSmLldC3R/AZzGRnGxQQKnKkoFVhFQhNUwEyJ" crossorigin="anonymous">"""

def printFileHeaders(filename):
    print "Content-Disposition: attachment; filename=\""+filename+"\""
    print "Content-Type:application/octet-stream; name=\""+filename+"\""
    print

# Flag can be 0, 1, or 2 depending on the desired output
# Flag = 1 will output all mixtures as "X"
# Flag = 2 will output all synonymous mixtures as they are and all non-synonymous mixtures as "X"
# Flag = 3 will output all mixtures in the format [A/B] if a mixture encodes for amino acid A or B
def translateDNA(sequence, resolvecharacter="X", flag=3):
    sequence = sequence.translate(None, ' \r\n').upper()
    aaseq = []
    i = 0
    while i < len(sequence):
        try:
            codon = Codon.resolveCodon(sequence[i:i+3])
        except IndexError:
            codon = Codon.resolveCodon('???')
        # If the codon has no mixture bases just add it to the amino acid chain
        if len(codon) <= 1:
            aaseq.append(Codon.codon_dict[codon[0]])
        # Codon contains mixture base
        else:
            # If flag is set to 1
            if (flag == 1):
                aaseq.append(resolvecharacter)
            # If flag is set to 2
            elif (flag == 2):
                unique = set([Codon.codon_dict[potential] for potential in codon])
                # If there is more than resolved one amino acid
                if (len(unique) > 1):
                    aaseq.append(resolvecharacter)
                else:
                    aaseq.append(unique.pop())
            # If flag is set to 3
            else:
                unique = set([Codon.codon_dict[potential] for potential in codon])
                # If there is more than resolved one amino acid
                if (len(unique) > 1):
                    aaseq.append('['+('/').join(unique)+']')
                else:
                    aaseq.append(unique.pop())
        i += 3
    return aaseq

def check_sequence(dna):
    valid_chars = [x for x in 'ACGTURYSWKMBDHVN.-']
    if any([x not in valid_chars for x in dna]):
        return False
    return True
    
def parse(input):
    val = [x.split('\t') for x in input.splitlines()]
    return val

def parseHLA(hla, res=4):
    rval = hla.strip().translate(None, "*:").upper()
    try:
        int(rval[-1])
    except (ValueError, IndexError) as e:
        rval = rval[:-1]
    return rval[:res+1]
    
def parseSeqs(sequences):
    return [translateDNA(x) for x in sequences.splitlines()]

def load_patients(patients):
    d = {}
    for i,patient in enumerate(patients):
        pid = patient[0]
        # Check if last column is DNA sequence
        if check_sequence(patient[-1]):
            d[pid] = {'hlas': set(), 'i': i, 'seq': patient[-1], 'aa': translateDNA(patient[-1])}
        else:
            d[pid] = {'hlas': set(), 'i': i, 'seq': None, 'aa': None}
        for hla in patient[1:-1]:
            hla = parseHLA(hla)
            if (hla == ""):
                continue
            d[pid]['hlas'].add(hla)
    return d

def parseEpitopes(epitopes_file):
    results = {}
    with open(epitopes_file, 'r') as f:
        lines = [x.strip().split('\t') for x in f.readlines()]
        results = lines
    return results
    
def analyze_patient(patient, patient_id, epitopes):
    logging.debug('Analyzing patient: {}'.format(patient))
    results = []
    for hla in patient['hlas']:
        result = {
            'pid': patient_id,
            'hla': hla,
            'epitopes': set(),
            'expanded': False
        }
        for epitope in epitopes:
            if (hla in epitope.hlas):
                logging.debug('patient {} with {} matches epitope: {}'.format(patient_id, hla, epitope.hlas))
                result['epitopes'].add(epitope)
            elif (hla in epitope.r2):
                logging.debug('patient {} with {} matches epitope: {}'.format(patient_id, hla, epitope.r2))
                result['epitopes'].add(epitope)
            elif (hla in epitope.r4):
                logging.debug('patient {} with {} matches epitope: {}'.format(patient_id, hla, epitope.r4))
                result['epitopes'].add(epitope)
                result['expanded'] = True
        if len(result['epitopes']) > 0:
            results.append(result)
    return results

def print_table_header(columns):
    return ('').join(['<th>{}</th>'.format(x) for x in columns])

def print_result(result, patients):
    expanded = '1' if result['expanded'] else '0'
    # print 'pid: {}'.format(result['pid'])
    # for ep in result['epitopes']:
        # print 'ep: {}'.format(ep)
        # print
    lines = [
        ('<td>{}</td>'*9).format(
            result['pid'],
            protein,
            result['hla'],
            (',').join([x for x in ep.epitope]),
            patients[result['pid']]['aa'][ep.start - 1 : ep.end] if patients[result['pid']]['aa'] else null_char,
            (',').join([hla for hla in ep.hlas]),
            '{}-{}'.format(ep.start, ep.end),
            ep.source,
            expanded
        ) for ep in result['epitopes']
    ]
    # print 'lines: {}'.format(lines)
    # print
    return ('').join(['<tr>{}</tr>'.format(x) for x in lines])
    
# def print_result(result, patients):
    # text0 = '<tr>' + \
    # ('<td>{}</td>' * 3).format(
        # result['pid'],
        # protein,
        # result['hla']
    # )
    # text1 = (',').join(['('+(',').join(x.epitope)+')' for x in result['epitopes']])
    # patient_seq = 'NAN'
    # if patients[result['pid']]['aa']:
        # patient_seq = (',').join(['({})'.format(('').join(patients[result['pid']]['aa'][x.start - 1 : x.end])) for x in result['epitopes']])
    # if text1:
        # text2 = ('<td>{}</td>' * 6).format(
            # text1,
            # patient_seq,
            # (',').join(['('+(',').join(x.hlas)+')' for x in result['epitopes']]),
            # (',').join(['({}-{})'.format(x.start, x.end) for x in result['epitopes']]),
            # (',').join(['({})'.format(x.source) for x in result['epitopes']]),
            # 'Y' if result['expanded'] else 'NAN'
        # ) + '</tr>'
    # else:
        # text2 = ('<td>NAN</td>' * 6) + '</tr>'
    # return text0 + text2

def display_results(results, patients, protein, output_columns):
    print '<table class="table table-bordered" id="output_table">'
    print print_table_header(output_columns)
    for result in results:
        print print_result(result, patients)
    print '</table>'
    
output_format = ('\t').join([
    'PATIENT_ID',
    'HIV_PROTEIN',
    'HLA_ALLELE',
    'CTL_EPITOPE',
    'PATIENT_CTL_SEQUENCE',
    'HLA_RESTRICTION',
    'EPITOPE_COORDINATES',
    'EPITOPE_SOURCE',
    'EXPANDED_HLA_DEFINITION'
])

epitopes = Epitope.parseEpitopes(epitopes_file, only_proteins=[protein])
for e in epitopes:
    for i,hla in enumerate(e.hlas):
        e.hlas[i] = parseHLA(hla)
    for i,hla in enumerate(e.r2):
        e.r2[i] = parseHLA(hla)
    for i,hla in enumerate(e.r4):
        e.r4[i] = parseHLA(hla)
    e.hlas = set(e.hlas)
    e.r2 = set(e.r2)
    e.r4 = set(e.r4)

#printHtmlHeaders(content_type='text/plain')
printHtmlHeaders()

patients = load_patients(parse(patients))

results = []
for patient_id in patients:
    hits = analyze_patient(patients[patient_id], patient_id, epitopes)
    for hit in hits:
        results.append(hit)

display_results(results, patients, protein, output_columns)