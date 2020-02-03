#!/usr/bin/env python
"""
Requires biopython and beautifulsoup4.
"""
from Bio import Entrez
from argparse import ArgumentParser
from bs4 import BeautifulSoup

from operator import itemgetter
import requests
import os
import time

def parse_meta(meta):
    """Parse out the N50, etc. from the XML in the "Meta" tag."""
    soup = BeautifulSoup(meta, 'html5lib')
    return dict([(s['category'], int(s.text)) for s in soup.find_all('stat') if s['sequence_tag'] == 'all'])

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--limit', type=int,
                        help='Only process this many new assemblies')
    return parser.parse_args()

def main():
    opts = parse_args()

    Entrez.email = 'joelcarmstrong@gmail.com'
    # Restrict search to vertebrate assemblies, and those that are
    # annotated as being whole-genome rather than single-chromosome or
    # mitochondrial.
    handle = Entrez.esearch(db='assembly', term='txid7742[Organism:exp] AND "full genome representation"[Properties]', retmax='10000')
    search_results = Entrez.read(handle)
    handle.close()
    # When this fails it's time to paginate
    assert int(search_results['Count']) < 10000

    search_results['IdList'].reverse()
    ids = map(int, search_results['IdList'])
    if opts.limit is not None:
        # Only do the first N new assemblies.
        ids = list(ids)[:opts.limit]

    print "\t".join(["AssemblyName", "Accession", "ScaffoldN50", "ContigN50", "TotalLength", "SubmissionDate", "Species", "SpeciesCommonName", "SpeciesTaxID"])
    for id in ids:
        handle = Entrez.efetch(db='assembly', id=id, rettype='docsum')
        record = Entrez.read(handle, validate=False)
        summary = record['DocumentSummarySet']['DocumentSummary'][0]
        stats = parse_meta(summary['Meta'])
        name = summary['AssemblyName']
        accession = summary['AssemblyAccession']
        scaffold_n50 = stats['scaffold_n50']
        contig_n50 = stats['contig_n50']
        size = stats['total_length']
        submission_date = summary['SubmissionDate']
        handle = Entrez.efetch(db='taxonomy', id=int(summary['SpeciesTaxid']), rettype='full', retmode='xml')
        record = Entrez.read(handle)[0]
        lineage = dict((i['Rank'], i['ScientificName']) for i in record['LineageEx'])
        species_common_name = record.get('CommonName') or record.get('OtherNames', {}).get('GenbankCommonName') or record['ScientificName']
        species_scientific=record['ScientificName']
        print "\t".join(map(str, [name, accession, scaffold_n50, contig_n50, size, submission_date, species_scientific, species_common_name, int(summary['SpeciesTaxid'])]))
        # Sleep to avoid Genbank IP-banning us.
        time.sleep(1)

if __name__ == '__main__':
    main()
