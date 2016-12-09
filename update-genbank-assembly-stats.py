#!/usr/bin/env python
"""

"""
from Bio import Entrez
from argparse import ArgumentParser
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine, Column, Integer, String, ForeignKey, Boolean
from sqlalchemy.orm import relationship, sessionmaker
from bs4 import BeautifulSoup
import inflect
from sonLib.bioio import popenCatch

from operator import itemgetter
import requests
import os
import time

inflect_eng = inflect.engine()

Base = declarative_base()

class Assembly(Base):
    __tablename__ = 'assemblies'
    id = Column(Integer, primary_key=True)
    species_id = Column(Integer, ForeignKey('species.id'))
    name = Column(String)
    accession = Column(String, unique=True)
    scaffold_n50 = Column(Integer)
    contig_n50 = Column(Integer)
    size = Column(Integer)
    num_ns = Column(Integer)
    num_masked = Column(Integer)
    ftp_url = Column(String)
    this_assembly_is_terrible = Column(Boolean)

    species = relationship('Species', back_populates='assemblies')

class Species(Base):
    __tablename__ = 'species'
    id = Column(Integer, primary_key=True)
    name = Column(String, unique=True)
    common_name = Column(String)
    genus = Column(String)
    family = Column(String, index=True)
    order = Column(String, index=True)
    class_ = Column(String)

    assemblies = relationship('Assembly', back_populates='species')

def get_position_in_rank(session, assembly, stat, rank, name):
    my_stat = assembly.__dict__[stat]
    assemblies = get_assemblies_in_rank(session, rank, name)
    stat_list = map(lambda x: x.__dict__[stat], assemblies)
    if assembly not in assemblies:
        stat_list.append(my_stat)
    ranked_list = sorted(stat_list, reverse=True)
    return ranked_list.index(my_stat)

def get_assemblies_in_rank(session, rank, name):
    if name is None:
        return []
    # Construct a filter-function for the rank and name.
    filters = {'species': lambda x: x.name == name,
               'genus': lambda x: x.genus == name,
               'family': lambda x: x.family == name,
               'order': lambda x: x.order == name,
               'class': lambda x: x.class_ == name}
    filter_fn = filters[rank]
    # Get a list-of-lists of the assemblies in the rank.
    assemblies_lists = map(lambda x: x.assemblies, session.query(Species).filter(filter_fn(Species)).all())
    # Flat-map
    return [i for list in assemblies_lists for i in list]

def get_ncbi_url(accession):
    return "https://www.ncbi.nlm.nih.gov/assembly/%s" % accession

def get_fasta_gz_url(ftp_url):
    paths = ftp_url.split('/')
    return ftp_url + '/' + paths[-1] + '_genomic.fna.gz'

def get_minimal_rank(session, species, at_least=1):
    nth = float('infinity')
    rank = None
    name = None
    for i in ['class', 'order', 'family', 'genus', 'species']:
        attr = 'class_' if i == 'class' else i
        attr = 'name' if attr == 'species' else attr
        if species.__dict__[attr] == None:
            continue
        n = len(get_assemblies_in_rank(session, i, species.__dict__[attr]))
        if at_least == 1 and nth == 1:
            # HACK: so the "first" is always the maximal rank, which
            # is typically the most interesting data point
            break
        if rank is None or (n <= nth and n >= at_least):
            nth = n
            rank = i
            name = species.__dict__[attr]
    return nth, rank, name

def percentile(session, assembly, stat, rank, name):
    pos = get_position_in_rank(session, assembly, stat, rank, name)
    num_assemblies = len(get_assemblies_in_rank(session, rank, name))
    return (num_assemblies - pos) * 100/num_assemblies

def median(session, assembly, stat, rank, name):
    assemblies = get_assemblies_in_rank(session, rank, name)
    stat_list = map(lambda x: x.__dict__[stat], assemblies)
    if assembly not in assemblies:
        stat_list.append(my_stat)
    ranked_list = sorted(stat_list, reverse=True)
    return ranked_list[len(ranked_list)/2]

def post_to_slack(hook_url, assembly, session):
    MIN_SPECIES_TO_COMPARE = 4

    species = assembly.species
    text = "New GenBank assembly <{url}|*{name}*> for *{common}* (_{species}_).\n".format(name=assembly.name, url=get_ncbi_url(assembly.accession), common=species.common_name, species=species.name)

    nth, rank, name = get_minimal_rank(session, species)
    text += "This is the {nth} assembly for the {rank} _{name}_.\n".format(nth=inflect_eng.ordinal(inflect_eng.number_to_words(nth)), rank=rank, name=name)

    text += "_{species}_ ({num_sp}) < _{genus}_ ({num_gn}) < _{family}_ ({num_fam}) < _{order}_ ({num_ord}) < _{class_}_ ({num_cl})\n".format(species=species.name, num_sp=len(get_assemblies_in_rank(session, "species", species.name)),
    genus=species.genus, num_gn=len(get_assemblies_in_rank(session, "genus", species.genus)),
    family=species.family, num_fam=len(get_assemblies_in_rank(session, "family", species.family)),
    order=species.order, num_ord=len(get_assemblies_in_rank(session, "order", species.order)),
    class_=species.class_, num_cl=len(get_assemblies_in_rank(session, "class", species.class_)))

    _, rank_to_compare, name_to_compare = get_minimal_rank(session, species, at_least=MIN_SPECIES_TO_COMPARE)

    num_assemblies_compared = len(get_assemblies_in_rank(session, rank_to_compare, name_to_compare))

    text += "*Scaffold N50*: {n50:,}".format(n50=assembly.scaffold_n50)
    if num_assemblies_compared >= MIN_SPECIES_TO_COMPARE:
        text += " ({pct} percentile within {rank} _{name}_)".format(
            pct=inflect_eng.ordinal(percentile(session, assembly, "scaffold_n50", rank_to_compare, name_to_compare)),
            rank=rank_to_compare,
            name=name_to_compare)
    text += "\n"
    text += "*Contig N50*: {n50:,}".format(n50=assembly.contig_n50)
    if num_assemblies_compared >= MIN_SPECIES_TO_COMPARE:
        text += " ({pct} percentile within {rank} _{name}_)".format(
            pct=inflect_eng.ordinal(percentile(session, assembly, "contig_n50", rank_to_compare, name_to_compare)),
            rank=rank_to_compare,
            name=name_to_compare)
    text += "\n"
    text += "*Size*: {size:,}".format(size=assembly.size)
    if num_assemblies_compared >= MIN_SPECIES_TO_COMPARE:
        text += " ({pct}% of median size within {rank} _{name}_)".format(
            pct=assembly.size*100/median(session, assembly, "size", rank_to_compare, name_to_compare),
            rank=rank_to_compare,
            name=name_to_compare)
    text += "\n"
    text += "*%N*: {n_pct:.2f}%\n".format(n_pct=float(assembly.num_ns)/assembly.size*100)
    text += "*%masked*: {pct:.2f}%\n".format(pct=float(assembly.num_masked)/assembly.size*100)
    if assembly.ftp_url != '':
        text += "<{fastagz}|gzipped FASTA> | <{ftp}|GenBank FTP>".format(ftp=assembly.ftp_url, fastagz=get_fasta_gz_url(assembly.ftp_url))
    message = {'text': text}
    print text
    requests.post(hook_url, json=message)

def get_masked_bases(ftp_url):
    if ftp_url == '':
        return 0
    paths = ftp_url.split('/')
    rm_out_url = ftp_url + '/' + paths[-1] + '_rm.out.gz'
    output = popenCatch("curl -s %s | gzip -d | sed '1,3d' | awk '{total += $7 - $6} END { print total }'" % rm_out_url)
    return int(output)

def parse_meta(meta):
    """Parse out the N50, etc. from the XML in the "Meta" tag."""
    soup = BeautifulSoup(meta, 'html5lib')
    return dict([(s['category'], int(s.text)) for s in soup.find_all('stat') if s['sequence_tag'] == 'all'])

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('db', help='Path to database')
    parser.add_argument('--slack-hook',
                        help='URL of slack webhook to post new releases to')
    parser.add_argument('--limit', type=int,
                        help='Only process this many new assemblies')
    return parser.parse_args()

def main():
    opts = parse_args()

    engine = create_engine('sqlite:///%s' % os.path.abspath(opts.db))
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)

    Entrez.email = 'joelcarmstrong@gmail.com'
    handle = Entrez.esearch(db='assembly', term='txid7742[Organism:exp] AND ("full genome representation"[Properties]) AND (latest[filter]))', retmax='10000')
    search_results = Entrez.read(handle)
    handle.close()
    # When this fails it's time to paginate
    assert int(search_results['Count']) < 10000

    search_results['IdList'].reverse()
    ids = map(int, search_results['IdList'])
    session = Session()
    ids_to_remove = map(itemgetter(0), session.query(Assembly.id).filter(Assembly.id.in_(ids)).all())
    session.close()
    ids = set(ids) - set(ids_to_remove)

    if opts.limit is not None:
        # Only do the first N new assemblies.
        ids = list(ids)[:opts.limit]

    for id in ids:
        handle = Entrez.efetch(db='assembly', id=id, rettype='docsum')
        record = Entrez.read(handle, validate=False)
        summary = record['DocumentSummarySet']['DocumentSummary'][0]
        stats = parse_meta(summary['Meta'])
        assembly = Assembly(id=id,
                            species_id=int(summary['SpeciesTaxid']),
                            this_assembly_is_terrible=False,
                            ftp_url=summary['FtpPath_GenBank'],
                            name = summary['AssemblyName'],
                            accession = summary['AssemblyAccession'],
                            scaffold_n50 = stats['scaffold_n50'],
                            contig_n50 = stats['contig_n50'],
                            size = stats['total_length'],
                            num_ns = stats['total_length'] - stats['ungapped_length'],
                            num_masked = get_masked_bases(summary['FtpPath_GenBank']))
        session = Session()
        if session.query(Species.id).filter(Species.id==assembly.species_id).scalar() is None:
            handle = Entrez.efetch(db='taxonomy', id=assembly.species_id, rettype='full', retmode='xml')
            record = Entrez.read(handle)[0]
            lineage = dict((i['Rank'], i['ScientificName']) for i in record['LineageEx'])
            common_name = record.get('CommonName') or record.get('OtherNames', {}).get('GenbankCommonName') or record['ScientificName']
            species = Species(id=assembly.species_id,
                              name=record['ScientificName'],
                              common_name=common_name,
                              genus=lineage['genus'],
                              family=lineage.get('family'),
                              order=lineage.get('order'),
                              class_=lineage.get('class'))
            session.add(species)
        session.add(assembly)
        session.commit()
        if opts.slack_hook is not None:
            post_to_slack(opts.slack_hook, assembly, session)
        time.sleep(1)

if __name__ == '__main__':
    main()
