#!/usr/bin/env python
import os
from argparse import ArgumentParser
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine, Column, Integer, String, ForeignKey, Boolean, or_
from sqlalchemy.orm import relationship, sessionmaker
from tabulate import tabulate
from update_genbank_assembly_stats import Species, Assembly

Base = declarative_base()

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('db', help='Path to database')
    parser.add_argument('query')
    return parser.parse_args()

def main():
    opts = parse_args()
    query = opts.query

    engine = create_engine('sqlite:///%s' % os.path.abspath(opts.db))
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()
    speciess = session.query(Species).filter(or_(Species.common_name == query, Species.name == query, Species.genus == query, Species.family == query, Species.order == query, Species.class_ == query)).all()
    table = []
    for species in speciess:
        subtable = []
        for assembly in species.assemblies:
            subtable.append(['%s (%s)' % (species.name, species.common_name), assembly.name, assembly.accession, float(assembly.scaffold_n50)/1000000, float(assembly.contig_n50)/1000000, float(assembly.num_ns)/assembly.size, float(assembly.num_masked)/assembly.size, assembly.size - assembly.num_ns])
        subtable.sort(key=lambda x: x[3], reverse=True)
        table += subtable
    session.close()
    headers = ["Species", "Name", "Accession", "Scaffold N50 (Mb)", "Contig N50 (Mb)", "%N", "%msk", "Effective size"]
    print tabulate(table, headers, tablefmt="fancy_grid")

if __name__ == '__main__':
    main()
