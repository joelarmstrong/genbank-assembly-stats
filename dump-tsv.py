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
    parser.add_argument('accessions', nargs='+')
    return parser.parse_args()

def main():
    opts = parse_args()

    engine = create_engine('sqlite:///%s' % os.path.abspath(opts.db))
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()
    print "\t".join(["Common Name", "Species Name", "Assembly Name", "Accession", "Scaffold N50", "Contig N50", "Family", "Order"])
    for accession in opts.accessions:
        assembly = session.query(Assembly).filter(Assembly.accession == accession).one()
        print "\t".join([assembly.species.common_name, assembly.species.name, assembly.name, assembly.accession, str(assembly.scaffold_n50), str(assembly.contig_n50), assembly.species.family, str(assembly.species.order)])
    session.close()

if __name__ == '__main__':
    main()
