#!/usr/bin/env python
import sys
from argparse import ArgumentParser
from Bio import Phylo
from StringIO import StringIO

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('tree', help='nexus tree file')
    parser.add_argument('species', nargs='+')
    return parser.parse_args()

def main():
    opts = parse_args()
    tree = Phylo.read(opts.tree, 'nexus')
    assert len(opts.species) >= 2
    arb = tree.find_any(name=opts.species[0])
    for species in opts.species:
        node = tree.find_any(name=species)
        if node is None:
            raise RuntimeError("Species %s not found in tree" % species)
        anc = tree.common_ancestor(node, arb)
        for n in anc.get_path(node):
            n.comment='[&!color=#ff3333]'
        for n in anc.get_path(arb):
            n.comment='[&!color=#ff3333]'

    s = StringIO()
    Phylo.write(tree, s, 'nexus')
    print s.getvalue().replace('\\[', '').replace('\\]', '').replace("\\'", '')

if __name__ == '__main__':
    main()
