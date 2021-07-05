#!/usr/bin/env sage-python

# Python imports
import os
import numpy as np
from tqdm import tqdm
from argparse import ArgumentParser
from itertools import permutations, combinations
from sympy.combinatorics.permutations import Permutation
from pymongo import MongoClient
from pymongo.database import Database

# Sage imports
from sage.misc.sage_eval import sage_eval
from sage.rings.rational_field import QQ
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.rings.polynomial.multi_polynomial_ideal import MPolynomialIdeal


# Database connection and formatting

class ToricCY(Database):
    """Abstract the MongoDB connection to the database"""
    def __init__(self, uri):
        client = MongoClient(uri)
        self.__dict__ = client.ToricCY.__dict__
    
    def sample_join(self, query):
        sample = self.INVOL.aggregate([{'$match': query}, {'$sample': {'size': 1}}])
        doc = next(sample)
        doc.update(
            self.TRIANG.find_one({k: doc[k] for k in ["POLYID", "GEOMN", "TRIANGN"]})
        )
        doc.update(
            self.GEOM.find_one({k: doc[k] for k in ["POLYID", "GEOMN"]})
        )
        doc.update(
            self.POLY.find_one({"POLYID": doc["POLYID"]})
        )
        return doc
    
    def full_join(self, query):
        full = self.INVOL.find(query)
        for doc in full:
            doc.update(
                self.TRIANG.find_one({k: doc[k] for k in ["POLYID", "GEOMN", "TRIANGN"]})
            )
            doc.update(
                self.GEOM.find_one({k: doc[k] for k in ["POLYID", "GEOMN"]})
            )
            doc.update(
                self.POLY.find_one({"POLYID": doc["POLYID"]})
            )
            yield doc

    def count(self, query):
        return self.INVOL.find(query).count()


def format_matrix(mat):
    """Format database matrix objects as numpy arrays"""
    return np.array(eval(mat.replace('{', '[').replace('}', ']')))


def format_invol(invol, k):
    """Format database version of involution as sympy Permutation"""
    pairs = []
    for sub_invol in invol.strip('{}').split(','):
        x, y = sub_invol.split('->')
        x = int(x.lstrip('D')) - 1
        y = int(y.lstrip('D')) - 1
        pair = sorted((x, y))
        if pair not in pairs:
            pairs.append(pair)
    perm = Permutation(range(k)) * Permutation(pairs)
    return perm


def format_inputs(doc):
    """Reformat database records to numpy"""
    doc['NVERTS'] = format_matrix(doc['NVERTS'])
    doc['DRESVERTS'] = format_matrix(doc['DRESVERTS'])
    doc['RESCWS'] = format_matrix(doc['RESCWS'])
    doc['TRIANG'] = format_matrix(doc['TRIANG'])
    
    k, n = doc['DRESVERTS'].shape
    doc['INVOL'] = format_invol(doc['INVOL'], k)
    
    return doc


def format_srideal(prev_R, curr_R, srideal):
    """Format database version of srideal to match srideal here"""
    terms = srideal.replace('{', '').replace('}', '').replace('D', 'x').split(',')
    new_srideal = []
    for t in terms:
        new_term = []
        for tt in t.split('*'):
            eval_tt = sage_eval(tt, locals=prev_R.gens_dict())
            new_term.append(eval_tt)
        new_term = prev_R.ideal(*new_term).change_ring(curr_R)
        new_srideal.append(new_term)
    return new_srideal


def format_oplanes(prev_R, curr_R, oplanes):
    """Format database version of oplanes to match fixed_loci here"""
    new_oplanes = []
    for op in oplanes:
        oideal = []
        for oi in op['OIDEAL']:
            eval_oi = sage_eval(oi, locals=prev_R.gens_dict())
            oideal.append(eval_oi)
        oideal = prev_R.ideal(*oideal).change_ring(curr_R)
        new_oplanes.append(oideal)
    return new_oplanes


# Orientifold computations

def build_ring(k, start=0):
    """Build polynomial ring in k variables over the rational numbers"""
    R = QQ[tuple(f'x{i + start}' for i in range(k))]
    return R


def compute_monomials(R, nverts, dresverts):
    """Eq. A.8 from arXiv:1411.1418"""
    X = np.array(R.gens())
    Delta = Polyhedron(vertices=nverts)
    npoints = np.array(Delta.integral_points())
    P_expon = (npoints @ dresverts.T) + 1
    P_monoms = np.power(X, P_expon).prod(axis=1)
    return P_monoms


def compute_srideal(R, triang):
    """
    Compute the Stanley-Reisner ideal, i.e. the sets of 
    divisors that never intersect
    """
    k = R.ngens()
    n = len(triang[0])
    srideal = []
    for i in range(n):
        for c in map(set, combinations(range(k), i + 1)):
            # SR ideal generators must be minimal
            for prev_c in srideal:
                if prev_c.issubset(c):
                    break
            else:
                for t in map(set, triang):
                    # 1-cones in a simplex correspond to intersecting divisors
                    if c.issubset(t):
                        break
                else:
                    srideal.append(c)
    X = R.gens()
    return [R.ideal(*[X[i] for i in s]) for s in srideal]


def compute_fixed_loci(R, P_monoms, srideal, invol):
    """Compute the loci on the CY3 that are fixed under the involution"""
    X = np.array(R.gens())
    X_sub = dict(zip(X, X[invol.array_form]))
    # Compare each monomial instead of the full polynomial,
    # so that coefficients remain arbitrary
    P_diff = [x - x.subs(X_sub) for x in P_monoms]
    P_diff_ideal = R.ideal(*P_diff).radical()  # P_i = \sigma(P_i)
    P_monoms_ideal = R.ideal(*P_monoms).radical()  # P_i = 0
    fixed_loci = []
    for loc in P_diff_ideal.minimal_associated_primes():
        # Union of ideals == intersection of varieties
        P_loc = P_monoms_ideal + loc
        for s in srideal:
            # Quotient of ideals == set difference of varieties
            if P_loc.quotient(s).dimension() < 0:  # total ring
                break
        else:
            fixed_loci.append(loc)
    return fixed_loci


def compute_odim(ideal):
    """Compute the orientifold plane dimension of an ideal"""
    codim = ideal.dimension()
    dim = ideal.ring().krull_dimension() - codim
    return 3 + 2 * (3 - dim)


def compute_volformparity(rescws, invol):
    k, m = rescws.shape
    n = k - m
    vform = np.empty(0)
    vform_invol = np.empty(0)
    for c in combinations(range(k), n):
        term = 0
        term_invol = 0
        nc = tuple(set(range(k)) - set(c))
        for p in permutations(nc, k - n):
            p_c = Permutation(p + c)
            term += p_c.signature() * rescws[p, range(k - n)].prod()
            term_invol += (invol * p_c).signature() * rescws[invol.array_form][p, range(k - n)].prod()
        vform = np.append(vform, [term])
        vform_invol = np.append(vform_invol, [term_invol])
    if (vform == vform_invol).all():
        return 1
    elif (vform == -vform_invol).all():
        return -1
    else:
        return 0

    
def compute_smoothness(R, P_monoms, fixed_loci, srideal):
    srideal_intersect = MPolynomialIdeal.intersection(*srideal)
    J_P_monoms = [x.jacobian_ideal() for x in P_monoms]
    trans_ideal = R.ideal(*P_monoms) + sum(J_P_monoms) + sum(fixed_loci)
    trans_ideal = trans_ideal.radical().quotient(srideal_intersect)
    return trans_ideal.dimension() < 0


def compute_all(doc):
    nverts = doc['NVERTS']
    dresverts = doc['DRESVERTS']
    rescws = doc['RESCWS']
    triang = doc['TRIANG']
    invol = doc['INVOL']
    oplanes = doc['OPLANES']
    volformparity = doc['VOLFORMPARITY']
    smooth = doc['SMOOTH']
    
    k = dresverts.shape[0]
    
    prev_R = build_ring(k, start=1)
    curr_R = build_ring(k, start=0)
    
    P_monoms = compute_monomials(curr_R, nverts, dresverts)
    srideal = compute_srideal(curr_R, triang)
    
    prev_fixed_loci = format_oplanes(prev_R, curr_R, oplanes)
    curr_fixed_loci = compute_fixed_loci(curr_R, P_monoms, srideal, invol)
    fixed_loci = {'prev': prev_fixed_loci, 'curr': curr_fixed_loci}
    
    prev_oplanes = [{'ODIM': compute_odim(x), 'OIDEAL': [str(x) for x in sorted(x.gens())]} for x in prev_fixed_loci]
    curr_oplanes = [{'ODIM': compute_odim(x), 'OIDEAL': [str(x) for x in sorted(x.gens())]} for x in curr_fixed_loci]
    oplanes = {'prev': sorted(prev_oplanes, key=lambda x: x['OIDEAL']), 'curr': sorted(curr_oplanes, key=lambda x: x['OIDEAL'])}
    
    prev_volformparity = volformparity
    curr_volformparity = compute_volformparity(rescws, invol)
    volformparity = {'prev': prev_volformparity, 'curr': curr_volformparity}
    
    prev_smooth = smooth
    curr_smooth = compute_smoothness(curr_R, P_monoms, curr_fixed_loci, srideal)
    smooth = {'prev': prev_smooth, 'curr': curr_smooth}
    
    return {'H11': doc['H11'],
            'POLYID': doc['POLYID'],
            'GEOMN': doc['GEOMN'],
            'TRIANGN': doc['TRIANGN'],
            'INVOLN': doc['INVOLN'],
            'INVOL': invol.cyclic_form,
            'SRIDEAL': [[str(y) for y in x.gens()] for x in sorted(srideal)],
            'OPLANES': oplanes,
            'VOLFORMPARITY': volformparity,
            'SMOOTH': smooth}


# CLI
def parse_args():
    """Build argument parser"""
    parent_parser = ArgumentParser(description='Comparison of new method to existing db')
    parent_parser.add_argument('--uri',
                               default=os.getenv('MONGO_URI'),
                               help='Mongo database URI')
    parent_parser.add_argument('--query', '-q', default=None,
                               help=f'Query [{__file__} --query "H11=int(3),VOLFORMPARITY=int(-1)"]')
    subparsers = parent_parser.add_subparsers(dest="subcommand")
    sample_parser = subparsers.add_parser('sample')
    full_parser = subparsers.add_parser('full')
    args = parent_parser.parse_args()
    if args.query:
        args.query = dict(x.split('=') for x in args.query.split(','))
        for k in args.query:
            args.query[k] = eval(args.query[k])
    else:
        args.query = {}
    return args


if __name__ == '__main__':
    args = parse_args()
    
    db = ToricCY(args.uri)
    
    if args.subcommand == 'sample':
        doc = db.sample_join(args.query)
        doc = format_inputs(doc)
        result = compute_all(doc)

        ids = ['H11', 'POLYID', 'GEOMN', 'TRIANGN', 'INVOLN']
        print('Query:', {x: result[x] for x in ids})
        print('')
        print('Involution:', result['INVOL'])
        print('')
        print('Stanley-Reisner ideal:', result['SRIDEAL'])
        print('')
        print('Fixed loci:')
        print('\t', 'Previous:', result['OPLANES']['prev'])
        print('\t', 'Current:', result['OPLANES']['curr'])
        print('')
        print('Volume Form Parity:')
        print('\t', 'Previous:', result['VOLFORMPARITY']['prev'])
        print('\t', 'Current:', result['VOLFORMPARITY']['curr'])
        print('')
        print('Smoothness:')
        print('\t', 'Previous:', result['SMOOTH']['prev'])
        print('\t', 'Current:', result['SMOOTH']['curr'])
    elif args.subcommand == 'full':
        for doc in tqdm(db.full_join(args.query), total=db.count(args.query)):
            doc = format_inputs(doc)
            result = compute_all(doc)
            print(result)
