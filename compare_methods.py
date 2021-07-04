#!/usr/bin/env sage-python

# Python imports
import numpy as np
from argparse import ArgumentParser
from itertools import permutations, combinations
from sympy.combinatorics.permutations import Permutation
from pymongo import MongoClient
from pymongo.database import Database

# Sage imports
from sage.rings.rational_field import QQ
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.rings.polynomial.multi_polynomial_ideal import MPolynomialIdeal


# Database connection and formatting

class ToricCY(Database):
    """Abstract the MongoDB connection to the database"""
    URI = 'mongodb://frontend:password@129.10.135.232:27017/ToricCY'
    def __init__(self, uri):
        client = MongoClient(uri)
        self.__dict__ = client.ToricCY.__dict__
    
    def find_one_join(self, query):
#         result = self.INVOL.find_one(query)
        result = list(self.INVOL.aggregate([{'$match': query}, {'$sample': {'size': 1}}]))[0]
        result.update(
            self.TRIANG.find_one({k: result[k] for k in ["POLYID", "GEOMN", "TRIANGN"]})
        )
        result.update(
            self.GEOM.find_one({k: result[k] for k in ["POLYID", "GEOMN"]})
        )
        result.update(
            self.POLY.find_one({"POLYID": result["POLYID"]})
        )
        return result


def format_matrix(mat):
    """Format database objects as matrices"""
    return np.array(eval(mat.replace('{', '[').replace('}', ']')))


def format_invol(invol, k):
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


# Orientifold computations

def build_ring(k):
    """Build polynomial ring in k variables over the rational numbers"""
    R = QQ[tuple(f'x{i + 1}' for i in range(k))]
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


def compute_dim(I):
    """Compute the orientifold plane dimension of an ideal"""
    codim = I.dimension()
    dim = I.ring().krull_dimension() - codim
    return 3 + 2 * (3 - dim)


def compute_fixed_loci(R, P_monoms, srideal, invol):
    """Compute the loci on the CY3 that are fixed under the involution"""
    X = np.array(R.gens())
    X_sub = dict(zip(X, X[invol.array_form]))
    # Compare each monomial, so that coefficients remain arbitrary
    P_diff = np.vectorize(lambda x: x - x.subs(X_sub))(P_monoms)
    fixed_loci = R.ideal(*P_diff).minimal_associated_primes()
    # Remove points in the SR ideal from the fixed loci
    # The intersection of ideals == union of varieties
    # This ideal sweeps out all points that must be removed from the ambient space
    srideal_intersection = MPolynomialIdeal.intersection(*srideal)
    fixed_loci_reduced = []
    dimensions = []
    for locus_ideal in fixed_loci:
        locus_ideal = locus_ideal.quotient(srideal_intersection)
        if locus_ideal.dimension() > -1: # Quotient is trivial if dimension = -1
            fixed_loci_reduced.append(locus_ideal.gens())
            dimensions.append(compute_dim(locus_ideal))
    return dimensions, fixed_loci_reduced


# CLI
def parse_args():
    parser = ArgumentParser(description='Comparison of new method to existing db')
    parser.add_argument('--uri',
                        default=ToricCY.URI,
                        help='Mongo database URI')
    parser.add_argument('--h11', type=int, default=None, help='H11')
    parser.add_argument('--polyid', type=int, default=None, help='POLYID')
    parser.add_argument('--geomn', type=int, default=None, help='GEOMN')
    parser.add_argument('--triangn', type=int, default=None, help='TRIANGN')
    parser.add_argument('--involn', type=int, default=None, help='INVOLN')
    args = parser.parse_args()
    args.query = {}
    if args.h11:
        args.query['H11'] = args.h11
    if args.polyid:
        args.query['POLYID'] = args.polyid
    if args.geomn:
        args.query['GEOMN'] = args.geomn
    if args.triangn:
        args.query['TRIANGN'] = args.triangn
    if args.involn:
        args.query['INVOLN'] = args.involn
    return args


def find_format_inputs(query, uri=ToricCY.URI):
    db = ToricCY(uri)
    result = db.find_one_join(query)
    
    result['NVERTS'] = format_matrix(result['NVERTS'])
    result['DRESVERTS'] = format_matrix(result['DRESVERTS'])
    result['RESCWS'] = format_matrix(result['RESCWS'])
    result['TRIANG'] = format_matrix(result['TRIANG'])
    
    k, n = result['DRESVERTS'].shape
    result['INVOL'] = format_invol(result['INVOL'], k)
    
    return result


if __name__ == '__main__':
    args = parse_args()
    result = find_format_inputs(args.query, args.uri)
    
    k = result['DRESVERTS'].shape[0]
    R = build_ring(k)
    P_monoms = compute_monomials(R, result['NVERTS'], result['DRESVERTS'])
    srideal = compute_srideal(R, result['TRIANG'])
    this_oplanes = list(zip(*compute_fixed_loci(R, P_monoms, srideal, result['INVOL'])))
    db_oplanes = [(x['ODIM'], x['OIDEAL']) for x in result['OPLANES']]
    
    ids = ['H11', 'POLYID', 'GEOMN', 'TRIANGN', 'INVOLN']
    print('Query:', {x: result[x] for x in ids})
    print('')
    print('Involution:', result['INVOL'])
    print('')
    print('Stanley-Reisner ideal:')
    print('\t', 'DB:', result['SRIDEAL'])
    print('\t', 'This:', [x.gens() for x in srideal])
    print('')
    print('Fixed loci:')
    print('\t', 'DB:', db_oplanes)
    print('\t', 'This:', this_oplanes)
