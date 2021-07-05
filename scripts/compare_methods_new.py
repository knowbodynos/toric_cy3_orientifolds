#!/usr/bin/env sage-python

# Python imports
import os
import numpy as np
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
        result = next(sample)
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
    hysurf_diff_ideal = R.ideal(*P_diff).radical()
    # The intersection of ideals == union of varieties
    # This ideal sweeps out all points that must be removed from the ambient space
    srideal_intersection = MPolynomialIdeal.intersection(*srideal).radical()
    # Quotient of ideals == set difference of varieties
    # Remove points in the SR ideal from the fixed loci
    hysurf_diff_ideal = hysurf_diff_ideal.quotient(srideal_intersection)
    fixed_loci = hysurf_diff_ideal.minimal_associated_primes()
    return fixed_loci


def compute_odim(ideal):
    """Compute the orientifold plane dimension of an ideal"""
    codim = ideal.dimension()
    dim = ideal.ring().krull_dimension() - codim
    return 3 + 2 * (3 - dim)


# CLI
def parse_args():
    """Build argument parser"""
    parser = ArgumentParser(description='Comparison of new method to existing db')
    parser.add_argument('--uri',
                        default=os.getenv('MONGO_URI'),
                        help='Mongo database URI')
    parser.add_argument('--query', '-q', default=None,
                        help=f'Query [{__file__} --query "H11=int(3),VOLFORMPARITY=int(-1)"]')
    args = parser.parse_args()
    if args.query:
        args.query = dict(x.split('=') for x in args.query.split(','))
        for k in args.query:
            args.query[k] = eval(args.query[k])
    else:
        args.query = {}
    return args


def sample_format_inputs(query, uri=os.getenv('MONGO_URI')):
    db = ToricCY(uri)
    result = db.sample_join(query)
    
    result['NVERTS'] = format_matrix(result['NVERTS'])
    result['DRESVERTS'] = format_matrix(result['DRESVERTS'])
    result['RESCWS'] = format_matrix(result['RESCWS'])
    result['TRIANG'] = format_matrix(result['TRIANG'])
    
    k, n = result['DRESVERTS'].shape
    result['INVOL'] = format_invol(result['INVOL'], k)
    
    return result


if __name__ == '__main__':
    args = parse_args()
    result = sample_format_inputs(args.query, args.uri)
    
    k = result['DRESVERTS'].shape[0]
    
    prev_R = build_ring(k, start=1)
    curr_R = build_ring(k, start=0)
    
    P_monoms = compute_monomials(curr_R, result['NVERTS'], result['DRESVERTS'])
    srideal = compute_srideal(curr_R, result['TRIANG'])
    
    prev_fixed_loci = format_oplanes(prev_R, curr_R, result['OPLANES'])
    curr_fixed_loci = compute_fixed_loci(curr_R, P_monoms, srideal, result['INVOL'])
    
    prev_oplanes = [(compute_odim(x), x.gens()) for x in prev_fixed_loci]
    curr_oplanes = [(compute_odim(x), x.gens()) for x in curr_fixed_loci]
    
    ids = ['H11', 'POLYID', 'GEOMN', 'TRIANGN', 'INVOLN']
    print('Query:', {x: result[x] for x in ids})
    print('')
    print('Involution:', result['INVOL'])
    print('')
    print('Stanley-Reisner ideal:')
    print('\t', 'Previous:', [x.gens() for x in format_srideal(prev_R, curr_R, result['SRIDEAL'])])
    print('\t', 'Current:', [x.gens() for x in srideal])
    print('')
    print('Fixed loci:')
    print('\t', 'Previous:', prev_oplanes)
    print('\t', 'Current:', curr_oplanes)
