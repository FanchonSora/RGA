import math, time, re, sympy, json, operator, os.path, yaml
from datetime import datetime
import numpy as np
from sympy import Matrix, symbols, linsolve
from numpy.linalg import matrix_rank
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.AtomPairs import Pairs, Torsions
from rdkit.Chem.Descriptors import NumValenceElectrons
from collections import defaultdict
from itertools import permutations
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from pulp import LpProblem, LpMinimize, LpVariable, LpInteger, LpStatus, PulpSolverError, PULP_CBC_CMD
from random import randrange, seed
from functools import partial
from multiprocessing import Pool, cpu_count, freeze_support
import traceback
from urllib.request import urlopen
from urllib.error import HTTPError, URLError
from sys import exit

def combine2Dicts(A, B):
    return {x: A.get(x, 0) + B.get(x, 0) for x in set(A).union(B)}

def multiplyCoef(L, coef):
    return {elem: L[elem] * coef for elem in L}

def compare2Dicts(L1, L2):
    return L1 == L2

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def adjustNoProcessors(no_processors):
    no_cpu = cpu_count()

    if no_processors > no_cpu or no_processors == -1:
        no_processors = no_cpu
    elif no_processors < 1:
        no_processors = 1

    return no_processors

def divide_work(f, args):
    return f(*args)

def printHeader(fout, species, executor, config, year=None):

    if year is None:
        year = datetime.now().year

    fout.write("*** ------------------------------------- ***\n")
    fout.write("*** RGA: Reaction Generator and Analyzer  ***\n")
    fout.write("*** Copyright @ Triet Le & Lam Huynh " + str(year) + " ***\n")
    fout.write("***         Version 2.0 (2026)            ***\n")
    fout.write("*** ------------------------------------- ***\n\n")
    fout.write("*** CALCULATION OF HEAT OF FORMATION at 298 K ***\n\n")
    fout.write("=========================================================================\n")

    leftWidth = 30
    fout.write("+ Species".ljust(leftWidth) + ": " + species + "\n")
    fout.write("+ Date".ljust(leftWidth) + ": " + datetime.now().strftime("%A, %Y-%m-%d, %H:%M:%S") + "\n")
    fout.write("+ Calculated by".ljust(leftWidth) + ": " + executor + "\n\n")
    fout.write("=========================================================================\n")

    fout.write("CONTRAINTS:\n")
    fout.write(yaml.dump(config.configFile, default_flow_style=False))
    fout.write("\n")
    fout.write("=========================================================================\n")

    return fout

def printExecutionDetails(fout, exe_time):

    fout.write("CALCULATION EXECUTION DETAILS:\n")
    fout.write("\t+ Execution time:".ljust(20) + str(exe_time) + " s.\n\n")
    fout.write("=========================================================================\n")
