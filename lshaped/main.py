#!/usr/bin/python

import numpy as np

import cplex
from cplex.exceptions import CplexSolverError

from lshaped import *

num_scen = 0
probs = []

data = np.loadtxt('./farmer/scenario.dat',delimiter="\n")
num_scen = int(data[0])
for i in range(1,num_scen+1):
    probs.append(data[i])


master = cplex.Cplex('./farmer/master.lp')

subproblems = []
subproblems.append(cplex.Cplex('./farmer/subproblem1.lp'))
subproblems.append(cplex.Cplex('./farmer/subproblem2.lp'))
subproblems.append(cplex.Cplex('./farmer/subproblem3.lp'))
# subproblems.append(cplex.Cplex('./subproblem4.lp'))
if len(subproblems) != num_scen:
    print "PANIC"

x = []#[0.0] * master.variables.get_num()
teta = [0.0]
obj_value = [-cplex.infinity]

lshaped(master, subproblems, probs, num_scen, x, teta, obj_value)
print " * Optimal solution:\t x0 = {} \t theta = {}".format(x, teta)
print " \t\t\t Obj. value = {}".format(obj_value)

print "---------------- OOP ----------------"
lshape = Lshaped(master, subproblems, probs, num_scen)
lshape.solve()
lshape.print_solution()
