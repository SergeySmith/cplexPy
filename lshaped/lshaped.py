# -*- coding: utf-8 -*-

import numpy as np

import cplex
from cplex.exceptions import CplexSolverError


class Lshaped(object):

    def __init__(self, *args, **kwargs):
        self.master = cplex.Cplex(args[0])
        self.master.set_results_stream(None)
        self.master.solve()
        self.xnames = self.master.variables.get_names()
        self.x0 = self.master.solution.get_values()

        self.names = None
        self.Er = None
        self.er = None

        self.subproblems = args[1]
        self.probs = args[2]
        self.num_scen = args[3]
        self.xvec = []
        self.theta = None
        self.teta = [0.0]
        self.obj_value = [-cplex.infinity]
        self.Stop_Feas = False


    def _add_feasibility_cut(self):
        for k in range(self.num_scen):
            subproblem = cplex.Cplex(self.subproblems[k])
            con_names = subproblem.linear_constraints.get_names()
            # T matrix:
            Tmp = []
            for con in con_names:
                Tmp.append([subproblem.linear_constraints.get_coefficients(con, xname) for xname in self.xnames])
            T = np.array(Tmp)
            h = np.array([subproblem.linear_constraints.get_rhs(con) for con in subproblem.linear_constraints.get_names()])

            for yname in subproblem.variables.get_names(): 
                subproblem.objective.set_linear( yname, 0.0)
            ## ad v+ and v-:
            for con in subproblem.linear_constraints.get_names():
                subproblem.variables.add( obj = [1.0], lb = [0.0], ub = [cplex.infinity],
                                          columns=[cplex.SparsePair(ind = [con], val = [1.0])] )
            for con in subproblem.linear_constraints.get_names():
                subproblem.variables.add( obj = [1.0], lb = [0.0], ub = [cplex.infinity],
                                          columns=[cplex.SparsePair(ind = [con], val = [-1.0])] )

            for (idx, xname) in enumerate(self.xnames):
                subproblem.linear_constraints.add( lin_expr = [cplex.SparsePair(ind = [str(xname)], val = [1.0])], senses=["E"], rhs=[self.x0[idx]] )
            subproblem.set_results_stream(None)
            subproblem.solve()
            duals = subproblem.solution.get_dual_values(con_names)

            if subproblem.solution.get_objective_value() > 0.0:
                self.Stop_Feas = False
                Dr = np.dot(T.transpose(), duals)
                dr = np.dot(h.transpose(), duals)
                self.master.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = self.xnames, val = Dr)], senses=["G"], rhs=[dr])
                print "{} >= {}".format(Dr, dr)
                self.master.solve()
                self.x0 = self.master.solution.get_values()
                print self.x0


    def _set_optimality_cut(self):
        self.Er = [0.0] * len(self.xnames)
        self.er = 0.0
        for k in range(self.num_scen):
            subproblem = cplex.Cplex(self.subproblems[k])
            con_names = subproblem.linear_constraints.get_names()
            # T matrix:
            Tmp = []
            for con in con_names:
                Tmp.append([subproblem.linear_constraints.get_coefficients(con, xname) for xname in self.xnames])
            T = np.array(Tmp)
            h = np.array([subproblem.linear_constraints.get_rhs(con) for con in subproblem.linear_constraints.get_names()])

            for i in range(len(self.xnames)):
                subproblem.linear_constraints.add( lin_expr = [cplex.SparsePair(ind = [str(self.xnames[i])], val = [1.0])], senses=["E"], rhs=[self.x0[i]] )
            subproblem.set_results_stream(None)
            subproblem.solve()
            duals = subproblem.solution.get_dual_values(con_names)

            self.Er += self.probs[k]*np.dot(T.transpose(), duals)
            self.er += self.probs[k]*np.dot(h.transpose(), duals)


    def _check_termination(self):
        if ( (self.er - sum(self.Er[:]*self.x0[:])) <= self.theta ):
            print "[STOP]"
            # print " * Optimal solution:\t x0 = {} \t theta = {}".format(x0, master.solution.get_values()[-1])
            # print " \t\t\t Obj. value = {}".format(master.solution.get_objective_value())
            self.teta[0] = self.master.solution.get_values()[-1]
            self.obj_value[0] = self.master.solution.get_objective_value()
            self.xvec += self.x0
            return True
        return False


    def _add_optimality_cut(self):
        self.Er = np.append( self.Er, 1.0 )
        self.master.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = self.names, val = self.Er)], senses=["G"], rhs=[self.er])
        print "{} >= {}".format(self.Er, self.er)
        self.master.solve()
        self.x0 = self.master.solution.get_values()
        del self.x0[-1]


    def solve(self):
        print ">> Feasibility cuts:"
        while (self.Stop_Feas == False):
            self.Stop_Feas = True
            self._add_feasibility_cut()
        # add theta:
        self.master.variables.add( obj=[1.0], lb=[-cplex.infinity], ub=[cplex.infinity], names=['theta'])
        self.names = self.master.variables.get_names()
        print ">> Optimality cut:"
        self._set_optimality_cut()
        self._add_optimality_cut()
        self.theta = self.x0[-1]

        while True:
            Stop_Feas = False
            print ">> Feasibility cuts:"
            while (self.Stop_Feas == False):
                self.Stop_Feas = True
                self._add_feasibility_cut()
            print ">> Optimality cut:"
            self._set_optimality_cut()
            if (self._check_termination()):
                break
            self._add_optimality_cut()


    def print_solution(self):
        print " * Optimal solution:\t x0 = {} \t theta = {}".format(self.xvec, self.teta)
        print " \t\t\t Obj. value = {}".format(self.obj_value)


def lshaped(master_arg, subproblems, probs, num_scen, x, teta, obj_value):

    master = cplex.Cplex(master_arg)
    master.set_results_stream(None)
    master.solve()
    xnames = master.variables.get_names()
    x0 = master.solution.get_values()

    Stop_Feas = False
    print ">> Feasibility cuts:"
    while (Stop_Feas == False):
        Stop_Feas = True
        for k in range(num_scen):
            subproblem = cplex.Cplex(subproblems[k])
            con_names = subproblem.linear_constraints.get_names()
            # T matrix:
            Tmp = []
            for con in con_names:
                Tmp.append([subproblem.linear_constraints.get_coefficients(con, xname) for xname in xnames])
            T = np.array(Tmp)
            h = np.array([subproblem.linear_constraints.get_rhs(con) for con in subproblem.linear_constraints.get_names()])

            for yname in subproblem.variables.get_names(): 
                subproblem.objective.set_linear( yname, 0.0)
            ## ad v+ and v-:
            for con in subproblem.linear_constraints.get_names():
                subproblem.variables.add( obj = [1.0], lb = [0.0], ub = [cplex.infinity],
                                          columns=[cplex.SparsePair(ind = [con], val = [1.0])] )
            for con in subproblem.linear_constraints.get_names():
                subproblem.variables.add( obj = [1.0], lb = [0.0], ub = [cplex.infinity],
                                          columns=[cplex.SparsePair(ind = [con], val = [-1.0])] )

            for i in range(len(xnames)):
                subproblem.linear_constraints.add( lin_expr = [cplex.SparsePair(ind = [str(xnames[i])], val = [1.0])], senses=["E"], rhs=[x0[i]] )
            subproblem.set_results_stream(None)
            subproblem.solve()
            duals = subproblem.solution.get_dual_values(con_names)

            if subproblem.solution.get_objective_value() > 0.0:
                Stop_Feas = False
                Dr = np.dot(T.transpose(), duals)
                dr = np.dot(h.transpose(), duals)
                master.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = xnames, val = Dr)], senses=["G"], rhs=[dr])
                print "{} >= {}".format(Dr, dr)
                master.solve()
                x0 = master.solution.get_values()
                print x0

    print ">> Optimality cut:"
    # add theta:
    master.variables.add( obj=[1.0], lb=[-cplex.infinity], ub=[cplex.infinity], names=['theta'])
    names = master.variables.get_names()
    Er = [0.0] * len(xnames)
    er = 0.0
    for k in range(num_scen):
        subproblem = cplex.Cplex(subproblems[k])
        con_names = subproblem.linear_constraints.get_names()
        # T matrix:
        Tmp = []
        for con in con_names:
            Tmp.append([subproblem.linear_constraints.get_coefficients(con, xname) for xname in xnames])
        T = np.array(Tmp)
        h = np.array([subproblem.linear_constraints.get_rhs(con) for con in subproblem.linear_constraints.get_names()])

        for i in range(len(xnames)):
            subproblem.linear_constraints.add( lin_expr = [cplex.SparsePair(ind = [str(xnames[i])], val = [1.0])], senses=["E"], rhs=[x0[i]] )
        subproblem.set_results_stream(None)
        subproblem.solve()
        duals = subproblem.solution.get_dual_values(con_names)

        Er += probs[k]*np.dot(T.transpose(), duals)
        er += probs[k]*np.dot(h.transpose(), duals)

    Er = np.append( Er, 1.0 )
    master.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = names, val = Er)], senses=["G"], rhs=[er])
    print "{} >= {}".format(Er, er)
    master.solve()
    x0 = master.solution.get_values()
    del x0[-1]
    theta = x0[-1]

    while True:
        Stop_Feas = False
        print ">> Feasibility cuts:"
        while (Stop_Feas == False):
            Stop_Feas = True
            for k in range(num_scen):
                subproblem = cplex.Cplex(subproblems[k])
                con_names = subproblem.linear_constraints.get_names()
                # T matrix:
                Tmp = []
                for con in con_names:
                    Tmp.append([subproblem.linear_constraints.get_coefficients(con, xname) for xname in xnames])
                T = np.array(Tmp)
                h = np.array([subproblem.linear_constraints.get_rhs(con) for con in subproblem.linear_constraints.get_names()])

                for yname in subproblem.variables.get_names(): 
                    subproblem.objective.set_linear( yname, 0.0)
                ## ad v+ and v-:
                for con in subproblem.linear_constraints.get_names():
                    subproblem.variables.add( obj = [1.0], lb = [0.0], ub = [cplex.infinity],
                                              columns=[cplex.SparsePair(ind = [con], val = [1.0])] )
                for con in subproblem.linear_constraints.get_names():
                    subproblem.variables.add( obj = [1.0], lb = [0.0], ub = [cplex.infinity],
                                              columns=[cplex.SparsePair(ind = [con], val = [-1.0])] )

                for i in range(len(xnames)):
                    subproblem.linear_constraints.add( lin_expr = [cplex.SparsePair(ind = [str(xnames[i])], val = [1.0])], senses=["E"], rhs=[x0[i]] )
                subproblem.set_results_stream(None)
                subproblem.solve()
                duals = subproblem.solution.get_dual_values(con_names)

                if subproblem.solution.get_objective_value() != 0.0:
                    Stop_Feas = False
                    Dr = np.dot(T.transpose(), duals)
                    dr = np.dot(h.transpose(), duals)
                    master.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = xnames, val = Dr)], senses=["G"], rhs=[dr])
                    print "{} >= {}".format(Dr, dr)
                    master.solve()
                    x0 = master.solution.get_values()
                    del x0[-1]

        print ">> Optimality cut:"
        Er = [0.0] * len(xnames)
        er = 0.0
        for k in range(num_scen):
            subproblem = cplex.Cplex(subproblems[k])
            con_names = subproblem.linear_constraints.get_names()
            # T matrix:
            Tmp = []
            for con in con_names:
                Tmp.append([subproblem.linear_constraints.get_coefficients(con, xname) for xname in xnames])
            T = np.array(Tmp)
            h = np.array([subproblem.linear_constraints.get_rhs(con) for con in subproblem.linear_constraints.get_names()])

            for i in range(len(xnames)):
                subproblem.linear_constraints.add( lin_expr = [cplex.SparsePair(ind = [str(xnames[i])], val = [1.0])], senses=["E"], rhs=[x0[i]] )
            subproblem.set_results_stream(None)
            subproblem.solve()
            duals = subproblem.solution.get_dual_values(con_names)

        Er += probs[k]*np.dot(T.transpose(), duals)
        er += probs[k]*np.dot(h.transpose(), duals)

        if ( (er - sum(Er[:]*x0[:])) <= theta ):
            print "[STOP]"
            # print " * Optimal solution:\t x0 = {} \t theta = {}".format(x0, master.solution.get_values()[-1])
            # print " \t\t\t Obj. value = {}".format(master.solution.get_objective_value())
            teta[0] = master.solution.get_values()[-1]
            obj_value[0] = master.solution.get_objective_value()
            x += x0
            break
        Er = np.append( Er, 1.0 )
        master.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = names, val = Er)], senses=["G"], rhs=[er])
        print "{} >= {}".format(Er, er)
        master.solve()
        x0 = master.solution.get_values()
        del x0[-1]

#-----------------------#

