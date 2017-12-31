import cplex
from cplex.exceptions import CplexSolverError
import numpy as np

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

