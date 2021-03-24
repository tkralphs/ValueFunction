'''
Created on Jan 18, 2014

@author: ted
'''

from pulp import LpVariable, lpSum, LpInteger, LpStatusOptimal, LpProblem
from pulp import LpMinimize, LpMaximize,LpStatus
import random
from pyomo.environ import *

def GenerateRandomMILP(VARIABLES, CONSTRAINTS, density = 0.2,
                       maxObjCoeff = 10, maxConsCoeff = 10,
                       tightness = 2, rand_seed = 2):
    random.seed(rand_seed)
    VARIABLES = VARIABLES
    CONSTRAINTS = CONSTRAINTS
    OBJ = dict((i, random.randint(1, maxObjCoeff)) for i in VARIABLES)
    MAT = dict(((i, j), random.randint(1, maxConsCoeff)
                        if random.random() <= density else 0)
                        for i in CONSTRAINTS for j in VARIABLES)
    RHS = dict((i, random.randint(int(len(VARIABLES)*density*maxConsCoeff/tightness),
                                  int(len(VARIABLES)*density*maxConsCoeff/1.5)))
                                  for i in CONSTRAINTS)
    
    return OBJ, MAT, RHS
    
#numVars = 40
#numIntVars = 20
#numCons = 20
numVars = 5
numIntVars = 3
numCons = 1
INTVARS = range(numIntVars)
CONVARS = range(numIntVars, numVars)
VARS = range(numVars)
CONS = range(numCons)
#CONS = ["C"+str(i) for i in range(numCons)]
intvars = LpVariable.dicts("x", INTVARS, 0, None, LpInteger)
convars = LpVariable.dicts("y", CONVARS, 0)
    
#Generate random MILP
#OBJ, MAT, RHS = GenerateRandomMILP(VARS, CONS, rand_seed = 3)
OBJ = [3, 3.5, 3, 6, 7]
MAT = {(0, 0):6, (0, 1):5, (0, 2):-4, (0, 3):2, (0, 4):-7}
RHS = [5]
LpRelaxation = LpProblem("relax", LpMinimize)
LpRelaxation += (lpSum(OBJ[j]*intvars[j] for j in INTVARS)
                     + lpSum(OBJ[j]*convars[j] for j in CONVARS)), "Objective"
for i in CONS:
    LpRelaxation += (lpSum(MAT[i, j]*intvars[j] for j in INTVARS)
                     + lpSum(MAT[i, j]*convars[j] for j in CONVARS) 
                     == RHS[i]), i
# Solve the LP relaxation
status = LpRelaxation.solve()
print(LpStatus[status])
for i in INTVARS:
    print(i, intvars[i].varValue)
for i in CONVARS:
    print(i, convars[i].varValue)
    
intvars = LpVariable.dicts("x", INTVARS, 0, 100, LpInteger)
IpRestriction = LpProblem("relax", LpMinimize)
IpRestriction += lpSum(OBJ[j]*intvars[j] for j in INTVARS), "Objective"
for i in CONS:
    IpRestriction += lpSum(MAT[i, j]*intvars[j] for j in INTVARS) == RHS[i], i
# Solve the LP relaxation
status = IpRestriction.solve()
print(LpStatus[status])
for i in INTVARS:
    print(i, intvars[i].varValue)

max_iters = 1000
listTemp = [] 
Master = AbstractModel()
Master.intIndices = Set(initialize=INTVARS)
Master.constraintSet = Set(initialize=CONS)
Master.conIndices = Set(initialize=CONVARS)
Master.intPartList = Set(initialize=listTemp)
Master.dualVarSet = Master.constraintSet * Master.intPartList
Master.theta = Var(domain=Reals, bounds = (None, None))
Master.intVars = Var(Master.intIndices, domain=NonNegativeIntegers, 
                     bounds=(0, 10))
Master.dualVars = Var(Master.dualVarSet, domain=Reals, bounds = (None, None))

def objective_rule(model):
    return model.theta
Master.objective = Objective(rule=objective_rule, sense=maximize)

def theta_constraint_rule(model, k):
    return (model.theta <=
            sum(OBJ[j]*model.int_part_list[k][j] for j in INTVARS)
            - sum(OBJ[j]*model.intVars[j] for j in INTVARS)
            + sum(MAT[(i, j)]*model.dualVars[(i, k)]*(model.intVars[j] - 
                                                    model.int_part_list[k][j]) 
                for j in INTVARS for i in CONS))
Master.theta_constraint = Constraint(Master.intPartList,
                                     rule=theta_constraint_rule)

def dual_constraint_rule(model, j, k):
    return (sum(MAT[(i, j)]*model.dualVars[(i, k)] for i in CONS) <= OBJ[j])
Master.dual_constraint = Constraint(Master.conIndices, Master.intPartList,
                                    rule=dual_constraint_rule)

Master.int_part_list = [dict((i, 0) for i in INTVARS)]

debug_print = False
opt = SolverFactory("couenne")
for i in range(max_iters):
    listTemp.append(i)
    instance = Master.create_instance()
    results = opt.solve(instance)
    instance.solutions.load_from(results)
    print('Solution in iteration', i) 
    for j in instance.intVars:
        print(j, instance.intVars[j].value)
    print('Theta:', instance.theta.value)
    if instance.theta.value < .01:
        print("Finished!")
        for int_part in Master.int_part_list:
            print('Solution:', int_part)
            print('Right-hand side: [',) 
            for k in CONS:
                print(sum(MAT[(k, l)]*int_part[l] 
                          for l in INTVARS),)
            print(']') 
        break
    if debug_print:
        for i in instance.dualVars:
            print(i, instance.dualVars[i].value)
        print(instance.dualVars[(0,0)].value)
        for k in range(len(Master.int_part_list)):
            for j in CONVARS:
                print(j, k, OBJ[j],) 
                print(sum(MAT[(i, j)]*instance.dualVars[(i, k)].value 
                          for i in CONS))
        for k in range(len(Master.int_part_list)):
            for j in INTVARS:
                print(k, sum((MAT[(i, j)])*instance.dualVars[(i, k)].value - 
                             OBJ[j] for i in CONS),)
    Master.int_part_list.append(dict((i, round(instance.intVars[i].value)) 
                                     for i in INTVARS))
