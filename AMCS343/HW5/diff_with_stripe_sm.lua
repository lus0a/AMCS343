-- Copyright (c) 2010-2016:  G-CSC, Goethe University Frankfurt
-- Authors: D. Logashenko, S. Höllbacher
-- Based on the poisson.lua script by Andreas Vogel and Sebastian Reiter
-- 
-- This file is part of UG4.
-- 
-- UG4 is free software: you can redistribute it and/or modify it under the
-- terms of the GNU Lesser General Public License version 3 (as published by the
-- Free Software Foundation) with the following additional attribution
-- requirements (according to LGPL/GPL v3 §7):
-- 
-- (1) The following notice must be displayed in the Appropriate Legal Notices
-- of covered and combined works: "Based on UG4 (https://urldefense.com/v3/__http://www.ug4.org/license__;!!Nmw4Hv0!glBe5elr0LPrhTSs9pWW4qDzIJO4LV-gZ_fXmoJEjwWecYMjUv0bAgPcnm5a2H_7wJBIgXRI$ )".
-- 
-- (2) The following notice must be displayed at a prominent place in the
-- terminal output of covered works: "Based on UG4 (https://urldefense.com/v3/__http://www.ug4.org/license__;!!Nmw4Hv0!glBe5elr0LPrhTSs9pWW4qDzIJO4LV-gZ_fXmoJEjwWecYMjUv0bAgPcnm5a2H_7wJBIgXRI$ )".
-- 
-- (3) The following bibliography is recommended for citation and must be
-- preserved in all covered files:
-- "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
--   parallel geometric multigrid solver on hierarchically distributed grids.
--   Computing and visualization in science 16, 4 (2013), 151-164"
-- "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
--   flexible software system for simulating pde based models on high performance
--   computers. Computing and visualization in science 16, 4 (2013), 165-179"
-- 
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
-- GNU Lesser General Public License for more details.


-- Load utility scripts (e.g. from from ugcore/scripts)
ug_load_script("ug_util.lua")

-- Parameters of the problem
dim = 2 -- Dimensionality of the problem
gridName = util.GetParam("-grid", "stripe.ugx", "stripe.ugx or stripe1.ugx") -- the grid
stripeD = util.GetParamNumber("-stripeD", 0.0001, "Diffusion coeff. in the stripe")
numRefs = util.GetParamNumber("-numRefs", 6, "Number of refinements")
numSmooth = util.GetParamNumber("-numSmooth", 2, "Number of pre- and postsmoothing steps")
smootherName = util.GetParam("-smoother", "ilu", "ilu, jac or gs: the smoother")
solverName = util.GetParam("-solver", "linear", "linear or cg: the iteration")

-- initialize ug with the world dimension and the algebra type
InitUG(dim, AlgebraType("CPU", 1));

-- Load a domain without initial refinements.
dom = util.CreateDomain(gridName, 0)

-- Refine the domain (redistribution is handled internally for parallel runs)
print("refining...")
util.refinement.CreateRegularHierarchy(dom, numRefs, true)

-- set up approximation space
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:init_levels()
approxSpace:init_top_surface()

print("approximation space:")
approxSpace:print_statistic()

-- set up the discretization

elemDisc_Bulk = ConvectionDiffusion("c", "Bulk", "fv1")
elemDisc_Bulk:set_diffusion(1.0)

elemDisc_Stripe = ConvectionDiffusion("c", "Stripe", "fv1")
elemDisc_Stripe:set_diffusion(stripeD)

dirichletBND = DirichletBoundary()
dirichletBND:add(1.0, "c", "LeftBnd")
dirichletBND:add(0.0, "c", "RightBnd")

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc_Bulk)
domainDisc:add(elemDisc_Stripe)
domainDisc:add(dirichletBND)

-- set up solver (using 'util/solver_util.lua')
solverDesc = {
	type = solverName,
	precond =
		{
			type = "gmg",
			approxSpace = approxSpace,
			smoother = smootherName,
			preSmooth = numSmooth,
			postSmooth = numSmooth,
			rap = true,
			baseSolver = "lu"
		},
	convCheck = {
		type		= "standard",
		iterations	= 64,
		absolute	= 1e-12,
		reduction	= 1e-10,
		verbose		= true
	}
}

solver = util.solver.CreateSolver(solverDesc)

print ("")
print ("Grid: "..gridName)
print ("Solver: GMG-preconditioned "..solverName.." iteration")
print ("GMG smoother: "..smootherName)
print ("Number of pre- and postsmooting steps: "..numSmooth)
print ("")

-- Assemble the discretized problem
print("solving...")
A = AssembledLinearOperator(domainDisc)
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)
u:set(0.0)
domainDisc:adjust_solution(u)
domainDisc:assemble_linear(A, b)

-- Solve the discretized linear problem
wall_clock_start = GetClockS()
solver:init(A, u)
solver:apply(u, b)
wall_clock_end = GetClockS()

print ("Elapsed time of the linear solver: "..(wall_clock_end-wall_clock_start).." [s]\n")

-- Write the numerical solution for visualization
solFileName = "stripe_diffusion_GL"..numRefs
print("writing solution to '" .. solFileName .. "'...")
WriteGridFunctionToVTK(u, solFileName)

print("done")
