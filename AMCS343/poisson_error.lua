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
-- of covered and combined works: "Based on UG4 (https://urldefense.com/v3/__http://www.ug4.org/license__;!!Nmw4Hv0!kWu2lyWgTASUdqIMloRBnfUZ45b-2R4Zi3fdFrcJBhlQ4QdGOWIifazW9Oyfq8557KZr87Zj$ )".
-- 
-- (2) The following notice must be displayed at a prominent place in the
-- terminal output of covered works: "Based on UG4 (https://urldefense.com/v3/__http://www.ug4.org/license__;!!Nmw4Hv0!kWu2lyWgTASUdqIMloRBnfUZ45b-2R4Zi3fdFrcJBhlQ4QdGOWIifazW9Oyfq8557KZr87Zj$ )".
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
gridName = "unit_square_01_tri_2x2.ugx" -- the grid
numRefs = util.GetParamNumber("-numRefs", 4, "Number of refinements")

function ExactSol (x, y)
	return math.sin (20 * x) + math.sin (20 * y)
end

function RHS (x, y)
	return 400 * (math.sin (20 * x) + math.sin (20 * y))
end

function Dirichlet (x, y)
	return true, ExactSol (x, y)
end

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

-- set up discretization

elemDisc = ConvectionDiffusion("c", "Inner", "fv1")
elemDisc:set_diffusion(1.0)
elemDisc:set_source("RHS")

dirichletBND = DirichletBoundary()
dirichletBND:add("Dirichlet", "c", "Boundary")

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
domainDisc:add(dirichletBND)

-- set up solver (using 'util/solver_util.lua')
solverDesc = {
	type = "linear",
	precond = { type = "gmg", approxSpace = approxSpace, smoother = "ilu", baseSolver = "lu" },
	convCheck = {
		type		= "standard",
		iterations	= 100,
		absolute	= 1e-12,
		reduction	= 1e-10,
		verbose		= true
	}
}

solver = util.solver.CreateSolver(solverDesc)

-- Assemble the discretized problem
print("\nsolving...")
A = AssembledLinearOperator(domainDisc)
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)
u:set(0.0)
domainDisc:adjust_solution(u)
domainDisc:assemble_linear(A, b)

-- Solve the discretized linear problem
solver:init(A, u)
solver:apply(u, b)

-- Compute the error
local L2_error = L2Error ("ExactSol", u, "c", 0, 4)
print("L2 error on grid level "..numRefs..": "..L2_error)

-- Write the numerical solution for visualization
solFileName = "poisson_sol_GL"..numRefs
print("writing solution to '" .. solFileName .. "'...")
WriteGridFunctionToVTK(u, solFileName)

print("done")
