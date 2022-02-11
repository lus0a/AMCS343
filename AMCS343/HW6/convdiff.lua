-- Convection-diffusion equation, 2d
--
-- Created for AMFS 2020
-- D. Logashenko
-- Based on Examples/conv_diff.lua

-- Load utility scripts (e.g. from from ugcore/scripts)
ug_load_script ("ug_util.lua")
ug_load_script ("util/load_balancing_util.lua")

------------------------------------------------------------------------------------------
-- Get command line parameters
------------------------------------------------------------------------------------------

numRefs 	= util.GetParamNumber("-numRefs", 5, "number of grid refinements")
numPreRefs 	= util.GetParamNumber("-numPreRefs", 0, "number of prerefinements (parallel)")
endTime		= util.GetParamNumber("-end", 0.04, "end time")
dt			= util.GetParamNumber("-dt", 0.001, "time step size")

------------------------------------------------------------------------------------------
-- Geometry data
------------------------------------------------------------------------------------------

dim 		= 2
gridName	= "unit_square_01_tri_2x2.ugx"

------------------------------------------------------------------------------------------
-- Print the parameters
------------------------------------------------------------------------------------------

print ("numRefs		= " .. numRefs)
print ("numPreRefs	= " .. numPreRefs)

------------------------------------------------------------------------------------------
-- Initialize UG4, load, refine and distribute the grid
------------------------------------------------------------------------------------------

InitUG (dim, AlgebraType("CPU", 1))

-- Create the domain, load the grid and refine it
dom = util.CreateDomain (gridName, numPreRefs, {"Inner", "Boundary"})
balancer.RefineAndRebalanceDomain (dom, numRefs - numPreRefs)

print ("Domain-info:")
print (dom:domain_info():to_string())

-- Create the vertex-centered approximation space
approxSpace = ApproximationSpace (dom)
approxSpace:add_fct ("c", "Lagrange", 1)

approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()

util.solver.defaults.approxSpace = approxSpace

------------------------------------------------------------------------------------------
-- Set up the coefficients
------------------------------------------------------------------------------------------

-- This scales the amount of diffusion of the problem
eps = 1e-1

-- The coordinates (cx, cy) specify the rotation center of the cone
cx = 0.5
cy = 0.5

-- The coordinates (ax, ay) specify the position of the highest point of the
-- cone at start time t=0.0
ax = 0.25
ay = 0.0

-- The parameter rot_vel specifies the rotation velocity
rot_vel = 100

-- The parameter delta is a scaling factor influencing the steepness of the cone
delta = 1e-1

-- This is the exact solution for our problem
function exactSolution(x, y, t)
	local xRot = math.cos(rot_vel*t) * (x-cx) - math.sin(rot_vel*t) * (y-cy) 
	local yRot = math.sin(rot_vel*t) * (x-cx) + math.cos(rot_vel*t) * (y-cy) 
	
	local expo = -((xRot - ax)*(xRot - ax) + (yRot - ay)*(yRot - ay)) / (delta + 4*eps*t)
	local scale = delta/(delta+4*eps*t)

	return scale * math.exp(expo)
end
	
-- The velocity field
function Velocity(x, y, t)
	return	rot_vel * (y - cx), rot_vel * (cy - x)
end
	
-- The dirichlet condition
function DirichletValue(x, y, t)
	return exactSolution(x, y, t)
end

------------------------------------------------------------------------------------------
-- Compose the discretization
------------------------------------------------------------------------------------------

-- upwind method
upwind = FullUpwind ()
-- upwind = PartialUpwind ()

-- local discretization of the convection-diffusion equation
elemDisc = ConvectionDiffusionFV1 ("c", "Inner")
elemDisc:set_upwind (upwind)
elemDisc:set_diffusion (eps)
elemDisc:set_velocity ("Velocity")

-- Dirichlet boundary conditions
dirichletBND = DirichletBoundary ()
dirichletBND:add ("DirichletValue", "c", "Boundary")

-- global discetization
domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
domainDisc:add(dirichletBND)

------------------------------------------------------------------------------------------
-- Set up the solver
------------------------------------------------------------------------------------------

solverDesc = {
	type = "newton",
	linSolver =
	{
--		type = "linear",
--		type = "cg",
		type = "bicgstab",
		precond = {
			type 		= "gmg",	-- preconditioner ["gmg", "ilu", "ilut", "jac", "gs", "sgs"]
			approxSpace = approxSpace,
			smoother 	= {
				type		= "jac",-- gmg-smoother ["ilu", "ilut", "jac", "gs", "sgs"]
				damping 	= 0.8
			},
			cycle		= "V",		-- gmg-cycle ["V", "F", "W"]
			preSmooth	= 2,		-- number presmoothing steps
			postSmooth 	= 2,		-- number postsmoothing steps
			baseSolver 	= "lu"
		},
		convCheck = {
			type		= "standard",
			iterations	= 100,
			absolute	= 1e-9,
			reduction	= 1e-12,
			verbose		= true
		}
	},
	convCheck = {
		type		= "standard",
		iterations	= 100,
		absolute	= 1e-9,
		reduction	= 1e-12
	}
}

solver = util.solver.CreateSolver (solverDesc)

------------------------------------------------------------------------------------------
-- Ploter for the results
------------------------------------------------------------------------------------------

out = VTKOutput()
out:clear_selection ()
out:select_nodal ("c", "C")

vtk_file_name = "scal-confdiff-lev" .. numRefs
print ("Output to file " .. vtk_file_name .. ".vtu")

vtkObserver = VTKOutputObserver (vtk_file_name, out)

------------------------------------------------------------------------------------------
-- Apply the solver
------------------------------------------------------------------------------------------

-- grid function for the solution
u = GridFunction (approxSpace)
Interpolate("exactSolution", u, "c", 0.0)

timeDisc = ThetaTimeStep (domainDisc, 1.0) -- Implicit Euler: 1.0
timeIntegrator = SimpleTimeIntegrator (timeDisc)
timeIntegrator:set_solver (solver)
timeIntegrator:attach_finalize_observer (vtkObserver)

timeIntegrator:set_time_step (dt)
timeIntegrator:apply (u, endTime, u, 0.0)

-- End of File
