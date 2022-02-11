-- Navier-Stokes equation, 2d
-- Discretization: Vertex-centered, stabilized
--
-- Created for AMFS 2020
-- D. Logashenko

-- Load utility scripts (e.g. from from ugcore/scripts)
ug_load_script ("ug_util.lua")
ug_load_script ("util/load_balancing_util.lua")

------------------------------------------------------------------------------------------
-- Get command line parameters
------------------------------------------------------------------------------------------

-- Physical parameters
viscosity 	= util.GetParamNumber("-visc", 1e-3, "kinematic viscosity")
inflow		= util.GetParamNumber("-inflow", 0.5, "max. inflow velocity")
bStokes 	= util.HasParamOption("-Stokes", "If defined, only Stokes Eq. computed")

-- Numerical parameters
numRefs 	= util.GetParamNumber("-numRefs", 2, "number of grid refinements")
numPreRefs 	= util.GetParamNumber("-numPreRefs", 0, "number of prerefinements (parallel)")
bPecletBlend= util.HasParamOption("-PecletBlend", "If defined, Peclet Blend used")
upwind      = util.GetParam("-upwind", "lps", "Upwind type")
bPac        = util.HasParamOption("-pac", "If defined, pac upwind used")
stab        = util.GetParam("-stab", "fields", "Stabilization type (fields or flow)")
diffLength  = util.GetParam("-difflength", "cor", "Diffusion length type (raw, fivepoint or cor)")

------------------------------------------------------------------------------------------
-- Geometry data
------------------------------------------------------------------------------------------

dim 		= 2
gridName = "cylinderp.ugx"
inflow_y_0	= 0 -- y coord. of the lower corner of the inflow
inflow_y_1	= 0.41 -- y coord. of the upper corner of the inflow

-- Grid file name

------------------------------------------------------------------------------------------
-- Print the parameters
------------------------------------------------------------------------------------------

print (" Geometry: " .. gridName .. ", dim = " .. dim)
print (" Physical parameter:")
print ("    visc		= " .. viscosity)
print ("    inflow		= " .. inflow)
print ("    Stokes		= " .. tostring (bStokes))
print (" Numerical parameter:")
print ("    numRefs		= " .. numRefs)
print ("    numPreRefs	= " .. numPreRefs)
print ("    PecletBlend = " .. tostring (bPecletBlend))
print ("    upwind		= " .. upwind)
print ("    pac			= " .. tostring (bPac))
print ("    stab		= " .. stab)
print ("    difflength	= " .. diffLength)

------------------------------------------------------------------------------------------
-- Initialize UG4, load, refine and distribute the grid
------------------------------------------------------------------------------------------

InitUG (dim, AlgebraType("CPU", dim + 1))

-- Create the domain, load the grid and refine it
dom = util.CreateDomain (gridName, numPreRefs, {"Inlet", "Outlet", "Inner", "UpperWall", "LowerWall", "CylinderWall"})
balancer.RefineAndRebalanceDomain (dom, numRefs - numPreRefs)

print ("Domain-info:")
print (dom:domain_info():to_string())

-- Create the vertex-centered approximation space
approxSpace = ApproximationSpace (dom)
approxSpace:add_fct ({"u", "v", "p"}, "Lagrange", 1)

approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()

util.solver.defaults.approxSpace = approxSpace

------------------------------------------------------------------------------------------
-- Compose the discretization
------------------------------------------------------------------------------------------

-- inner space
NavierStokesDisc = NavierStokesFV1 ({"u", "v", "p"}, {"Inner"})
NavierStokesDisc:set_stokes (bStokes)
NavierStokesDisc:set_kinematic_viscosity (viscosity)
NavierStokesDisc:set_upwind (upwind)
NavierStokesDisc:set_peclet_blend (bPecletBlend)
NavierStokesDisc:set_stabilization (stab, diffLength)
NavierStokesDisc:set_pac_upwind (bPac)

-- boundary condition at the inlet
inflow_R = (inflow_y_1 - inflow_y_0) / 2
function inflowVel2d (x, y, t)
	return inflow * (inflow_y_1 - y) * (y - inflow_y_0) / (inflow_R * inflow_R), 0.0
end
InletDisc = NavierStokesInflow (NavierStokesDisc)
InletDisc:add ("inflowVel2d", "Inlet")

-- boundary condition at the outlet
OutletDisc = NavierStokesNoNormalStressOutflow (NavierStokesDisc)
OutletDisc:add ("Outlet")

-- boundary condition at the impermeable walls
WallDisc = NavierStokesWall (NavierStokesDisc)
WallDisc:add ("UpperWall,LowerWall,CylinderWall")

-- the global discretization
domainDisc = DomainDiscretization (approxSpace)
domainDisc:add (NavierStokesDisc)
domainDisc:add (InletDisc)
domainDisc:add (OutletDisc)
domainDisc:add (WallDisc)

-- create the assembled operator for the solver
assembledOp = AssembledOperator (domainDisc)

------------------------------------------------------------------------------------------
-- Set up the solver
------------------------------------------------------------------------------------------

solverDesc =
{
	type = "newton",
	linSolver =
	{
		type = "bicgstab",
		precond =
		{
			type = "gmg",
			rap = false,
			smoother =
			{
				type = "ilu",
				beta = -0.2
			},
			preSmooth = 1,
			postSmooth = 1,
			baseSolver = "lu",
			baseLevel = numPreRefs
		},
		convCheck =
		{
			type		= "standard",
			iterations	= 128,
			absolute	= 1e-12,
			reduction	= 1e-10,
			verbose		= true
		}
	},
	lineSearch =
	{
		type			= "standard",
		maxSteps		= 10,
		lambdaStart		= 1,
		lambdaReduce	= 0.5,
		acceptBest 		= true,
		checkAll		= false
	},
	convCheck =
	{
		type		= "standard",
		iterations	= 64,
		absolute	= 2.5e-12,
		reduction	= 2.5e-10,
		verbose		= true
	}
}

solver = util.solver.CreateSolver(solverDesc)

------------------------------------------------------------------------------------------
-- Apply the solver
------------------------------------------------------------------------------------------

-- grid function for the solution
u = GridFunction (approxSpace)
u:set (0)

-- Order the DoFs:
OrderLex (approxSpace,  "x")

-- initialize the solver
solver:init (assembledOp)
solver:prepare (u)

-- apply the solver
if not solver:apply (u) then
	print ("The solver failed.")
	exit ()
end

------------------------------------------------------------------------------------------
-- "Plot" the results
------------------------------------------------------------------------------------------

out = VTKOutput()
out:clear_selection ()
out:select_nodal ({"u", "v"}, "vel")
out:select_nodal ("u", "vel_u")
out:select_nodal ("v", "vel_v")
out:select_nodal ("p", "p")

vtk_file_name = "Flow-lev" .. numRefs
if bStokes then
	vtk_file_name = vtk_file_name .. "-Stokes"
end
print ("Output to file " .. vtk_file_name .. ".vtu")
out:print (vtk_file_name, u)

-- End of File
