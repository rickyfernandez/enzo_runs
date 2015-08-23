import numpy as np

#def cloud_grids(fp, grid_level, DomainLeftEdge, DomainRightEdge, TopGridDimensions, cloud_center, cloud_radius):
def cloud_grids(DomainLeftEdge, DomainRightEdge, TopGridDimensions, cloud_center, cloud_radius):

    edge_l = []; edge_r = []
    new_center = []
    dQ = []

    for i in range(3):

        # grid spacing
        dq = (DomainRightEdge[i] - DomainLeftEdge[i])/np.float(TopGridDimensions[i])
        dQ.append(dq)
        print "dq:", dq

        # create left, right edge positions
        ql = DomainLeftEdge[i] + np.arange(TopGridDimensions[i])*dq

        # cell center position 
        qc = ql + 0.5*dq

        # find what cell the cloud lives in
        idq = (np.abs(qc-cloud_center[i])).argmin()
        new_center.append(qc[idq])
        print "cloud left position:", ql[idq]

        # how many cells resolved the radius 
        radius_cells = np.int(np.ceil(cloud_radius/dq))
        print "radius resolve:", radius_cells

        edge_l.append(ql[idq - radius_cells])
        edge_r.append(ql[idq + radius_cells + 1])
        print "left position:",   ql[idq - radius_cells]
        print "right position:",  ql[idq + radius_cells]
        print 

    pl1 = [edge_l[0]-2*dQ[0], edge_l[1]-2*dQ[1], edge_l[2]-2*dQ[2]]
    pr1 = [edge_r[0]+2*dQ[0], edge_r[1]+2*dQ[1], edge_r[2]+2*dQ[2]]

    pl2 = [edge_l[0], edge_l[1], edge_l[2]]
    pr2 = [edge_r[0], edge_r[1], edge_r[2]]

#    # write out the grids - this is the first grid
#    fp.write("CloudCollisionGridLeftEdge[%d]  = %g %g %g\n"\
#    "CloudCollisionGridRightEdge[%d] = %g %g %g\n"\
#    "CloudCollisionGridLevel[%d]     = 1\n\n" %\
#    (grid_level, edge_l[0]-2*dQ[0], edge_l[1]-2*dQ[1], edge_l[2]-2*dQ[2],\
#    grid_level,  edge_r[0]+2*dQ[0], edge_r[1]+2*dQ[1], edge_r[2]+2*dQ[2],\
#    grid_level))
#
#    grid_level += 1
#
#    # second grid that lies inside the first grid
#    fp.write("CloudCollisionGridLeftEdge[%d]  = %g %g %g\n"\
#    "CloudCollisionGridRightEdge[%d] = %g %g %g\n"\
#    "CloudCollisionGridLevel[%d]     = 2\n\n" %\
#    (grid_level, edge_l[0], edge_l[1], edge_l[2],\
#    grid_level,  edge_r[0], edge_r[1], edge_r[2],\
#    grid_level))

    return new_center, pl1, pr1, pl2, pr2

# cloud collision problem
ProblemType = 213

# domain properties 
DomainLeftEdge = [0.0, 0.0, 0.0]
DomainRightEdge = [150, 150, 150]
LeftFaceBoundaryCondition = [3, 3, 3]    # periodic domain
RightFaceBoundaryCondition = [3, 3, 3]
TopGridDimensions = [128, 128, 128]
TopGridGravityBoundary = 1               # Isolated gravity
TopGridRank = 3

# gravitational properties
SelfGravity = 1
# spacing of 15pc between spheres
# problem specific parameters
CloudCollisionCloud1Filename = "be_sphere"
CloudCollisionCloud2Filename = "be_sphere"
CloudCollisionCloud1Temperature	= 6000               # [K]
CloudCollisionCloud2Temperature	= 6000               # [K]
CloudCollisionBackgroundDensity = 0.0027             # [density in Msun/pc^3]
CloudCollisionBackgroundTemperature = 1E6            # [K]
CloudCollisionCloud1Position = [75.0,  75.0, 75.0]   # [pc]
CloudCollisionCloud2Position = [1000.0, 1000.0, 1000.0] # [pc]
CloudCollisionNumberOfInitialGrids = 5               # must be set to 5: top grid puls two nested for each cloud
cloud1_radius = 48.3  # in pc
cloud2_radius = 48.3  

# turbulence settings
CloudCollisionAddTurbulence = 1       # 1 - yes, add turbulence
CloudCollisionCloud1K1 = 9            # turbulence k-mode range
CloudCollisionCloud1K2 = 19           # for cloud 1
CloudCollisionCloud2K1 = 9            #
CloudCollisionCloud2K2 = 19           # for cloud 2
CloudCollisionDK = 1                  # k-mode step size 
CloudCollisionCloudMachNumber = 1     # mach number: amplitude of turbulence 

# build initial nested grids for each cloud
# need to add code here

# output properties
StopTime = 5.0
dtDataDump = 0.25
DataDumpName = "CC_singletrans128amr2_"

# units
LengthUnits = 3.08568E18   # [pc in cm]
MassUnits = 1.98892E33     # [Msolar in g]
TimeUnits = 3.08568E13     # [velocity in km/s ~ 1Myr]

# hydro properties
Gamma = 1.67
Mu = 2.0
CourantSafetyNumber = 0.3
RadiativeCooling = 1
MultiSpecies = 0
FluxCorrection = 1
HydroMethod = 2
UseMinimumPressureSupport = 2
MinimumPressureSupportParameter = 9.7344E-4 # 4*Jeans length squared in code units
PhotoelectricHeating = 0
PhotoelectricHeatingRate = 6.07E-26

# amr
MinimumEfficiency = 0.4
StaticHierarchy = 0
MaximumRefinementLevel = 2
RefineBy = 2
CellFlaggingMethod = 2
MinimumMassForRefinement = 0.00005

# constants
#Msun = 1.989e33;
#Grav = 6.673e-8;
#mh   = 1.673e-24;
#pc   = 3.086e18; kpc  = 3.086e21; Mpc  = 3.086e24;
#yr   = 3.1557e7; kyr  = 1e3*yr;   Myr  = 1e6*yr;
#kms  = 1e5;
#k_b  = 1.38e-16;

cloud1_c, cloud1_pl1, cloud1_pr1, cloud1_pl2, cloud1_pr2 = cloud_grids(DomainLeftEdge, DomainRightEdge, TopGridDimensions, CloudCollisionCloud1Position, cloud1_radius)
#cloud2_c, cloud2_pl1, cloud2_pr1, cloud2_pl2, cloud2_pr2 = cloud_grids(DomainLeftEdge, DomainRightEdge, TopGridDimensions, CloudCollisionCloud2Position, cloud2_radius)
cloud2_c = cloud2_pl1 = cloud2_pr1 = cloud2_pl2 = cloud2_pr2 = np.zeros(3)

ParamFile = "Test"

# generate parameter file
fptr = open(ParamFile, "w");
#fptr.write("# HALO COLLAPSE TEST\n" \
#           "# -- M_vir = %g Msun\n" \
#           "# -- z     = %f\n" \
#           "# -- r_vir = %g kpc\n" \
#           "# -- T_vir = %g K\n" \
#           "# -- V_c   = %g km/s\n" \
#           "#\n" \
#           "# -- f_b   = %f\n" \
#           "# -- spin  = %f\n" \
#           "# -- turb  = %f\n" \
#           "# -- f_star = %f\n" \
#           "#\n" % \
#           (mvir*mass_correction, redshift, rvir*rvir_correction/1e3, \
#            tvir*tvir_correction, vc*vc_correction, BaryonFraction, \
#            KeplerianFraction, OriginalTurbulentMach, SFEfficiency));

fptr.write("# cloud collision problem\n"\
           "ProblemType                = %d\n\n" %\
           (ProblemType))

fptr.write("# domain properties\n"\
           "TopGridRank                = %d\n"\
           "TopGridGravityBoundary     = %d\n"\
           "TopGridDimensions          = %d %d %d\n"\
           "LeftFaceBoundaryCondition  = %d %d %d\n"\
           "RightFaceBoundaryCondition = %d %d %d\n"\
           "DomainLeftEdge             = %d %d %d\n"\
           "DomainRightEdge            = %d %d %d\n\n" %\
           (TopGridRank, TopGridGravityBoundary, TopGridDimensions[0], TopGridDimensions[1], TopGridDimensions[2],\
           LeftFaceBoundaryCondition[0], LeftFaceBoundaryCondition[1], LeftFaceBoundaryCondition[2],\
           RightFaceBoundaryCondition[0], RightFaceBoundaryCondition[1], RightFaceBoundaryCondition[2],\
           DomainLeftEdge[0], DomainLeftEdge[1], DomainLeftEdge[2],\
           DomainRightEdge[0], DomainRightEdge[1], DomainRightEdge[2]))

fptr.write("# gravitational properties\n"\
           "SelfGravity = %d\n\n" %\
           (SelfGravity))

fptr.write("# problem specific parameters\n"\
           "CloudCollisionCloud1Filename        = %s\n"\
           "CloudCollisionCloud2Filename        = %s\n"\
           "CloudCollisionCloud1Temperature	    = %g\n"\
           "CloudCollisionCloud2Temperature	    = %g\n"\
           "CloudCollisionBackgroundDensity     = %g\n"\
           "CloudCollisionBackgroundTemperature = %g\n"\
           "CloudCollisionCloud1Position        = %g %g %g\n"\
           "CloudCollisionCloud2Position        = %g %g %g\n\n" %\
           (CloudCollisionCloud1Filename, CloudCollisionCloud2Filename,\
           CloudCollisionCloud1Temperature, CloudCollisionCloud2Temperature,\
           CloudCollisionBackgroundDensity, CloudCollisionBackgroundTemperature,
           cloud1_c[0], cloud1_c[1], cloud1_c[2],\
           cloud2_c[0], cloud2_c[1], cloud2_c[2]))

fptr.write("# turbulence settings\n"\
           "CloudCollisionAddTurbulence   = %d\n"\
           "CloudCollisionCloud1K1        = %d\n"\
           "CloudCollisionCloud1K2        = %d\n"\
           "CloudCollisionCloud2K1        = %d\n"\
           "CloudCollisionCloud2K2        = %d\n"\
           "CloudCollisionDK              = %d\n"\
           "CloudCollisionCloudMachNumber = %d\n\n" %\
           (CloudCollisionAddTurbulence, CloudCollisionCloud1K1, CloudCollisionCloud1K2,\
           CloudCollisionCloud2K1, CloudCollisionCloud2K2, CloudCollisionDK,\
           CloudCollisionCloudMachNumber))

# build initial nested grids for each cloud
# need to add code here
# initial grids
# **grids around clouds must not overlap**
#
fptr.write("# initial nested grids\n"\
        "CloudCollisionNumberOfInitialGrids = %d\n\n" %\
        (CloudCollisionNumberOfInitialGrids))

#cloud_grids(fptr, 1, DomainLeftEdge, DomainRightEdge, TopGridDimensions, CloudCollisionCloud1Position, cloud1_radius)
#cloud_grids(fptr, 3, DomainLeftEdge, DomainRightEdge, TopGridDimensions, CloudCollisionCloud2Position, cloud2_radius)

# two nested grids for cloud 1
fptr.write("CloudCollisionGridLeftEdge[%d]  = %g %g %g\n"\
           "CloudCollisionGridRightEdge[%d] = %g %g %g\n"\
           "CloudCollisionGridLevel[%d]     = 1\n\n" %\
           (1, cloud1_pl1[0], cloud1_pl1[2], cloud1_pl1[2], 1,  cloud1_pr1[0], cloud1_pr1[2], cloud1_pr1[2], 1))

fptr.write("CloudCollisionGridLeftEdge[%d]  = %g %g %g\n"\
           "CloudCollisionGridRightEdge[%d] = %g %g %g\n"\
           "CloudCollisionGridLevel[%d]     = 2\n\n" %\
           (2, cloud1_pl2[0], cloud1_pl2[2], cloud1_pl2[2], 2,  cloud1_pr2[0], cloud1_pr2[2], cloud1_pr2[2], 2))

# two nested grids for cloud 2
fptr.write("CloudCollisionGridLeftEdge[%d]  = %g %g %g\n"\
           "CloudCollisionGridRightEdge[%d] = %g %g %g\n"\
           "CloudCollisionGridLevel[%d]     = 1\n\n" %\
           (3, cloud2_pl1[0], cloud2_pl1[2], cloud2_pl1[2], 3,  cloud2_pr1[0], cloud2_pr1[2], cloud2_pr1[2], 3))

fptr.write("CloudCollisionGridLeftEdge[%d]  = %g %g %g\n"\
           "CloudCollisionGridRightEdge[%d] = %g %g %g\n"\
           "CloudCollisionGridLevel[%d]     = 2\n\n" %\
           (4, cloud2_pl2[0], cloud2_pl2[2], cloud2_pl2[2], 4,  cloud2_pr2[0], cloud2_pr2[2], cloud2_pr2[2], 4))

fptr.write("# output properties\n"\
           "StopTime     = %d\n"\
           "dtDataDump   = %g\n"\
           "DataDumpName = %s\n\n" %\
           (StopTime, dtDataDump, DataDumpName))

fptr.write("# units\n"\
           "LengthUnits = %g\n"\
           "MassUnits   = %g\n"\
           "TimeUnits   = %g\n\n" %\
           (LengthUnits, MassUnits, TimeUnits))

fptr.write("# hydro properties\n"\
           "Gamma                           = %g\n"\
           "Mu                              = %g\n"\
           "CourantSafetyNumber             = %g\n"\
           "RadiativeCooling                = %g\n"\
           "MultiSpecies                    = %d\n"\
           "FluxCorrection                  = %d\n"\
           "HydroMethod                     = %d\n"\
           "UseMinimumPressureSupport       = %d\n"\
           "MinimumPressureSupportParameter = %g\n"\
           "PhotoelectricHeating            = %d\n"\
           "PhotoelectricHeatingRate        = %g\n\n" %\
           (Gamma, Mu, CourantSafetyNumber, RadiativeCooling, MultiSpecies, FluxCorrection,\
           HydroMethod, UseMinimumPressureSupport, MinimumPressureSupportParameter,\
           PhotoelectricHeating, PhotoelectricHeatingRate))

fptr.write("# amr\n"\
           "MinimumEfficiency        = %g\n"\
           "StaticHierarchy          = %d\n"\
           "MaximumRefinementLevel   = %d\n"\
           "RefineBy                 = %d\n"\
           "CellFlaggingMethod       = %d\n"\
           "MinimumMassForRefinement = %g\n" %\
           (MinimumEfficiency, StaticHierarchy, MaximumRefinementLevel,\
           RefineBy, CellFlaggingMethod, MinimumMassForRefinement))
fptr.close()

print "*** Wrote parameter file %s" % (ParamFile);
