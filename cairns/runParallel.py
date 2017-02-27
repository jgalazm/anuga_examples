
import os
import time
import sys

# Related major packages
import anuga

# The parallel interface
from anuga import distribute, myid, numprocs, finalize, barrier
if myid==0:
    base_scale = 100000 # 162170 triangles # 45sec fit
    #base_scale = 400000 # 42093
    default_res = 100 * base_scale   # Background resolution
    islands_res = base_scale
    cairns_res = base_scale
    shallow_res = 5 * base_scale

    bounding_polygon = anuga.read_polygon('extent.csv')
    poly_cairns = anuga.read_polygon('cairns.csv')
    poly_island0 = anuga.read_polygon('islands.csv')
    poly_island1 = anuga.read_polygon('islands1.csv')
    poly_island2 = anuga.read_polygon('islands2.csv')
    poly_island3 = anuga.read_polygon('islands3.csv')
    poly_shallow = anuga.read_polygon('shallow.csv')

    # Define list of interior regions with associated rezsolutions
    interior_regions = [[poly_cairns,  cairns_res],
                        [poly_island0, islands_res],
                        [poly_island1, islands_res],
                        [poly_island2, islands_res],
                        [poly_island3, islands_res],
                        [poly_shallow, shallow_res]]

    domain = anuga.create_domain_from_regions(bounding_polygon,
                                            boundary_tags={'top': [0],
                                                           'ocean_east': [1],
                                                           'bottom': [2],
                                                           'onshore': [3]},
                                            maximum_triangle_area=default_res,
                                            mesh_filename = 'cairns.msh',
                                            interior_regions = interior_regions,
                                            use_cache = True,
                                            verbose = True)

    # Print some stats about mesh and domain
    print 'Number of triangles = ', len(domain)
    print 'The extent is ', domain.get_extent()
    print domain.statistics()
    

    #------------------------------------------------------------------------------
    # Setup parameters of computational domain
    #------------------------------------------------------------------------------
    domain.set_name('cairns_fixed_wave') # Name of sww file
    domain.set_datadir('.')                       # Store sww output here
    domain.set_minimum_storable_height(0.01)      # Store only depth > 1cm
    domain.set_flow_algorithm('DE2')

    #------------------------------------------------------------------------------
    # Setup initial conditions
    #------------------------------------------------------------------------------
    tide = 0.0
    domain.set_quantity('stage', tide)
    domain.set_quantity('friction', 0.0)
    domain.set_quantity('elevation',
                        filename='cairns.pts',
                        use_cache = True,
                        verbose = True,
                        alpha=0.1)    
else:
    domain = None
#------------------------------------------------------------------------------
# Now produce parallel domain
#------------------------------------------------------------------------------
domain = distribute(domain,verbose=True)

domain.set_store_vertices_uniquely(False)
#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
print 'Available boundary tags', domain.get_boundary_tags()
tide = 0.0
Bd = anuga.Dirichlet_boundary([tide, 0, 0]) # Mean water level
Bs = anuga.Transmissive_stage_zero_momentum_boundary(domain) # Neutral boundary

# Huge 50m wave starting after 60 seconds and lasting 1 hour.
Bw = anuga.Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(
                    domain=domain, 
                    function=lambda t: [(60<t<3660)*10.0, 0, 0])

domain.set_boundary({'ocean_east': Bw,
                     'bottom': Bs,
                     'onshore': Bd,
                     'top': Bs})
#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
import time
t0 = time.time()

# Save every two mins leading up to wave approaching land
for t in domain.evolve(yieldstep=2*60, finaltime=5000):
    if myid == 0:
        print domain.timestepping_statistics()
        #print domain.boundary_statistics(tags='ocean_east')

# Save every 30 secs as wave starts inundating ashore
for t in domain.evolve(yieldstep=60*0.5, finaltime=10000, 
                       skip_initial_step=True):
    if myid == 0:
        print domain.timestepping_statistics()
        #print domain.boundary_statistics(tags='ocean_east')
domain.sww_merge(delete_old=True,verbose=True)
if myid == 0:
    print 'That took %.2f seconds' %(time.time()-t0)
finalize()