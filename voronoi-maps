from shapely.geometry import MultiPoint, Point, Polygon
from scipy.spatial import Voronoi, voronoi_plot_2d
import operator
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter


#Voronoi-delaunay tesselation (see an example in Shi+2021: https://arxiv.org/abs/2102.06499)

def voronoi_polygons(vor):
    '''
    Create voronoi regions in a 2D (RA, Dec) map of galaxies.
    
    As input, the function requires: 
    vor = Voronoi (your input 2D region, in voronoi terms)
    
    The function returns the regions (indices of vertices; tuple) and vertices (coordinates 
    of vertices; tuple) of the voronoi tesselation
    '''
    voronoi_regions = []
    voronoi_vertices = vor.vertices.tolist()

    centre = vor.points.mean(axis=0)
    radius = vor.points.ptp().max()*2

    #Create map containing all ridges for a given galaxy
    all_ridges = {}
    for (point1, point2), (vrtc1, vrtc2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(point1, []).append((point2, vrtc1, vrtc2))
        all_ridges.setdefault(point2, []).append((point1, vrtc1, vrtc2))

    #Construct regions
    for point1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(vrtc >= 0 for vrtc in vertices):
            #finite region
            voronoi_regions.append(vertices)
            continue

        #Construct a non-finite region
        ridges = all_ridges[point1]
        new_vor_region = [vrtc for vrtc in vertices if vrtc >= 0]

        for point2, vrtc1, vrtc2 in ridges:
            if vrtc2 < 0:
                vrtc1, vrtc2 = vrtc2, vrtc1
            if vrtc1 >= 0:
                continue

            #Compute the missing endpoint of an infinite ridge
            tangent = vor.points[point2] - vor.points[point1] 
            tangent /= np.linalg.norm(tangent)
            normal = np.array([-tangent[1], tangent[0]])  

            midpoint = vor.points[[point1, point2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - centre, normal)) * normal
            point_far = vor.vertices[vrtc2] + direction * radius

            new_vor_region.append(len(voronoi_vertices))
            voronoi_vertices.append(point_far.tolist())

        #Sort region counterclockwise
        vrtcs = np.asarray([voronoi_vertices[v] for v in new_vor_region])
        vrtcs_mean = vrtcs.mean(axis=0)
        angle = np.arctan2(vrtcs[:,1] - vrtcs_mean[1], vrtcs[:,0] - vrtcs_mean[0])
        new_vor_region = np.array(new_vor_region)[np.argsort(angle)]

        voronoi_regions.append(new_vor_region.tolist())

    return voronoi_regions, np.asarray(voronoi_vertices)


####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################

#Example of how to use the code

#your 3D data
Zf_mxdf = np.array([...]) #redshift
RAf_mxdf = np.array([...]) #right ascension
Decf_mxdf = np.array([...]) #declination

#Redshift slices to create the 2D maps of galaxies
z_step = 0.1
z_grid = np.arange(3,6.1,z_step) #from min to max redshift of your sample in steps of 0.1

#run Voronoi function to get voronoi maps in the 2D galaxy maps
for k, zk in enumerate(z_grid):
    if k<len(z_grid)-1:
        #define your MUSE footprint edges in RA and Dec (this coordinates are for MUSE-Deep)
        footprint_edges_ra = np.array([53.1238,53.163,53.20,53.161])
        footprint_edges_dec = np.array([-27.788,-27.829,-27.787,-27.75])
        points_footprint_edges = sorted(zip(footprint_edges_ra,footprint_edges_dec), key = operator.itemgetter(0))
        
        #select galaxies within the redshift slice
        select = (Zf_mxdf>=zk) & (Zf_mxdf<z_grid[k+1]) 
        ra = RAf_mxdf[select]
        dec = Decf_mxdf[select]
        z = Zf_mxdf[select]
        
        points=sorted(zip(ra,dec), key = operator.itemgetter(0))

        # compute Voronoi tesselation
        vor = Voronoi(sorted(points))
        regions, vertices = voronoi_polygons(vor)
        pts = MultiPoint([Point(i) for i in points_footprint_edges]) 
        mask = pts.convex_hull  
        
        voronoi_vertices = []
        densities = np.array([])
        density_contrast = np.array([])
        fig = plt.figure().add_subplot(111)
        for region in regions:
            polygon = vertices[region]
            shape = list(polygon.shape)
            shape[0] += 1
            p = Polygon(np.append(polygon, polygon[0]).reshape(*shape)).intersection(mask) 
            try:
                poly = np.array(list(zip(p.boundary.coords.xy[0][:-1], p.boundary.coords.xy[1][:-1])))
                voronoi_vertices.append(poly)
                plt.fill(*zip(*poly),'None',  edgecolor = 'black', zorder=2)
            except:
                print('Warning: The redshift',zk,'crashed the Voronoi tesselation!')

            #calculate area of polygons
            area = Polygon(list(zip(poly[:,0],poly[:,1]))).area  
            densities = np.append(densities, 1/area)

        density_mean = np.mean(densities)
        density_contrast = np.append(density_contrast, densities/density_mean)
        
        #plot the voronoi tesselation along with your z sliced galaxies
        cm = plt.cm.get_cmap('jet')
        plt.scatter(np.asarray(points).T[0], np.asarray(points).T[1], marker='o',c=density_contrast, zorder=3, cmap=cm) 
        cbar=plt.colorbar()
        cbar.set_label('Density contrast', fontsize = 14)
        plt.title(str(round(zk,2))+'<z<'+str(round(z_grid[k+1],2)))
        cbar.ax.tick_params( direction='in')
        fig.xaxis.set_ticks_position('both')
        fig.yaxis.set_ticks_position('both')
        fig.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        fig.xaxis.set_tick_params(direction='in', which='both')
        fig.yaxis.set_tick_params(direction='in', which='both')
        #your survey limits (this are for MUSE-Deep)
        plt.xlim(53.20,53.12)
        plt.ylim(-27.83,-27.75)
        #plt.savefig('Voronoi-Delaunay tesselation for MXDF z_slice'+str(round(zk,1))+'.png',dpi=200)
        plt.show()
