# -*- coding: utf-8 -*-
"""

Code to plot two images as a vector plot. Takes in two images and creates a
normalized 2d vector at each point. From this we can calculate the angle of
magnetization at this point and create a hsl colour plot.

In theory this technique is extendable to a 3D vector (i.e. transverse,
longitudinal and polar images), but this is not implemented.

Dependencies:
        skimage
        matplotlib
        numpy
        json

Defined also for the scale bars is a json file with information about
common lenses used with these microscopes. If you wish to change these,
please change this file first.

@author: Robbie G. Hunt for the Stoner group @ Leeds.

To Do:
    improve quiver functionality, try to automate point choice
    add a scale bar dependent upon the lens choice
    
    big job: extract arbitrary line profiles. another bit of code maybe.

"""

from skimage import io
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import matplotlib as mpl
import json
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from skimage.restoration import (denoise_tv_chambolle, denoise_bilateral,
                                 denoise_wavelet, estimate_sigma)

def norm(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def VecMap(x, y, denoise=True,):
    '''
    Read in x and y data and turn into a vector map.

    Parameters
    ----------
    x : String
        filename (including .png) in direction of x-data.
    y : String
        filename (including .png) in direction of x-data.
    denoise : bool, optional
        Toggle use of wavelet denoising. The default is True.

    Returns
    -------
    vector : np.array
        Array containing 2D data at each point which represents
        the magnetization vector at that pixel.

    '''
    
    xdata = io.imread(x)
    ydata = io.imread(y)
    
    # Crop out square region of data. Specify as 500 pixels, makes life easy.
    x_comp = xdata[:500, :500]
    y_comp = ydata[:500, :500]
    
    ## Mid here is the midpoint of the data range. The denoising adjusts this
    ## So must be accounted for. 
    
    if denoise==True:
        x_comp = denoise_wavelet(x_comp)
        y_comp = denoise_wavelet(y_comp)
        mid = 0.5
    else:
        mid = 32768
    # Create an empty vector to overwite - xpos, ypos, 2d vector
    vector = np.zeros( (500, 500, 2) )
    
    xindex = np.arange( np.shape(x_comp)[0] )
    yindex = np.arange( np.shape(x_comp)[1] )
    
    for i in xindex:
        for j in yindex:
            
            vec = np.array([(x_comp[i,j] - mid), (y_comp[i,j] - mid)])
    #        vec = vec - mid     
            #normalise the new vector so that mx^2 + my^2 = 1
            vec = norm(vec)
        
            vector[i,j] = vec
            
    return vector
    
def angles(vectors):
    '''
    Parameters
    ----------
    vectors : numpy array
        3D numpy array. (x,y,2) or (x,y,3) type array where x and y are
        positions and vector information is contained as either 2d or 3d data.

    Returns
    -------
    angle : numpy array
        (x,y) array contain angular output.

    '''
    shape = np.shape(vectors)
    angle = np.zeros( (shape[0], shape[1] ))
    
    xidx, yidx = np.arange(np.shape(angle)[0]) , np.arange(np.shape(angle)[1])
    
    for i in xidx:
        for j in yidx:
            vec = vectors[i,j,:]
            
            #arctan2 specifies y then x - weird function.
            theta = np.arctan2( vec[1], vec[0] )
            
            angle[i,j] = theta
    return angle

def vec_cmap(vector, arrow=True, 
             wheel=False, sep=False, scale=True,
             brightness = 0.8,
             arrow_settings=[61, 0.006, 5.0, 4.5, 3],
             skip=15,
             lens='100',
             ):
    '''
    

    Parameters
    ----------
    vector : numpy array
        Vector data. (500,500,2) numpy 3D array by default.
    arrow : bool, optional
        Define if you want vector arrows. The default is True.
    wheel : bool, optional
        Include a color wheel indicating direction. The default is False.
    sep : bool, optional
        Decide if you want the wheel separate to the main plot. The default is False.
    brightness : float, optional
        Adjust brightness of vector map. The default is 0.8.
    arrow_settings : list, optional
        Settings for the quiver plot. In order, these are:
            scale
            width
            headlength
            headaxislength
            headwidth
        The default is [61, 0.006, 5.0, 4.5, 3].
    skip : int, optional
        control the amount of points to skip when plotting the quiver
    lens : string
        Lens used for imaging - include only if you display the scale bar.
        Pixels per micrometer to convert to scale bar. Default based on 100x lens.
        Defined in a json as magnification : bar_length (px), bar_size (um)
        Common px/um for different lenses are as follows:
            100x - 5.32 px/um
            60x - 3.16 px/um
            50x - 3.78 px/um
            20x - 1.5 px/um
        
    Returns
    -------
    f : matplotlib.figure.Figure
        Vector map output as a figure.

    '''
    ##Convert vectors into angles
    theta = angles(vector)
    
    ## Use hsv style as a basis
    cmap = cm.get_cmap('hsv')

    # Define a new colormap with reduced brightness
    cmap_bright = cmap(np.linspace(0, 1, 256))
    cmap_bright[:, :3] *= brightness

    # Create the new colormap
    bright_cmap = ListedColormap(cmap_bright)
    
    f, ax = plt.subplots(dpi=200)
    f.tight_layout()
    ax.frame_alpha=0
    ax.axis('off')
    ax.imshow(theta,
              vmin=(-np.pi),
              vmax=np.pi,
              cmap=bright_cmap,
              )
    
    #Plot local magnetization directions as arrows.
    if arrow==True:
        pq = ax.quiver( np.arange(500)[::skip],
                        np.arange(500)[::skip],
                        vector[::skip, ::skip, 0],
                        vector[::skip, ::skip, 1],
                        
                        pivot='mid',
                        scale = arrow_settings[0],
                        width=arrow_settings[1],
                        headlength=arrow_settings[2],
                        headaxislength=arrow_settings[3],
                        headwidth=arrow_settings[4],
                      )    

    ## This could be converted to a json dictionary lookup, I'm sure.
    if scale==True:

        try:
            ### Lookup lens information from json file. Do nothing if this is
            ### empty.
            info = json.load( open( "lens.json" ))
            bar_length, bar_size = info[lens]
            
            bar_x = 20
            bar_y = vector.shape[0] - 20
            ax.plot([bar_x, bar_x + bar_length], [bar_y, bar_y], color='w', linewidth=4)
            ax.text(bar_x + bar_length/2, bar_y - 15, str(bar_size) + '$\mu$m',
                    ha='center', va='bottom', color='w', fontsize=12, family='sans-serif', ) 
            
        except(KeyError):
            print("This lens is not defined.")
            
        except(FileNotFoundError):
            print("Lens information file missing - please create new json.")
            
        
    #Colour wheel of angle
    if wheel==True:
        if sep==False:
            ax2 = f.add_subplot(488, polar=True,)
            # Define colormap normalization for 0 to 2*pi
            norm = mpl.colors.Normalize(-np.pi, np.pi) 
            n = 40000  #the number of secants for the mesh
            phi = np.linspace(-np.pi,np.pi,n)   #theta values
            radius = np.linspace(.7,1,2)        #radius values change 0.6 to 0 for full circle
            rg, tg = np.meshgrid(radius, phi)      #create a radius,phi meshgrid
            c = tg                         #define color values as theta value
            ax2.grid(False)
            im = ax2.pcolormesh(phi, radius, c.T, norm=norm, cmap=bright_cmap,)  #plot the colormesh on axis with colormap
            ax2.set_yticklabels([])                   #turn of radial tick labels (yticks)
            ax2.tick_params(pad=3,labelsize=8)      #cosmetic changes to tick labels
            ax2.spines['polar'].set_visible(False)    #turn off the axis spine.
        
        
        ## I think this works but i dont remember - check
        if sep==True:
            f2= plt.figure(dpi=200)
            ax2 = f2.add_axes([0.1,0.1,0.8,0.8], projection='polar')
            # Define colormap normalization for 0 to 2*pi
            norm = mpl.colors.Normalize(-np.pi, np.pi) 
            n = 40000  #the number of secants for the mesh
            phi = np.linspace(-np.pi,np.pi,n)   #theta values
            radius = np.linspace(.7,1,2)        #radius values change 0.6 to 0 for full circle
            rg, tg = np.meshgrid(radius, phi)      #create a radius,phi meshgrid
            c = tg                         #define color values as theta value
            ax2.grid(False)
            im = ax2.pcolormesh(phi, radius, c.T, norm=norm, cmap=bright_cmap,)  #plot the colormesh on axis with colormap
            ax2.set_yticklabels([])                   #turn of radial tick labels (yticks)
            ax2.tick_params(pad=3,labelsize=8)      #cosmetic changes to tick labels
            ax2.spines['polar'].set_visible(False)    #turn off the axis spine.
            

            
    return f


def line_prof(theta, axis, px, graph=True):
    '''
    Take an x- or y- profile through the image. This will return angular data
    that may be useful.
    
    Parameters
    ----------
    theta : 2D np.array
        Angular data in a pixel-wise fashion.
    axis : 0, 1
        choose whether to scan along an x or y line.
        0 = x
        1 = y
    px : int
        integer to index within the length of theta.
    graph : bool, optional
        Choose whether to plot data. The default is True.

    Returns
    -------
    line_prof : np.array
        slice of the angular data.
    fig : matplotlib.figure.Figure
        Figure, if plotted.

    '''
    if axis == 0:
        ## x-line scan
        line_prof = np.degrees( theta[px, :])
        
    if axis == 1:
        line_prof = np.degrees ( theta[:, px])
        
    if graph==True:
        line_x = np.arange(len(line_prof))
        fig, ax_p = plt.subplots(dpi=200)
        ax_p.plot(line_x, line_prof, '--',)
        ax_p.scatter(line_x,line_prof, s=8, marker='x',)
        ax_p.set_xlabel('Pixel')
        ax_p.set_ylabel('Magnetization Angle ($^\circ$)')
        
    return line_prof, fig


          
if __name__=="__main__":
    vector = VecMap('x.png', 'y.png',)
    plot = vec_cmap(vector, lens='100',)
    plot.savefig('demo_output.png', transparent=True, dpi=200)
          
          
          
              
          
          
          
          
          
          
    