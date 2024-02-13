#!/usr/bin/python
#
#
# Alexis Rohou, September 2015
#
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1 as axes_grid1
import numpy as np
import sys

from matplotlib import rcParams, rc
rcParams['mathtext.default'] = 'regular'

num_bins_theta=18
num_bins_phi=90

# Read command line argument
if (len(sys.argv) != 3):
    print "Wrong number of arguments"
    sys.exit(2)
else:
    text_filename=sys.argv[1]
    figure_filename=sys.argv[2]

# Load data from the text file
theta, phi = np.loadtxt(text_filename,unpack=True)

# Make sure phi & theta are both between 0 and 360
for i in range(phi.size):
    theta[i] = theta[i] % 360.0
    if theta[i] < 0.0:
        theta[i] += 360.0
    phi[i] = phi[i] % 360.0
    if phi[i] < 0.0:
        phi[i] += 360.0


# Fold phi theta back into our range
for i in range(phi.size):
    if theta[i] > 270.0:
        theta[i] = theta[i] - 360.0
    elif theta[i] > 90.0:
        theta[i] = theta[i] - 180.0
        if phi[i] > 180.0:
            phi[i] = phi[i] - 180.0
        else:
            phi[i] = phi[i] + 180.0
    phi[i] = np.radians(phi[i])
    theta[i] = abs(theta[i])



# Compute histogram
hist, theta_edges, phi_edges = np.histogram2d(theta,phi,bins=[num_bins_theta,num_bins_phi],normed=False)

if (False):
    # For the normalisation to work nicely, let's treat the North pole differently and make sure that all cells in the first theta bin have been averaged
    average_count = np.zeros(4)
    for phi_index in range(hist.shape[1]):
        average_count[0] += hist[0,phi_index]
        average_count[1] += hist[1,phi_index]
        average_count[2] += hist[2,phi_index]
        average_count[3] += hist[3,phi_index]

    average_count[0] /= hist.shape[1]
    average_count[1] /= hist.shape[1]
    average_count[2] /= hist.shape[1]
    average_count[3] /= hist.shape[1]
    for phi_index in range(hist.shape[1]):
        hist[0,phi_index] = average_count[0]
        hist[1,phi_index] = average_count[1]
        hist[2,phi_index] = average_count[2]
        hist[3,phi_index] = average_count[3]



#print hist[0:4,0]
#print np.sum(hist,axis=1)

# Let's renormalise the histogram so that we have units of counts per steradians (the surface are of a sphere is 4*pi*R^2 steradians)
# To compute the area on the sphere bound by two latitudes and two longitudes, see:
# http://mathforum.org/library/drmath/view/63767.html
# Matlab has a built-in function for this, called areaquad: http://www.mathworks.com/help/map/ref/areaquad.html
#print hist.shape
for theta_index in range(hist.shape[0]):
    #print theta_index
    #print theta_edges[theta_index+1]
    for phi_index in range(hist.shape[1]):
        area = (np.cos(np.radians(theta_edges[theta_index]))-np.cos(np.radians(theta_edges[theta_index+1]))) * (phi_edges[phi_index+1]-phi_edges[phi_index])
        #hist[theta_index,phi_index] = hist[theta_index,phi_index] / np.sin(np.radians(theta_edges[theta_index]+0.5*(theta_edges[theta_index+1]-theta_edges[theta_index])))
        hist[theta_index,phi_index] = hist[theta_index,phi_index] / area  / 1000.0

#print np.sum(hist,axis=1)

# Remove outliers
mean = np.mean(hist)
sig  = np.std(hist)
np.clip(hist,a_min=hist.min(),a_max=np.percentile(hist,99.0),out=hist)


#print theta_edges
#print phi_edges

# Plot
gs = mpl.gridspec.GridSpec(1, 2,
                       width_ratios=[10,1],
                       )

ax1 = plt.subplot(gs[0],projection="polar",aspect=1.0)
ax2 = plt.subplot(gs[1])
ax1.axes.set_rgrids([30.0,60.0],labels=[r"$\theta$"u' = 30\N{DEGREE SIGN}',r"$\theta$"u'= 60\N{DEGREE SIGN}'],angle=90.0,color='white',horizontalalignment='center')
ax1.axes.set_thetagrids([0.0,45.0,90.0,135.0,180.0,225.0,270.0,315.0],frac=1.17,labels=[r"$\phi$"u' = 0\N{DEGREE SIGN}',r"$\phi$"u' = 45\N{DEGREE SIGN}',r"$\theta$"u' = 'r"$\phi$"u' = 90\N{DEGREE SIGN}',r"$\phi$"u' = 135\N{DEGREE SIGN}',r"$\phi$"u' = 180\N{DEGREE SIGN}',r"$\phi$"u' = 225\N{DEGREE SIGN}',r"$\phi$"u' = 270\N{DEGREE SIGN}',r"$\phi$"u' = 315\N{DEGREE SIGN}',r"$\phi$"u' = 0\N{DEGREE SIGN}'])
ax1.grid(True,which='both',linewidth=1,color='white')


# Draw circular grid lines manually
theta = np.arange(0,2*np.pi,0.01)
r = np.zeros(theta.shape)
r[:] = 30.0
ax1.plot(theta,r,color='w',linewidth=1,linestyle='--')
r[:] = 60.0
ax1.plot(theta,r,color='w',linewidth=1,linestyle='--')

# Do the color plotting
X, Y = np.meshgrid(phi_edges, theta_edges)
im = ax1.pcolormesh(X, Y, hist)
color_bar = plt.colorbar(im, cax=ax2)
color_bar.solids.set_edgecolor("face")
plt.ylabel("Projection density ($10^3/sr$)")

# Save the figure
plt.savefig(figure_filename)
#plt.show()
