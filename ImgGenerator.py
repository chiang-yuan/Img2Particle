"""
Created on Fri Jul  3 11:27:35 2020
Author: HsienChun Chan, Yuan Chiang

MIT License

Copyright (c) 6th/June/2020

"""
import sys
import argparse
import glob
import numpy as np
import math
from scipy.spatial import Voronoi, voronoi_plot_2d
from matplotlib import pyplot as plt


type = "s"
n = 21

if type == "v":
    points = np.random.rand(n,2)   # 0703,0706 using n=128

elif type == "s":
    points = np.zeros(((n+1)*(n+1),2))
    for i in range(n+1):
        for j in range(n+1):
            points[(n+1)*i+j,0] = i
            points[(n+1)*i+j,1] = j

elif type == "h":
     points = np.zeros(((n+1)*(n+2),2))
     for i in range(n+1):
         for j in range(n+1):
             points[(n+1)*i+j,0] = j
     for i in range(n+1):
         for j in range(0,n+1,2):
             points[j+(n+1)*i,1] = 2*(3**0.5)/3 * i
         for j in range(0,n,2):
             points[j+(n+1)*i+1,1] =  2*(3**0.5)/3 * (i+0.5)
points = np.append(points, [[999,999], [-999,999], [999,-999], [-999,-999]], axis = 0)

vor = Voronoi(points,furthest_site=False)
fig = voronoi_plot_2d(vor,line_width=5,show_vertices=False,show_points=False)


plt.axis("square")

if type == "v":
    plt.xlim([0,1]), plt.ylim([0,1])
elif type == "s":
    plt.xlim([0,n]), plt.ylim([0,n])
elif type == "h":
    plt.xlim([0,n]), plt.ylim([0,n])


plt.gca().set_axis_off()

#plt.gca().xaxis.set_major_locator(plt.NullLocator())
#plt.gca().yaxis.set_major_locator(plt.NullLocator())

#plt.margins(0,0)
plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
plt.savefig('test.png', bbox_inches='tight', pad_inches=0)
plt.show()
'''
if type == "v":
    plt.savefig("/home/swclab/Documents/lammps/picture/Voronoi.jpg",dpi=330, bbox_inches = 'tight',pad_inches = 0)
elif type == "s":
    plt.savefig("/home/swclab/Documents/lammps/picture/Square.jpg",dpi=330, bbox_inches = 'tight',pad_inches = 0)
elif type == "h":
    plt.savefig("/home/swclab/Documents/lammps/picture/Hexagon.jpg",dpi=330, bbox_inches = 'tight',pad_inches = 0)
'''
