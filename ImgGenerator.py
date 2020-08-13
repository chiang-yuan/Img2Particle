"""
Created on Fri Jul  3 11:27:35 2020
Author: HsienChun Chen, Yuan Chiang



MIT License
Copyright (c) 6th/June/2020

======================================================================
version 1.1
    Add Fibonacci style
28th/July/2020 by HsienChun Chen
======================================================================

"""
import os
import argparse
import glob
import numpy as np
import math
from scipy.spatial import Voronoi, voronoi_plot_2d
from matplotlib import pyplot as plt
from PIL import Image

def import_img(infile):
    img = Image.open(infile).convert('L')
    img.load()
    array = np.asarray(img, dtype="int32")
    return array

version = '1.1'

parser = argparse.ArgumentParser(prog='ImgGenerator',
                                 description='Generate images with different pattern')

#parser.add_argument('style', metavar='style', type=str,
#                    help='available styles: s (square), h (hexagon), r (random)')
parser.add_argument('style', metavar='style', type=str,
                    help='available styles: s (square), h (hexagon), r (random) , f (fibonacci)')
parser.add_argument('-o', dest='outfile', metavar='image file', type=str)
parser.add_argument('-lx', dest='lx', metavar='Lx', type=float, default=256,
                    help='image size in x dimension\t(default = 256)')
parser.add_argument('-ly', dest='ly', metavar='Ly',type=float, default=256,
                    help='image size in y dimension\t(default = 256)')
parser.add_argument('-n', dest='num', metavar='num', type=int, default=256,
                    help='number of particles\t(default = 256)')
parser.add_argument('-p', dest='porosity', metavar='porosity',type=float, default=0.6,
                    help='porosity\t(default = 0.6)')
parser.add_argument('--version', action='version', version='%(prog)s {:s}'.format(version))

args = parser.parse_args()

style = args.style
num = args.num

lx = args.lx
ly = args.ly

if style == 'r':
    points = np.concatenate((np.random.uniform(0,lx,size=(num,1)),np.random.uniform(0,ly,size=(num,1))),
                            axis=1)
    #points = np.random.rand(num,2)   # 0703,0706 using n=128
elif style == 's':

    s = (lx+ly)/(2*math.sqrt(num))
    unit = s*np.array([[1,0],[0,1]], dtype=float) # distant unit

    nx = int(round(lx/unit[0,0]))
    ny = int(round(ly/unit[1,1]))

    mx = 2
    my = 2

    nx = nx + 2*mx
    ny = ny + 2*my

    points = np.zeros((2*nx*ny,2))
    id_ = 0
    for i in range(nx):
        for j in range(ny):
            points[id_,:] = np.dot(np.transpose(unit),[[i-mx],[j-my]]).reshape(-1)
            id_ = id_ + 1
    points[:,0] = points[:,0] - ((np.max(points[:,0]) + np.min(points[:,0]))/2.0 - lx/2.0) - np.dot(np.transpose(unit),[[0.5],[0.5]])[0]
    points[:,1] = points[:,1] - ((np.max(points[:,1]) + np.min(points[:,1]))/2.0 - ly/2.0) - np.dot(np.transpose(unit),[[0.5],[0.5]])[1]

elif style == 'h':

    s = (lx+ly)/(2*math.sqrt(num))
    unit = s*np.array([[1,0],[0,math.sqrt(3.)]], dtype=float) # distant unit

    nx = int(round(lx/unit[0,0]))
    ny = int(round(ly/unit[1,1]))

    mx = 2
    my = 2

    nx = nx + 2*mx
    ny = ny + 2*my

    points = np.zeros((2*nx*ny,2))
    id_ = 0
    for i in range(nx):
        for j in range(ny):
            points[id_,:] = np.dot(np.transpose(unit),[[i-mx],[j-my]]).reshape(-1)
            id_ = id_ + 1
            points[id_,:] = np.dot(np.transpose(unit),[[i+0.5-mx],[j+0.5-my]]).reshape(-1)
            id_ = id_ + 1
    points[:,0] = points[:,0] - ((np.max(points[:,0]) + np.min(points[:,0]))/2.0 - lx/2.0) - np.dot(np.transpose(unit),[[0.25],[0.25]])[0]
    points[:,1] = points[:,1] - ((np.max(points[:,1]) + np.min(points[:,1]))/2.0 - ly/2.0) - np.dot(np.transpose(unit),[[0.25],[0.25]])[1]


elif style == 'f':

    golden_angle = math.pi * (3 - math.sqrt(5))

    def fibonacci(n):
        point = np.zeros((n,2))
        count = 0
        for i in range(n):
            theta = i * golden_angle
            r = np.sqrt(lx*ly) * math.sqrt(i) / math.sqrt(n)
            point[i,0]=(r * math.cos(theta)+0.5*lx)
            point[i,1]=(r * math.sin(theta)+0.5*ly)

        for i in range(n):
            if point[i,0]<=lx and point[i,0]>=0 and point[i,1]<=ly and point[i,1]>=0:
                count = count + 1
        return point , count


    f_iter = 0
    d_num = 0

    while True:

        dn = math.floor(num - d_num * num)
        if abs ((fibonacci(dn)[1]-num)//num) < 1e-2 or f_iter >= 200 :
            points = fibonacci(dn)[0]
            break
        else:
            if f_iter == 0:
                d_num = (fibonacci(dn)[1]-num)//num
            else:
                d_num = d_num + (fibonacci(dn)[1]-num)/num

            f_iter = f_iter + 1


else:
    print('Error: Unreconized Style')
    quit()



points = np.append(points, [[1000*lx,1000*ly], [-999*lx,1000*ly], [1000*lx,-999*ly], [-999*lx,-999*ly]], axis = 0)

vor = Voronoi(points,furthest_site=False)

if args.outfile is None:
    i = 0
    while True:
        outfile = '{:s}_{:3d}_{:.1f}_{:02d}.png'.format(args.style.upper(),args.num,args.porosity,i)
        if os.path.isfile(outfile):
            i = i + 1
        else:
            break
else:
    if args.outfile.lower().endswith('png'):
        outfile = args.outfile
    else:
        outfile = args.outfile + '.png'

iter = 0
width = 5
rate = 0.5

while True:
    fig, ax = plt.subplots(facecolor=(0, 0, 0))
    voronoi_plot_2d(vor, ax=ax,
                    line_width=width, line_colors='white',
                    show_vertices=False, show_points=False)
    ax.set(xlim=[0,lx], ylim=[0,ly], aspect='equal')
    ax.set_axis_off()
    plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    plt.savefig(outfile, facecolor='black', bbox_inches='tight', pad_inches=0, dpi=110)
    #plt.show()
    img = import_img(outfile)
    hist, bin_edges = np.histogram(img.reshape(-1),bins=3)
    porosity = float(hist[0])/float(np.sum(hist))
    if abs((porosity - args.porosity)/args.porosity) < 1e-2 or iter >= 50 :
        break
    else:
        if iter == 0:
            descent = ((porosity - args.porosity)/args.porosity)
        else:
            descent = 0.5*descent + 0.5*((porosity - args.porosity)/args.porosity)

        width = width + descent*rate*width

        iter = iter + 1
        print('Iter {} \t Porosity = {:f} \t Try width {:f}'.format(iter,porosity,width),end='\r')

        plt.close()

print('Iter {} \t Porosity = {:f} \t Final wall width {:f}'.format(iter,porosity,width))
