from ROOT import *
import root_numpy

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

print 1
import numpy

import sys

file=TFile(sys.argv[1])
tree=file.Get("clusterInfo/tree")
tree=root_numpy.tree2array(tree, selection="hitFound&&broken_cluster&&abs(hit_localPixel_y-cluster_center_y)<10")
 
print 1 
seq1=numpy.stack([tree["track_pt"],tree["track_eta"],tree["track_phi"]])
print seq1.shape
 
 
seq1=numpy.stack([tree["track_pt"],tree["track_eta"],tree["track_phi"]])
print seq1.shape
track=numpy.swapaxes(seq1,0,1)
print track.shape
 
seq1=numpy.stack([tree["track_global_z"],tree["track_global_phi"],tree["track_exp_sizeX"],tree["track_exp_sizeY"] ,tree["track_exp_charge"], tree["track_alpha"], tree["track_beta"]])
print seq1.shape
track1=numpy.swapaxes(seq1,0,1)
print track1.shape

Max = 41530#numpy.max(tree["cluster_chargeBroken_in_hits"])
print Max, "maxxxx"
 
seq1=numpy.stack([((tree["cluster_chargeBroken_in_hits"])/(1.0*Max)).reshape(len(tree["track_pt"]),21,7),tree["cluster_columnBroken_ON"].reshape(len(tree["track_pt"]),21,7)])
#print seq1.shape
seq1=numpy.swapaxes(seq1,0,1)
seq1=numpy.swapaxes(seq1,1,2)
seq1=numpy.swapaxes(seq1,2,3)
print seq1.shape
#image=seq1.reshape(len(seq1), 21, 7, 2)
image=seq1


seq1=numpy.stack([tree["hit_localPixel_x"]-tree["cluster_center_x"],tree["hit_localPixel_y"]-tree["cluster_center_y"]])
#seq1=numpy.stack([tree["hit_localPixel_x"]-tree["track_localPixel_x"]+3,tree["hit_localPixel_y"]-tree["track_localPixel_y"]+10])
print seq1.shape
hit_pos=numpy.swapaxes(seq1,0,1)
print hit_pos.shape

seq1=numpy.stack([tree["hit_localPixel_x"]-tree["track_localPixel_x"],tree["hit_localPixel_y"]-tree["track_localPixel_y"]])
print seq1.shape
delta_TH=numpy.swapaxes(seq1,0,1)
print delta_TH.shape


seq1=numpy.stack([-tree["cluster_center_x"]+tree["track_localPixel_x"],-tree["cluster_center_y"]+tree["track_localPixel_y"]])
print seq1.shape
track_pos=numpy.swapaxes(seq1,0,1)
print track_pos.shape


seq1=numpy.stack([tree["cluster_localPixel_x"]-tree["cluster_center_x"],tree["cluster_localPixel_y"]-tree["cluster_center_y"]])
print seq1.shape
local_cluster=numpy.swapaxes(seq1,0,1)
print local_cluster.shape

seq1=numpy.stack([tree["brokenCluster_localPixel_x"]-tree["cluster_center_x"],tree["brokenCluster_localPixel_y"]-tree["cluster_center_y"]])
print seq1.shape
local_brokencluster=numpy.swapaxes(seq1,0,1)
print local_brokencluster.shape


seq1=numpy.stack([tree["cluster_center_x"],tree["cluster_center_y"]])
print seq1.shape
theCenter=numpy.swapaxes(seq1,0,1)
print theCenter.shape

seq1=numpy.stack([((tree["cluster_charge_in_hits"])/(1.0*Max)).reshape(len(tree["track_pt"]),21,7),tree["cluster_column_ON"].reshape(len(tree["track_pt"]),21,7)])
#print seq1.shape
seq1=numpy.swapaxes(seq1,0,1)
seq1=numpy.swapaxes(seq1,1,2)
seq1=numpy.swapaxes(seq1,2,3)
print seq1.shape
#image=seq1.reshape(len(seq1), 21, 7, 2)
Truth_image=seq1

from tempfile import TemporaryFile
outfile = TemporaryFile()

numpy.savez(sys.argv[1].split(".")[0]+"broken_cluster", track=track, track1=track1, image=image, hit_pos=hit_pos, delta_TH=delta_TH, track_pos=track_pos, Truth_image=Truth_image,
local_brokencluster=local_brokencluster, local_cluster=local_cluster, theCenter=theCenter)