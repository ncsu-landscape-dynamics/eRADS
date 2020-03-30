#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ee
#ee.Authenticate()
ee.Initialize()


# In[2]:


import folium
import geemap.eefolium as emap
import subprocess
#import geemap as emap
from IPython.display import Image
import pandas as pd


# In[3]:



collection = ee.ImageCollection('USDA/NAIP/DOQQ')
aoi = ee.Geometry.Polygon([
    [-74.04,40.55],
          [-74.04,40.90],
          [-73.80,40.90],
          [-73.80,40.55]
])
centroid = aoi.centroid()
long, lat = centroid.getInfo()['coordinates']
print("long = {}, lat = {}".format(long,lat))


# In[4]:


collection = ee.ImageCollection('USDA/NAIP/DOQQ')
aoi = ee.Geometry.Polygon([
    [-73.98,40.65],
          [-73.98,40.70],
          [-74.03,40.70],
          [-74.03,40.65]
])
centroid = aoi.centroid()
long, lat = centroid.getInfo()['coordinates']
print("long = {}, lat = {}".format(long,lat))


# In[5]:


long_lat = ee.Geometry.Point(long, lat)
naip = collection.filterBounds(aoi)
naip15 = collection.filterDate('2015-05-01','2015-10-30')
imgs = naip15.mosaic().clip(aoi)
count = naip15.size().getInfo()
print('Count:', count)


# In[6]:


Map = emap.Map(center=[lat,long], zoom=14)

Map.add_basemap('SATELLITE') 
#vis = {'bands': ['N', 'R', 'G']}

vis = {'bands': ['R', 'G', 'B']}
imgs = naip15.mosaic().clip(aoi)
#Map.addLayer(aoi)
Map.addLayer(imgs,vis)
Map


# In[7]:


#nir, r = imgs.select('N'), imgs.select('R')
ndvi = imgs.normalizedDifference(["N", "R"])
ndvi_vis = {'min': -1, 'max': 1, 'palette':['red',  'yellow', 'green']}

Map.addLayer(ndvi,ndvi_vis)
Map


# In[8]:


veg_mask = ndvi.updateMask(ndvi.gte(0.1))
veg_vis = {'min': 0, 'max': 1, 'palette': ['blue']}
Map.addLayer(veg_mask,veg_vis)
Map


# In[9]:


seed = ee.Algorithms.Image.Segmentation.seedGrid(6)
#seg = ee.Algorithms.Image.Segmentation.GMeans(image=imgs,numIterations=100,pValue=50,neighborhoodSize=500)
#seg = ee.Algorithms.Image.Segmentation.SNIC(image=imgs, size=10,compactness= 0, neighborhoodSize=500,connectivity= 8, seeds=seed).select(['R_mean', 'G_mean', 'B_mean', 'N_mean', 'clusters'], ['R', 'G', 'B', 'N', 'clusters'])
seg = ee.Algorithms.Image.Segmentation.KMeans(imgs, 6, 50, 50)
clusters = seg.select('clusters')


# In[10]:


seg_vis = {'bands': ['R', 'G', 'B'], 'min':0, 'max':1, 'gamma':0.8}
Map.addLayer(clusters.randomVisualizer(), {}, 'clusters',opacity=0.5)
Map


# In[11]:


## ndvi
seg_ndvi = ndvi.addBands(clusters).reduceConnectedComponents(ee.Reducer.mean(),'clusters').rename('seg_ndvi')
Map.addLayer(seg_ndvi,{},'seg_ndvi')

## standard-deviation
std = ndvi.addBands(clusters).reduceConnectedComponents(ee.Reducer.stdDev(),'clusters').rename('std')
Map.addLayer(std,{},'StdDev')

## area
area = ee.Image.pixelArea().addBands(clusters).reduceConnectedComponents(ee.Reducer.sum(), 'clusters')
Map.addLayer(area,{}, 'Area')

## perimeter
minMax = clusters.reduceNeighborhood(ee.Reducer.minMax(), ee.Kernel.square(1))
perimeterPixels = minMax.select(0).neq(minMax.select(1)).rename('perimeter')
Map.addLayer(perimeterPixels,{},'perimeterPixels')
perimeter = perimeterPixels.addBands(clusters).reduceConnectedComponents(ee.Reducer.sum(), 'clusters')
Map.addLayer(perimeter, {}, 'Perimeter')

## width and height
sizes = ee.Image.pixelLonLat().addBands(clusters).reduceConnectedComponents(ee.Reducer.minMax(), 'clusters')
width = sizes.select('longitude_max').subtract(sizes.select('longitude_min')).rename('width')
height = sizes.select('latitude_max').subtract(sizes.select('latitude_min')).rename('height')
Map.addLayer(width, {},'Width')
Map.addLayer(height, {}, 'Height')


# In[17]:


seg_veg = clusters.updateMask(seg_ndvi.gt(0.1))
seg_veg2 = seg_veg.updateMask(area.gt(6))
seg_veg3 = seg_veg2.updateMask(area.lt(400))
Map2 = emap.Map(center=[lat,long], zoom=16)
Map2.add_basemap('SATELLITE') 
Map2.addLayer(seg_veg3.randomVisualizer(), {})
Map2


# In[18]:


vector = seg_veg3.reduceToVectors(scale=1, maxPixels=2000000000,geometryType= 'polygon',labelProperty='label')
#vector =clusters.reduceToVectors(scale=1, maxPixels=1540930650,geometryType= 'polygon',labelProperty='label')

vector= ee.FeatureCollection(vector)


# In[19]:


Map3 = emap.Map(center=[lat,long], zoom=16)
Map3.add_basemap('SATELLITE') 
Map3.addLayer(vector, {'color': 'FF00FF'},opacity=0.5)
Map3


# In[20]:


# read cords
import os
os.chdir("C:/users/liang/Desktop/host_mapping")
crds=pd.read_csv("crds.csv")
crds.head(5)


# In[21]:


crdls = crds.values.tolist()
crdls[1:10]


# In[22]:


pts = ee.List(crdls)
pts2 = ee.Geometry.MultiPoint(pts,proj=aoi.projection())
pts3 = ee.FeatureCollection(pts2)
Map2.addLayer(pts3,{'fillColor':'#4285F4','color':'#4285F4', 'width':1},opacity=0.3)
Map2


# In[26]:


join_filter = ee.Filter.withinDistance(20, '.geo', None, '.geo')
close_veg = ee.Join.simple().apply(vector,  pts3, join_filter)


# In[27]:


Map5 = emap.Map(center=[lat,long], zoom=16)
Map5.add_basemap('SATELLITE') 
Map5.addLayer(close_veg, {'color': 'FF0000'},opacity=0.9)
Map5.addLayer(pts3,{'fillColor':'#4285F4','color':'#4285F4', 'width':1},opacity=0.5)
Map5


# In[20]:


naip13 = collection.filterDate('2013-03-01','2013-12-30')
naip17 = collection.filterDate('2017-03-01','2017-12-30')
#naip19 = collection.filterDate('2019-03-01','2019-12-30')

np13 = naip13.mosaic().clip(aoi)
np17 = naip17.mosaic().clip(aoi)
#np19 = naip19.mosaic().clip(aoi)

print('Count 13:', naip13.size().getInfo())
print('Count 17:', naip17.size().getInfo())
#print('Count 19:', naip19.size().getInfo())


# In[21]:


Map3 = emap.Map(center=[lat,long], zoom=14)

Map3.add_basemap('SATELLITE') 
vis = {'bands': ['N', 'R', 'G']}
#vis = {'bands': ['R', 'G', 'B']}
#Map.addLayer(aoi)
Map3.addLayer(np13,vis)
Map3

