#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ee
#ee.Authenticate()
ee.Initialize()

import folium
import geemap.eefolium as emap
import subprocess
#import geemap as emap
from IPython.display import Image
import pandas as pd


# In[2]:


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


# In[3]:


collection = ee.ImageCollection('USDA/NAIP/DOQQ')
aoi = ee.Geometry.Polygon([
    [-73.99,40.67],
          [-73.99,40.68],
          [-74.00,40.68],
          [-74.00,40.67]
])
centroid = aoi.centroid()
long, lat = centroid.getInfo()['coordinates']
print("long = {}, lat = {}".format(long,lat))


# In[4]:


long_lat = ee.Geometry.Point(long, lat)
naip = collection.filterBounds(aoi)
naip15 = collection.filterDate('2015-05-01','2015-10-30')
np15 = naip15.mosaic().clip(aoi)
count = naip15.size().getInfo()
print('Count:', count)


# In[5]:


Map = emap.Map(center=[lat,long], zoom=14)

Map.add_basemap('SATELLITE') 
#vis = {'bands': ['N', 'R', 'G']}

vis = {'bands': ['R', 'G', 'B']}
#Map.addLayer(aoi)
Map.addLayer(np15,vis)
Map


# In[6]:


#nir, r = imgs.select('N'), imgs.select('R')
ndvi = np15.normalizedDifference(["N", "R"])
ndvi_vis = {'min': -1, 'max': 1, 'palette':['red',  'yellow', 'green']}

Map.addLayer(ndvi,ndvi_vis)
Map


# In[7]:


veg_mask = ndvi.updateMask(ndvi.gte(0.1))
veg_vis = {'min': 0, 'max': 1, 'palette': ['blue']}
Map.addLayer(veg_mask,veg_vis)
Map


# In[8]:


seed = ee.Algorithms.Image.Segmentation.seedGrid(6)
#seg = ee.Algorithms.Image.Segmentation.GMeans(image=imgs,numIterations=100,pValue=50,neighborhoodSize=500)
#seg = ee.Algorithms.Image.Segmentation.SNIC(image=imgs, size=10,compactness= 0, neighborhoodSize=500,connectivity= 8, seeds=seed).select(['R_mean', 'G_mean', 'B_mean', 'N_mean', 'clusters'], ['R', 'G', 'B', 'N', 'clusters'])
seg = ee.Algorithms.Image.Segmentation.KMeans(np15, 6, 50, 50)
clusters = seg.select('clusters')


# In[9]:


seg_vis = {'bands': ['R', 'G', 'B'], 'min':0, 'max':1, 'gamma':0.8}
Map.addLayer(clusters.randomVisualizer(), {}, 'clusters',opacity=0.5)
Map


# In[10]:


## ndvi
seg_ndvi = ndvi.addBands(clusters).reduceConnectedComponents(ee.Reducer.mean(),'clusters').rename('seg_ndvi')

## area
area = ee.Image.pixelArea().addBands(clusters).reduceConnectedComponents(ee.Reducer.sum(), 'clusters')


# In[11]:


seg_veg = clusters.updateMask(seg_ndvi.gt(0.1))
seg_veg2 = seg_veg.updateMask(area.gt(6))
seg_veg3 = seg_veg2.updateMask(area.lt(400))
Map2 = emap.Map(center=[lat,long], zoom=16)
Map2.add_basemap('SATELLITE') 
Map2.addLayer(seg_veg3.randomVisualizer(), {})
Map2


# In[12]:


vector = seg_veg3.reduceToVectors(scale=1, maxPixels=2000000000,geometryType= 'polygon',labelProperty='label')
#vector =clusters.reduceToVectors(scale=1, maxPixels=1540930650,geometryType= 'polygon',labelProperty='label')

vector= ee.FeatureCollection(vector)


# In[13]:


Map3 = emap.Map(center=[lat,long], zoom=16)
Map3.add_basemap('SATELLITE') 
Map3.addLayer(vector, {'color': 'FF00FF'},opacity=0.5)
Map3


# In[14]:


# read cords
import os
os.chdir("C:/users/liang/Desktop/host_mapping")
crds=pd.read_csv("crds.csv")
crds.head(5)


# In[15]:


crdls = crds.values.tolist()
crdls[1:10]


# In[16]:


pts = ee.List(crdls)
pts2 = ee.Geometry.MultiPoint(pts,proj=aoi.projection())
pts3 = ee.FeatureCollection(pts2)
Map2.addLayer(pts3,{'fillColor':'#4285F4','color':'#4285F4', 'width':1},opacity=0.3)
Map2


# In[17]:


join_filter = ee.Filter.withinDistance(10, '.geo', None, '.geo')
close_veg = ee.Join.simple().apply(vector,  pts3, join_filter)


# In[18]:


Map5 = emap.Map(center=[lat,long], zoom=14)
Map5.add_basemap('SATELLITE') 
Map5.addLayer(close_veg, {'color': 'FF0000'},opacity=0.9)
Map5.addLayer(pts3,{'fillColor':'#4285F4','color':'#4285F4', 'width':1},opacity=0.5)
Map5


# In[19]:


seg_veg4 = seg_veg3.clipToCollection(close_veg)
Map5 = emap.Map(center=[lat,long], zoom=16)

Map5.add_basemap('SATELLITE') 
Map5.addLayer(seg_veg4.randomVisualizer(), {})
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


#nir, r = imgs.select('N'), imgs.select('R')
ndvi15 = np15.normalizedDifference(["N", "R"])
ndvi13 = np13.normalizedDifference(["N", "R"])
ndvi17 = np17.normalizedDifference(["N", "R"])
fea_imgs = ee.Image([np13, np15, np17,ndvi13,ndvi15,ndvi17])
fea_imgs2 = ee.Image([np13.select("N"),np15.select("N"),np17.select("N"),ndvi13,ndvi15, ndvi17])


# In[22]:


cluster2=clusters.clipToCollection(close_veg)

## mean of all bands and ndvi
npmeans = fea_imgs.addBands(cluster2).reduceConnectedComponents(ee.Reducer.mean(),'clusters')

##  standard deviation of nir and ndvi
npstds =  fea_imgs2.addBands(cluster2).reduceConnectedComponents(ee.Reducer.stdDev(),'clusters')

## geometric features
## area
area = ee.Image.pixelArea().addBands(cluster2).reduceConnectedComponents(ee.Reducer.sum(), 'clusters').rename('area')
## perimeter
minMax = cluster2.reduceNeighborhood(ee.Reducer.minMax(), ee.Kernel.square(1))
perimeterPixels = minMax.select(0).neq(minMax.select(1)).rename('perimeter')
perimeter = perimeterPixels.addBands(cluster2).reduceConnectedComponents(ee.Reducer.sum(), 'clusters').rename('perimeter')
## width and height
sizes = ee.Image.pixelLonLat().addBands(cluster2).reduceConnectedComponents(ee.Reducer.minMax(), 'clusters')
width = sizes.select('longitude_max').subtract(sizes.select('longitude_min')).rename('width')
height = sizes.select('latitude_max').subtract(sizes.select('latitude_min')).rename('height')

## textural features
glcm13=np13.select('N').glcmTexture(size= 3)
glcm15=np15.select('N').glcmTexture(size= 3)
glcm17=np17.select('N').glcmTexture(size= 3)
glcm_sele=['N_contrast','N_var','N_ent','N_savg','N_diss']
txtr_sele13=glcm13.select(glcm_sele)
txtr_sele15=glcm15.select(glcm_sele)
txtr_sele17=glcm17.select(glcm_sele)
texture_fea = ee.Image([txtr_sele13,txtr_sele15,txtr_sele17])

texture = texture_fea.addBands(cluster2).reduceConnectedComponents(ee.Reducer.mean(),'clusters')

## compile all feature images
all_fea = ee.Image.cat([npmeans, npstds, area, perimeter, width, height, texture])
all_fea.bandNames().getInfo()


# In[23]:


training = all_fea.sampleRegions(close_veg,scale=1)
clusterer = ee.Clusterer.wekaCascadeKMeans(2, 10).train(training)
result = all_fea.cluster(clusterer,"cluster")
#result.getInfo()


# In[24]:


Map6 = emap.Map(center=[lat,long], zoom=12)
Map6.add_basemap('SATELLITE') 
Map6.addLayer(pts3,{'fillColor':'#4285F4','color':'#4285F4', 'width':1},opacity=0.5)
Map6.addLayer(result.randomVisualizer(), {}, 'clusters')
Map6

