#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ee
#ee.Authenticate()
ee.Initialize()

import folium
import geemap.eefolium as emap
#import subprocess
#import geemap as emap
#from IPython.display import Image
import pandas as pd


# # <font color='blue'> Table of Contents </font><br/>
# 
# 
# <font color='blue' size=4>
# 
# 1. [Define area of interest (AOI) and load NAIP images](#step1)
# 2. [Load presences data and convert to feature_collection](#step2)
# 3. [Image segmentation and select vegatation objects with area between 5-500 m$^2$](#step3)
# 4. [Convert remaining objects  into feature_collection and only keep objects within a given distance of any presence points](#step4)
# 5. [Clustering the remaining vegetation objects into 2 clusters using spectral, phenological, and textural features](#step5)
# 6. [Export major vegetation objects and evaluate by zoomming in Google Earth ](#step6)
# 7. [Summary](#summary)
#     
# </font>

# ## <font color="blue">Step 1-- Define area of interest (AOI) and load NAIP images</font>
# <a id='step1'></a>

# ### <font color="blue"> AOI is defined based on presences density (used a very small AOI here, could be a large bigger AOI)</font>

# In[2]:


collection = ee.ImageCollection('USDA/NAIP/DOQQ')
aoi = ee.Geometry.Polygon([
    [-73.95,40.775],
          [-73.95,40.768],
          [-73.958,40.768],
          [-73.958,40.775]
])
centroid = aoi.centroid()
long, lat = centroid.getInfo()['coordinates']
print("long = {}, lat = {}".format(long,lat))


# ### <font color="blue"> NAIP 2013, 2015, and 2017 are used here. NAIP 2019 is not available yet </font>

# In[3]:


long_lat = ee.Geometry.Point(long, lat)
naip = collection.filterBounds(aoi)
naip15 = collection.filterDate('2015-05-01','2015-10-30')
np15 = naip15.mosaic().clip(aoi)
count15 = naip15.size().getInfo()
naip13 = collection.filterDate('2013-03-01','2013-12-30')
naip17 = collection.filterDate('2017-03-01','2017-12-30')
#naip19 = collection.filterDate('2019-03-01','2019-12-30')  # not available yet

np13 = naip13.mosaic().clip(aoi)
np17 = naip17.mosaic().clip(aoi)
#np19 = naip19.mosaic().clip(aoi)

print('Count of NAIP 13:', naip13.size().getInfo())
print('Count of NAIP 15:', count15)
print('Count of NAIP 17:', naip17.size().getInfo())
#print('Count 19:', naip19.size().getInfo())


# ### <font color="blue"> Calculate NDVI for all three years of NAIP</font>

# In[4]:


#nir, r = imgs.select('N'), imgs.select('R')
ndvi15 = np15.normalizedDifference(["N", "R"])
ndvi13 = np13.normalizedDifference(["N", "R"])
ndvi17 = np17.normalizedDifference(["N", "R"])

ndvi = ndvi15


# ### <font color="blue"> Visualize image and NDVI</font>

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
ndvi_vis = {'min': -1, 'max': 1, 'palette':['black','red','orange' ,'yellow', 'green','blue','black']}

Map.addLayer(ndvi,ndvi_vis)
Map


# ### <font color="blue"> Show areas with NDVI > 0.1 in blue</font>

# In[7]:


veg_mask = ndvi.updateMask(ndvi.gte(0.1))
veg_vis = {'min': 0, 'max': 1, 'palette': ['blue']}
Map.addLayer(veg_mask,veg_vis)
Map


# ## <font color="blue"> Step 2 -- Load presences data and convert to feature_collection</font>
# <a id='step2'></a>

# In[8]:


# read cords
import os
os.chdir("C:/users/liang/Desktop/host_mapping")
crds=pd.read_csv("crds.csv")
print(crds.head(5))
print(crds.shape)


# In[9]:


crdls = crds.values.tolist()
print(crdls[1:10])

pts = ee.List(crdls)
pts2 = ee.Geometry.MultiPoint(pts,proj=aoi.projection())
pts3 = ee.FeatureCollection(pts2)


# ### <font color="blue"> Visualize points data </font>

# In[10]:


Map = emap.Map(center=[lat,long], zoom=12)

Map.add_basemap('SATELLITE') 
Map.addLayer(pts3,{'fillColor':'#4285F4','color':'#4285F4', 'width':1},opacity=0.8)
Map


# ## <font color="blue"> Step 3 -- Image segmentation and select vegatation objects with area between 5-500 m$^{2}$</font>
# <a id='step3'></a>

# ### <font color="blue"> Merge NAIP 15 and NAIP 17 for image segemention</font>

# In[11]:


np_all = np15.addBands(np17.select('R','G','B','N'))
np_all.bandNames().getInfo()


# ### <font color="blue"> Image segmentation and visualize segments </font>

# In[12]:


seed = ee.Algorithms.Image.Segmentation.seedGrid(6)
#seg = ee.Algorithms.Image.Segmentation.GMeans(image=imgs,numIterations=100,pValue=50,neighborhoodSize=500)
#seg = ee.Algorithms.Image.Segmentation.SNIC(image=np15, size=10,compactness= 0, neighborhoodSize=500,connectivity= 8, seeds=seed).select(['R_mean', 'G_mean', 'B_mean', 'N_mean', 'clusters'], ['R', 'G', 'B', 'N', 'clusters'])
seg = ee.Algorithms.Image.Segmentation.KMeans(np_all, 6, 100, 10)
clusters = seg.select('clusters')

Map.addLayer(clusters.randomVisualizer(), {}, 'clusters')
Map


# ### <font color="blue"> Calculate mean ndvi and area of each segment</font>

# In[13]:


## calculate mean ndvi of each segment
seg_ndvi15 = ndvi15.addBands(clusters).reduceConnectedComponents(ee.Reducer.mean(),'clusters').rename('seg_ndvi')
seg_ndvi17 = ndvi17.addBands(clusters).reduceConnectedComponents(ee.Reducer.mean(),'clusters').rename('seg_ndvi')

## calcualte area of each cluster
area = ee.Image.pixelArea().addBands(clusters).reduceConnectedComponents(ee.Reducer.sum(), 'clusters')


# ### <font color="blue"> Filter objects by NDVI and area, keep vegetation objects, and visualize remaining vegetation objects</font>

# In[14]:


seg_veg = clusters.updateMask(seg_ndvi15.gte(0.20))
seg_veg = seg_veg.updateMask(seg_ndvi17.gte(0.20))
seg_veg2 = seg_veg.updateMask(area.gt(5))
seg_veg3 = seg_veg2.updateMask(area.lt(500))
Map2 = emap.Map(center=[lat,long], zoom=16)
Map2.add_basemap('SATELLITE') 
Map2.addLayer(seg_veg3.randomVisualizer(), {})
Map2


# ## <font color="blue"> Step 4 -- Convert remaining objects  into feature_collection and only keep objects within a given distance of any presence points </font>
# <a id='step4'></a>

# ### <font color="blue"> Tried 3 different distance: 15-meter, 10-meter, and 8-meter </font>

# ### <font color="blue"> Convert remaining vegetation objects to features </font>

# In[15]:


vector = seg_veg3.reduceToVectors(scale=1,bestEffort= True, geometryType= 'polygon',labelProperty='label')
vector= ee.FeatureCollection(vector)


# ### <font color="blue"> Only keep objects within 15/10-meter of any presence points </font>

# In[16]:


#join_filter = ee.Filter.withinDistance(15, '.geo', None, '.geo')
#join_filter = ee.Filter.withinDistance(10, '.geo', None, '.geo')
join_filter = ee.Filter.withinDistance(8, '.geo', None, '.geo')
close_veg = ee.Join.simple().apply(vector,  pts3, join_filter)
seg_veg4 = seg_veg3.clipToCollection(close_veg)


# ### <font color="blue"> Visualize remaining vegetation objects </font>

# In[17]:


Map5 = emap.Map(center=[lat,long], zoom=16)
Map5.add_basemap('SATELLITE') 
Map5.addLayer(seg_veg4.randomVisualizer(),{})
Map5.addLayer(pts3,{'fillColor':'#4285F4','color':'#4285F4', 'width':1})
Map5


# ## <font color="blue"> Step 5 -- Clustering the remaining vegetation objects into 2 clusters using spectral, phenological, and textural features  </font>
# <a id='step5'></a>

# ### <font color="blue"> Calculate all features  </font>

# In[18]:


fea_imgs = ee.Image([np13, np15, np17,ndvi13,ndvi15,ndvi17])
fea_imgs2 = ee.Image([np13.select("N"),np15.select("N"),np17.select("N"),ndvi13,ndvi15, ndvi17])

cluster2=clusters.clipToCollection(close_veg)

## mean of all bands and ndvi
npmeans = fea_imgs.addBands(cluster2).reduceConnectedComponents(ee.Reducer.mean(),'clusters')

##  standard deviation of nir and ndvi
npstds =  fea_imgs2.addBands(cluster2).reduceConnectedComponents(ee.Reducer.stdDev(),'clusters')

## geometric features
## area
#area = ee.Image.pixelArea().addBands(cluster2).reduceConnectedComponents(ee.Reducer.sum(), 'clusters').rename('area')

## perimeter
#minMax = cluster2.reduceNeighborhood(ee.Reducer.minMax(), ee.Kernel.square(1))
#perimeterPixels = minMax.select(0).neq(minMax.select(1)).rename('perimeter')
#perimeter = perimeterPixels.addBands(cluster2).reduceConnectedComponents(ee.Reducer.sum(), 'clusters').rename('perimeter')

## width and height
#sizes = ee.Image.pixelLonLat().addBands(cluster2).reduceConnectedComponents(ee.Reducer.minMax(), 'clusters')
#width = sizes.select('longitude_max').subtract(sizes.select('longitude_min')).rename('width')
#height = sizes.select('latitude_max').subtract(sizes.select('latitude_min')).rename('height')


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
#all_fea = ee.Image.cat([npmeans, npstds, area, perimeter, width, height, texture])
all_fea = ee.Image.cat([npmeans, npstds, area, texture])
all_fea.bandNames().getInfo()


# ### <font color="blue"> Cluster all remaining vegetation objects into 2 clusters  </font>

# In[19]:


training = all_fea.sampleRegions(close_veg,scale=1)
#clusterer = ee.Clusterer.wekaCascadeKMeans(2, 10,distanceFunction= "Manhattan").train(training)
#clusterer = ee.Clusterer.wekaCascadeKMeans(minClusters=2, maxClusters=15, maxIterations=100).train(training)
clusterer = ee.Clusterer.wekaKMeans(2,2).train(training)
result = all_fea.cluster(clusterer,"cluster")


# ### <font color="blue"> Check number of features classified in each cluster </font>

# In[20]:


for i in range(2):
    major_cluster = result.updateMask(result.eq(i))
    feaCo = major_cluster.reduceToVectors(scale=1,bestEffort= True)
    n=feaCo.size().getInfo()
    print(n)


# ### <font color="blue"> Extract the major cluster and visualize </font>

# In[21]:


#clu1 = unsuper_clusters .filter(ee.Filter.eq('cluster', 0))
major_cluster = result.updateMask(result.eq(0))
#all_cluster = result.clip(close_veg)
feaCo = major_cluster.reduceToVectors(scale=1,bestEffort= True)
feaCo
#clu1.size().getInfo()


# In[22]:


Map6 = emap.Map(center=[lat,long], zoom=16)
Map6.add_basemap('SATELLITE') 
Map6.addLayer(pts3,{'fillColor':'#4285F4','color':'#4285F4', 'width':1})
Map6.addLayer(major_cluster.randomVisualizer(), {})
#Map6.addLayer(result.randomVisualizer(), {})
Map6


# ###<font color="blue"> Export the major vegation cluster into google drive </font>

# In[23]:


task = ee.batch.Export.table.toDrive(
            collection=feaCo,
            description="TOH_majorCluster_8m",
            folder="myFolder",
            fileFormat="KML",
        )
task.start()


# ## <font color="blue">Step 6 -- Export major vegetation objects and evaluate by zoomming in Google Earth </font>
# <a id='step6'></a>

# #### With the decrease of distance from presences (within which the vegetation objects are passed to the final unsupervised clustering), the number of remaining vegetation objects also decreases.There are situations when the Google Earth Street View is not available, therefore the accuracy metrics is define as TOH / Non-TOH.
# 
# #### True/False = 24/8 for 15-m distance, 24/7 for 10-m, and 15/0 for 8-m. 

# ## <font color="blue"> Summary </font>
# <a id='summary'></a>

# ### To summarize, selecting vegetation objects within 8 meter of any presences points, then clustering these selected objects based on spectral, phenological (NDVI), and textural features, and finally picking up the largest cluster as the TOH. This workflow could pick up TOH at satisfacotry accuracy.  Meanwhile, based on this first quick research, I think there is a very promising potential of using Google Street View images and these open-source presences to map certain trees. 
