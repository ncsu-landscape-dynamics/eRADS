
import pandas as pd
import numpy as np
import os
import requests
from bs4 import BeautifulSoup
import multiprocessing
from joblib import Parallel, delayed
import time

os.chdir("C:/Users/wliang5/Desktop/cabi")

dt=pd.read_csv('urlCABIsheet.csv',sep=',',header=0,encoding = "ISO-8859-1")
dt.shape
dt.head()
dt["DistributionMap"]=0
dt.head()


def mapCheck2(i):
    
    import requests
    from bs4 import BeautifulSoup

    url=dt["URL"][i]
    page = requests.get(url)
    soup = BeautifulSoup(page.content, 'html.parser')
    maps = soup.find(id='toDistributionMaps')

    try:
        n=len(maps)
    except :
        n=-1

    if n!= -1:
        dt.loc[i,"DistributionMap"]=1
    return dt
        

#pool=Pool(8)
#pool.map(mapCheck2(dt,np.arange(200)))

le=dt.shape[0]

ncores=multiprocessing.cpu_count()-1

start=time.time()

re=Parallel(n_jobs=ncores)(delayed(mapCheck2)(i) for i in np.arange(le))

end=time.time()
end-start

np.savetxt("CabiDistributionMap.csv",re,  delimiter=', ',header="'Scientific name','Common name','Coverage','URL','DistributionMap'")

    
