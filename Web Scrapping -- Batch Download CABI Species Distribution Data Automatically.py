#!/usr/bin/env python
# coding: utf-8

# ###### This code shows how to download species distribution data from cabi.org automatically using python. The key python library used here is Selenium. Details on how to use selenium with python can be found here: https://selenium-python.readthedocs.io/

# In[2]:


from selenium import webdriver
import time
import os
import numpy as np
import pandas as pd


# ###### Another notebook documented the code to check whether cabi has the distribution map for species or not (can be found here https://github.com/wanwanliang/For-Open-Science/blob/master/Web%20Scrapping%20CABI%20to%20find%20information%20available%20for%20Species%20with%20Parallel.ipynb). This code was used to find all the species on cabi with distribution map/data.

# In[6]:


os.chdir("C:/Users/liang/Desktop/cabi")

dt=pd.read_csv("CABIwithmap.csv",sep=",",header=0,encoding = "ISO-8859-1")
display(dt.shape)
display(dt.head())
dt["URL"][0]


# ###### The URLs in the dataframe will be the websites to do the scrapping.  Here I'll first show an example, and then use a loop to batch download data.

# In[8]:


# start webdriver and load the website to do the scrapping
driver=webdriver.Chrome()
url=dt["URL"][0]
driver.get(url)
time.sleep(2)

# find the button where we need to click in order to download data. the class of that button in cabi named "Product_export-sml"
button=driver.find_elements_by_class_name("Product_export-sml")
# the function of the below code is to click that button
button[1].click() # there are two buttons with class_name "Product_export-sml", use button[0] to download kml file, and button[1] to download csv file
time.sleep(5) 

# in cabi, a new window pop up and need another click action. the class of this button is named “btn-primary”
driver.find_element_by_class_name("btn-primary").click()
driver.quit()


# ###### And that's it. This simple code do the magic and check your folder for the downloaded zip file. I guess the above code doesn't work well in the jupyter notebook. The loop below works well to batch download the data. 

# In[ ]:


driver=webdriver.Chrome()

for i in np.arange(10):

    url=dt["URL"][i]

    driver.get(url)
    time.sleep(2)

    # find and click 1st button
    driver.find_elements_by_class_name("Product_export-sml")[1].click()
    time.sleep(5)

    # find and click 2nd button
    driver.find_element_by_class_name("btn-primary").click()
   

driver.quit()

