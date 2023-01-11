# -*- coding: utf-8 -*-
from PyQt5.QtCore import (QSettings, QTranslator, qVersion, QCoreApplication,
                          QObject, QThread, pyqtRemoveInputHook, pyqtSignal,
                          QVariant)
from PyQt5.QtGui import QIcon, QDoubleValidator, QColor
from PyQt5.QtWidgets import QAction

from qgis.core import *
from qgis.gui import *

import os.path
import processing
from osgeo import ogr, osr
from math import isnan

import sys
tooldir = r"./Documents"
sys.path.append(tooldir)
import os

#from icworker import ICWorker
import icsimple as ic

import time
import datetime

try: 
    import pandas as pd
    pandasimport = True
except: 
    print("Pandas not installed, using csv module instead.")
    pandasimport = False
    import csv

#fp = r"./Documents/icsimple.py"


iter_start = time.time()


# INPUT DATA
blocks_fp = r"./Documents/IC_inputs/blocks-3067.gpkg"

#filename = "kantakaupunki-points-3067"
filename = "kaarela-TESTAUGUST2"
points_fp = fr"./Documents/IC_inputs/{filename}.gpkg"
print("points fp:\n", points_fp)

# OUTPUT SETTINGS
out_folder = fr"./Documents/IC_results/{filename}"
outcsvfp = os.path.join(out_folder, f"IC_results_{filename}.csv")

if os.path.exists(out_folder) == False:
    os.makedirs(out_folder)
    print("created folder", out_folder)
    
    
#Empty dataframe or dict for results
if pandasimport == True: 
    results = pd.DataFrame()
    results.index.name = 'xyind'

else: 
    results = {} 

# Create layers via QGIS tools to work with IC plugin script
blocks = QgsVectorLayer(blocks_fp, "blocks layer", "ogr")
starting_point_layer_all = QgsVectorLayer(points_fp, "start points layer", "ogr")

print("#############################################")
print("---------------IC ITERATOR ------------------")
print("#############################################")
print("Number of start points :",starting_point_layer_all.featureCount())

project_crs = blocks.crs().toWkt()

counter = 0

for feature in starting_point_layer_all.getFeatures():
    
    counter += 1
    
    print("---------------------------------")
    print("processing point", f"{counter}/{starting_point_layer_all.featureCount()}")
    
    #print(feature.geometry().asPoint().x())
    #print(feature.geometry().asPoint().y())
    #print("xyind:", feature["xyind"])
    
    # Grab xyind at this stage..
    feature_id = feature["xyind"]
    
    # Create new dummy layer to be fed into IC run
    temp_point_layer =  QgsVectorLayer('Point?crs='+project_crs,
                                                'temp_point' , "memory")
                                                
    # get the data provider of the layer
    provider = temp_point_layer.dataProvider()
    
    #add the feature to this temp layer
    with edit(temp_point_layer): 
        provider.addFeatures([feature])
        
    temp_point_layer.updateFields()
    type(temp_point_layer)
    
    print("---------------------------------------------------")
    print("temp layer points :",temp_point_layer.featureCount())

    print("starting IC runner..")
    IC = ic.run(blocks, 
                temp_point_layer, 
                xyind=feature_id, 
                out_folder=out_folder)
    
    print("Finished IC run for", feature_id, "result:", IC, "meters")
    
    
    if pandasimport == True: 
        results.loc[feature_id, "IC"] = IC
        results.to_csv(os.path.join(outcsvfp))
        
    else: 
        results[feature_id] = IC

        with open(outcsvfp, 'w') as csvfile:
            writer = csv.writer(csvfile)
            for key, value in results.items():
                writer.writerow([key, value])
        
print(outcsvfp)

meantime = time.time()
end_time = str(datetime.timedelta(seconds=round(meantime - iter_start, 0)))
print("All iterations, run time: ", end_time)

