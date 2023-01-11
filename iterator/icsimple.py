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

import time
import datetime

#import pandas as pd

MESSAGE_CATEGORY = 'InterfaceCatchment'

def run(blocks_layer, 
        starting_point_layer, 
        xyind, 
        walking_distance=300, 
        deadend_solution=False, 
        out_folder = r""):
            
    """Here you implement your heavy lifting
    Should periodically test for isCanceled() to gracefully
    abort.
    This method MUST return True or False.
    Raising exceptions will crash QGIS, so we handle them
    internally and raise them in self.finished
    """
    
    start = time.time()

    icWorkerContext = QgsProcessingContext()
    
    ################ MANUAL INPUT
    blocks_layer = blocks_layer
    starting_point_layer = starting_point_layer
    
    print("blocks, n:", blocks_layer.featureCount())
    print("points, n:", starting_point_layer.featureCount())
    
    if starting_point_layer.featureCount() == 1:
        
        for feature in starting_point_layer.getFeatures():
            geom = feature.geometry().asPoint()
            x = geom.x()
            y = geom.y()
            #xyind = feature["xyind"]
    else:
        print("invalid input for this test run")
        
    print("X (EPSG:3067):", x)
    print("Y (EPSG:3067):",y)
    print("xyind:", xyind) # from iterator..
    
    walking_distance = walking_distance
    deadend_solution = deadend_solution
    buffering_distance = 40/2 # This is max street width divided by 2. Not used if deadend_solution == False
    x_coordinate = x
    y_coordinate = y
    #project_crs = 'EPSG:3067'
    project_crs = blocks_layer.crs().toWkt()
    
    try:
        
        # read the GUI parameters at run time
        #blocks_layer = self.parent.dlg.mMapLayerComboBox.currentLayer()
        #starting_point_layer = self.parent.dlg.mMapLayerComboBox_2.currentLayer()
        #walking_distance = self.parent.dlg.mQgsDoubleSpinBox.value()
        #deadend_solution = self.parent.dlg.checkBox.checkState()
        #buffering_distance = self.parent.dlg.mQgsDoubleSpinBox_2.value()/2
        #x_coordinate = float(self.parent.dlg.lineEdit.text())
        #y_coordinate = float(self.parent.dlg.lineEdit_2.text())
        #project_crs = self.parent.canvas.mapSettings().destinationCrs().authid()

        #self.log('Started task "%s"' %self.description())
        # abort parameter for stopping function in loop execution

        # debug parameter for enabling or disabling the debugging with vscode
        __debug = False

        starttime = time.time()

        # Preliminary fix the blocks layer
        blocks_layer = processing.run(
            'qgis:fixgeometries',
            {
                'INPUT': blocks_layer,
                'OUTPUT': 'memory:',
            }
        )['OUTPUT']

        # REPROJECT both blocks and sp layers
        blocks_layer = processing.run(
            'qgis:reprojectlayer',
            {
                'INPUT': blocks_layer,
                'TARGET_CRS': project_crs,
                'OUTPUT': 'memory:',
            }
        )['OUTPUT']

        if starting_point_layer:
            starting_point_layer = processing.run(
                'qgis:reprojectlayer',
                {
                    'INPUT': starting_point_layer,
                    'TARGET_CRS': project_crs,
                    'OUTPUT': 'memory:',
                }
            )['OUTPUT']

    #############################################################################################
        # HANDLING THE STARTING POINT COORDINATES
        # Create a new layer that will hold the starting point
        starting_point_layer =  QgsVectorLayer('Point?crs='+project_crs,
                                                    'IC_starting_point' , "memory")
        # get the data provider of the starting point layer
        provider = starting_point_layer.dataProvider()
        with edit(starting_point_layer):
            # create the geometry of the starting point
            starting_pointXY = QgsPointXY(x_coordinate, y_coordinate)
            point_geom = QgsGeometry.fromPointXY(starting_pointXY)
            # Create a new feature with the geometry of the
            # block and the id, perimeter and area attributes
            feature = QgsFeature()
            feature.setGeometry(point_geom)
            feature.setId(1)
            provider.addFeatures( [feature] )
        starting_point_layer.updateFields()

        # This is where I make a buffer around the starting point, the size of
        # a given walking distance
        walking_buffer_layer = processing.run(
            'qgis:buffer',
            {
                'INPUT' : starting_point_layer,
                'END_CAP_STYLE' : 0,
                'OUTPUT' : 'memory:buffer',
                'SEGMENTS' : 360,
                'MITER_LIMIT' : 2,
                'DISTANCE' : walking_distance,
                'JOIN_STYLE' : 0,
                'DISSOLVE' : False
            }
        )['OUTPUT']

        # this is where I need to select from the blocks layer only those features which intersect with the buffer created above
        blocks_layer = processing.run(
            'qgis:extractbylocation',
            {
                'INPUT' : blocks_layer,
                'PREDICATE' : [0],
                'INTERSECT' : walking_buffer_layer,
                'OUTPUT' : 'memory:',
            }
        )['OUTPUT']

        #if self.isCanceled():
        #    return False


    ###################################################################
    ####### HANDLING THE BLOCKS LAYER

        # check if the blocks layer == lines
        if blocks_layer.geometryType() == 1:
            blocks_layer = processing.run(
                'qgis:linestopolygons',
                {
                    'INPUT': blocks_layer,
                    'OUTPUT': 'memory:',
                }
            )['OUTPUT']

        #if self.isCanceled():
        #    return False

        # fix the polygons by using the fix geometries processing
        # algorithm
        fixedgeometries = processing.run(
            'qgis:fixgeometries',
            {
                'INPUT': blocks_layer,
                'OUTPUT': 'memory:',
            }
        )['OUTPUT']

        #if self.isCanceled():
         #   return False

        # dissolve the created polygons in order to make touching polygons
        # into one block
        pathdissolve = processing.run(
            'qgis:dissolve',
            {
                'INPUT': fixedgeometries,
                'OUTPUT': 'memory:',
            }
        )['OUTPUT']

        #if self.isCanceled():
        #    return False

        # convert dissolve result from multipart to singleparts
        blocks_layer = processing.run(
            'qgis:multiparttosingleparts',
            {
                'INPUT': pathdissolve,
                'OUTPUT': 'memory:',
            }
        )['OUTPUT']

        #if self.isCanceled():
        #    return False

        # see if the deadend solution has been selected when the plugin was run
        # if true, run the whole buffering solution for removing deadends from
        # the results
        if deadend_solution:
            # first I need to buffer out the blocks by the given distance
            # amount
            buffer_out_layer = processing.run(
                'qgis:buffer',
                {
                    'INPUT' : blocks_layer,
                    'END_CAP_STYLE' : 1,
                    'OUTPUT' : 'memory:',
                    'SEGMENTS' : 5,
                    'MITER_LIMIT' : 2,
                    'DISTANCE' : buffering_distance,
                    'JOIN_STYLE' : 1,
                    'DISSOLVE' : False
                }
            )['OUTPUT']

            #if self.isCanceled():
            #    return False

            # then I need to buffer back in by the same distance amount
            buffer_in_layer = processing.run(
                'qgis:buffer',
                {
                    'INPUT' : buffer_out_layer,
                    'END_CAP_STYLE' : 1,
                    'OUTPUT' : 'memory:',
                    'SEGMENTS' : 5,
                    'MITER_LIMIT' : 2,
                    'DISTANCE' : -buffering_distance,
                    'JOIN_STYLE' : 1,
                    'DISSOLVE' : False
                }
            )['OUTPUT']

            # finally, fill in any wholes that are left in the blocks
            blocks_layer = processing.run(
                'qgis:deleteholes',
                {
                    'INPUT' : buffer_in_layer,
                    'MIN_AREA' : 0,
                    'OUTPUT' : 'memory:',
                }
            )['OUTPUT']


        # Create a new layer that will hold all the blocks
        final_blocks_layer =  QgsVectorLayer('Polygon?crs='+project_crs, 'IC_blocks' , "memory")
        # get the data provider of the blocks layer
        provider = final_blocks_layer.dataProvider()
        i = 1

        with edit(final_blocks_layer):
            for block in blocks_layer.getFeatures():
                block_geom = block.geometry()
                # Create a new feature with the geometry of the
                # block and the id, perimeter and area attributes
                feature = QgsFeature()
                feature.setGeometry(block_geom)
                feature.setId(i)
                provider.addFeatures( [feature] )
                i += 1
        final_blocks_layer.updateFields()

        #if self.isCanceled():
        #    return False

    ##########################################
        # This is where I extract the boundary of the blocks layer in order to
        # work out the vertices
        blocks_boundary_layer = processing.run(
            'qgis:boundary',
            {
                'INPUT' : blocks_layer,
                'OUTPUT' : 'memory:boundary',
            }
        )['OUTPUT']
        
        #if self.isCanceled():
        #    return False

        # This is where I clip the boundary layer with the walking distance
        # buffer
        clipped_boundary_layer = processing.run(
            'qgis:clip',
            {
                'INPUT' : blocks_boundary_layer,
                'OVERLAY' : walking_buffer_layer, 
                'OUTPUT' : 'memory:clip1',
            }
        )['OUTPUT']
        
        #if self.isCanceled():
        #    return False

        clipped_boundary_layer = processing.run(
            'qgis:advancedpythonfieldcalculator',
            {
                'INPUT' : clipped_boundary_layer,
                'OUTPUT' : 'memory:clip2',
                'FIELD_LENGTH' : 10,
                'GLOBAL' : '',
                'FIELD_TYPE' : 0,
                'FIELD_NAME' : 'ic_boundary_id',
                'FORMULA' : 'value = $id',
                'FIELD_PRECISION' : 1,
            }
        )['OUTPUT']

        #if self.isCanceled():
        #    return False

        # This is where I dissolve the blocks layer 
        dissolved_blocks_layer = processing.run(
            'qgis:dissolve',
            {
                'INPUT' : blocks_layer,
                'OUTPUT' : 'memory:dissolve',
            }
        )['OUTPUT']
        
        print("dissolved_blocks_layer, n:", dissolved_blocks_layer.featureCount())

    ############################################################################
    ######## This is where I will test if the points are visible from the starting

        # Create a new layer that will hold all the lines representing the
        # walkable portions of the boundaries
        lines_layer = QgsVectorLayer('Polygon?crs='+project_crs,
                                                    'lines' , "memory")
        # get the data provider of the lines layer
        lines_provider = lines_layer.dataProvider()
        lines_provider.addAttributes(
            [
                QgsField('boundary_id', QVariant.Int),
            ]
        )
        lines_layer.updateFields()
        lines_boundary_id = lines_provider.fieldNameIndex('boundary_id')


        # Create a new layer that will hold all the resulting points
        new_vertices_layer = QgsVectorLayer('Point?crs='+project_crs,
                                                    'IC_vertices' , "memory")
        # get the data provider of the new vertices layer
        provider = new_vertices_layer.dataProvider()
        provider.addAttributes(
            [
                QgsField('iteration', QVariant.Int),
                QgsField('prev_id', QVariant.Int),
                QgsField('distance', QVariant.Double),
                QgsField('boundary_id', QVariant.Int),
            ]
        )
        # Tell vector layer to update fields in order to get the new layers
        # from data provider
        new_vertices_layer.updateFields()
        distance_id = provider.fieldNameIndex('distance')
        iteration_id = provider.fieldNameIndex('iteration')
        prev_id_id = provider.fieldNameIndex('prev_id')
        boundary_id_id = provider.fieldNameIndex('boundary_id')
        fields = new_vertices_layer.fields()

        iteration = 1

        # Get the feature where all the blocks are dissolved into a single
        # feature
        for f in dissolved_blocks_layer.getFeatures():
            # I have to buffer in the blocks by a small amount because of
            # the floating point error on the newly created points
            blocks_geom = f.geometry().buffer(-0.05, 360)

        # First I need to take in the layer with existing points
        for starting_point in starting_point_layer.getFeatures():
            # take the geometry of the starting point
            sp_geom = starting_point.geometry()
            sp_aspoint = sp_geom.asPoint()

            for boundary in clipped_boundary_layer.getFeatures():
                boundary_geom = boundary.geometry()
                boundary_id = boundary['ic_boundary_id']
                boundary_geom.convertToMultiType()
                intersected_boundaries = boundary_geom.asMultiPolyline()

                writepoints_list = []
                boundary_points = [sp_aspoint]

                for boundary_line in intersected_boundaries:
                    for p_geom in boundary_line:
                        # create a line between the starting point and the vertice
                        # point
                        #if self.isCanceled():
                        #    return False
                        line = QgsLineString(
                            [p_geom.x(), sp_aspoint.x()],
                            [p_geom.y(), sp_aspoint.y()]
                        )
                        line_geom = QgsGeometry(line)
                        if not line_geom.crosses(blocks_geom):
                            distance = walking_distance - line_geom.length()
                            pointarea = QgsRectangle(p_geom.x() - 0.005,p_geom.y() - 0.005,
                                                    p_geom.x() + 0.005,p_geom.y() + 0.005)
                            writepoint = True
                            for point in new_vertices_layer.getFeatures(
                                QgsFeatureRequest().setFilterRect(pointarea)):
                                if distance > point['distance']:
                                    provider.deleteFeatures( [point.id()] )
                                    new_vertices_layer.updateFields()
                                else:
                                    writepoint = False
                            if writepoint:
                                feature = QgsFeature(fields)
                                feature.setGeometry(QgsGeometry(QgsPoint(p_geom)))
                                feature[distance_id] = distance
                                feature[iteration_id] = iteration
                                feature[boundary_id_id] = boundary_id
                                provider.addFeatures( [feature] )
                                new_vertices_layer.updateFields()

                            boundary_points.append(p_geom)

                newlines = []
                for j, p1_geom in enumerate(boundary_points):
                    for k, p2_geom in enumerate(boundary_points[j+1:]):
                        #if self.isCanceled():
                        #    return False
                        line = QgsGeometry(
                            QgsLineString(
                                [p1_geom.x(), p2_geom.x()],
                                [p1_geom.y(), p2_geom.y()]
                            )
                                                        )
                        if line.within(boundary_geom.buffer(0.1, 5)):
                            #newlines.append(line.buffer(0.01, 5, 2, 2, 2))
                            newlines.append(line.buffer(0.01, 5)) # TESTI!!!

                
                if newlines:
                    feature = QgsFeature(lines_layer.fields())
                    walkable_line_geometry = QgsGeometry().unaryUnion(newlines)
                    boundary_geom = boundary_geom.difference(walkable_line_geometry)
                    feature.setGeometry(walkable_line_geometry)
                    feature.setId(boundary_id)
                    feature[lines_boundary_id] = boundary_id
                    lines_provider.addFeatures( [feature] )

        lines_layer.updateFields()
        
        print("lines_layer after first iteration, n:", lines_layer.featureCount())

        iteration += 1
        #if __debug:
        #    self.log('itaration = %s' %iteration)

        while any( new_vertices_layer.getFeatures( '"iteration" = %s AND "distance" > 0.001' % (iteration-1))):
            for boundary in clipped_boundary_layer.getFeatures():
                
                #if self.isCanceled():
                #    return False
                
                boundary_geom = boundary.geometry()
                boundary_id = boundary['ic_boundary_id']
                writepoints_list = []
                # newlines = []
                newlines = QgsGeometry()
                existing_line_features = [l for l in lines_layer.getFeatures('"boundary_id"=%s' %boundary_id)]
                            
                
                if len(existing_line_features) == 1:
                    existing_line = existing_line_features[0]
                    existing_lines = existing_line.geometry()
                    # boundary_geom = boundary_geom.difference(existing_lines)
                else:
                    existing_lines = QgsGeometry()

                for sp in new_vertices_layer.getFeatures(
                    '"iteration" = %s AND "distance" > 0.001' % (iteration-1)
                ):
                    sp_boundary_id = sp['boundary_id']
                    sp_distance = sp['distance']
                    sp_geom = sp.geometry()
                    sp_buffer = sp_geom.buffer(sp_distance, 180)
                    sp_aspoint = sp_geom.asPoint()

                    if boundary_geom.intersects(sp_buffer):

                        intersected_boundaries = boundary_geom.intersection(sp_buffer)
                        intersected_boundaries.convertToMultiType()

                        boundary_points = [sp_aspoint]

                        for boundary_line in intersected_boundaries.asMultiPolyline():
                            for p_geom in boundary_line:
                                # create a line between the starting point and the vertice
                                # point
                                #if self.isCanceled():
                                #    return False
                                line = QgsLineString(
                                    [p_geom.x(), sp_aspoint.x()],
                                    [p_geom.y(), sp_aspoint.y()]
                                )
                                line_geom = QgsGeometry(line)
                                if not line_geom.crosses(blocks_geom):
                                    distance = sp_distance - line_geom.length()
                                    pointarea = QgsRectangle(p_geom.x() - 0.005,p_geom.y() - 0.005,
                                                            p_geom.x() + 0.005,p_geom.y() + 0.005)
                                    writepoint = True
                                    for point in new_vertices_layer.getFeatures(
                                        QgsFeatureRequest().setFilterRect(pointarea)):
                                        if distance > point['distance']:
                                            provider.deleteFeatures( [point.id()] )
                                            new_vertices_layer.updateFields()
                                        else:
                                            writepoint = False
                                    if writepoint:
                                        feature = QgsFeature(fields)
                                        feature.setGeometry(QgsGeometry(QgsPoint(p_geom)))
                                        feature[distance_id] = distance
                                        feature[iteration_id] = iteration
                                        feature[prev_id_id] = sp.id()
                                        feature[boundary_id_id] = boundary_id
                                        provider.addFeatures( [feature] )
                                        new_vertices_layer.updateFields()

                                    boundary_points.append(p_geom)
                                    
                        for j, p1_geom in enumerate(boundary_points):
                            for k, p2_geom in enumerate(boundary_points[j+1:]):
                                #if self.isCanceled():
                                #    return False
                                line = QgsGeometry(
                                    QgsLineString(
                                        [p1_geom.x(), p2_geom.x()],
                                        [p1_geom.y(), p2_geom.y()]
                                    )
                                    )
                                if line.within(boundary_geom.buffer(0.1, 5)):
                                    # newlines.append(line.buffer(0.01, 5, 2, 2, 2))
                                    #newlines = QgsGeometry().unaryUnion([newlines, line.buffer(0.01, 5, 2, 2, 2)])
                                    newlines = QgsGeometry().unaryUnion([newlines, line.buffer(0.01, 5)])

                walkable_line_geometry = QgsGeometry().unaryUnion( [newlines, existing_lines] )
                if newlines:
                    if not existing_lines.isEmpty():
                        with edit(lines_layer):
                            lines_layer.changeGeometry(
                                existing_line.id(),
                                walkable_line_geometry
                            )
                    else:
                        feature = QgsFeature(lines_layer.fields())
                        feature.setGeometry(walkable_line_geometry)
                        feature.setId(boundary_id)
                        feature[lines_boundary_id] = boundary_id
                        lines_provider.addFeatures( [feature] )

            iteration += 1
            #print("round", iteration)
            
            #if __debug:
            #    self.log('iteration = %s' %iteration)

            new_vertices_layer.updateFields()
            lines_layer.updateFields()

        # Fix any invalid geometries made with the lines layer 
        lines_layer = processing.run(
            'qgis:fixgeometries',
            {
                'INPUT': lines_layer,
                'OUTPUT': 'memory:',
            }
        )['OUTPUT']
        
        print("lines_layer after 2nd iteration, n:", lines_layer.featureCount())

        
        fp =os.path.join(out_folder, rf"lines_layer_temp_{xyind}.gpkg")
        _writer = QgsVectorFileWriter.writeAsVectorFormat(lines_layer, fp, 'utf-8',driverName='GPKG')

        # This is where I clip the boundary layer with the walking distance
        # buffer
        walkable_lines_layer = processing.run(
            'qgis:clip',
            {
                'INPUT' : clipped_boundary_layer,
                'OVERLAY' : lines_layer, 
                'OUTPUT' : 'memory: IC_reachable',
            }
        )['OUTPUT']

        # first calculate the lengths of all lines representing IC result
        walkable_lines_layer = processing.run(
            'qgis:fieldcalculator',
            {
                'INPUT' : walkable_lines_layer,
                'OUTPUT' : 'memory:IC_reachable',
                'FIELD_LENGTH' : 20,
                'GLOBAL' : '',
                'FIELD_TYPE' : 0,
                'FIELD_NAME' : 'length',
                'FORMULA' : '$length',
                'FIELD_PRECISION' : 3,
            }
        )['OUTPUT']

        # check the calculated lengths for "nan" values and replace them
        layer_provider = walkable_lines_layer.dataProvider()
        for ic_feature in walkable_lines_layer.getFeatures():
            length = ic_feature['length']
            if isnan(length):
                # replace any nan lengths with calculated ones
                fid = ic_feature.id()
                idx = walkable_lines_layer.fields().names().index('length')
                calculated_length = ic_feature.geometry().length()
                attr_value={idx:calculated_length}
                layer_provider.changeAttributeValues({fid : attr_value})

        ## DOES THIS WORK MODIFIED?
        # reproject start_point_layer with processing context
        starting_point_layer = processing.runAndLoadResults(
            'qgis:reprojectlayer',
            {
                'INPUT': starting_point_layer,
                'TARGET_CRS': project_crs,
                'OUTPUT': 'memory:',
            },
            context=icWorkerContext
        )['OUTPUT']

        walkable_lines_layer = processing.run(
            'qgis:fieldcalculator',
            {
                'INPUT' : walkable_lines_layer,
                'OUTPUT' : 'memory:IC_reachable',
                'FIELD_LENGTH' : 20,
                'GLOBAL' : '',
                'FIELD_TYPE' : 1,
                'FIELD_NAME' : 'IC',
                'FORMULA' : 'sum("length")',
                'FIELD_PRECISION' : 1,
            },
        )['OUTPUT']
        
        print(" walkable_lines_layer, final:",walkable_lines_layer.featureCount())
        
        #fp =os.path.join(out_folder, rf"walkable-lines_final_{xyind}.gpkg")
        #_writer = QgsVectorFileWriter.writeAsVectorFormat(walkable_lines_layer, fp, 'utf-8',driverName='GPKG')

        IC = round(next(walkable_lines_layer.getFeatures())['IC'])
        temp = next(walkable_lines_layer.getFeatures())['IC']

        
        print("IC:", IC)
        print("IC, not rounded:", temp)

        meantime = time.time()
        current_time = str(datetime.timedelta(seconds=round(meantime - start, 0)))
        print("Run time: ", current_time)
        
        # self.log('total running time: %s s' %(endtime-starttime))

        #walkable_lines_layer.setName('IC_' + str(IC))
        #self.walkable_lines_layer = walkable_lines_layer

        #self.starting_point_layer = starting_point_layer

        #return True
        # Return the IC value
        return IC
        
    except:
        print("ERROR, maybe point is not in built up area..")
        f = open('./Documents/errorlog.txt', 'w')
        f.write(f'{xyind}\n')
        
        return None
    
    
def iterate():
    """ This is where we iterate the run-function over all points in a layer, 
    and get the IC values as pandas dataframe --> save as csv"""
    
    iter_start = time.time()
    
    print("Hello World!")
    
    #Filepaths
    blocks_fp = r"W:\STYLE\walkability\IC\blocks-3067.gpkg"
    points_fp = r"W:\STYLE\walkability\IC\test-point-3067.gpkg"
    
    #points_fp = r"W:\STYLE\walkability\IC\kaarela-points-3067.gpkg"
    
    out_folder = r"W:\STYLE\walkability\IC_results"
    
    results = pd.DataFrame()

    # Create layers
    blocks_layer = QgsVectorLayer(blocks_fp, "blocks layer", "ogr")
    starting_point_layer_all = QgsVectorLayer(points_fp, "start points layer", "ogr")
    
    print("#############################################")
    print("---------------IC ITERATOR ------------------")
    print("#############################################")
    print("Number of start points :",starting_point_layer_all.featureCount())

    project_crs = blocks_layer.crs().toWkt()
    
    counter = 0

    for feature in starting_point_layer_all.getFeatures():
        
        counter += 1
        
        print("---------------------------------")
        print("processing point", f"{counter}/{starting_point_layer_all.featureCount()}")
        
        #print(feature.geometry().asPoint().x())
        #print(feature.geometry().asPoint().y())
        #print("xyind:", feature["xyind"])
        
        # Grab xyind at this stage..
        xyind = feature["xyind"]
        
        # Create new dummy layer to be fed into IC run
        temp_point_layer =  QgsVectorLayer('Point?crs='+project_crs,
                                                    'temp_point' , "memory")
                                                    
        # get the data provider of the blocks layer
        provider = temp_point_layer.dataProvider()
        
        #add the feature to this temp layer
        with edit(temp_point_layer): 
            provider.addFeatures([feature])
            
        temp_point_layer.updateFields()
        type(temp_point_layer)
        
        
        
        print("---------------------------------------------------")
        print("temp layer points :",temp_point_layer.featureCount())

        for feat in temp_point_layer.getFeatures():
            #print("Temp points layer check:")
            #print(feat.geometry().asPoint().x())
            #print(feat.geometry().asPoint().y())
            #print("xyind:", feature["xyind"])
            
            #Pass this layer on to run-function
            
            print("starting IC runner..")
            
            IC = run(blocks_layer, temp_point_layer, xyind = xyind)
            
            print("Finished IC run for", xyind, "result:", IC, "meters")
            results.loc[xyind, "IC"] = IC
            
    results.to_file(os.path.join(out_folder, "IC_results.csv"))
    
    meantime = time.time()
    end_time = str(datetime.timedelta(seconds=round(meantime - iter_start, 0)))
    print("All iterations, run time: ", end_time)
    
if __name__ == "__main__":
    
    print("Starting IC iteration..")
    
     #Run the whole thing by running this script...
    iterate()