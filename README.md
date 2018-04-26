# CSFD Tools
A tool to conduct crater size-frequency distribution (CSFD) measurements on planetary bodies from polygon shapefiles.
Currently in beta version.

The compiled version of this tool can be downloaded here:  
http://www.geo.fu-berlin.de/en/geol/fachrichtungen/planet/software/index.html#faq_csfdtools.

# Requirements

The tool uses the following external libraries: 

GDAL and Shapely for geospatial applications and map projections  
PyQT4, Matplotlib Basemap, and descartes for the user interface layout  
Numpy for array management  

All libraries are included in the compiled version of CSFD Tools. No further software is required.  
The compiled version was tested on Windows 7. 

# Description

The tool uses two polygon shapefiles for CSFD measurements on planetary bodies. One shapefile contains impact craters (circular polygons), the other shapefile contains reference area(s) to be investigated. Any spherical or ellipsoidal reference body and any map projection can be assigned to the shapefiles. However, the same spatial reference has to be applied for both, the crater and the area shapefile to avoid inconsistencies from map distortion effects between the two datasets. Measurements are conducted according to the shapefile geometries and are independent from the attributes of the given shapefiles. Hence, any Geographic Information System which supports the generation of circular polygon features can be used to digitize the data. 

CSFD Tools supports "Traditional Crater Counting", "Buffered Crater Counting", "Non-sparseness Correction", and "Buffered Non-sparseness Correction" measurement techniques (for details, see https://doi.org/10.1016/j.icarus.2014.12.008 and https://doi.org/10.1016/j.icarus.2016.05.015). The modified reference areas can be saved as shapefile geometries. 

CSFD measurements are saved in an SCC or DIAM textfile for further statistical analysis with Craterstats (available at http://www.geo.fu-berlin.de/en/geol/fachrichtungen/planet/software/index.html#faq_craterstats). 

CSFD Tools supports measurements on spheres and biaxial ellipsoids. This includes areas which intersect the Date Line. The implemented techniques are valid for CSFD measurements on planetary bodies that can be approximated by a biaxial ellipsoid with a flattening of 0.3 or lower. Data processing can be parallelized for multi-core support to increase performance. 

# Limitations

CSFD measurements on reference areas which intersect the poles of a planetary body are not supported at this point. We suggest to shift the datasets towards the equator to avoid polar intersections. We will consider this issue for our future work. 

Reference areas which intersect the Date Line are not displayed correctly in the user interface. However, this does not affect the accuracy of the CSFD measurements. 

Data processing in multi-core mode is usually faster when compared to single-core. However, since it takes longer to start a multi-core process, the increase in performance is not significant when the number of impact craters is relatively low (< ~150 craers). In such cases, single-core data processing may be faster. 
