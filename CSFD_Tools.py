""" This script includes methods for geodesic measurements with GDAL/OGR as well as workarounds to ensure that GDAL libraries 
are correctly passed to multicore processes and functions. This strongly affects how functions are defined and called in this script. """

import gdal, ogr, osr, datetime, time, numpy, math, os, multiprocessing, multiprocessing.forking, shapely, shapely.ops, re, sys, ctypes, operator, pyproj, images_qr, matplotlib.pyplot as plt
from collections import defaultdict, Counter
from PyQt4.QtGui import *
from PyQt4.QtCore import *
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from descartes import PolygonPatch
from mpl_toolkits.basemap import Basemap
from shapely.wkt import loads 

""" When compiling multiprocessing scripts with Pyinstaller on Windows, the multiprocessing classes Popen and Process need to be redefined. See 
https://github.com/pyinstaller/pyinstaller/wiki/Recipe-Multiprocessing for details. """

class _Popen(multiprocessing.forking.Popen):
	def __init__(self, *args, **kw):
		if hasattr(sys, 'frozen'):
			os.putenv('_MEIPASS2', sys._MEIPASS)
		try:
			super(_Popen, self).__init__(*args, **kw)
		finally:
			if hasattr(sys, 'frozen'):
				if hasattr(os, 'unsetenv'):
					os.unsetenv('_MEIPASS2')
				else:
					os.putenv('_MEIPASS2', '')

class Process(multiprocessing.Process):
	_Popen = _Popen

""" Define user interface layout and functions """

class GUI_window(QWidget):
	def __init__(self):
		super(GUI_window, self).__init__()
		
		self.setFixedSize(750, 625)
		self.setWindowTitle("CSFD Tools")
		self.setWindowIcon(QIcon(":/CSFD_ICON6.ico"))
		self.file_dialog_area()
		
		self.show()
	
	""" ----- Define Buttons ----- """
	
	def file_dialog_area(self):
		global map_params, draw_meridians, patches, path_to_area_file, path_to_crater_file, path_to_craterstats_outfile
		global approach, bufferfactor, bufferfactor_crater, generate_polygon_file, multicore_operation, write_logfile, path_to_outfile
		
		""" Label where area file name is shown when selected """
		
		self.area_label = QLabel("Select Area Shapefile		  ",self)
		self.area_label.move(15, 20)
		
		""" Label where crater file name is shown when selected """
		
		self.crater_label = QLabel("Select Crater Shapefile		  ",self)
		self.crater_label.move(15, 60)
		
		""" Label to describe area field combo box """
		
		self.combo_box_label = QLabel("Select Area by Field",self)
		self.combo_box_label.move(15, 100)
		
		""" Select area shapefile dialog """
		
		path_to_area_file = ""
		button_area_file = QPushButton("...", self)
		button_area_file.resize(30, 25)
		button_area_file.move(220, 15)
		button_area_file.clicked.connect(self.get_area_file)
		
		""" Select crater shapefile dialog """
		
		path_to_crater_file = ""
		button_crater_file = QPushButton("...", self)
		button_crater_file.resize(30, 25)
		button_crater_file.move(220, 55)
		button_crater_file.clicked.connect(self.get_crater_file)
		
		""" Select area by field combo box """
		
		self.select_area_field_combo_box = QComboBox(self)
		self.select_area_field_combo_box.resize(90, 25)
		self.select_area_field_combo_box.move(160, 95) 
		self.connect(self.select_area_field_combo_box, SIGNAL('currentIndexChanged(QString)'), self.get_field_value)
		
		""" Select area by attribute table """
		
		patches = {}
		self.select_area_list = QListWidget(self)
		self.select_area_list.resize(235, 150)
		self.select_area_list.move(15, 135)
		self.select_area_list.itemChanged.connect(self.get_checked_areas) 
		
		""" Map of selected reference areas - draw initial map """
		
		map_params = [1000000, 0, 0, 200000, 200000, 5]
		
		self.overview_map = QWidget(self)
		self.map_layout = QVBoxLayout(self.overview_map)
		self.fig = Figure()
		self.canvas = FigureCanvas(self.fig)
		self.axes = self.fig.add_subplot(111)
		
		draw_meridians = True
		self.update_map(map_params, draw_meridians)##
		#self.map_layout.addWidget(self.canvas)
		
		self.overview_map.resize(465, 555)
		self.overview_map.move(275, 10)
		
		""" Select crater counting technique combo box """
		
		self.select_crater_counting_technique_combo_box = QComboBox(self)
		self.select_crater_counting_technique_combo_box.resize(235, 25)
		self.select_crater_counting_technique_combo_box.move(15, 300) 
		self.select_crater_counting_technique_combo_box.addItems(['Traditional Crater Counting','Buffered Crater Counting','Non-sparseness Correction','Buffered Non-sparseness Correction'])
		approach = "TRAD" # Set initial approach to traditional crater counting if no changes in combo box selection are made
		self.connect(self.select_crater_counting_technique_combo_box, SIGNAL('currentIndexChanged(QString)'), self.get_crater_counting_technique)
		
		""" Label to describe buffer factor """
		
		self.buffer_factor_label = QLabel("Buffer Factor",self)
		self.buffer_factor_label.move(15, 345)
		self.buffer_factor_label.hide() # hide buffer factor field - only show when BCC or BNSC approach is selected
		
		""" Double spin box to describe buffer factor """
		
		self.buff_factor_spin_box = QDoubleSpinBox(self)
		self.buff_factor_spin_box.setRange(1, 4)
		self.buff_factor_spin_box.setValue(1)
		self.buff_factor_spin_box.setDecimals(1)
		self.buff_factor_spin_box.setSingleStep(0.1)
		self.buff_factor_spin_box.resize(90, 25)
		self.buff_factor_spin_box.move(160, 340)
		bufferfactor = 1 # Set initial buffer factor to 1 if no changes in combo box selection are made
		self.buff_factor_spin_box.valueChanged.connect(self.get_buffer_factor)
		self.buff_factor_spin_box.hide() # hide buffer factor field - only show when BCC or BNSC approach is selected
		
		""" Label to describe obliteration factor """
		
		self.obliteration_factor_label = QLabel("Obliteration Factor", self)
		self.obliteration_factor_label.move(15, 385)
		self.obliteration_factor_label.hide() # hide obliteration factor field - only show when NSC or BNSC approach is selected
		
		""" Double spin box to describe obliteration factor """
		
		self.excl_factor_spin_box = QDoubleSpinBox(self)
		self.excl_factor_spin_box.setRange(1, 4)
		self.excl_factor_spin_box.setValue(2)
		self.excl_factor_spin_box.setDecimals(1)
		self.excl_factor_spin_box.setSingleStep(0.1)
		self.excl_factor_spin_box.resize(90, 25)
		self.excl_factor_spin_box.move(160, 380)
		bufferfactor_crater = 2 # Set initial obliteration factor to 2 if no changes in combo box selection are made
		self.excl_factor_spin_box.valueChanged.connect(self.get_obliteration_factor)
		self.excl_factor_spin_box.hide() # hide obliteration factor field - only show when NSC or BNSC approach is selected
		
		""" SCC outfile label """
		
		self.SCC_label = QLabel("Save SCC/DIAM Outfile		  ", self)
		self.SCC_label.move(15, 425)
		
		""" SCC outfile dialog """
		
		path_to_craterstats_outfile = ""
		self.button_SCC_file = QPushButton("...", self)
		self.button_SCC_file.resize(30, 25)
		self.button_SCC_file.move(220, 420)
		self.button_SCC_file.clicked.connect(self.get_SCC_outfile)
		
		""" Shapefile outfile label """
		
		self.out_shapefile_label = QLabel("Save Modified Areas (optional)	  ", self)
		self.out_shapefile_label.move(15, 465)
		self.out_shapefile_label.hide() # hide label when traditional crater counting is selected
		
		""" Shapefile outfile dialog """
		
		self.button_out_shapefile = QPushButton("...", self)
		self.button_out_shapefile.resize(30, 25)
		self.button_out_shapefile.move(220, 460)
		self.button_out_shapefile.clicked.connect(self.get_shapefile_outfile)
		path_to_outfile = ""
		generate_polygon_file = False
		self.button_out_shapefile.hide() # hide button when traditional crater counting is selected
		
		""" Single-core checkbox """
		
		self.singlecore_checkbox = QCheckBox("Multi-core Mode", self)
		self.singlecore_checkbox.move(15, 505)
		self.singlecore_checkbox.stateChanged.connect(self.single_core_mode)
		multicore_operation = False
		
		""" Logfile checkbox """
		
		self.write_logfile_checkbox = QCheckBox("Write Logfile", self)
		self.write_logfile_checkbox.move(15, 545)
		self.write_logfile_checkbox.stateChanged.connect(self.write_logfile_mode)
		write_logfile = False
		
		""" Start button """
		
		button_start_CSFD = QPushButton("Start Crater Counting", self)
		button_start_CSFD.resize(140, 25)
		button_start_CSFD.move(15, 585)
		button_start_CSFD.clicked.connect(self.start_CSFD)
		
		""" Processes labels """
		
#		self.Processes_label = QLabel("Processes:", self)
#		self.Processes_label.move(285, 590)
#		self.no_of_cores = multiprocessing.cpu_count()
#		self.no_of_processes = 1
#		self.process_no_label = {}
#		
#		process_no = 0
#		if self.no_of_cores <= 5:
#			start_placing_y = 590
#		if self.no_of_cores > 5:
#			start_placing_y = 575
#		if self.no_of_processes > 10:
#			ctypes.windll.user32.MessageBoxA(0, "All " + str(self.no_of_processes) + " cores are used for multi-core data processing.\nThe status is only shown for the first ten processes in the user interface.", "Info", 0)
#		while process_no < self.no_of_cores:
#			process_no += 1
#			if process_no <= 5:
#				self.process_no_label[process_no] = QLabel("     ", self)
#				self.process_no_label[process_no].move(315 + (process_no * 70), start_placing_y)
#			if process_no > 5 and process_no <= 10:
#				self.process_no_label[process_no] = QLabel("     ", self)
#				self.process_no_label[process_no].move(315 + ((process_no-5) * 70), start_placing_y + 30)
#		self.process_no_label[1].move(385, 590)
#		self.process_no_label[1].setText("1:")
		
	""" ----- Define Button Functions ----- """
	
	""" Select area shapefile dialog """
	
	def get_area_file(self):
		global path_to_area_file, no_of_area_features_shapefile, extent_areas_xmin, extent_areas_xmax, extent_areas_ymin, extent_areas_ymax, no_of_area_features_shapefile
		global map_params, proj_to_geog_GUI, draw_meridians, sr_area
		path_to_area_file = str(QFileDialog.getOpenFileName(self, 'Select Area Shapefile', 'C:\\', "Shapefiles (*.shp)"))
		area_file_name = str(os.path.basename(path_to_area_file))
		
		""" Show area file name in label """
		
		self.area_label.setText(area_file_name)
		
		""" Read area file """
		
		driver = ogr.GetDriverByName('ESRI Shapefile')
		area_file = driver.Open(path_to_area_file)
		layer = area_file.GetLayer()
		no_of_area_features_shapefile = layer.GetFeatureCount()
		layerdefinition = layer.GetLayerDefn()
		no_of_area_file_attribute_fields = layerdefinition.GetFieldCount()
		
		""" Check if area file is polygon shapefile """
		
		if layerdefinition.GetGeomType() != 3:
			ctypes.windll.user32.MessageBoxA(0, str(area_file_name) + " is not a polygon shapefile. Please select different reference area file. ", "Error", 0)
			path_to_area_file = ""
			area_file_name = "Select Area Shapefile"
			self.area_label.setText(area_file_name)
			return
		
		""" Get field names """
		
		area_file_fields = []
		for field_count in range(no_of_area_file_attribute_fields):
			field_name = layerdefinition.GetFieldDefn(field_count).name
			area_file_fields.append(field_name)
		
		
		""" Update fields in combo box """
		
		self.select_area_field_combo_box.clear()
		self.select_area_field_combo_box.addItems(area_file_fields)
		
		""" Get extent of areas """
		
		extent_areas_xmin, extent_areas_xmax, extent_areas_ymin, extent_areas_ymax = layer.GetExtent()
		
		""" Make polygon from extent and use it for map parameters """
		
		extent_polygon_ring = ogr.Geometry(ogr.wkbLinearRing)
		extent_polygon_ring.AddPoint(extent_areas_xmin, extent_areas_ymin)
		extent_polygon_ring.AddPoint(extent_areas_xmax, extent_areas_ymin)
		extent_polygon_ring.AddPoint(extent_areas_xmax, extent_areas_ymax)
		extent_polygon_ring.AddPoint(extent_areas_xmin, extent_areas_ymax)
		extent_polygon_ring.AddPoint(extent_areas_xmin, extent_areas_ymin)
		
		extent_polygon = ogr.Geometry(ogr.wkbPolygon)
		extent_polygon.AddGeometry(extent_polygon_ring)
		
		""" Get new map parameters (major axis; lat,lon centroid & geogr. extent for meridian distance; height, width in LAEA) and draw map """
		
		sr = layer.GetSpatialRef()
		sr_area = sr
		major_axis = sr.GetSemiMajor()
		
		geogr_sr = sr.CloneGeogCS()
		proj_to_geog_GUI = osr.CoordinateTransformation(sr, geogr_sr)
		extent_polygon.Transform(proj_to_geog_GUI)
		extent_polygon_centroid = extent_polygon.Centroid()
		centroid_X, centroid_Y, centroid_Z = extent_polygon_centroid.GetPoint()
		extent_polygon_geogr = extent_polygon.GetEnvelope()
		xmin_geogr = extent_polygon_geogr[0]
		ymin_geogr = extent_polygon_geogr[2]
		xmax_geogr = extent_polygon_geogr[1]
		ymax_geogr = extent_polygon_geogr[3]
		geogr_width = abs(xmax_geogr - xmin_geogr)
		geogr_height = abs(ymax_geogr - ymin_geogr)
		meridiandistance = round(max(geogr_width, geogr_height)/6)
		
		""" Round meridian distance to 0.1, 0.01, or 0.001 when reference areas are small or round meridian distance to nearest 10 when areas are large """
		
		if meridiandistance >= 10:
			meridiandistance = round(max(geogr_width, geogr_height)/6, -1) 
			if meridiandistance >= 20:
				meridiandistance = round(max(geogr_width, geogr_height)/5, -1) 
				if meridiandistance >= 30:
					meridiandistance = round(max(geogr_width, geogr_height)/4, -1) 
		if meridiandistance < 1: 
			meridiandistance = round(max(geogr_width, geogr_height)/5, 1)
			if meridiandistance < 0.1:
				meridiandistance = round(max(geogr_width, geogr_height)/4, 2)
				if meridiandistance < 0.01:
					meridiandistance = round(max(geogr_width, geogr_height)/3, 3)
					if meridiandistance < 0.001:
						meridiandistance = 1
		
		sr_laea_text = 'PROJCS["PROJECTED_LAMBERT_AEA",'+str(geogr_sr)+',PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["central_meridian",'+str(centroid_X)+'],PARAMETER["latitude_of_origin",'+str(centroid_Y)+'],UNIT["Meter",1.0]]'
		sr_laea_projection = osr.SpatialReference()
		sr_laea_projection.ImportFromWkt(sr_laea_text)
		geog_to_laea = osr.CoordinateTransformation(geogr_sr, sr_laea_projection)
		extent_polygon.Transform(geog_to_laea)
		extent_polygon_laea = extent_polygon.GetEnvelope()
		xmin_LAEA = extent_polygon_laea[0]
		ymin_LAEA = extent_polygon_laea[2]
		xmax_LAEA = extent_polygon_laea[1]
		ymax_LAEA = extent_polygon_laea[3]
		LAEA_width = abs(xmax_LAEA - xmin_LAEA)
		LAEA_height = abs(ymax_LAEA - ymin_LAEA)
		
		map_params = [major_axis, centroid_X , centroid_Y, LAEA_height*1.1, LAEA_width*1.1, meridiandistance]
		draw_meridians = True
		
		self.axes.clear() # clear old plot and add new one with new file selection
		self.update_map(map_params, draw_meridians)
		
	""" Select crater shapefile dialog """
	
	def get_crater_file(self):
		global path_to_crater_file, sr_crater
		
		path_to_crater_file = str(QFileDialog.getOpenFileName(self, 'Select Crater Shapefile', 'C:\\', "Shapefiles (*.shp)"))
		crater_file_name = str(os.path.basename(path_to_crater_file))
		
		""" Read crater file """
		
		driver = ogr.GetDriverByName('ESRI Shapefile')
		crater_file = driver.Open(path_to_crater_file)
		layer = crater_file.GetLayer()
		layerdefinition = layer.GetLayerDefn()
		sr_crater = layer.GetSpatialRef()
		
		""" Check if area file is polygon shapefile """
		
		if layerdefinition.GetGeomType() != 3:
			ctypes.windll.user32.MessageBoxA(0, str(crater_file_name) + " is not a polygon shapefile. Please select different crater file. ", "Error", 0)
			path_to_crater_file = ""
			crater_file_name = "Select Crater Shapefile"
			self.crater_label.setText(crater_file_name)
			return
		
		""" Show crater file name in label """
		
		self.crater_label.setText(crater_file_name)
	
	""" Select area by field combo box """
	
	def get_field_value(self):
		global area_file_field, draw_meridians
		
		""" Read field for area shapefile from combo box """
		
		area_file_field = str(self.select_area_field_combo_box.currentText())
		
		""" read attribute values """
		
		driver = ogr.GetDriverByName('ESRI Shapefile')
		area_file = driver.Open(path_to_area_file)
		layer = area_file.GetLayer()
		
		self.select_area_list.clear()
		
		for area_feature_count in range(no_of_area_features_shapefile):
			feature = layer.GetFeature(area_feature_count)
			attribute_value = str(feature.GetField(area_file_field))
			
			""" Pass attribute values to list and make items checkable """
			
			item = QListWidgetItem(attribute_value)
			item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
			item.setCheckState(Qt.Unchecked)
			self.select_area_list.addItem(item)
		
		""" Clear map when new field is selected """
		
		draw_meridians = True
		self.axes.clear()
		self.update_map(map_params, draw_meridians)
	
	""" Detect when area is checked and read outlines for overview map """
	
	def get_checked_areas(self, item):
		global polygon, patches
		
		state = ['UNCHECKED', 'TRISTATE',  'CHECKED'][item.checkState()]
		
		if state == 'CHECKED':
			
			""" Read polygon geometries """
			
			driver = ogr.GetDriverByName('ESRI Shapefile')
			area_file = driver.Open(path_to_area_file)
			layer = area_file.GetLayer()
			
			for area_feature_count in range(no_of_area_features_shapefile):
				area_feature = layer.GetFeature(area_feature_count)
				field_value = area_feature.GetField(area_file_field)
				
				if str(field_value) == str(item.text()):
					area_shape = area_feature.GetGeometryRef()
					area_shape.Segmentize(20000)
					area_shape.Transform(proj_to_geog_GUI) # needs to be transformed as Basemap can only use geographic / its own projection input 
					
					""" Distinguish between simple polygon and polygon with holes """
					
					if area_shape.GetGeometryCount() == 1: # Polygon without hole
						number_of_inner_polygons = 0
						area_shape_geometry = area_shape.GetGeometryRef(0)
						no_of_polygon_vertices = area_shape_geometry.GetPointCount()
						
						""" Transform coordinates to Basemap and add polygon to plot """
						
						area_shape_wkt = area_shape.ExportToWkt()
						area_shape_shapely = loads(area_shape_wkt) # must be in shapely format due to dependencies of PolygonPatch
						area_shape_basemap = shapely.ops.transform(self.mm, area_shape_shapely)
						
						""" Add polygon as plottable patch and update map """
						
						patches[item.text()] = PolygonPatch(area_shape_basemap, facecolor= 'gray', edgecolor='k')
						self.axes.add_patch(patches[item.text()])
						draw_meridians = False # Skip drawing meridians when plot is refreshed
						
						""" Notify if multipart feature is present. They could be displayed in the map but they can't be analyzed with the tool. That's why they are not shown in the map. """
						
						if no_of_polygon_vertices == 0:
							ctypes.windll.user32.MessageBoxA(0, "No vertices were identified. Is this a multipart feature? If so, please explode multipart feature and try again.", "Error", 0)
							item.setCheckState(Qt.Unchecked)
							return
						
						self.update_map(map_params, draw_meridians) 
					
					if area_shape.GetGeometryCount() > 1: # Polygon with hole 
						number_of_inner_polygons = area_shape.GetGeometryCount() - 1
						area_shape_geometry = area_shape.GetGeometryRef(0)
						no_of_polygon_vertices = area_shape_geometry.GetPointCount()
						
						#print no_of_polygon_vertices, "vertices,", number_of_inner_polygons, "holes"
						
						""" Transform coordinates to Basemap and add polygon to plot """
						
						area_shape_wkt = area_shape.ExportToWkt()
						area_shape_shapely = loads(area_shape_wkt) # must be in shapely format due to dependencies of PolygonPatch
						area_shape_basemap = shapely.ops.transform(self.mm, area_shape_shapely)
						
						""" Add polygon as plottable patch and update map """
						
						patches[item.text()] = PolygonPatch(area_shape_basemap, facecolor= 'gray', edgecolor='k')
						self.axes.add_patch(patches[item.text()])
						draw_meridians = False # Skip drawing meridians when plot is refreshed
						
						""" Notify if multipart feature is present. They could be displayed in the map but they can't be analyzed with the tool. That's why they are not shown in the map. """
						
						if no_of_polygon_vertices == 0:
							ctypes.windll.user32.MessageBoxA(0, "No vertices were identified. Is this a multipart feature? If so, please explode multipart feature and try again.", "Error", 0)
							item.setCheckState(Qt.Unchecked)
							patches[item.text()] = PolygonPatch(area_shape_basemap, facecolor= 'gray', edgecolor='k')
							return
						
						self.update_map(map_params, draw_meridians) 
		
		if state == 'UNCHECKED':
			
			""" Remove area from basemap """
			
			patches[item.text()].remove()
			draw_meridians = False
			self.update_map(map_params, draw_meridians) 
			
	""" Draw map """
	
	def update_map(self, map_params, draw_meridians):
		
		rsphere = map_params[0]
		lon_centr = map_params[1]
		lat_centr = map_params[2]
		LAEA_height = map_params[3]
		LAEA_width = map_params[4]
		meridiandistance = map_params[5]
		
		self.mm = Basemap( 
		lon_0 = lon_centr, lat_0 = lat_centr, 
		height = LAEA_height, width = LAEA_width,
		resolution=None,
		projection='laea',
		rsphere = rsphere,
		ax=self.axes)
		
		if draw_meridians == True:
			merid = self.mm.drawmeridians(numpy.arange(-360, 360, meridiandistance), labels=[False, False, False, True])
			parall = self.mm.drawparallels(numpy.arange(-90, 90, meridiandistance), labels=[True, True, False, False])
		self.fig.canvas.draw() 
		self.map_layout.addWidget(self.canvas)
	
	""" Select crater counting technique combo box """
	
	def get_crater_counting_technique(self):
		global approach, generate_polygon_file, path_to_outfile
		
		""" Get selected crater counting technique """
		
		crater_counting_technique = str(self.select_crater_counting_technique_combo_box.currentText())
		
		if crater_counting_technique == 'Traditional Crater Counting':
			approach = "TRAD"
			#hide buffer factor and obliteration factor fields and shapefile output 
			self.buffer_factor_label.hide()
			self.obliteration_factor_label.hide()
			self.excl_factor_spin_box.hide()
			self.buff_factor_spin_box.hide()
			self.button_out_shapefile.hide()
			self.out_shapefile_label.hide()
			generate_polygon_file = False 
			path_to_outfile = ""
			self.out_shapefile_label.setText("Save Modified Areas (optional)	  ")
		if crater_counting_technique == 'Buffered Crater Counting':
			approach = "BCC"
			#show buffer factor and shapefile output, hide obliteration factor fields 
			self.buffer_factor_label.show()
			self.obliteration_factor_label.hide()
			self.excl_factor_spin_box.hide()
			self.buff_factor_spin_box.show()
			self.button_out_shapefile.show()
			self.out_shapefile_label.show()
		if crater_counting_technique == 'Non-sparseness Correction':
			approach = "NSC"
			#show obliteration factor and shapefile output, hide buffer factor fields 
			self.buffer_factor_label.hide()
			self.obliteration_factor_label.show()
			self.excl_factor_spin_box.show()
			self.buff_factor_spin_box.hide()
			self.button_out_shapefile.show()
			self.out_shapefile_label.show()
		if crater_counting_technique == 'Buffered Non-sparseness Correction':
			approach = "BNSC"
			#show buffer factor and obliteration factor fields and shapefile output 
			self.buffer_factor_label.show()
			self.obliteration_factor_label.show()
			self.excl_factor_spin_box.show()
			self.buff_factor_spin_box.show()
			self.button_out_shapefile.show()
			self.out_shapefile_label.show()
	
	""" Get buffer factor """
	
	def get_buffer_factor(self):
		global bufferfactor
		bufferfactor = float(self.buff_factor_spin_box.value())
	
	""" Get obliteration factor """
	
	def get_obliteration_factor(self):
		global bufferfactor_crater
		bufferfactor_crater = float(self.excl_factor_spin_box.value())
	
	""" Save Craterstats outfile dialog """
	
	def get_SCC_outfile(self):
		global path_to_craterstats_outfile, craterstats_file_name, outfile_type
		
		path_to_craterstats_outfile = str(QFileDialog.getSaveFileName(self, 'Save Outfile for Craterstats Analysis', 'C:\\', "SCC Files (*.scc);;DIAM Files (*.diam)"))
		craterstats_file_name = str(os.path.basename(path_to_craterstats_outfile))
		craterstats_file_name_2, craterstats_file_extension = os.path.splitext(craterstats_file_name)
		if craterstats_file_extension == ".scc":
			outfile_type = "SCC"
		if craterstats_file_extension == ".diam":
			outfile_type = "DIAM"
		
		""" Show Craterstats file name in label """
		
		self.SCC_label.setText(craterstats_file_name)
	
	""" Save shapefile dialog """
	
	def get_shapefile_outfile(self):
		global path_to_outfile, generate_polygon_file
		
		path_to_outfile = str(QFileDialog.getSaveFileName(self, 'Save Shapefile (optional)', 'C:\\', "Shapefiles (*.shp)"))
		outfile_file_name = str(os.path.basename(path_to_outfile))
		
		""" Show shapefile name in label """
		
		self.out_shapefile_label.setText(outfile_file_name)
		
		if path_to_outfile != "":
			generate_polygon_file = True
		else: 
			generate_polygon_file = False
	
	""" Set single-core mode """
	
	def single_core_mode(self):
		global multicore_operation
		
		if self.singlecore_checkbox.isChecked() == True:
			multicore_operation = True
			self.no_of_cores = multiprocessing.cpu_count()#
			self.no_of_processes = self.no_of_cores
#			self.draw_processes(self.no_of_cores, self.no_of_processes)
		else:
			multicore_operation = False
			self.no_of_processes = 1
#			self.draw_processes(self.no_of_cores, self.no_of_processes)
	
	""" Write logfile """
	
	def write_logfile_mode(self):
		global write_logfile
		
		if self.write_logfile_checkbox.isChecked() == True:
			write_logfile = True
		else:
			write_logfile = False
	
	""" Change process labels with single-core / multi-core selection """
	
#	def draw_processes(self, no_of_cores, no_of_processes):
#		
#		if no_of_processes <= 5:
#			start_placing_y = 590
#		if no_of_processes > 5:
#			start_placing_y = 575
#		
#		process_no = 0
#		
#		""" Only display and align process 1 and hide the other ones in the user interface """
#		
#		if no_of_processes == 1:
#			while process_no < self.no_of_cores:
#				process_no += 1
#				if process_no == 1:
#					self.process_no_label[process_no].move(315 + (process_no * 70), start_placing_y)
#				else:
#					self.process_no_label[process_no].setText("")
#		
#		""" Display and allign all processes in the user interface """
#		
#		if no_of_processes > 1:
#			if self.no_of_cores > 10:
#				ctypes.windll.user32.MessageBoxA(0, "All " + str(self.no_of_cores) + " cores are used for multi-core data processing.\nThe status is only shown for the first ten processes in the user interface.", "Info", 0)
#			while process_no < self.no_of_cores:
#				process_no += 1
#				if process_no <= 5:
#					if process_no == 1:
#						self.process_no_label[process_no].move(315 + (process_no * 70), start_placing_y)
#					else:
#						self.process_no_label[process_no].setText(str(process_no) + ":")
#				if process_no > 5 and process_no <= 10:
#					if process_no == 1:
#						self.process_no_label[process_no].move(315 + (process_no * 70), start_placing_y)
#					else:
#						self.process_no_label[process_no].setText(str(process_no) + ":")
	
	""" Start button - Get attributes for crater counting """
	
	def start_CSFD(self):
		
		""" Get checked field attributes from list """
		
		Area_Names = []
		
		for index in range(self.select_area_list.count()):
			if self.select_area_list.item(index).checkState() == Qt.Checked:
				Area_Names.append(str(self.select_area_list.item(index).text()))
		
		""" Show warnings when something is missing """
		
		if path_to_area_file == "":
			ctypes.windll.user32.MessageBoxA(0, "Please select reference area file. ", "Error", 0)
		
		if path_to_crater_file == "":
			ctypes.windll.user32.MessageBoxA(0, "Please select crater file. ", "Error", 0)
		
		if len(Area_Names) == 0:
			ctypes.windll.user32.MessageBoxA(0, "Please select reference areas. ", "Error", 0)
		
		if path_to_craterstats_outfile == "":
			ctypes.windll.user32.MessageBoxA(0, "Please select SCC/DIAM outfile. ", "Error", 0)
		
		if path_to_area_file != "" and path_to_crater_file != "" and len(Area_Names) != 0 and path_to_craterstats_outfile != "":
			if str(sr_area) != str(sr_crater):
				ctypes.windll.user32.MessageBoxA(0, "Spatial reference of area shapefile and crater shapefile doesn't match. Please use the same spatial reference on both files.", "Error", 0)
				return
			print "Selected Parameters: "
			print "Field:", area_file_field, "\n", "Areas:", Area_Names, "\n", "Approach:", approach, "\n", "BF:", bufferfactor, "\n", "OF:", bufferfactor_crater, "\n", "Outfile Type:", outfile_type, "\n", "Outfile Path:", path_to_craterstats_outfile,"\n", "Multicore:", multicore_operation, "\n", "Logfile:", write_logfile, "\n", "Polygon File:", generate_polygon_file, "\n", "Path to Polygon File:", path_to_outfile, "\n_____\n"
			
			CSFD_measurement_main(generate_polygon_file, path_to_outfile, outfile_type, path_to_craterstats_outfile, area_file_field, Area_Names, approach, bufferfactor, bufferfactor_crater,  multicore_operation, write_logfile)

#########################################
#										#
# FUNCTIONS FOR CSFD MEASUREMENTS BELOW	#
#										#
#########################################

""" Main function for processing CSFD measurement """

def CSFD_measurement_main(generate_polygon_file, path_to_outfile, outfile_type, path_to_craterstats_outfile, area_file_field, Area_Names, approach, bufferfactor, bufferfactor_crater,  multicore_operation, write_logfile):
	global union_polygon_centroid, total_area_size, ring_index, generate_point_file, generate_connectors_crater_polygon, layer_polygon, layer_crater, craters_for_counting_list, craters_within_range, craters_inside_area, craters_outside_area, craters_excluded_from_list, buffered_craters_wkt_list, crater_area_list, logfile
	
	""" Write Logfile """
	
	if write_logfile == True:
		logfile_file_path = str(os.path.dirname(path_to_craterstats_outfile))
		craterstats_file_name = str(os.path.basename(path_to_craterstats_outfile))
		craterstats_file_name_2, craterstats_file_extension = os.path.splitext(craterstats_file_name)
		logfile_name = str(craterstats_file_name_2) + ".log"
		logfile_path = logfile_file_path + "/" + logfile_name
		logfile = open(logfile_path, "w") 
		logfile = open(logfile_path, "a") 
		logfile.write("CSFD measurement Logfile" + "\n_____\n\n" + "Field: " + str(area_file_field) + "\n" + "Areas: " + str(Area_Names) + "\n" + "Approach: " + str(approach) + "\n" + "BF: " + str(bufferfactor) + "\n" + "OF: " + str(bufferfactor_crater) + "\n" + "Outfile Type: " + str(outfile_type) + "\n" + "Outfile Path: " + str(path_to_craterstats_outfile) + "\n" + "Multicore: " + str(multicore_operation) + "\n" + "Logfile: " + str(write_logfile) + "\n" + "Polygon File: " + str(generate_polygon_file) + "\n" + "Path to Polygon File: " + str(path_to_outfile) + "\n_____\n\n")
	if write_logfile == False:
		logfile = ""
	
	start_time = time.time()
	start_time2 = time.time()
	
	""" Parameters to generate further output shapefiles for tesing geodesic polygon buffer and crater identification. Not connected to user interface. """
	
	generate_point_file = False # Used to identify possible errors during modification - Not implemented in user interface - Singlecore only
	generate_connectors_crater_polygon = False # Used to identify possible errors during crater identification - Not implemented in user interface - Singlecore only
	
	""" Define number of multicore processes according to number of CPU cores. One process is subtracted since splitting of lists (for multicore process definition) 
	would otherwise result in one list (one process) extra.  """
	
	if multicore_operation == True:
		no_of_processes = multiprocessing.cpu_count() - 1 
	
	""" Open files """
	
	driver = ogr.GetDriverByName('ESRI Shapefile')
	area_file = driver.Open(path_to_area_file)
	crater_file = driver.Open(path_to_crater_file)
	layer = area_file.GetLayer()
	layer_crater = crater_file.GetLayer()
	
	print "Number of Areas to be examined:", len(Area_Names)
	if write_logfile == True:
		logfile.write("Number of Areas to be examined: " + str(len(Area_Names)) + "\n")
	
	get_spheroid_information(layer, layer_crater)
	
	""" Single outfiles as well as point and connector shapefiles are only generated using define_outfile_parameters during singlecore computation """
	if multicore_operation == False:
		define_outfile_parameters(driver, path_to_outfile, sr, geogr_sr, sr_equal_area_projection)
	
	define_coordinate_transformations(sr, geogr_sr)
	
	""" Write SCC / DIAM Header """
	
	write_crater_stats_file_header()
	read_crater_features(layer_crater, proj_to_geog, geog_to_eq_area_proj, write_logfile, logfile)
	
	no_of_area_features_shapefile = layer.GetFeatureCount()
	crater_area_list = []
	area_sizes = []
	Area_IDs = []
	all_craters = []
	n_area = 1
	ring_index = 0
	
	""" Spatial_reference needs to be transferred as WKT during multiprocessing. GDAL objects such as geometries, spatial references etc. can't be passed to multicore process ('can't pickle object'). """
	
	sr_wkt = sr.ExportToWkt()
	
	""" Iterate area features and merge research areas if more than one research area is selected (as fractions are calculated from one coherent area which needs to be buffered). 
	Also, craters within research area (and buffered research area) are determined """
	
	union_polygon = ogr.Geometry(ogr.wkbPolygon)
	for an in Area_Names:
		i = 0
		while i < no_of_area_features_shapefile:
			area_feature = layer.GetFeature(i)
			area_feature2 = layer.GetFeature(i)
			area_feature_SCC_File = layer.GetFeature(i)
			Area_ID = area_feature.GetField(area_file_field)
			if str(Area_ID) == str(an):
				geom = area_feature.GetGeometryRef()
				union_polygon = union_polygon.Union(geom)
				union_polygon_centroid = union_polygon.Centroid()
				
				""" Write SCC/DIAM area information. geom is reprojected to geographic coordinates here. """
				
				write_crater_stats_file_area(n_area, area_feature, area_feature2, area_feature_SCC_File, union_polygon, Area_Names, Area_ID, an, write_logfile, logfile)
				
				#####################################################################################
				#																					#
				#			Approaches for crater count (Step 1: find relevant craters):			#
				#																					#
				#####################################################################################
				
				########################
				# Traditional Approach #
				########################
				
				if approach == "TRAD":
					
					print "Determining craters within research area..."
					if write_logfile == True:
						logfile.write("\nDetermining craters within research area...\n")
					st = time.time()
					
					""" Find relevant craters (get_craters_for_trad function): """
					
					out_q = multiprocessing.Queue()
					lock = multiprocessing.Lock()
					get_craters_for_trad(out_q, crater_features_list, geom, generate_connectors_crater_polygon, lock, sr_wkt)
					
					if len(all_craters) == 0:
						all_craters = craters_for_counting_list
					else:
						all_craters = all_craters + craters_for_counting_list
					
					print "Done.\n-\n", len(all_craters), "craters inside area\n-\nElapsed time for crater detection:", round(time.time()-st,2), "sec.\n\n-----\n"
					if write_logfile == True:
						logfile.write("Done.\n-\n" + str(len(all_craters)) + " craters inside area\n-\nElapsed time for crater detection: " + str(round(time.time()-st,2)) + " sec.\n\n-----\n")
				
				############################
				# Buffered Crater Counting #
				############################
				
				if approach == "BCC":
					
					print "Determining craters within research area and buffered research area..."
					if write_logfile == True:
						logfile.write("Determining craters within research area and buffered research area...\n")
					st = time.time()
					
					""" Find relevant craters (in get_craters_for_buffering_BCC function): """
					
					out_q = multiprocessing.Queue()
					out_q_2 = multiprocessing.Queue()
					out_q_3 = multiprocessing.Queue()
					lock = multiprocessing.Lock()
					get_craters_for_buffering_BCC(out_q, out_q_2, out_q_3, crater_features_list, geom, vertices_list, flattening, major_axis, bufferfactor, generate_connectors_crater_polygon, lock, sr_wkt)
					if len(all_craters) == 0:
						all_craters = craters_for_counting_list
					else:
						all_craters = all_craters + craters_for_counting_list
					""" get total number of craters within buffer range """
					total_no_of_craters_within_range = len(craters_within_range)
					""" get total number of craters inside area """
					total_no_of_craters_within_area = len(craters_inside_area)
						
					print "Done.\n-\n", total_no_of_craters_within_area, "craters inside area, ", total_no_of_craters_within_range, "craters inside buffered area, ", no_of_crater_features-len(all_craters), "craters outside area.\n-\nElapsed time for crater detection:", round(time.time()-st,2), "sec.\n\n-----\n"
					if write_logfile == True:
						logfile.write("Done.\n-\n" + str(total_no_of_craters_within_area) + " craters inside area, " + str(total_no_of_craters_within_range) + " craters inside buffered area, " + str(no_of_crater_features-len(all_craters)) + " craters outside area.\n-\nElapsed time for crater detection: " + str(round(time.time()-st,2)) + " sec.\n\n-----\n")
				
				##################################################################
				# Non-sparseness correction & Buffered non-sparseness correction #
				##################################################################
				
				if approach == "NSC" or approach == "BNSC":
					
					print "Determining craters within research area and buffered research area..."
					if write_logfile == True:
						logfile.write("\nDetermining craters within research area and buffered research area...\n")
					st = time.time()
					
					""" Find relevant craters (in get_craters_for_buffering_NSC_BNSC function): """
					
					out_q = multiprocessing.Queue()
					out_q_2 = multiprocessing.Queue()
					out_q_3 = multiprocessing.Queue()
					out_q_4 = multiprocessing.Queue()
					lock = multiprocessing.Lock()
					
					""" Only one process is started during multi-core computation. It consumes more time to split the lists and determine craters 
					of interest than to just have a singlecore process. The reason we do a multicore operation anyways is that GDAL is not always  
					passed to a multicore process on Windows (forking issue). This effect doesn't occur with Linux. If we just do a single-core process here and have the crater buffering as a MC process later, 
					GDAL will not be passed to the MC process at a later stage. This is really a pain. In order to make multi-core avaliable at a later
					stage of the script, we have to start the process here even though there is no benefit time-wise. """
					
					if multicore_operation == True:
						
						pro = Process(target = get_craters_for_buffering_NSC_BNSC, args=(approach, out_q, out_q_2, out_q_3, out_q_4, crater_features_list, geom, vertices_list, flattening, major_axis, bufferfactor, bufferfactor_crater, generate_connectors_crater_polygon, lock, sr_wkt))
						
						pro.start()
						
						craters_from_multiprocessing = out_q.get()
						craters_within_range = out_q_2.get()
						craters_inside_area = out_q_3.get()
						craters_excluded_from_list = out_q_4.get()
						
						pro.join()
						pro.terminate()
						
						""" If at least one crater relevant for age determination was determined during multiprocess operation, append list to all_craters """
						
						if len(craters_from_multiprocessing) > 0:
							if len(all_craters) == 0:
								all_craters = craters_from_multiprocessing
							else:
								all_craters = all_craters + craters_from_multiprocessing
					
					if multicore_operation == False:
						get_craters_for_buffering_NSC_BNSC(approach, out_q, out_q_2, out_q_3, out_q_4, crater_features_list, geom, vertices_list, flattening, major_axis, bufferfactor, bufferfactor_crater, generate_connectors_crater_polygon, lock, sr_wkt)
						if len(all_craters) == 0:
							all_craters = craters_for_counting_list
						else:
							all_craters = all_craters + craters_for_counting_list
						
					""" get total number of craters within buffer range """
					total_no_of_craters_within_range = len(craters_within_range)
					""" get total number of craters inside area """
					total_no_of_craters_within_area = len(craters_inside_area)

					print "\n", total_no_of_craters_within_area, "craters in area, ", total_no_of_craters_within_range, "craters in buffered area, ", no_of_crater_features - total_no_of_craters_within_area - total_no_of_craters_within_range, "craters outside area. \n", len(craters_excluded_from_list), "/", total_no_of_craters_within_area + total_no_of_craters_within_range, "craters excluded from CSFD analysis due to obliteration effects. \n\nElapsed time for crater detection:", round(time.time()-st,2), "sec.\n\n-----\n"
					if write_logfile == True:
						logfile.write("\n" + str(total_no_of_craters_within_area) + " craters in area, " + str(total_no_of_craters_within_range) + " craters in buffered area, " + str(no_of_crater_features - total_no_of_craters_within_area - total_no_of_craters_within_range) + " craters outside area. \n" + str(len(craters_excluded_from_list)) + " / " + str(total_no_of_craters_within_area + total_no_of_craters_within_range) + " craters excluded from CSFD analysis due to obliteration effects. \n\nElapsed time for crater detection: " + str(round(time.time()-st,2)) + " sec.\n\n-----\n")
					del total_no_of_craters_within_range
					del total_no_of_craters_within_area
					
					""" Add craters that are located on an ejecta blanket. They are used to buffer (to consider their ejecta blanket for the 
					reference area modification) but not in the actual CSFD measurement (since they don't obliterate the actual count area but an area
					which was affected by resurfacing). """
					
					all_craters = all_craters + craters_excluded_from_list
				
				area_sizes.append(Area_Size)
				Area_IDs.append(Area_ID)
				n_area += 1
				i += 1
			else:
				i += 1
	
	""" Get total area of reference areas. """
	
	total_area_size = sum(area_sizes)
	
	""" Clean crater lists from duplicates! This needs to be done as craters might falsly be added multiple times. If area 1 overlaps crater A and 
	area 2 overlaps crater A as well, area 1 and area 2 would be buffered two times. During BNSC, duplicate craters are deleted according to 
	their distance to the reference area. """
	
	if approach == "BCC":
		
		all_craters_sorted = []
		
		for orig_crater in all_craters:
			if orig_crater not in all_craters_sorted:
				all_craters_sorted.append(orig_crater)
		
		all_craters = all_craters_sorted
	
	if approach == "NSC" or approach == "BNSC":
		
		""" Identify crater IDs which occur multiple times. """
		
		all_craters_id_count = Counter(crt[0] for crt in all_craters)
		all_craters_duplicate_ids = [k for k, v in all_craters_id_count.iteritems() if v > 1]
		
		""" If a crater ID occurs multiple times (crater overlaps more than one area), compare distance to reference areas and keep entry for 
		the crater that is closest. """
		
		if len(all_craters_duplicate_ids) > 0:
			
			for duplicate_id in all_craters_duplicate_ids:
				
				""" Save all duplicate craters is a list. """
				
				overlapping_craters = []
				
				for overlapping_crater in all_craters:
					overlapping_crater_id = overlapping_crater[0]
					
					if overlapping_crater_id == duplicate_id:
						overlapping_craters.append(overlapping_crater)
						
						""" Remove all overlapping craters from all_craters """
						
						all_craters.remove(overlapping_crater)
					
				""" Identify overlapping crater with minimum distance value from overlapping_craters list and add it back to all_craters """
				
				duplicate_crater_with_minimum_distance = min(overlapping_craters, key=lambda item: item[4])
				all_craters.append(duplicate_crater_with_minimum_distance)
		
	""" Sort craters according to crater diameter (descending). """
	
	all_craters = sorted(all_craters, key = operator.itemgetter(1), reverse = True)
	
	#####################################################################################
	#																					#
	#			Approaches for crater count (Step 2: modify reference area):			#
	#																					#
	#####################################################################################
	
	st_buffer = time.time()
	
	if approach == "TRAD":
		
		write_crater_stats_file_craters_TRAD(all_craters)
	
	if approach == "BCC":
		
		""" Buffer research areas """
		
		if multicore_operation == True:
			
			lock = multiprocessing.Lock()
			
			""" It may happen that the global variable layer_polygon is not correctly passed to the buffer_area function. 
			To avoid that, the variable is directly passed to the function. During multi-core, the layer_polygon variable is generated within the 
			buffer_area function (since data is not merged but multiple files are generated). Accordingly, the layer_polygon variable definition 
			here is a placeholder. """
			
			layer_polygon = 0
			
			""" Adapt no of processes to no of craters: if no of craters is lower than no of cores, no of processes is equal to no of craters """
			
			if len(all_craters) < no_of_processes: 
				no_of_processes = len(all_craters)
			
			""" Rearange list of craters for optimized multicore operation. Large craters lead to simple polygon modifications; small crater lead to complex polygon modifications.
			Craters are sorted in alternating crater diametesr so that the workload is distributed equally among the multi-core processes. """
			
			# sort list according to element 1 (crater diameter) in list's list
			all_craters = sorted(all_craters, key = operator.itemgetter(1))
			
			# rearange alternating
			all_craters = sum(zip(reversed(all_craters), all_craters), ())[:len(all_craters)]
			
			""" split all_craters list according to no of processes """
			
			all_craters_splitted = [all_craters[w:w + (len(all_craters)/no_of_processes)] for w in xrange(0, len(all_craters), (len(all_craters)/no_of_processes))]
			
			""" Queue for passing list results """
			
			crater_area_out_q = multiprocessing.Queue()
			
			processes = dict()
			process_count = 0
			
			""" Process definition """
			
			for all_craters_splitted_part in all_craters_splitted:
				processes[process_count] = Process(target = buffer_area, args=(union_polygon, crater_area_list, all_craters_splitted_part, sr_wkt, generate_point_file, generate_polygon_file, generate_connectors_crater_polygon, flattening, major_axis, bufferfactor, crater_area_out_q, Area_IDs, multicore_operation, path_to_outfile, process_count, lock, layer_polygon, write_logfile, logfile))
				process_count += 1
			process_count2 = 0
			
			""" Process start """
			
			for all_craters_splitted_part in all_craters_splitted:
				processes[process_count2].start()
				process_count2 += 1
				
			""" Get crater area lists from multiprocessing operations.  """
			
			for all_craters_splitted_part in all_craters_splitted:
				crater_area_from_multiprocessing = crater_area_out_q.get()
				
				""" Merge lists from multiprocessing to one list. """
				
				if len(crater_area_list) == 0:
					crater_area_list = crater_area_from_multiprocessing
				else:
					crater_area_list = crater_area_list + crater_area_from_multiprocessing
				
			""" Terminate processes. """
			
			process_count3 = 0
			for all_craters_splitted_part in all_craters_splitted:
				processes[process_count3].join()
				processes[process_count3].terminate()
				process_count3 += 1
			
		if multicore_operation == False:
			if generate_polygon_file == False:
				layer_polygon = 0
			crater_area_out_q = multiprocessing.Queue()
			process_count = 0
			lock = multiprocessing.Lock()
			buffer_area(union_polygon, crater_area_list, all_craters, sr, generate_point_file, generate_polygon_file, generate_connectors_crater_polygon, flattening, major_axis, bufferfactor, crater_area_out_q, Area_IDs, multicore_operation, path_to_outfile, process_count, lock, layer_polygon, write_logfile, logfile)
		
		print "Done.\n-\n", len(all_craters), "Buffers created. \n-\nElapsed time for buffering:", str(round(time.time() - st_buffer, 2)), "sec.\n\n_____\n" 
		if write_logfile == True:
			logfile.write("\nDone.\n-\n" + str(len(all_craters)) + " Buffers created. \n-\nElapsed time for buffering: " + str(round(time.time() - st_buffer, 2)) + " sec.\n\n_____\n" )
		create_crater_fractions_list(crater_area_list)
		write_crater_stats_file_craters(crater_fractions_list)
	
	if approach == "NSC" or approach == "BNSC":
		
		craters_for_counting_list = all_craters 
		
		print "Buffering " + str(len(craters_for_counting_list)) + " Craters for non-sparseness correction..."
		if write_logfile == True:
			logfile.write("\nBuffering " + str(len(craters_for_counting_list)) + " Craters for non-sparseness correction...\n")
		
		if multicore_operation == True:
			
			lock = multiprocessing.Lock()
			
			""" Aadapt no of processes to no of craters: if no of craters is lower than no of cores, no of processes is equal to no of craters """
			
			if len(craters_for_counting_list) < no_of_processes: 
				no_of_processes = len(craters_for_counting_list)
			
			""" split all_craters list according to no of processes """
			
			craters_for_counting_list_splitted = [craters_for_counting_list[w:w + (len(craters_for_counting_list)/no_of_processes)] for w in xrange(0, len(craters_for_counting_list), (len(craters_for_counting_list)/no_of_processes))]
			
			""" Queue for passing list results """
			buffered_craters_out_q = multiprocessing.Queue()
			
			""" Process definition """
			
			processes = dict()
			process_count = 0
			for craters_for_counting_list_splitted_part in craters_for_counting_list_splitted:
				processes[process_count] = Process(target = NSC_BNSC_buffer_craters, args=(buffered_craters_out_q, craters_for_counting_list_splitted[process_count], flattening, major_axis, sr_wkt, lock, bufferfactor_crater))
				process_count += 1
			process_count2 = 0
			
			""" Process start """
			
			for craters_for_counting_list_splitted_part in craters_for_counting_list_splitted:
				processes[process_count2].start()
				process_count2 += 1
			
			""" Get buffered crater WKT and crater diameter lists from multiprocessing operations  """
			
			for craters_for_counting_list_splitted_part in craters_for_counting_list_splitted:
				buffered_craters_from_multiprocessing = buffered_craters_out_q.get()
				
				""" Merge lists from multiprocessing to one list """
				
				if len(crater_area_list) == 0:
					crater_area_list = buffered_craters_from_multiprocessing
				else:
					crater_area_list = crater_area_list + buffered_craters_from_multiprocessing
			
			buffered_craters_wkt_list = crater_area_list
			
			""" A new crater area list is used at a later point to calculate fractions. Thus, crater_area_list can be deleted. """
			
			del crater_area_list
			crater_area_list = []
			
			""" Terminate processes. """
			
			process_count3 = 0
			for craters_for_counting_list_splitted_part in craters_for_counting_list_splitted:
				processes[process_count3].join()
				processes[process_count3].terminate()
				process_count3 += 1
			
			""" Sort buffered crater wkt list according to descending crater diameter. Since multi-core processes might end at a different time, the merged crater lists 
			are not in order anymore. """
			
			buffered_craters_wkt_list_sorted = sorted(buffered_craters_wkt_list, key = lambda x:x[3], reverse=True)
			buffered_craters_wkt_list = buffered_craters_wkt_list_sorted
		
		if multicore_operation == False:
			buffered_craters_out_q = multiprocessing.Queue()
			NSC_BNSC_buffer_craters(buffered_craters_out_q, craters_for_counting_list, flattening, major_axis, sr_wkt, lock, bufferfactor_crater)
		
		if approach == "BNSC":
			print "\nApplying buffer to reference area..."
			if write_logfile == True:
				logfile.write("\nApplying buffer to reference area...\n")
		
		""" It may happen that the global variable layer_polygon is not correctly passed to the buffer_area function. 
		To avoid that, the variable is directly passed to the function. During multi-core, the layer_polygon variable is generated within the 
		buffer_area function (since data is not merged but multiple files are generated). Accordingly, the layer_polygon variable definition 
		here is a placeholder. """
		
		if multicore_operation == True or generate_polygon_file == False:
			layer_polygon = 0
		
		""" This list is used during BNSC to buffer the reference area. Since during multi-core a full list with crater information is also needed 
		(craters_for_counting_list will not be splitted), a separate list with crater information is required. """
		
		craters_for_counting_list_BNSC = craters_for_counting_list
		
		""" Use geometries in buffered_craters_wkt_list to erase from initial reference area. For every crater, all larger craters are excluded from the count area. """
		
		if multicore_operation == True and approach == "BNSC":
			
			""" Rearange list of craters for optimized multicore operation: sort alternating accordig to crater diameter so that the 'workload' is distributed 
			equally among the processes (large craters require simple polygon modifications; small craters require complex polygon modifications). """
			
			# sort list according to element 1 (crater diameter) in list's list
			craters_for_counting_list_BNSC = sorted(craters_for_counting_list_BNSC, key = operator.itemgetter(1))
			
			# rearange alternating
			craters_for_counting_list_BNSC = sum(zip(reversed(craters_for_counting_list_BNSC), craters_for_counting_list_BNSC), ())[:len(craters_for_counting_list_BNSC)]
			
			""" split craters_for_counting list according to no of processes """
			
			craters_for_counting_list_BNSC_splitted = [craters_for_counting_list_BNSC[w:w + (len(craters_for_counting_list_BNSC)/no_of_processes)] for w in xrange(0, len(craters_for_counting_list_BNSC), (len(craters_for_counting_list_BNSC)/no_of_processes))]
			
			""" Queue for passing list results """
			
			crater_area_out_q = multiprocessing.Queue()
			
			processes = dict()
			process_count = 0
			
			""" Process definition """
			
			for craters_for_counting_list_BNSC_splitted_part in craters_for_counting_list_BNSC_splitted:
				
				processes[process_count] = Process(target = NSC_BNSC_exclude_craters, args=(approach, buffered_craters_wkt_list, union_polygon, sr_wkt, generate_polygon_file, path_to_outfile, craters_for_counting_list, craters_for_counting_list_BNSC_splitted_part, multicore_operation, layer_polygon, bufferfactor, bufferfactor_crater, process_count, crater_area_out_q, flattening, major_axis, Area_IDs, write_logfile, logfile))
				process_count += 1
			process_count2 = 0
			
			""" Process start """
			
			for craters_for_counting_list_BNSC_splitted_part in craters_for_counting_list_BNSC_splitted:
				processes[process_count2].start()
				process_count2 += 1
			
			""" Get crater area lists from multiprocessing operations  """
			
			for craters_for_counting_list_BNSC_splitted_part in craters_for_counting_list_BNSC_splitted:
				crater_area_from_multiprocessing = crater_area_out_q.get()

				""" merge lists from multiprocessing to one list """
				
				if len(crater_area_list) == 0:
					crater_area_list = crater_area_from_multiprocessing
				else:
					crater_area_list = crater_area_list + crater_area_from_multiprocessing

			""" Terminate Processes. """
			
			process_count3 = 0
			for craters_for_counting_list_BNSC_splitted_part in craters_for_counting_list_BNSC_splitted:
				processes[process_count3].join()
				processes[process_count3].terminate()
				process_count3 += 1
			
		""" During NSC, craters are removed from the reference area in single-core mode. Multi-core is not efficient at this point. """
		
		if multicore_operation == False or approach == "NSC":
			process_count = 0 # placeholder
			crater_area_out_q = multiprocessing.Queue() # placeholder
			NSC_BNSC_exclude_craters(approach, buffered_craters_wkt_list, union_polygon, sr_wkt, generate_polygon_file, path_to_outfile, craters_for_counting_list, craters_for_counting_list_BNSC, multicore_operation, layer_polygon, bufferfactor, bufferfactor_crater, process_count, crater_area_out_q, flattening, major_axis, Area_IDs, write_logfile, logfile)
		
		if approach == "NSC":
			print "\n-----\n\n", len(craters_for_counting_list) - len(crater_area_list), "craters removed due to location outside reference area.\n", len(crater_area_list), "areas created. \n\nElapsed time for area modification:", str(round(time.time() - st_buffer, 2)), "sec.\n\n_____\n" 
			if write_logfile == True:
				logfile.write("\n-----\n\n" + str(len(craters_for_counting_list) - len(crater_area_list)) + " craters removed due to location outside reference area.\n" + str(len(crater_area_list)) + " areas created. \n\nElapsed time for area modification: " + str(round(time.time() - st_buffer, 2)) + " sec.\n\n_____\n")
		if approach == "BNSC":
			print "\n-----\n\n", len(crater_area_list), "areas created. \n\nElapsed time for area modification:", str(round(time.time() - st_buffer, 2)), "sec.\n\n_____\n" 
			if write_logfile == True:
				logfile.write("\n-----\n\n" + str(len(crater_area_list)) + " areas created. \n\nElapsed time for area modification: " + str(round(time.time() - st_buffer, 2)) + " sec.\n\n_____\n")
		
		create_crater_fractions_list(crater_area_list)
		write_crater_stats_file_craters(crater_fractions_list)
	
	elapsed_time = time.time() - start_time
	
	print "Done. Elapsed time in total: " + str(round(elapsed_time, 2)) + " sec."
	if write_logfile == True:
		logfile.write("\nDone. Elapsed time in total: " + str(round(elapsed_time, 2)) + " sec.")
		logfile.close()
	
	ctypes.windll.user32.MessageBoxA(0, "Done. Elapsed time in total: " + str(round(elapsed_time, 2)) + " sec." , "Done", 0)

""" Get information on reference body from Shapefile metadata. """

def get_spheroid_information(layer, layer_crater): 
	global geogr_sr, sr, sr2, sr_equal_area_projection, major_axis, minor_axis, flattening, sr_name
	sr = layer.GetSpatialRef()
	sr2 = layer_crater.GetSpatialRef()
	
	""" Check if spatial reference of area and crater files match. """
	
	if sr.IsProjected:
		sr_name = sr.GetAttrValue('projcs')
	else:
		sr_name = sr.GetAttrValue('geogcs')
	if sr2.IsProjected:
		sr2_name = sr2.GetAttrValue('projcs')
	else:
		sr2_name = sr2.GetAttrValue('geogcs')
	if str(sr) != str(sr2):
		print "Spatial references", sr_name, "and", sr2_name, "don't match. Please use the same coordinate system for Area and Crater Shapefile. Please note that this may also happen due to different coordinate system parameters (Projection, Center of Projection, etc.). "
		exit()
	
	major_axis = sr.GetSemiMajor()
	minor_axis = sr.GetSemiMinor()
	flattening = (major_axis-minor_axis) / major_axis 
	geogr_sr = sr.CloneGeogCS()
	
	""" Define equal area projected coordinate system (LAEA) from spatial reference. """
	
	sr_eq_area_text = 'PROJCS["PROJECTED_LAMBERT_AEA",'+str(geogr_sr)+',PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["central_meridian",0.0],PARAMETER["latitude_of_origin",0.0],UNIT["Meter",1.0]]'
	sr_equal_area_projection = osr.SpatialReference()
	sr_equal_area_projection.ImportFromWkt(sr_eq_area_text)


""" Define parameters for output files. """

def define_outfile_parameters(driver, path_to_outfile, sr, geogr_sr, sr_equal_area_projection):
	global point_data_source, polygon_data_source, line_crater_polygon_data_source, layer_points, layer_polygon, layer_line_crater_polygon, outfile_path
	
	outfile_split = os.path.split(path_to_outfile)
	outfile_path = str(outfile_split[0]) + "\\"
	outfile_name = str(outfile_split[1])
	final_outfile = str(outfile_path) + "Geodesic_Buffer_" + str(outfile_name)
	
	if generate_polygon_file == True:
		polygon_data_source = driver.CreateDataSource(outfile_path + outfile_name)
		layer_polygon = polygon_data_source.CreateLayer('Buffer_Polygon', sr, geom_type = ogr.wkbPolygon)
		
	if generate_connectors_crater_polygon == True:
		line_crater_polygon_data_source = driver.CreateDataSource(outfile_path + 'Line_Crater_Polygon.shp')
		layer_line_crater_polygon = line_crater_polygon_data_source.CreateLayer('line_crater_polygon', geogr_sr, geom_type = ogr.wkbLineString)
		idField = ogr.FieldDefn('Area', ogr.OFTString)
		layer_line_crater_polygon.CreateField(idField)
		idField = ogr.FieldDefn('Dist2Dgeo', ogr.OFTReal)
		layer_line_crater_polygon.CreateField(idField)
		idField = ogr.FieldDefn('Dist3Dmet', ogr.OFTReal)
		layer_line_crater_polygon.CreateField(idField)
	
	""" Define fields. """
	
	if generate_point_file == True:
		idField = ogr.FieldDefn('prev', ogr.OFTString) 
		layer_points.CreateField(idField)
		idField = ogr.FieldDefn('next', ogr.OFTString)
		layer_points.CreateField(idField)
		idField = ogr.FieldDefn('angle', ogr.OFTString)
		layer_points.CreateField(idField)
	
	if generate_polygon_file == True:
		if approach == "BCC":
			idField = ogr.FieldDefn('Crater_ID', ogr.OFTString)
			layer_polygon.CreateField(idField)
			idField = ogr.FieldDefn('Size_sq_km', ogr.OFTReal)
			layer_polygon.CreateField(idField)
			idField = ogr.FieldDefn('Buffer_m', ogr.OFTReal)
			layer_polygon.CreateField(idField)
		
		if approach == "NSC" or approach == "BNSC":
			idField = ogr.FieldDefn('Size_sq_km', ogr.OFTReal)
			layer_polygon.CreateField(idField)
			idField = ogr.FieldDefn('crater_ID', ogr.OFTInteger)
			layer_polygon.CreateField(idField)

""" Definition for reference system transformation from PROJCS to GEOGCS and vice versa (definition used from projected reference in .prj file) 
in order to get geographic coordinates (mandatory for geodesic buffer calculation). """

def define_coordinate_transformations(sr, geogr_sr):
	global proj_to_geog, geog_to_proj, proj_to_proj, proj_to_eq_area_proj, geog_to_eq_area_proj 
	
	proj_to_geog = osr.CoordinateTransformation(sr, geogr_sr)
	geog_to_proj = osr.CoordinateTransformation(geogr_sr, sr)
	
	proj_to_eq_area_proj = osr.CoordinateTransformation(sr, sr_equal_area_projection)
	geog_to_eq_area_proj = osr.CoordinateTransformation(geogr_sr, sr_equal_area_projection)
	

""" Calculate geodesic buffer points around Polygons according to Vincenty's direct formula: Calculate coordinates of point 2 from point 1, 
azimuth and distance Point1-Point2. """

def direct_vincenty(flattening, major_axis, vertices_angle_list, dist_buffer, bufferfactor): 
	global buffer_vertices_list
	buffer_vertices_list = []
	for pair in vertices_angle_list:
		
		""" Calculation of Point 1 (phi) on auxiliary sphere, azimuth of geodesic at equator and length of geodesic between 
		equator and Point 1. """
		
		f = flattening
		a = major_axis
		phi1 = pair[1]
		lambda1 = pair[0]
		alpha12 = pair[2]
		s = bufferfactor * dist_buffer
		piD4 = math.atan(1.0) 
		two_pi = piD4 * 8.0 
		phi1 = phi1 * piD4 / 45.0 
		lambda1 = lambda1 * piD4 / 45.0 
		alpha12 = alpha12 * piD4 / 45.0 
		if alpha12 < 0.0: 
			alpha12 = alpha12 + two_pi 
		if alpha12 > two_pi: 
			alpha12 = alpha12 - two_pi
		b = a * (1.0 - f) 
		tanU1 = (1-f) * math.tan(phi1) 
		U1 = math.atan(tanU1) 
		sigma1 = math.atan2(tanU1, math.cos(alpha12)) 
		Sinalpha0 = math.cos(U1) * math.sin(alpha12) 
		cosalpha0_sq = 1.0 - Sinalpha0 * Sinalpha0 
		u_sq = cosalpha0_sq * (a * a - b * b ) / (b * b) 
		A = 1.0 + (u_sq / 16384) * (4096 + u_sq * (-768 + u_sq * \
			(320 - 175 * u_sq))) 
		B = (u_sq / 1024) * (256 + u_sq * (-128 + u_sq * (74 - 47 * u_sq))) 
		sigma = (s / (b * A)) 
		last_sigma = -9999999.9
		
		""" Approximation of sigma (distance Point1-Point2 on auxiliary sphere). """
		
		while abs(last_sigma - sigma) > 1.0e-12:
			two_sigma_m = 2 * sigma1 + sigma 
			delta_sigma = B * math.sin(sigma) * (math.cos(two_sigma_m) + (B/4) * (math.cos(sigma) * \
				(-1 + 2 * math.pow(math.cos(two_sigma_m), 2) - (B/6) * math.cos(two_sigma_m) * \
				(-3 + 4 * math.pow(math.sin(sigma), 2 )) * (-3 + 4 * math.pow(math.cos(two_sigma_m), 2)))))
			last_sigma = sigma 
			sigma = (s / (b * A)) + delta_sigma 
		
		""" Calculation of Point 2 coordinates on ellipsoid. """
		
		phi2 = math.atan2 ((math.sin(U1) * math.cos(sigma) + math.cos(U1) * math.sin(sigma) * math.cos(alpha12)), \
			((1-f) * math.sqrt(math.pow(Sinalpha0, 2) + pow(math.sin(U1) * math.sin(sigma) - math.cos(U1) * \
			math.cos(sigma) * math.cos(alpha12), 2))))
		lambda_new = math.atan2((math.sin(sigma) * math.sin(alpha12)), (math.cos(U1) * math.cos(sigma) -  \
			math.sin(U1) *  math.sin(sigma) * math.cos(alpha12))) 
		C = (f/16) * cosalpha0_sq * (4 + f * (4 - 3 * cosalpha0_sq)) 
		L = lambda_new - (1-C) * f * Sinalpha0 * (sigma + C * math.sin(sigma) * (math.cos(two_sigma_m) + \
			C * math.cos(sigma) * (-1 + 2 * math.pow(math.cos(two_sigma_m), 2)))) 
		lambda2 = lambda1 + L 
		alpha21 = math.atan2 (Sinalpha0, (-math.sin(U1) * math.sin(sigma) + math.cos(U1) * math.cos(sigma) * math.cos(alpha12))) 
		alpha21 = alpha21 + two_pi / 2.0 # backwards azimuth
		
		if alpha21 < 0.0: 
			alpha21 = alpha21 + two_pi 
		if alpha21 > two_pi: 
			alpha21 = alpha21 - two_pi 
		
		phi2 = phi2 * 45.0 / piD4 
		lambda2 = lambda2 * 45.0 / piD4 
		alpha21 = alpha21 * 45.0 / piD4 
		
		buffer_vertices_list.append([lambda2, phi2, pair[3], pair[4]])
		
""" Calculate geodesic distance and azimuth between two points (Vincenty's inverse formula). """
	
def inverse_vincenty(flattening, major_axis, phi1, lambda1, phi2, lambda2 ):
	global geodesic_distance_crater_area, geodesic_distance_crater_area2, direction12
	
	""" Calculation of Points 1 and 2 (phi) on auxiliary sphere and difference in latitude. """
	
	a = major_axis
	f = flattening
	piD4 = math.atan(1.0)
	two_pi = piD4 * 8.0
	phi1 = phi1 * piD4 / 45.0
	lambda1 = lambda1 * piD4 / 45.0
	phi2 = phi2 * piD4 / 45.0
	lambda2 = lambda2 * piD4 / 45.0
	
	b = a * (1.0 - f)
	TanU1 = (1-f) * math.tan(phi1)
	TanU2 = (1-f) * math.tan(phi2)
	U1 = math.atan(TanU1)
	U2 = math.atan(TanU2)
	L = lambda2 - lambda1
	lambda_new = L
	last_lembda = -9999999.9

	""" Approximation of lambda_new (difference in longitude on auxiliary sphere, sigma (distance Point1Point2 on auxiliary sphere) and alpha0 
	(azimuth of geodesic at the equator). """
	
	while abs(last_lembda - lambda_new) > 1.0e-12:
		sqr_sin_sigma = pow(math.cos(U2) * math.sin(lambda_new), 2) + pow((math.cos(U1) * math.sin(U2) - \
			math.sin(U1) *  math.cos(U2) * math.cos(lambda_new)), 2)
		Sin_sigma = math.sqrt(sqr_sin_sigma)
		Cos_sigma = math.sin(U1) * math.sin(U2) + math.cos(U1) * math.cos(U2) * math.cos(lambda_new)
		sigma = math.atan2(Sin_sigma, Cos_sigma)
		if sigma == 0:
			sigma += 0.0000000001 # to avoid division by zero with Sin_alpha0 computation
		Sin_alpha0 = math.cos(U1) * math.cos(U2) * math.sin(lambda_new) / math.sin(sigma)
		alpha0 = math.asin(Sin_alpha0)
		Cos2sigma_m = math.cos(sigma) - (2 * math.sin(U1) * math.sin(U2) / pow(math.cos(alpha0), 2))
		C = (f/16) * pow(math.cos(alpha0), 2) * (4 + f * (4 - 3 * pow(math.cos(alpha0), 2)))
		last_lembda = lambda_new
		lambda_new = L + (1-C) * f * math.sin(alpha0) * (sigma + C * math.sin(sigma) * \
			(Cos2sigma_m + C * math.cos(sigma) * (-1 + 2 * pow(Cos2sigma_m, 2))))
	
	u2 = pow(math.cos(alpha0), 2) * (a * a - b * b) / (b * b)
	A = 1 + (u2/16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
	B = (u2/1024) * (256 + u2 * (-128+ u2 * (74 - 47 * u2)))
	delta_sigma = B * Sin_sigma * (Cos2sigma_m + (B/4) * (Cos_sigma * (-1 + 2 * pow(Cos2sigma_m, 2)) - \
		(B/6) * Cos2sigma_m * (-3 + 4 * sqr_sin_sigma) * (-3 + 4 * pow(Cos2sigma_m, 2))))
	
	""" Calculation of distance and azimuth on ellipsoid. """
	
	s = b * A * (sigma - delta_sigma)
	alpha12 = math.atan2((math.cos(U2) * math.sin(lambda_new)), \
		(math.cos(U1) * math.sin(U2) - math.sin(U1) * math.cos(U2) * math.cos(lambda_new)))
	alpha21 = math.atan2((math.cos(U1) * math.sin(lambda_new)), \
		(-math.sin(U1) * math.cos(U2) + math.cos(U1) * math.sin(U2) * math.cos(lambda_new)))

	if alpha12 < 0.0: 
		alpha12 =  alpha12 + two_pi
	if alpha12 > two_pi: 
		alpha12 = alpha12 - two_pi
		
	alpha21 = alpha21 + two_pi / 2.0 # backwards azimuth
	if alpha21 < 0.0: 
		alpha21 = alpha21 + two_pi
	if alpha21 > two_pi: 
		alpha21 = alpha21 - two_pi

	alpha12 = alpha12 * 45.0 / piD4
	alpha21 = alpha21 * 45.0 / piD4
	geodesic_distance_crater_area = s
	geodesic_distance_crater_area2 = s
	direction12 = alpha12
	
""" Iterate crater features: add centroid coordinates and Diameter(km) to list. """

def read_crater_features(layer_crater, proj_to_geog, geog_to_eq_area_proj, write_logfile, logfile):
	global crater_features_list, no_of_crater_features
	
	n = 0
	no_of_crater_features =  layer_crater.GetFeatureCount()
	crater_features_list = []
	print "Number of digitized craters: ", no_of_crater_features, "\n_____\n"
	if write_logfile == True:
		logfile.write("Number of digitized craters: " + str(no_of_crater_features) + "\n_____\n\n")
	
	for crater in layer_crater:
		crater_feature = layer_crater.GetFeature(n)
		
		""" Get information from attribute table when digitization was done in CraterTools. """
		
		#Crater_Diam = crater_feature.GetField("Diam_km") 
		#Crater_X = crater_feature.GetField("x_coord")
		#Crater_Y = crater_feature.GetField("y_coord")
		
		""" Autodetect crater centroid coordinates and crater diameter from crater polygon size (equal area projection). """
		
		crater2 = crater.GetGeometryRef()
		
		if sr2.IsProjected:
			crater2.Transform(proj_to_geog)
			crater2_centroid = crater2.Centroid()
			crater2_centroid_X, crater2_centroid_Y, crater2_centroid_Z = crater2_centroid.GetPoint()
			Crater_X = crater2_centroid_X
			Crater_Y = crater2_centroid_Y
			crater2.Transform(geog_to_eq_area_proj)
		else:
			crater2_centroid = crater2.Centroid()
			crater2_centroid_X, crater2_centroid_Y, crater2_centroid_Z = crater2_centroid.GetPoint()
			Crater_X = crater2_centroid_X
			Crater_Y = crater2_centroid_Y
			crater2.Transform(geog_to_eq_area_proj)
		
		Crater_Diam = (2*(math.sqrt(crater2.GetArea()/math.pi)))/1000
		
		crater_features_list.append([n, Crater_Diam, Crater_X, Crater_Y])
		n += 1
	
""" Get craters inside reference area for traditional crater counting. """

def get_craters_for_trad(out_q, crater_features_list_part, geom, generate_connectors_crater_polygon, lock, sr_wkt): 
	from shapely.geometry import Point, LineString
	global craters_for_counting_list
	
	craters_for_counting_list = []
	
	""" Convert spatial reference from wkt to spatial reference type (during multiprocessing). Also, Transformations etc. need to be defined since spatial reference data type can't
	be handed to function during multiprocessing on Windows OS (forking GDAL objects). """
	
	if isinstance(sr_wkt, basestring) == True:
		sr = osr.SpatialReference()
		sr.ImportFromWkt(sr_wkt)
		
	if isinstance(sr_wkt, basestring) == False:
		sr = sr_wkt

	""" Transform vertices to geographic coordinates if spatial reference is projected. This doesn't apply when area_polygon is already geographic. 
	Reprojecting would lead to false coordinates otherwise. Grographic area_polygon and projected sr can occur when (only) one research area is in 
	union_polygon and that research area has at least one hole. The error occurs with iterating the first inner ring. 
	
	Also, a geographic and a projected reference system from the input spatial reference is defined. """
	
	if sr.IsProjected():
		geogr_sr = sr.CloneGeogCS()
		proj_to_geog = osr.CoordinateTransformation(sr, geogr_sr)
		geog_to_proj = osr.CoordinateTransformation(geogr_sr, sr)
		geom_sr = geom.GetSpatialReference()
	
	if sr.IsGeographic():
		geogr_sr = sr
		proj_to_geog = osr.CoordinateTransformation(sr, geogr_sr)
		geog_to_proj = osr.CoordinateTransformation(geogr_sr, sr)
	
	""" Get center of research area to define projection center (if geom directly intersects the dateline, the center is not in the dateline area, 
	but somewhere between min/max longitude on the map. It's not very nice but as long as the central meridian is not at 0 deg lon in that case, 
	everything is fine. """
	
	geom_centroid = geom.Centroid()
	geom_centroid_X, geom_centroid_Y, geom_centroid_Z = geom_centroid.GetPoint()
	projection_center_X = round(geom_centroid_X, 1)
	projection_center_Y = round(geom_centroid_Y, 1)

	""" Change reference meridian in geographic coordinate system to center of input polygon to avoid problems during dateline intersection. This new geographic 
	coordinate system is used for an equal area projection which consists of the geographic reference with new reference meridian and the parameters from the sr 
	projected reference system. """
	
	geogr_sr_reprojection_text = re.sub('(PRIMEM)(.*)(,)', r'\1["Reference_Meridian",' + str(projection_center_X) + '],', str(geogr_sr))
	geogr_sr_reprojection = osr.SpatialReference(geogr_sr_reprojection_text)
	
	sr_reprojection_text = 'PROJCS["MOLLWEIDE_REPROJECTION_EQ_AREA",'+str(geogr_sr_reprojection_text)+',PROJECTION["Mollweide"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["central_meridian",' + str(projection_center_X) + '],PARAMETER["latitude_of_origin",0.0],UNIT["Meter",1.0]]'
	sr_reprojection = osr.SpatialReference(sr_reprojection_text)
	
	""" Define reprojections """
	
	geogr_to_geogr_reprojection = osr.CoordinateTransformation(geogr_sr, geogr_sr_reprojection)
	geogr_to_proj_reprojection = osr.CoordinateTransformation(geogr_sr, sr_reprojection)
	proj_reprojection_to_geogr_reprojection = osr.CoordinateTransformation(sr_reprojection, geogr_sr_reprojection)
	geogr_reprojection_to_geogr = osr.CoordinateTransformation(geogr_sr_reprojection, geogr_sr)
	geogr_reprojection_to_eq_area = osr.CoordinateTransformation(geogr_sr_reprojection, sr_reprojection)
	eq_area_proj_to_proj = osr.CoordinateTransformation(sr_reprojection, sr)
	geogr_reprojection_to_geogr = osr.CoordinateTransformation(geogr_sr_reprojection, geogr_sr)
	geogr_reprojection_to_proj_reprojection = osr.CoordinateTransformation(geogr_sr_reprojection, sr_reprojection)
	
	""" Reproject input geometry """
	
	geom.Transform(geogr_to_proj_reprojection)
	geom.Transform(proj_reprojection_to_geogr_reprojection)
	
	if sr.IsProjected():
		geom.Transform(geogr_reprojection_to_proj_reprojection)
		geom.Transform(eq_area_proj_to_proj)
	
	for crater_centroid in crater_features_list_part:
		
		crater_centroid_geometry = Point(crater_centroid[2], crater_centroid[3])
		crater_centroid_geometry_2 = ogr.CreateGeometryFromWkt(crater_centroid_geometry.wkt)
		
		""" Reproject craters. First to projected, then to geographic system. This is mandatory, as the generation of geographic coordinates from projected coordinates 
		considers the pure distance from the central meridian. This avoids misinterpretations of nearest neighbors by (lon 175deg - lon -175deg) distances. 
		Accordingly, dateline intersections will be considered during processing. """
		
		crater_centroid_geometry_2.Transform(geogr_to_proj_reprojection)
		crater_centroid_geometry_2.Transform(proj_reprojection_to_geogr_reprojection)
		
		if sr.IsProjected():
			crater_centroid_geometry_2.Transform(geogr_reprojection_to_proj_reprojection) 
			crater_centroid_geometry_2.Transform(eq_area_proj_to_proj)
		
		distance_crater_polygon = crater_centroid_geometry_2.Distance(geom)
		
		""" Get craters inside research area """
		
		if distance_crater_polygon <= 0:
			craters_for_counting_list.append(crater_centroid)
		
	""" Pass list """
	
	out_q.put(craters_for_counting_list)

""" Get distance between crater centroid and area polygon and determine which craters are inside bufferfactor * crater radius range without using 
ogr.intersection method (time consuming) """

def get_craters_for_buffering_BCC(out_q, out_q_2, out_q_3, crater_features_list_part, geom, vertices_list, flattening, major_axis, bufferfactor, generate_connectors_crater_polygon, lock, sr_wkt): 
	from shapely.geometry import Point, LineString
	global craters_inside_area, craters_outside_area, craters_within_range, craters_for_counting_list
	
	craters_inside_area = []
	craters_outside_area = []
	craters_within_range = []
	st = time.time()
	
	""" Convert spatial reference from wkt to spatial reference type (during multiprocessing). Also, Transformations etc. need to be defined as Spatial reference data type can't
	be handed to function during multiprocessing. """
	
	if isinstance(sr_wkt, basestring) == True:
		sr = osr.SpatialReference()
		sr.ImportFromWkt(sr_wkt)
	
	if isinstance(sr_wkt, basestring) == False:
		sr = sr_wkt
	
	""" Transform vertices to geographic coordinates if spatial reference is projected. This doesn't apply when area_polygon is already geographic. 
	Reprojecting would lead to false coordinates otherwise. Grographic area_polygon and Projected sr can occur when (only) one research area is in 
	union_polygon and that research area has at least one hole. The error occurs with iterating the first inner ring. 
	
	Also, a geographic and a projected reference system from the input spatial reference is defined. """
	
	if sr.IsProjected():
		geogr_sr = sr.CloneGeogCS()
		proj_to_geog = osr.CoordinateTransformation(sr, geogr_sr)
		geog_to_proj = osr.CoordinateTransformation(geogr_sr, sr)
		geom_sr = geom.GetSpatialReference()

	if sr.IsGeographic():
		geogr_sr = sr
		proj_to_geog = osr.CoordinateTransformation(sr, geogr_sr)
		geog_to_proj = osr.CoordinateTransformation(geogr_sr, sr)
	
	""" Get center of research area to define projection center (if geom directly intersects the dateline, the center is not in the dateline area, 
	but somewhere between min/max longitude on the map. It's not very elegant but as long as the central meridian is not at 0 deg lon in that case, 
	it's all good. """
	
	geom_centroid = geom.Centroid()
	geom_centroid_X, geom_centroid_Y, geom_centroid_Z = geom_centroid.GetPoint()
	projection_center_X = round(geom_centroid_X, 1)
	projection_center_Y = round(geom_centroid_Y, 1)

	""" Change reference meridian in geographic coordinate system to center of input polygon to avoid problems during dateline intersection. This new geographic 
	coordinate system is used for an equal area projection which consists of the geographic reference with new reference meridian and the parameters from the sr 
	projected reference system. """
	
	geogr_sr_reprojection_text = re.sub('(PRIMEM)(.*)(,)', r'\1["Reference_Meridian",' + str(projection_center_X) + '],', str(geogr_sr))
	geogr_sr_reprojection = osr.SpatialReference(geogr_sr_reprojection_text)
	sr_reprojection_text = 'PROJCS["MOLLWEIDE_REPROJECTION_EQ_AREA",'+str(geogr_sr_reprojection_text)+',PROJECTION["Mollweide"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["central_meridian",' + str(projection_center_X) + '],PARAMETER["latitude_of_origin",0.0],UNIT["Meter",1.0]]'
	sr_reprojection = osr.SpatialReference(sr_reprojection_text)
	
	""" Define reprojections """
	
	geogr_to_geogr_reprojection = osr.CoordinateTransformation(geogr_sr, geogr_sr_reprojection)
	geogr_to_proj_reprojection = osr.CoordinateTransformation(geogr_sr, sr_reprojection)
	proj_reprojection_to_geogr_reprojection = osr.CoordinateTransformation(sr_reprojection, geogr_sr_reprojection)
	geogr_reprojection_to_geogr = osr.CoordinateTransformation(geogr_sr_reprojection, geogr_sr)
	geogr_reprojection_to_eq_area = osr.CoordinateTransformation(geogr_sr_reprojection, sr_reprojection)
	eq_area_proj_to_proj = osr.CoordinateTransformation(sr_reprojection, sr)
	geogr_reprojection_to_geogr = osr.CoordinateTransformation(geogr_sr_reprojection, geogr_sr)
	geogr_reprojection_to_input_sr = osr.CoordinateTransformation(geogr_sr_reprojection, sr)
	input_sr_to_geogr_reprojection = osr.CoordinateTransformation(sr, geogr_sr_reprojection)
	input_sr_to_sr_reprojection = osr.CoordinateTransformation(sr, sr_reprojection)
	geogr_reprojection_to_proj_reprojection = osr.CoordinateTransformation(geogr_sr_reprojection, sr_reprojection)
	
	geogr_to_sr = osr.CoordinateTransformation(geogr_sr, sr)
	
	geom.Transform(geogr_to_proj_reprojection)
	geom.Transform(proj_reprojection_to_geogr_reprojection)
	
	if sr.IsProjected():
		geom.Transform(geogr_reprojection_to_proj_reprojection)
		geom.Transform(eq_area_proj_to_proj)
	
	#################
	#	OUTER RING	#
	#################
	
	""" Generate list of reprojected geometry outline. This one is used to find closest point to reprojected crater centroids. """
	
	vertices_list_reprojected = []
	geom_reprojected = geom.GetGeometryRef(0)
	no_of_polygon_vertices = geom_reprojected.GetPointCount()
	for vertex in xrange(no_of_polygon_vertices):
		current_vertex_X, current_vertex_Y, z = geom_reprojected.GetPoint(vertex)
		vertices_list_reprojected.append([current_vertex_X, current_vertex_Y])
	area_linestring = LineString(list(vertices_list_reprojected))
	
	#################
	#	INNER RING	#
	#################
	
	""" Generate lists of inner ring geometry lines. They are used to find the closest point between crater centroids and polygons / inner rings.  """
	
	no_of_inner_rings = geom.GetGeometryCount() - 1
	if no_of_inner_rings > 0:
		inner_rings = dict()
		for inner_polygon_count in range(no_of_inner_rings):
			inner_polygon_count += 1
			vertices_inner_rings = []
			inner_rings[inner_polygon_count] = []
			inner_ring = geom.GetGeometryRef(inner_polygon_count)
			no_of_ring_vertices = inner_ring.GetPointCount()
			for vertex in xrange(no_of_ring_vertices):
				current_vertex_X, current_vertex_Y, z = inner_ring.GetPoint(vertex)
				vertices_inner_rings.append([current_vertex_X, current_vertex_Y])
			inner_linestring = LineString(list(vertices_inner_rings))
			inner_rings[inner_polygon_count].append([inner_linestring])
	
	for crater_centroid in crater_features_list_part:
		
		""" Exclude craters within polygon for distance determination (geodesic distance between area and crater from two point coordinates) 
		(as closest point between crater centroid and area can later only be conducted using a line, not a polygon feature 
		- difficult du destinguish between inside and outside craters) """
		
		""" crater_centroid_geometry and crater_centroid_geometry_2 are the same reprojected point geometries. crater_centroid_geometry is in 
		Shapely, crater_centroid_geometry_2 in OGR format. Both libraries can apply their functions only to their own geometry format. """
		
		crater_centroid_geometry = Point(crater_centroid[2], crater_centroid[3])
		crater_centroid_geometry_2 = ogr.CreateGeometryFromWkt(crater_centroid_geometry.wkt)
		
		""" Reproject craters. First to projected, then to geographic system. This is mandatory, as the generation of geographic coordinates from projected coordinates 
		considers the pure distance from the central meridian. This avoids misinterpretations of nearest neighbors by (lon 175deg - lon -175deg) distances. 
		Accordingly, dateline intersections will be considered during processing. """
		
		crater_centroid_geometry_2.Transform(geogr_to_proj_reprojection)
		crater_centroid_geometry_2.Transform(proj_reprojection_to_geogr_reprojection)
		
		""" Reproject data back to original spatial reference when original spatial reference is projected. The interpolate function to 
		identify the nearest point on the reference area identifies the coordinates with respect to its spatial reference. A straight line 
		between two vertices looks different in every projection. In order to identify craters with respect to the polygon boundaries as they were 
		digitized, this step needs to be conducted. """
		
		if sr.IsProjected():
			crater_centroid_geometry_2.Transform(geogr_reprojection_to_proj_reprojection)
			crater_centroid_geometry_2.Transform(eq_area_proj_to_proj)
		
		crater_centroid_geometry_reprojected = Point(crater_centroid_geometry_2.GetX(), crater_centroid_geometry_2.GetY())
		
		distance_crater_polygon = crater_centroid_geometry_2.Distance(geom)
		crater_centroid_diameter_km = crater_centroid[1]
		
		""" Distinguish between craters inside and outside research area. """
		
		""" Investigate craters outside research area. Distance > 0 would lead to craters with centroid exactly on polygon boundary to always 
		be considered closer to the inner ring than the outer boundary. That's why a little extra distance is considered (> 0.001). """

		if distance_crater_polygon > 0.001:  
			
			#################
			#	OUTER RING	#
			#################
			
			intersect_point = area_linestring.interpolate(area_linestring.project(crater_centroid_geometry_reprojected)) 
			craters_outside_area.append(crater_centroid)
			 
			""" Determine distance of intersecting point to crater using geographic coordinates - Vincenty needs them. """
			 
			intersect_point_2 = ogr.CreateGeometryFromWkt(intersect_point.wkt)

			if sr.IsProjected():
				intersect_point_2.Transform(proj_to_geog)
				intersect_point_2.Transform(geogr_to_proj_reprojection)
				intersect_point_2.Transform(proj_reprojection_to_geogr_reprojection)
				crater_centroid_geometry_2.Transform(proj_to_geog)
				crater_centroid_geometry_2.Transform(geogr_to_proj_reprojection)
				crater_centroid_geometry_2.Transform(proj_reprojection_to_geogr_reprojection)
			
			intersect_point = Point(intersect_point_2.GetX(),intersect_point_2.GetY()) 
			
			""" Calculate geodesic distance between crater centroid and intersection with research area. """
			
			lembda1 = crater_centroid_geometry_2.GetX() 
			phi1 = crater_centroid_geometry_2.GetY() 
			lembda2 = intersect_point.x
			phi2 = intersect_point.y
			
			inverse_vincenty(flattening, major_axis, phi1, lembda1, phi2, lembda2)
			geodesic_distance_crater_area = geodesic_distance_crater_area2
			
			#################
			#	INNER RING	#
			#################
			
			if no_of_inner_rings > 0:
				
				""" Copy values of distance to outer polygon boundary analysis to compare distances to the inner ring. The shortest distance is used for 
				the CSFD analysis. """
				
				distance_value_closest = geodesic_distance_crater_area
				intersect_point_closest = intersect_point
				
				for inner_polygon_count in range(no_of_inner_rings):
					inner_polygon_count += 1
					inner_ring = inner_rings[inner_polygon_count][0][0]
					
					intersect_point = inner_ring.interpolate(inner_ring.project(crater_centroid_geometry_reprojected))
					
					intersect_point_2 = ogr.CreateGeometryFromWkt(intersect_point.wkt)
					
					if sr.IsProjected():
						intersect_point_2.Transform(proj_to_geog)
					
						intersect_point_2.Transform(geogr_to_proj_reprojection) 
						intersect_point_2.Transform(proj_reprojection_to_geogr_reprojection)
						intersect_point = Point(intersect_point_2.GetX(),intersect_point_2.GetY())
					
					lembda2 = intersect_point.x
					phi2 = intersect_point.y
					
					inverse_vincenty(flattening, major_axis, phi1, lembda1, phi2, lembda2)
					
					""" Save intersect point for geometry line generation and overwrite distance value if crater is closer to the inner ring than the outer ring. """
					
					if geodesic_distance_crater_area2 < distance_value_closest:
						intersect_point_closest = intersect_point
						distance_value_closest = geodesic_distance_crater_area2
					
				""" Get final intersect point and distance """
				
				intersect_point = intersect_point_closest
				geodesic_distance_crater_area = distance_value_closest
			
			""" Check if geodesic distance between crater centroid and area is smaller than / equal to craterradius * bufferfactor. 
			If so, craters are considered relevant for buffered crater counting. """
			
			if geodesic_distance_crater_area <= crater_centroid_diameter_km / 2 * 1000 * bufferfactor:
				craters_within_range.append(crater_centroid)
			
			""" Append coordinates to craters within range list """
			
			if generate_connectors_crater_polygon == True:
				
				""" Transform intersect_point to OGR format - which is intersect_point_2 - to reproject data back to original geographic coordinates. 
				Transform intersect_point_2 back to Shapely format - which is intersect_point_orig_geogr. Transformation can only be done in OGR.  """
				
				intersect_point_2 = ogr.CreateGeometryFromWkt(intersect_point.wkt)
				intersect_point_2.Transform(geogr_reprojection_to_geogr) 
				intersect_point_orig_geogr = Point(intersect_point_2.GetX(), intersect_point_2.GetY())
				
				line_crater_polygon = LineString([Point(crater_centroid[2], crater_centroid[3]), intersect_point_orig_geogr])
				line_crater_polygon_layer = ogr.CreateGeometryFromWkt(line_crater_polygon.wkt)
				ll_featureDefn = layer_line_crater_polygon.GetLayerDefn()
				ll_feature = ogr.Feature(ll_featureDefn)
				ll_feature.SetGeometry(line_crater_polygon_layer) 
				ll_feature.SetField('Area', Area_ID)
				ll_feature.SetField('Dist2Dgeo', distance_crater_polygon)
				ll_feature.SetField('Dist3Dmet', geodesic_distance_crater_area)
				layer_line_crater_polygon.CreateFeature(ll_feature)
		
		else: # craters inside research area
			craters_inside_area.append(crater_centroid)
	
	""" Merge lists of craters within research area and craters within buffer range of research area. """
	
	craters_for_counting_list = craters_inside_area + craters_within_range
	
	""" Pass list """
	
	out_q.put(craters_for_counting_list)
	out_q_2.put(len(craters_within_range))
	out_q_3.put(len(craters_inside_area))

def get_craters_for_buffering_NSC_BNSC(approach, out_q, out_q_2, out_q_3, out_q_4, crater_features_list_part, geom, vertices_list, flattening, major_axis, bufferfactor, bufferfactor_crater, generate_connectors_crater_polygon, lock, sr_wkt): 
	from shapely.geometry import Point, LineString
	global craters_inside_area, craters_outside_area, craters_within_range, craters_for_counting_list, craters_excluded_from_list
	
	craters_inside_area = []
	craters_outside_area = []
	craters_within_range = []
	st = time.time()
	
	""" Reprojection parameters to consider dateline and polar intersections """
	
	""" Convert spatial reference from wkt to spatial reference type (during multiprocessing). Also, Transformations etc. need to be defined as Spatial reference data type can't
	be handed to function during multiprocessing """
	
	if isinstance(sr_wkt, basestring) == True:
		sr = osr.SpatialReference()
		sr.ImportFromWkt(sr_wkt)
		
	if isinstance(sr_wkt, basestring) == False:
		sr = sr_wkt
	
	""" Transform vertices to geographic coordinates if spatial reference is projected. This doesn't apply when area_polygon is already geographic. 
	Reprojecting would lead to false coordinates otherwise. Grographic area_polygon and Projected sr can occur when (only) one research area is in 
	union_polygon and that research area has at least one hole. The error occurs with iterating the first inner ring. 
	
	Also, a geographic and a projected reference system from the input spatial reference is defined. """
	
	if sr.IsProjected():
		geogr_sr = sr.CloneGeogCS()
		proj_to_geog = osr.CoordinateTransformation(sr, geogr_sr)
		geog_to_proj = osr.CoordinateTransformation(geogr_sr, sr)
		geom_sr = geom.GetSpatialReference()
	if sr.IsGeographic():
		geogr_sr = sr
		proj_to_geog = osr.CoordinateTransformation(sr, geogr_sr)
		geog_to_proj = osr.CoordinateTransformation(geogr_sr, sr)
	
	""" Get center of research area to define projection center (if geom directly intersects the dateline, the center is not in the dateline area, 
	but somewhere between min/max longitude on the map. It's not very elegant but as long as the central meridian is not at 0 deg lon in that case, 
	it's all good. """
	
	geom_centroid = geom.Centroid()
	geom_centroid_X, geom_centroid_Y, geom_centroid_Z = geom_centroid.GetPoint()
	projection_center_X = round(geom_centroid_X, 1)
	projection_center_Y = round(geom_centroid_Y, 1)

	""" Change reference meridian in geographic coordinate system to center of input polygon to avoid problems during dateline intersection. This new geographic 
	coordinate system is used for an equal area projection which consists of the geographic reference with new reference meridian and the parameters from the sr 
	projected reference system. """
	
	geogr_sr_reprojection_text = re.sub('(PRIMEM)(.*)(,)', r'\1["Reference_Meridian",' + str(projection_center_X) + '],', str(geogr_sr))
	geogr_sr_reprojection = osr.SpatialReference(geogr_sr_reprojection_text)

	sr_reprojection_text = 'PROJCS["MOLLWEIDE_REPROJECTION_EQ_AREA",'+str(geogr_sr_reprojection_text)+',PROJECTION["Mollweide"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["central_meridian",' + str(projection_center_X) + '],PARAMETER["latitude_of_origin",0.0],UNIT["Meter",1.0]]'
	sr_reprojection = osr.SpatialReference(sr_reprojection_text)
	
	""" define reprojections """
	
	geogr_to_geogr_reprojection = osr.CoordinateTransformation(geogr_sr, geogr_sr_reprojection)
	geogr_to_proj_reprojection = osr.CoordinateTransformation(geogr_sr, sr_reprojection)
	proj_reprojection_to_geogr_reprojection = osr.CoordinateTransformation(sr_reprojection, geogr_sr_reprojection)
	geogr_reprojection_to_geogr = osr.CoordinateTransformation(geogr_sr_reprojection, geogr_sr)
	geogr_reprojection_to_eq_area = osr.CoordinateTransformation(geogr_sr_reprojection, sr_reprojection)
	eq_area_proj_to_proj = osr.CoordinateTransformation(sr_reprojection, sr)
	geogr_reprojection_to_geogr = osr.CoordinateTransformation(geogr_sr_reprojection, geogr_sr)
	geogr_reprojection_to_input_sr = osr.CoordinateTransformation(geogr_sr_reprojection, sr)
	input_sr_to_geogr_reprojection = osr.CoordinateTransformation(sr, geogr_sr_reprojection)
	input_sr_to_sr_reprojection = osr.CoordinateTransformation(sr, sr_reprojection)
	geogr_reprojection_to_proj_reprojection = osr.CoordinateTransformation(geogr_sr_reprojection, sr_reprojection)
	geogr_to_sr = osr.CoordinateTransformation(geogr_sr, sr)
	
	""" Reproject input geometry """
	
	geom.Transform(geogr_to_proj_reprojection)
	geom.Transform(proj_reprojection_to_geogr_reprojection)
	
	if sr.IsProjected():
		geom.Transform(geogr_reprojection_to_proj_reprojection) 
		geom.Transform(eq_area_proj_to_proj)

	#################
	#	OUTER RING	#
	#################
	
	""" Generate list of reprojected geometry outline. This one is used to find closest point to reprojected crater centroids. """
	
	vertices_list_reprojected = []
	geom_reprojected = geom.GetGeometryRef(0)
	no_of_polygon_vertices = geom_reprojected.GetPointCount()
	for vertex in xrange(no_of_polygon_vertices):
		current_vertex_X, current_vertex_Y, z = geom_reprojected.GetPoint(vertex)
		vertices_list_reprojected.append([current_vertex_X, current_vertex_Y])
	area_linestring = LineString(list(vertices_list_reprojected))
	
	#################
	#	INNER RING	#
	#################
	
	""" Generate lists of inner ring geometry lines. They are used to find the closest point between crater centroids and polygons / inner rings. """
	
	no_of_inner_rings = geom.GetGeometryCount() - 1
	if no_of_inner_rings > 0:
		inner_rings = dict()
		for inner_polygon_count in range(no_of_inner_rings):
			inner_polygon_count += 1
			vertices_inner_rings = []
			inner_rings[inner_polygon_count] = []
			inner_ring = geom.GetGeometryRef(inner_polygon_count)
			no_of_ring_vertices = inner_ring.GetPointCount()
			for vertex in xrange(no_of_ring_vertices):
				current_vertex_X, current_vertex_Y, z = inner_ring.GetPoint(vertex)
				vertices_inner_rings.append([current_vertex_X, current_vertex_Y])
			inner_linestring = LineString(list(vertices_inner_rings))
			inner_rings[inner_polygon_count].append([inner_linestring])
	
	for crater_centroid in crater_features_list_part:
		
		""" Create copy of crater centroid for later use. Distance between crater and area is added. If more than one reference area is investigated, 
		distances would sum up otherwise. """
		
		crater_centroid_2 = []
		crater_centroid_2.append(crater_centroid[0])
		crater_centroid_2.append(crater_centroid[1])
		crater_centroid_2.append(crater_centroid[2])
		crater_centroid_2.append(crater_centroid[3])
		
		""" Exclude craters within polygon for distance determination (geodesic distance between area and crater from two point coordinates) 
		(as closest point between crater centroid and area can later only be conducted using a line, not a polygon feature 
		- difficult du destinguish between inside and outside craters) """
		
		""" crater_centroid_geometry and crater_centroid_geometry_2 are the same reprojected point geometries. crater_centroid_geometry is in 
		Shapely, crater_centroid_geometry_2 in OGR format. Both libraries can apply their functions only to their own geometry format. """
		
		crater_centroid_geometry = Point(crater_centroid[2], crater_centroid[3])
		crater_centroid_geometry_2 = ogr.CreateGeometryFromWkt(crater_centroid_geometry.wkt)
		
		""" Reproject craters. First to projected, then to geographic. This is mandatory as the generation of geographic coordinates from projected coordinates 
		considers the pure distance from the central meridian. This avoids misinterpretations of nearest neighbors by (lon 175deg - lon -175deg) distances. 
		Accordingly, dateline intersections will be considered during processing. """
		
		crater_centroid_geometry_2.Transform(geogr_to_proj_reprojection)
		crater_centroid_geometry_2.Transform(proj_reprojection_to_geogr_reprojection)
		
		""" Reproject data back to original spatial reference when original spatial reference is projected. The interpolate function to 
		identify the nearest point on the reference area identifies the coordinates with respect to its spatial reference. A straight line 
		between two vertices looks different in every projection. In order to identify craters with respect to the polygon boundaries as they were 
		digitized, this step needs to be conducted. """
		
		if sr.IsProjected():
			crater_centroid_geometry_2.Transform(geogr_reprojection_to_proj_reprojection)
			crater_centroid_geometry_2.Transform(eq_area_proj_to_proj)
		
		crater_centroid_geometry_reprojected = Point(crater_centroid_geometry_2.GetX(), crater_centroid_geometry_2.GetY())
		
		distance_crater_polygon = crater_centroid_geometry_2.Distance(geom)
		crater_centroid_diameter_km = crater_centroid[1]
		
		""" Distinguish between craters inside and outside research area. """
		
		""" Investigate craters outside research area. Distance > 0 would lead to craters with centroid exactly on polygon boundary to always 
		be considered closer to inner ring than outer boundary. That's why a little extra distance is considered (> 0.001). """

		if distance_crater_polygon > 0.001:  
			
			#################
			#	OUTER RING	#
			#################
			
			intersect_point = area_linestring.interpolate(area_linestring.project(crater_centroid_geometry_reprojected)) 
			craters_outside_area.append(crater_centroid)
			 
			""" Determine distance of intersecting point to crater using geographic coordinates - Vincenty needs them. """
			 
			intersect_point_2 = ogr.CreateGeometryFromWkt(intersect_point.wkt)
			
			if sr.IsProjected():
				intersect_point_2.Transform(proj_to_geog)
				
				intersect_point_2.Transform(geogr_to_proj_reprojection)
				intersect_point_2.Transform(proj_reprojection_to_geogr_reprojection)
				
				crater_centroid_geometry_2.Transform(proj_to_geog)
				crater_centroid_geometry_2.Transform(geogr_to_proj_reprojection)
				crater_centroid_geometry_2.Transform(proj_reprojection_to_geogr_reprojection)
			
				intersect_point = Point(intersect_point_2.GetX(),intersect_point_2.GetY())
			
			""" Calculate geodesic distance between crater centroid and intersection with research area. """
			
			lembda1 = crater_centroid_geometry_2.GetX() 
			phi1 = crater_centroid_geometry_2.GetY() 
			lembda2 = intersect_point.x
			phi2 = intersect_point.y
			
			inverse_vincenty(flattening, major_axis, phi1, lembda1, phi2, lembda2)
			geodesic_distance_crater_area = geodesic_distance_crater_area2
			
			#################
			#	INNER RING	#
			#################
			
			if no_of_inner_rings > 0:
				
				""" Copy values of distance to outer polygon boundary analysis to compare distances to inner ring. The shortest distance is used for the
				CSFD measurement. """
				
				distance_value_closest = geodesic_distance_crater_area
				intersect_point_closest = intersect_point
				
				for inner_polygon_count in range(no_of_inner_rings):
					inner_polygon_count += 1
					inner_ring = inner_rings[inner_polygon_count][0][0]
					
					intersect_point = inner_ring.interpolate(inner_ring.project(crater_centroid_geometry_reprojected))
					intersect_point_2 = ogr.CreateGeometryFromWkt(intersect_point.wkt)
					
					if sr.IsProjected():
						intersect_point_2.Transform(proj_to_geog)
						intersect_point_2.Transform(geogr_to_proj_reprojection)
						intersect_point_2.Transform(proj_reprojection_to_geogr_reprojection)
						intersect_point = Point(intersect_point_2.GetX(),intersect_point_2.GetY())
					
					lembda2 = intersect_point.x
					phi2 = intersect_point.y
					
					inverse_vincenty(flattening, major_axis, phi1, lembda1, phi2, lembda2)
					
					""" Save intersect point for geometry line generation and overwrite distance value if crater is closer to the inner ring than the outer ring. """
					
					if geodesic_distance_crater_area2 < distance_value_closest:
						intersect_point_closest = intersect_point
						distance_value_closest = geodesic_distance_crater_area2
					
				""" Get final intersect point and distance. """
				
				intersect_point = intersect_point_closest
				geodesic_distance_crater_area = distance_value_closest
			
			crater_centroid_2.append(geodesic_distance_crater_area)
			
			""" Check if geodesic distance between crater centroid and area is smaller than / equal to craterradius * bufferfactor (BNSC) 
			or craterradius * ( bufferfactor + bufferfactor_crater ) (NSC). If so, craters are considered relevant for buffered crater counting. 
			bufferfactor_crater is considererd at a later stage for BNSC. During NSC, it is needed to buffer craters outside of the 
			reference area and to erase parts of the area even though such craters are not included in the actual counting. """
			
			if approach == "NSC":
				
				if geodesic_distance_crater_area <= crater_centroid_diameter_km / 2 * 1000 * (bufferfactor_crater + 1):
					craters_within_range.append(crater_centroid_2)
			
			if approach == "BNSC":
				if geodesic_distance_crater_area <= crater_centroid_diameter_km / 2 * 1000 * (bufferfactor):
					craters_within_range.append(crater_centroid_2)
			
			""" Append coordinates to craters within range list. """
			
			if generate_connectors_crater_polygon == True:
				
				""" Transform intersect_point to OGR format - which is intersect_point_2 - to reproject data back to original geographic coordinates. 
				Transform intersect_point_2 back to Shapely format - which is intersect_point_orig_geogr. Transformation can only be done in OGR.  """
				
				intersect_point_2 = ogr.CreateGeometryFromWkt(intersect_point.wkt)
				intersect_point_2.Transform(geogr_reprojection_to_geogr) 
				intersect_point_orig_geogr = Point(intersect_point_2.GetX(), intersect_point_2.GetY())
				
				line_crater_polygon = LineString([Point(crater_centroid[2], crater_centroid[3]), intersect_point_orig_geogr])
				line_crater_polygon_layer = ogr.CreateGeometryFromWkt(line_crater_polygon.wkt)
				
				ll_featureDefn = layer_line_crater_polygon.GetLayerDefn()
				ll_feature = ogr.Feature(ll_featureDefn)
				ll_feature.SetGeometry(line_crater_polygon_layer) 
				ll_feature.SetField('Area', Area_ID)
				ll_feature.SetField('Dist2Dgeo', distance_crater_polygon)
				ll_feature.SetField('Dist3Dmet', geodesic_distance_crater_area)
				layer_line_crater_polygon.CreateFeature(ll_feature)
			
		else: # craters inside research area
			
			""" Modified value of -1 for distance to polygon. """
			
			crater_centroid_2.append(-1)
			craters_inside_area.append(crater_centroid_2)
	
	""" Merge lists of craters within research area and craters within buffer range of research area. """
	
	craters_for_counting_list = craters_inside_area + craters_within_range
	
	#############################################################################################################
	#																											#
	# (B)NSC: Find craters on top of larger craters and exclude them from craters inside area and within range.	#
	#																											#
	#############################################################################################################
	
	""" Sort list. Largest crater on top. """
	
	craters_for_counting_list = sorted(craters_for_counting_list, key = operator.itemgetter(1), reverse = True)
	craters_on_resurfaced_area = []
	
	""" For crater inside area: Determine width of buffered crater. We start with the largest and exclude smaller craters from the list as they are situated on a 
	resurfaced area. """
	
	for crater in craters_for_counting_list:
		
		crater_id = crater[0]
		crater_diameter = crater[1]
		crater_centroid_X = crater[2]
		crater_centroid_Y = crater[3]
		
		if approach == "NSC":
			dist_buffer_meter = ((crater_diameter * 1000) / 2) + ((crater_diameter * 1000)/2 * (bufferfactor_crater - 1))
		
		""" For other craters in list (smaller craters, since we sorted the list before): check distance to buffered crater. """
		
		for other_crater in craters_for_counting_list[1:]: # skip the first one because we don't need the distance between the same crater
			
			other_crater_diameter = other_crater[1]
			other_crater_centroid_X = other_crater[2]
			other_crater_centroid_Y = other_crater[3]
			
			if approach == "BNSC":
				dist_buffer_meter = ((crater_diameter * 1000) / 2) + ((crater_diameter * 1000)/2 * (bufferfactor_crater - 1)) - ((other_crater_diameter * 1000)/2 * bufferfactor)
			
			""" Check if crater diameter is smaller (only smaller crater can obliterate larger crater). """
			
			if other_crater_diameter < crater_diameter:
			
				""" Get distance between crater and other crater. """
				
				inverse_vincenty(flattening, major_axis, crater_centroid_Y, crater_centroid_X, other_crater_centroid_Y, other_crater_centroid_X)
				
				""" Find craters obliterating larger crater. If geodesic distance between crater centroids is smaller than original crater diameter 
				+ surrounding buffer of bufferfactor * crater radius, then other crater is obliterating the original crater. The craters are added to the 
				craters_on_resurfaced_area list. """
				
				if geodesic_distance_crater_area2 < dist_buffer_meter:
					#print "Crater", other_crater, "excluded since it obliterates", crater
					craters_on_resurfaced_area.append(other_crater)
	
	""" Eliminate duplicates in craters_on_resurfaced_area list (craters are on top of more than one ejecta blanket) and delete crater from 
	craters_for_counting_list. Such craters are located on top of larger craters and will not be considered for CSFD analysis. """
	
	craters_on_resurfaced_area_clean = []
	for resurfaced_crater in craters_on_resurfaced_area:
		
		if resurfaced_crater in craters_for_counting_list:
			craters_for_counting_list.remove(resurfaced_crater)
		
		if resurfaced_crater not in craters_on_resurfaced_area_clean:
			
			""" Inclute identifier 'obliterates' for obliterating crater. Used later to include obliterating craters in the buffering (ejecta blanket 
			still affects reference area) but not in the CSFD measurement (crater is located on a resurfaced area). """
			
			resurfaced_crater.append("obliterates")
			craters_on_resurfaced_area_clean.append(resurfaced_crater)
			
	craters_excluded_from_list = craters_on_resurfaced_area_clean
	
	""" Pass list """
	
	out_q.put(craters_for_counting_list)
	out_q_2.put(craters_within_range)
	out_q_3.put(craters_inside_area)
	out_q_4.put(craters_excluded_from_list)

""" Get area size from Lambert azimuthal equal area projection. """

def get_area_size(union_polygon):
	global Area_Size
	
	geog_to_proj = osr.CoordinateTransformation(geogr_sr, sr)
	proj_to_geogr = osr.CoordinateTransformation(sr, geogr_sr)
	
	""" Transform vertices to geographic coordinates if spatial reference is projected. """
	
	if sr.IsProjected():
		union_polygon_centroid.Transform(proj_to_geogr)
		union_polygon_sr = union_polygon.GetSpatialReference()
		if union_polygon_sr:
			if union_polygon_sr.IsProjected():
				union_polygon.Transform(proj_to_geogr)
		if not union_polygon_sr:
			union_polygon.Transform(proj_to_geogr)
	
	""" Get center of research area to define projection center.  """
	
	union_polygon_centroid_X, union_polygon_centroid_Y, union_polygon_centroid_Z = union_polygon_centroid.GetPoint()
	projection_center_X = round(union_polygon_centroid_X, 1)
	projection_center_Y = round(union_polygon_centroid_Y, 1)
	
	""" If a polygon intersects a dateline, it may happen that the centroid is wrongly identified (Polygon intersects dateline vs. 
	Polygon spans globe with a hole over dateline). If polygon cuts lines (which we assume to be datelines) at lon -179.8 and lon 179.8, we assume that a dateline 
	intersection is present. In this case, we set the projection center to (false projection center (if closer than 80 deg lon to central meridian) + 100 deg 
	lon). It is quite experimental to set a fixed lon value for such cases but it worked in all cases we tested. """
	
	""" Draw Date Lines """
	
	dateline_1 = ogr.Geometry(ogr.wkbLineString)
	dateline_1.AddPoint(-179, 90)
	dateline_1.AddPoint(-179, -90)
	
	dateline_2 = ogr.Geometry(ogr.wkbLineString)
	dateline_2.AddPoint(179, 90)
	dateline_2.AddPoint(179, -90)
	
	if union_polygon.Intersects(dateline_1) == True and union_polygon.Intersects(dateline_2) == True:
		if projection_center_X < 80 and projection_center_X > -80:
			projection_center_X = projection_center_X + 100
	
	sr_text = 'PROJCS["PROJECTED_LAMBERT_AEA",'+str(geogr_sr)+',PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["central_meridian",0.0],PARAMETER["latitude_of_origin",0.0],UNIT["Meter",1.0]]'
	
	if sr.IsProjected():
		union_polygon.Transform(geog_to_proj)
	
	""" Change reference meridian in geographic coordinate system to center of input polygon to avoid problems during dateline intersection. This new geographic 
	coordinate system is used for an equal area projection which consists of the geographic reference with a new reference meridian. """
	
	geogr_sr_reprojection_text = re.sub('(PRIMEM)(.*)(,)', r'\1["Reference_Meridian",' + str(projection_center_X) + '],', str(geogr_sr))
	geogr_sr_reprojection = osr.SpatialReference(geogr_sr_reprojection_text)
	
	proj_cs_section_text = sr_text.partition('PROJECTION')
	
	sr_eq_area_reprojection_text = 'PROJCS["LAEA_REPROJECTION_EQ_AREA",'+str(geogr_sr_reprojection_text)+',PROJECTION' + proj_cs_section_text[2]		
	sr_eq_area_reprojection = osr.SpatialReference(sr_eq_area_reprojection_text)
	
	sr_to_eq_area = osr.CoordinateTransformation(sr, sr_eq_area_reprojection)
	eq_area_to_sr = osr.CoordinateTransformation(sr_eq_area_reprojection, sr)
	
	union_polygon.Transform(sr_to_eq_area)
	
	Area_Size = union_polygon.GetArea()/1000000

""" Buffer research areas for buffered crater counting (geodesic buffer). """ 

def buffer_area(union_polygon, crater_area_list, all_craters, sr_wkt, generate_point_file, generate_polygon_file, generate_connectors_crater_polygon, flattening, major_axis, bufferfactor, crater_area_out_q, Area_IDs, multicore_operation, path_to_outfile, process_count, lock, layer_polygon, write_logfile, logfile):
	from shapely.geometry import Point
	global dist_buffer, vertices_angle_list, vertices_list, polygon_features_pass
	
	polygon_features_pass = []
	
	""" convert spatial reference from wkt to spatial reference type (during multiprocessing). Also, Transformations etc. need to be defined as Spatial reference data type can't
	be handed to function during multiprocessing """
	
	if isinstance(sr_wkt, basestring) == True:
		sr = osr.SpatialReference()
		sr.ImportFromWkt(sr_wkt)
		
	if isinstance(sr_wkt, basestring) == False:
		sr = sr_wkt
	
	""" GDAL objects such as layers, spatial references etc. can't be passed to a multi-core process (Swig objects can't 
	be pickled). This is why we declare a new output shapefile layer here. There are as many output files generated as there are multi-core processes. """
	
	if multicore_operation == True:
		if generate_polygon_file == True:
			outfile_split = os.path.split(path_to_outfile)
			outfile_path = str(outfile_split[0]) + "\\"
			outfile_name = str(outfile_split[1])
			outfile_name_no_extension = os.path.splitext(outfile_name)[0]
			
			driver = ogr.GetDriverByName('ESRI Shapefile')
			polygon_data_source = driver.CreateDataSource(outfile_path + outfile_name_no_extension + "_" + str(process_count) + ".shp")
			layer_polygon = polygon_data_source.CreateLayer('Buffer_Polygon', sr, geom_type = ogr.wkbPolygon)
			
			idField = ogr.FieldDefn('Crater_ID', ogr.OFTString)
			layer_polygon.CreateField(idField)
			idField = ogr.FieldDefn('Size_sq_km', ogr.OFTReal)
			layer_polygon.CreateField(idField)
			idField = ogr.FieldDefn('Buffer_m', ogr.OFTReal)
			layer_polygon.CreateField(idField)
	
	""" Transform vertices to geographic coordinates if spatial reference is projected. This doesn't apply when area_polygon is already geographic. 
	Reprojecting would lead to false coordinates otherwise. Grographic area_polygon and Projected sr can occur when (only) one research area is in 
	union_polygon and that research area has at least one hole. The error occurs with iterating the first inner ring. 
	
	Also, a geographic and a projected reference system from the input spatial reference is defined. """
	
	if sr.IsProjected():
		union_polygon.Segmentize(10000)
	if sr.IsGeographic():
		union_polygon.Segmentize(1)
	
	if sr.IsProjected():
		geogr_sr = sr.CloneGeogCS()
		sr_text = 'PROJCS["PROJECTED_LAMBERT_AEA",'+str(geogr_sr)+',PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["central_meridian",0.0],PARAMETER["latitude_of_origin",0.0],UNIT["Meter",1.0]]'
		proj_to_geog = osr.CoordinateTransformation(sr, geogr_sr)
		geog_to_proj = osr.CoordinateTransformation(geogr_sr, sr)
		union_polygon_sr = union_polygon.GetSpatialReference()
		if union_polygon_sr:
			if union_polygon_sr.IsProjected():
				union_polygon.Transform(proj_to_geog)
		if not union_polygon_sr:
			union_polygon.Transform(proj_to_geog)
	if sr.IsGeographic():
		geogr_sr = sr
		proj_to_geog = osr.CoordinateTransformation(sr, geogr_sr)
		geog_to_proj = osr.CoordinateTransformation(geogr_sr, sr)
		sr_text = 'PROJCS["PROJECTED_LAMBERT_AEA",'+str(geogr_sr)+',PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["central_meridian",0.0],PARAMETER["latitude_of_origin",0.0],UNIT["Meter",1.0]]'
	
	""" Get center of research area to define projection center  """
	
	union_polygon_centroid = union_polygon.Centroid()
	union_polygon_centroid_X, union_polygon_centroid_Y, union_polygon_centroid_Z = union_polygon_centroid.GetPoint()
	
	projection_center_X = round(union_polygon_centroid_X, 1)
	projection_center_Y = round(union_polygon_centroid_Y, 1)
	
	""" If a polygon intersects a dateline, it may happen that the centroid is wrongly identified (Polygon intersects dateline vs. 
	Polygon spans globe with a hole over dateline). If polygon cuts lines (which we assume to be datelines) at lon -179.8 and lon 179.8, we assume that a dateline 
	intersection is present. In this case, we set the projection center to (false projection center (if closer than 80 deg lon to central meridian) + 100 deg 
	lon). It is quite experimental to set a fixed lon value for such cases but it worked in all cases we tested. """
	
	dateline_1 = ogr.Geometry(ogr.wkbLineString)
	dateline_1.AddPoint(-179, 90)
	dateline_1.AddPoint(-179, -90)

	dateline_2 = ogr.Geometry(ogr.wkbLineString)
	dateline_2.AddPoint(179, 90)
	dateline_2.AddPoint(179, -90)
	
	""" If both lines (lon +179.8 & -179.8) are intersecting the polygon we assume a dateline intersection rather than a global polygon. """
	
	if union_polygon.Intersects(dateline_1) == True and union_polygon.Intersects(dateline_2) == True:
		if projection_center_X < 80 and projection_center_X > -80:
			projection_center_X = projection_center_X + 100
	
	""" Change reference meridian in geographic coordinate system to center of input polygon to avoid problems during dateline intersection. This new geographic 
	coordinate system is used for an equal area projection which consists of the geographic reference with new reference meridian and a Mollweide projected system
	with a central meridian at 0deg. """
	
	geogr_sr_reprojection_text = re.sub('(PRIMEM)(.*)(,)', r'\1["Reference_Meridian",' + str(projection_center_X) + '],', str(geogr_sr))
	geogr_sr_reprojection = osr.SpatialReference(geogr_sr_reprojection_text)
	
	proj_cs_section_text = sr_text.partition('PROJECTION')
	
	sr_eq_area_reprojection_text = 'PROJCS["LAEA_REPROJECTION_EQ_AREA",'+str(geogr_sr_reprojection_text)+',PROJECTION' + proj_cs_section_text[2]		
	sr_eq_area_reprojection = osr.SpatialReference(sr_eq_area_reprojection_text)
	
	sr_reprojection_text = 'PROJCS["LAEA_REPROJECTION_EQ_AREA",'+str(geogr_sr_reprojection_text)+',PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["central_meridian",' + str(projection_center_X) + '],PARAMETER["latitude_of_origin",' + str(projection_center_Y) + '],UNIT["Meter",1.0]]'
	sr_reprojection = osr.SpatialReference(sr_reprojection_text)
	
	""" Define reprojections. """
	
	geogr_to_geogr_reprojection = osr.CoordinateTransformation(geogr_sr, geogr_sr_reprojection)
	geogr_to_proj_reprojection = osr.CoordinateTransformation(geogr_sr, sr_reprojection)
	proj_reprojection_to_geogr_reprojection = osr.CoordinateTransformation(sr_reprojection, geogr_sr_reprojection)
	geogr_reprojection_to_proj_reprojection = osr.CoordinateTransformation(geogr_sr_reprojection, sr_reprojection)
	proj_reprojection_to_geogr = osr.CoordinateTransformation(sr_reprojection, geogr_sr)
	geogr_reprojection_to_geogr = osr.CoordinateTransformation(geogr_sr_reprojection, geogr_sr)
	geogr_reprojection_to_eq_area = osr.CoordinateTransformation(geogr_sr_reprojection, sr_eq_area_reprojection)
	eq_area_proj_to_proj = osr.CoordinateTransformation(sr_eq_area_reprojection, sr)
	geogr_to_eq_area = osr.CoordinateTransformation(geogr_sr, sr_eq_area_reprojection)
	sr_to_eq_area = osr.CoordinateTransformation(sr, sr_eq_area_reprojection)
	
	""" Reproject input data to new geographic projection in order to avoid problems during dateline intersections. """
	
	union_polygon.Transform(geogr_to_proj_reprojection)
	union_polygon.Transform(proj_reprojection_to_geogr_reprojection)
	
	""" Detailled logfile is only written in single-core mode """
	
	print "Buffering Area Polygons (", len(all_craters), ")..."
	if write_logfile == True and multicore_operation == False:
		logfile.write("\nBuffering Area Polygons ( " + str(len(all_craters)) + " )...\n")
		logfile.flush()
		
	st = time.time()
	
	counter = 0
	linear_ring_count = 0
	vertices_angle_list = dict()
	
	for area in union_polygon:
		
		number_of_inner_polygons = 0
		number_of_holes = 0
		
		#############################
		#	Step 1: CUT POLYGON		#
		#############################
		
		""" Determine number of polygons and number of holes to correctly split and buffer the data. """
		
		if area.GetGeometryName() == "LINEARRING":
			area_polygon = ogr.Geometry(ogr.wkbPolygon)
			area_polygon.AddGeometry(area)
		else:
			area_polygon = ogr.Geometry(ogr.wkbPolygon)
			area_polygon.AddGeometry(area.GetGeometryRef(0)) 
		
		if area.GetGeometryName() == "LINEARRING":
			number_of_holes = union_polygon.GetGeometryCount() - 1
			#print "There is one inner polygon and", number_of_holes, "holes in this polygon."
			
		if area.GetGeometryName() == "POLYGON":
			number_of_inner_polygons = union_polygon.GetGeometryCount()
			number_of_holes = area.GetGeometryCount() - 1 # When union_polygon becomes MULTIPOLYGON due to clip and buffer (formation of islands) 
			#print "There are", number_of_inner_polygons, "inner polygons and", area.GetGeometryCount() - 1, "holes in this polygon."
		
		""" Get inner rings: Get centroid of hole and split union_polygon or area into multiple parts according to holes. This way, only outlines need to be buffered (buffering 
		outlines is faster than buffering an inner ring). """
		
		if number_of_holes > 0:
			if area.GetGeometryName() == "LINEARRING":
				union_polygon_2 = ogr.CreateGeometryFromWkt(union_polygon.ExportToWkt())
			if area.GetGeometryName() == "POLYGON":
				area_2 = ogr.CreateGeometryFromWkt(area.ExportToWkt())
			
			for hole_count in range(number_of_holes):
				if number_of_inner_polygons == 0:
					inner_ring_geometry = union_polygon.GetGeometryRef(hole_count + 1)
				if number_of_inner_polygons > 0:
					inner_ring_geometry = area.GetGeometryRef(hole_count + 1)
				
				""" Get center of inner ring. Inner ring geometry has to be transformed from LINEARRING to POLYGON geometry. """
				
				inner_ring_polygon_geometry = ogr.Geometry(ogr.wkbPolygon)
				inner_ring_polygon_geometry.AddGeometry(inner_ring_geometry)
				
				inner_ring_centroid = inner_ring_polygon_geometry.Centroid()
				inner_ring_centroid_X, inner_ring_centroid_Y, inner_ring_centroid_Z = inner_ring_centroid.GetPoint()
				
				""" Generate lines which intersect centroid of inner rings.  """
				
				cut_line = ogr.Geometry(ogr.wkbLineString)
				cut_line.AddPoint(0, 90)
				cut_line.AddPoint(inner_ring_centroid_X, inner_ring_centroid_Y)
				cut_line.AddPoint(0, -90)
				cut_line.Segmentize(1)
				
				""" Generate splitted area from reference area and (buffered) split lines. OGR doesn't support polygon splitting 
				by lines. """
				
				buffered_cut_line = cut_line.Buffer(0.00000000001) 
				if area.GetGeometryName() == "LINEARRING":
					union_polygon_2 = union_polygon_2.Difference(buffered_cut_line)
					area_polygon = union_polygon_2
				if area.GetGeometryName() == "POLYGON":
					area_2 = area_2.Difference(buffered_cut_line)
					area_polygon = area_2
		
		#################################
		#	Step 2: READ VERTICES		#
		#################################
		
		""" Get each linear ring in the area polygon (polygon outlines or splitted polygon parts when polygon has holes) and buffer outlines. 
		Holes are not present anymore. """
		
		for linear_ring in area_polygon:
			vertices_angle_list[linear_ring_count] = []
			
			""" area_polygon is MULTIPOLYGON when polygon with holes (splitted polygon) is present. area_polygon is LINEARRING when no holes 
			are present (no splitted polygon). Using GetGeometryRef(0) we get the linear ring from the polygon when splitting was conducted. """
			
			if linear_ring.GetGeometryName() == "LINEARRING":
				linear_ring = linear_ring
			else:
				linear_ring = linear_ring.GetGeometryRef(0)
			
			no_of_polygon_vertices = linear_ring.GetPointCount()
			
			""" Get Buffer Points for outer polygon boundary. """
	
			for vertex in xrange(no_of_polygon_vertices):
				current_vertex_X, current_vertex_Y, z = linear_ring.GetPoint(vertex)
				
				""" Check if polygon is closed (first and last vertex share same coordinates). If not, change neighboring vertices accordingly. """
				
				if linear_ring.GetPoint(0) == linear_ring.GetPoint(no_of_polygon_vertices - 1):
					if vertex == 0:
						vertex_position = "first"
						previous_vertex_X, previous_vertex_Y, z = linear_ring.GetPoint(no_of_polygon_vertices - 2)
						next_vertex_X, next_vertex_Y, z = linear_ring.GetPoint(vertex + 1)
					if vertex == no_of_polygon_vertices - 1:
						vertex_position = "last"
						previous_vertex_X, previous_vertex_Y, z = linear_ring.GetPoint(vertex - 1)
						next_vertex_X, next_vertex_Y, z = linear_ring.GetPoint(1)
					if vertex > 0 and vertex < no_of_polygon_vertices - 1:
						vertex_position = "middle" 
						previous_vertex_X, previous_vertex_Y, z = linear_ring.GetPoint(vertex - 1)
						next_vertex_X, next_vertex_Y, z = linear_ring.GetPoint(vertex + 1)
				else:
					if vertex == 0:
						vertex_position = "first"
						previous_vertex_X, previous_vertex_Y, z = linear_ring.GetPoint(no_of_polygon_vertices - 1)
						next_vertex_X, next_vertex_Y, z = linear_ring.GetPoint(vertex + 1)
					if vertex == no_of_polygon_vertices - 1:
						vertex_position = "last"
						previous_vertex_X, previous_vertex_Y, z = linear_ring.GetPoint(vertex - 1)
						next_vertex_X, next_vertex_Y, z = linear_ring.GetPoint(0)
					if vertex > 0 and vertex < no_of_polygon_vertices - 1:
						vertex_position = "middle" 
						previous_vertex_X, previous_vertex_Y, z = linear_ring.GetPoint(vertex - 1)
						next_vertex_X, next_vertex_Y, z = linear_ring.GetPoint(vertex + 1)
				
				""" Calculate angle between previous vertex - current vertex and current vertex - next vertex from vincenty's inverse formula 
				(calculation on a spheroid). """
				
				inverse_vincenty(flattening, major_axis, previous_vertex_Y, previous_vertex_X, current_vertex_Y, current_vertex_X )
				angle_prev = direction12
				
				inverse_vincenty(flattening, major_axis, current_vertex_Y, current_vertex_X, next_vertex_Y, next_vertex_X )
				angle_next = direction12
				
				""" Ensure that angles remain within 0-360 deg range. """
				
				if angle_prev < 0:
					angle_prev = angle_prev + 360
				if angle_next < 0:
					angle_next = angle_next + 360
				
				angle_prev_BP = angle_prev - 90 
				angle_next_BP = angle_next - 90 
				
				""" Ensure that buffer points are perpendicular to reference area and remain within 0-360 deg range. """
				
				if angle_prev_BP < 0:
					angle_prev_BP = angle_prev_BP + 360  
				if angle_next_BP < 0:
					angle_next_BP = angle_next_BP + 360 
				
				""" Add to list in which coordinates and angles for buffer vertices calculation are stored. """
				
				vertices_angle_list[linear_ring_count].append([current_vertex_X, current_vertex_Y, angle_prev_BP, angle_prev, angle_next])
				
				""" Remember coordinates of the buffer polygon's first vertex - used later to close polygon. """
				
				if vertex == 0:
					X_0 = current_vertex_X
					Y_0 = current_vertex_Y
					angle_0 = angle_prev_BP
					angle_prev_0 = angle_prev
					angle_next_0 = angle_next
				
				""" Calculate buffer points between angle_prev-90 and angle_next-90 (used for round buffer edges). """
				
				if angle_next_BP > angle_prev_BP:
					if (angle_next_BP) - (angle_prev_BP) <= 180:
						
						""" Angles between previous and next buffer point - generate buffer points which are outside the original polygon. """
						
						angles = numpy.arange(angle_prev_BP, angle_next_BP, 8)
						for angle in angles:
							vertices_angle_list[linear_ring_count].append([current_vertex_X, current_vertex_Y, angle, angle_prev, angle_next])
					if (angle_next_BP) - (angle_prev_BP) > 180:
						
						""" Scissor intersection: angles between next buffer point and 360 degrees and between 0 degrees and previous buffer point - generate buffer points which are outside the original polygon. """
						
						angles = numpy.arange(angle_next_BP, 360, 8)
						for angle in angles:
							vertices_angle_list[linear_ring_count].append([current_vertex_X, current_vertex_Y, angle, angle_prev, angle_next])
						angles = numpy.arange(0, angle_prev_BP, 8)
						for angle in angles:
							vertices_angle_list[linear_ring_count].append([current_vertex_X, current_vertex_Y, angle, angle_prev, angle_next])
				
				if angle_next_BP < angle_prev_BP:
					if angle_next_BP - angle_prev_BP > -180:
						
						""" Angles between next and previous buffer point - generate buffer points which are outside the original polygon. """
						
						angles = numpy.arange(angle_next_BP, angle_prev_BP, 8)
						for angle in angles:
							vertices_angle_list[linear_ring_count].append([current_vertex_X, current_vertex_Y, angle, angle_prev, angle_next])
					if angle_next_BP - angle_prev_BP <= -180:
						
						""" Scissor intersection - angles between previous buffer point and 360 degrees and between 0 degrees and next buffer point - generate buffer points which are outside the original polygon. """
						
						angles = numpy.arange(angle_prev_BP, 360, 8)
						for angle in angles:
							vertices_angle_list[linear_ring_count].append([current_vertex_X, current_vertex_Y, angle, angle_prev, angle_next])
						angles = numpy.arange(0, angle_next_BP, 8)
						for angle in angles:
							vertices_angle_list[linear_ring_count].append([current_vertex_X, current_vertex_Y, angle, angle_prev, angle_next])
				vertices_angle_list[linear_ring_count].append([current_vertex_X, current_vertex_Y, angle_next_BP, angle_prev, angle_next])
				
			""" Close polygon using the first vertex. """
			
			vertices_angle_list[linear_ring_count].append([X_0, Y_0, angle_0, angle_prev_0, angle_next_0])
			
			linear_ring_count += 1
		
		""" Special case: If only one research area with hole(s) is investigated, iteration in union_polygon would not consider polygon1, polygon2, 
		polygon3, etc., but ring1, ring2, ring3, etc. As the inner ring is already considered during iteration (because it is assumed that in 'for area in union_polygon' 
		area is a research area and not a linear ring), the function has to be stopped here. Otherwise, every further iteration in union polygon would consider 
		the linear rings which have already been considered in the function. This would lead to too many resulting polygons (rings*no_of_craters and not no_of_craters). """
		
		if len(Area_IDs) == 1 and number_of_holes >= 1:
			break
		
	#################################
	#	Step 3: BUFFER POLYGON		#
	#################################
	
	""" Iterate crater features """
	
	cr_cnt = 0
	for crater in all_craters:
		
		if len(all_craters) > 0:
			print "Process", process_count, ": Processing crater", crater[0], ":", round(((float(cr_cnt) / float(len(all_craters)))*100), 1), "%"
		
		""" Detailled logfile is only written in single-core mode """
		
		if write_logfile == True and multicore_operation == False:
			logfile.write("Process " + str(process_count) +  ": Processing crater " + str(crater[0]) + ": " + str(round(((float(cr_cnt) / float(len(all_craters)))*100), 1)) + " %\n")
			logfile.flush()
		
		""" Define geometries. """
		
		BCC_union_polygon = ogr.Geometry(ogr.wkbPolygon)
		point_geometry = ogr.Geometry(ogr.wkbPoint)
	
		if generate_point_file == True:
			point_featureDefn = layer_points.GetLayerDefn()
			point_feature = ogr.Feature(point_featureDefn)
			point_feature.SetGeometry(point_geometry)
			
		if generate_polygon_file == True:
			polygon_featureDefn = layer_polygon.GetLayerDefn()
			polygon_feature = ogr.Feature(polygon_featureDefn)
			polygon_feature.SetGeometry(BCC_union_polygon)
			
		""" Get buffer distance """
		
		dist_buffer = crater[1]*1000/2
		
		""" Calculate coordinates of buffer points for outer polygon """
		
		geometry_dict_count = 0
		len_vertices_angle_list = len(vertices_angle_list)
		
		for geometry_dict_count in range(len_vertices_angle_list):
			
			splitted_buffered_ring = ogr.Geometry(ogr.wkbLinearRing)
			splitted_buffered_polygon = ogr.Geometry(ogr.wkbPolygon)
			
			direct_vincenty(flattening, major_axis, vertices_angle_list[geometry_dict_count], dist_buffer, bufferfactor)
			
			for buffer_vertex in buffer_vertices_list:
				
				""" Add point to ring geometry """
				
				splitted_buffered_ring.AddPoint(buffer_vertex[0], buffer_vertex[1])
				
				if generate_point_file == True:
					
					point_geometry.AddPoint(buffer_vertex[0], buffer_vertex[1])
					point_geometry.Transform(geogr_reprojection_to_proj_reprojection)
					point_geometry.Transform(proj_reprojection_to_geogr)
					point_feature.SetGeometry(point_geometry)
					point_feature.SetField('prev', buffer_vertex[2])
					point_feature.SetField('next', buffer_vertex[3])
					point_feature.SetField('angle', angle_prev_BP)
					layer_points.CreateFeature(point_feature)
			
			""" Add ring to polygon geometry """
			
			splitted_buffered_polygon.AddGeometry(splitted_buffered_ring)
			
			""" eliminate unwanted holes due to self-intersections on outer boundary using a planar 0 buffer """
			
			splitted_buffered_polygon = splitted_buffered_polygon.Buffer(0)
			
			""" Errors may occur during Buffer(0) so that two polygons are formed from one polygon due to severe self-intersections. 
			This would lead to an invalid geometry which could not be added to the BNSC_union_polygon. To prevent this, all geometries 
			in the splitted_buffered_polygon are investigated, Buffered (0) again and then added to the BNSC_union_polygon. """
			
			if splitted_buffered_polygon.IsValid() == False:
			
				for linear_ring_splitted_buffered_polygon in splitted_buffered_polygon:
					
					""" Add linear_ring_splitted_buffered_polygon to new polygon. """
					
					polygon_part_splitted_buffered_polygon = ogr.Geometry(ogr.wkbPolygon)
					polygon_part_splitted_buffered_polygon.AddGeometry(linear_ring_splitted_buffered_polygon)
					polygon_part_splitted_buffered_polygon = polygon_part_splitted_buffered_polygon.Buffer(0)
					
					if polygon_part_splitted_buffered_polygon.IsValid() == False:
						print "Error due to severe self-intersection during buffering. Please generate and check the output shapefile for errors."
						if write_logfile == True and multicore_operation == False:
							logfile.write("Error due to severe self-intersection during buffering. Please generate and check the output shapefile for errors.")
							logfile.flush()
						exit()
					
					BCC_union_polygon = BCC_union_polygon.Union(polygon_part_splitted_buffered_polygon)
			
			if splitted_buffered_polygon.IsValid() == True:
				
				BCC_union_polygon = BCC_union_polygon.Union(splitted_buffered_polygon)
			
			geometry_dict_count += 1
		
		BCC_union_polygon.Transform(geogr_reprojection_to_proj_reprojection)
		BCC_union_polygon.Transform(proj_reprojection_to_geogr)
		BCC_union_polygon.Transform(geogr_to_eq_area)
		
		geodesic_area = BCC_union_polygon.GetArea()/1000000
		
		if generate_polygon_file == True:
			
			""" Project back to original spatial reference. """
			
			BCC_union_polygon.Transform(eq_area_proj_to_proj) 
			polygon_feature.SetGeometry(BCC_union_polygon) 
			polygon_feature.SetField('Crater_ID', crater[0])
			polygon_feature.SetField('Size_sq_km', geodesic_area)
			polygon_feature.SetField('Buffer_m', dist_buffer * bufferfactor)
			layer_polygon.CreateFeature(polygon_feature)
			polygon_features_pass.append(polygon_feature)
		
		ejecta_diam = crater[1] * bufferfactor
		crater_area_list.append([crater[1], crater[2], crater[3], ejecta_diam, geodesic_area])
		
		""" Delete geometries """
		
		del splitted_buffered_ring
		del splitted_buffered_polygon
		del BCC_union_polygon
		del point_geometry
		if generate_point_file == True:
			del point_feature
		if generate_polygon_file == True:
			del polygon_feature
		
		cr_cnt += 1
		counter += 1
	
	crater_area_out_q.put(crater_area_list)
	
	""" Buffer craters to define the size of ejecta blankets. """
	
def NSC_BNSC_buffer_craters(buffered_craters_out_q, craters_for_counting_list, flattening, major_axis, sr_wkt, lock, bufferfactor_crater):
	global buffered_craters_wkt_list
	
	""" Convert spatial reference from wkt to spatial reference type (during multiprocessing). Also, Transformations etc. need to be defined as Spatial reference data type can't
	be handed to function during multiprocessing. """
	
	if isinstance(sr_wkt, basestring) == True:
		sr = osr.SpatialReference()
		sr.ImportFromWkt(sr_wkt)
		
	if isinstance(sr_wkt, basestring) == False:
		sr = sr_wkt
	
	""" Get geographic reference system. """
	
	if sr.IsProjected():
		geogr_sr = sr.CloneGeogCS()
	if sr.IsGeographic():
		geogr_sr = sr
	
	buffered_craters_wkt_list = []
	
	""" Create list of angles used for polar coordinates (equivalent to the number of geodesic buffer polygon's vertices). """
	
	no_of_buffer_vertices = 180.0
	divisor = 360.0/no_of_buffer_vertices
	ejecta_angles = numpy.arange(0.0, 360.0, divisor)
	ejecta_angles = numpy.append(ejecta_angles, 0.0)
	
	""" Get crater information """
	
	for crater in craters_for_counting_list:
		
		vertices_angle_list = []
		
		crater_id = crater[0]
		crater_diameter = crater[1]
		crater_centroid_X = crater[2]
		crater_centroid_Y = crater[3]
		distance_crater_polygon = crater[4]
		
		""" Check if crater obliterates larger crater from previously added field. """
		
		if len(crater) == 6:
			obliterates = "obliterates"
		if len(crater) == 5:
			obliterates = "included"
		
		""" Determine buffer distance in meter """
		
		dist_buffer = ((crater_diameter * 1000) / 2) + ((crater_diameter * 1000)/2 * (bufferfactor_crater - 1))
		
		""" Calculate vertices of buffered crater polygon """
		
		for ejecta_angle in ejecta_angles:
			vertices_angle_list.append([crater_centroid_X, crater_centroid_Y, ejecta_angle, 0, 0])
		
		""" bufferfactor_crater = 1 in this case since buffer distance has already been determined in this function """
		
		direct_vincenty(flattening, major_axis, vertices_angle_list, dist_buffer, 1) 
		
		""" Define buffered crater polygon """
		
		crater_ring_geometry = ogr.Geometry(ogr.wkbLinearRing)
		crater_polygon_geometry = ogr.Geometry(ogr.wkbPolygon)
		
		crater_ring_geometry.AssignSpatialReference(geogr_sr)
		crater_polygon_geometry.AssignSpatialReference(geogr_sr)
		
		""" Create buffered crater polygon """
		
		for crater_buffer_vertex in buffer_vertices_list:
			crater_ring_geometry.AddPoint(crater_buffer_vertex[0], crater_buffer_vertex[1])
		
		crater_polygon_geometry.AddGeometry(crater_ring_geometry)
		
		""" Convert polygon to WKT to pass it to other functions. Multicore doesn't support GDAL objects to pass. Text is OK. """
		
		buffered_craters_wkt_list.append([crater_polygon_geometry.ExportToWkt(), distance_crater_polygon, crater_id, crater_diameter, obliterates])
		
	buffered_craters_out_q.put(buffered_craters_wkt_list)

""" Modify initial reference areas for NSC and BNSC. """

def NSC_BNSC_exclude_craters(approach, buffered_craters_wkt_list, union_polygon, sr_wkt, generate_polygon_file, path_to_outfile, craters_for_counting_list, craters_for_counting_list_BNSC, multicore_operation, layer_polygon, bufferfactor, bufferfactor_crater, process_count, crater_area_out_q, flattening, major_axis, Area_IDs, write_logfile, logfile):
	from shapely.geometry import Point
	global crater_area_list
	
	crater_area_list = []
	
	""" Convert spatial reference from wkt to spatial reference type (during multiprocessing). Also, Transformations etc. need to be defined as Spatial reference data type can't
	be handed to function during multiprocessing """
	
	if isinstance(sr_wkt, basestring) == True:
		sr = osr.SpatialReference()
		sr.ImportFromWkt(sr_wkt)
		
	if isinstance(sr_wkt, basestring) == False:
		sr = sr_wkt
	
	geogr_sr = sr.CloneGeogCS()
	
	if multicore_operation == True:
		if generate_polygon_file == True:
			outfile_split = os.path.split(path_to_outfile)
			outfile_path = str(outfile_split[0]) + "\\"
			outfile_name = str(outfile_split[1])
			outfile_name_no_extension = os.path.splitext(outfile_name)[0]
			
			driver = ogr.GetDriverByName('ESRI Shapefile') 
			polygon_data_source = driver.CreateDataSource(outfile_path + outfile_name_no_extension + "_" + str(process_count) + ".shp") 
			
			layer_polygon = polygon_data_source.CreateLayer('Buffer_Polygon', sr, geom_type = ogr.wkbPolygon)
			
			idField = ogr.FieldDefn('Size_sq_km', ogr.OFTReal)
			layer_polygon.CreateField(idField)
			idField = ogr.FieldDefn('crater_ID', ogr.OFTInteger)
			layer_polygon.CreateField(idField)
	
	if generate_polygon_file == True:
		polygon_featureDefn = layer_polygon.GetLayerDefn()
		polygon_feature = ogr.Feature(polygon_featureDefn)
	
	buffered_craters_geometry_list = []
	
	""" Define coordinate reprojections: craters are in geogr sr and have to be projected to proj sr if input reference system is projected. """
	
	geog_to_proj = osr.CoordinateTransformation(geogr_sr, sr)
	proj_to_geogr = osr.CoordinateTransformation(sr, geogr_sr) 
	
	if sr.IsProjected():
		union_polygon_sr = union_polygon.GetSpatialReference()
		if union_polygon_sr:
			if union_polygon_sr.IsProjected():
				union_polygon.Transform(proj_to_geogr)
		if not union_polygon_sr:
			union_polygon.Transform(proj_to_geogr)
	
	""" Get center of research area to define projection center. """
	
	union_polygon_centroid = union_polygon.Centroid()
	union_polygon_centroid_X, union_polygon_centroid_Y, union_polygon_centroid_Z = union_polygon_centroid.GetPoint()
	
	projection_center_X = round(union_polygon_centroid_X, 1)
	projection_center_Y = round(union_polygon_centroid_Y, 1)
	
	
	""" If a polygon intersects a dateline, it may happen that the centroid is wrongly identified (Polygon intersects dateline vs. 
	Polygon spans globe with a hole over dateline). If polygon cuts lines (which we assume to be datelines) at lon -179.8 and lon 179.8, we assume that a dateline 
	intersection is present. In this case, we set the projection center to (false projection center (if closer than 80 deg lon to central meridian) + 100 deg 
	lon). It is quite experimental to set a fixed lon value for such cases but it worked in all cases we tested. """
	
	dateline_1 = ogr.Geometry(ogr.wkbLineString)
	dateline_1.AddPoint(-179, 90)
	dateline_1.AddPoint(-179, -90)

	dateline_2 = ogr.Geometry(ogr.wkbLineString)
	dateline_2.AddPoint(179, 90)
	dateline_2.AddPoint(179, -90)
	
	""" If both lines (lon +179.8 & -179.8) are intersecting the polygon, we assume a dateline intersection rather than a global polygon. """
	
	if union_polygon.Intersects(dateline_1) == True and union_polygon.Intersects(dateline_2) == True:
		if projection_center_X < 80 and projection_center_X > -80:
			projection_center_X = projection_center_X + 100
	
	sr_text = 'PROJCS["PROJECTED_LAMBERT_AEA",'+str(geogr_sr)+',PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["central_meridian",0.0],PARAMETER["latitude_of_origin",0.0],UNIT["Meter",1.0]]'
	
	if sr.IsProjected():
		union_polygon.Transform(geog_to_proj)

	""" Change reference meridian in geographic coordinate system to center of input polygon to avoid problems during dateline intersection. """
	
	geogr_sr_reprojection_text = re.sub('(PRIMEM)(.*)(,)', r'\1["Reference_Meridian",' + str(projection_center_X) + '],', str(geogr_sr))
	geogr_sr_reprojection = osr.SpatialReference(geogr_sr_reprojection_text)
	
	proj_cs_section_text = sr_text.partition('PROJECTION')
	
	sr_eq_area_reprojection_text = 'PROJCS["LAEA_REPROJECTION_EQ_AREA",'+str(geogr_sr_reprojection_text)+',PROJECTION' + proj_cs_section_text[2]		
	sr_eq_area_reprojection = osr.SpatialReference(sr_eq_area_reprojection_text)
	
	sr_to_eq_area = osr.CoordinateTransformation(sr, sr_eq_area_reprojection)
	eq_area_to_sr = osr.CoordinateTransformation(sr_eq_area_reprojection, sr)
	
	""" Get Geometry from WKT, reproject buffered craters to input spatial reference and store geometry in buffered_craters_geometry list. """
	
	for buffered_crater in buffered_craters_wkt_list:
		buffered_crater_polygon = ogr.CreateGeometryFromWkt(buffered_crater[0])
		buffered_crater_polygon.Transform(geog_to_proj)
		buffered_craters_geometry_list.append(buffered_crater_polygon)
	
	buffered_craters_geometry_list_2 = list(buffered_craters_geometry_list)
	crater_index = 0
	
	for buffered_crater in buffered_craters_geometry_list:
		
		""" Get distance to crater from previous list (determined in get_craters_for_buffering_NSC_BNSC function via Vincenty's formulae) """
		
		distance_crater_area = buffered_craters_wkt_list[crater_index][1]
		original_crater_id = buffered_craters_wkt_list[crater_index][2]
		obliteration_status = buffered_craters_wkt_list[crater_index][4]
		
		if len(craters_for_counting_list) > 0:
			print "Process", process_count, ": Processing crater", original_crater_id, ":", round(((float(crater_index) / float(len(craters_for_counting_list)))*100), 1), "%"
			
			""" Detailled logfile is only written in single-core mode """
			
			if write_logfile == True and multicore_operation == False:
				logfile.write("Process " + str(process_count) + ": Processing crater " + str(original_crater_id) + ": " + str(round(((float(crater_index) / float(len(craters_for_counting_list)))*100), 1)) + " %\n")
				logfile.flush()
		
		""" Start with largest crater, erase larger craters (craters with index -1) from initial reference area """
		
		if crater_index > 0:
			buffered_crater_geometry_larger = buffered_craters_geometry_list[crater_index - 1]
			union_polygon = union_polygon.Difference(buffered_crater_geometry_larger)
			
			""" Ignore craters that are located on top of an ejecta blanket for CSFD measurement. """
			
			if obliteration_status == "obliterates":
				crater_index += 1
				continue
			
		""" Buffer reference area during BNSC after a crater + ejecta blanket is removed. """
			
		if approach == "BNSC":
			
			""" The reference area with removed craters and two lists are available at this point: buffered_craters_geometry_list (all the buffered crater geometries with 
			respective IDs and distance from initial reference area) and craters_for_counting_list_BNSC (information on craters which are included in the counting process but 
			splitted for multicore processing). During multiprocessing, craters_for_counting_list_BNSC is used to buffer the reference area according to the respective 
			crater diameters. When a crater is removed from the reference area, the script checks if that particular crater is included in craters_for_counting_list_BNSC. 
			If it is included, the reference area is buffered. This procedure is conducted to optimize the BNSC buffering for multicore computation. In every process, all craters in 
			buffered_craters_geometry_list are removed from the reference area but the area is only buffered if the crater is included in the splitted 
			craters_for_counting_list_BNSC list. The process of removing craters is conducted multiple times (one time in each process - extra effort but not very time consuming) 
			while the areas for buffering are unique for each process. """
			
			""" Check if crater which was currently removed from the reference area is in craters_for_counting_list_BNSC """
			
			for BNSC_crater in craters_for_counting_list_BNSC:
				BNSC_crater_id = BNSC_crater[0]
				
				if BNSC_crater_id == original_crater_id:
					BNSC_crater_diameter = BNSC_crater[1]
					buffer_distance_polygon_BNSC = ((BNSC_crater_diameter * 1000)/2) * bufferfactor
					
					#####################################################################
					#																	#
					# BNSC BUFFER (not in a separate function due to multicore support)	#
					#																	#
					#####################################################################
					
					BNSC_union_polygon = ogr.Geometry(ogr.wkbPolygon)
					
					""" Segmentize polygon to decreacse the effect of angular distortion in polar regions. """
					
					if sr.IsProjected():
						union_polygon.Segmentize(10000)
						union_polygon.Transform(proj_to_geogr)
					if sr.IsGeographic():
						union_polygon.Segmentize(1)
					
					""" Consider dateline intersections: Get center of research area to define projection center. """
					
					union_polygon_centroid = union_polygon.Centroid()
					union_polygon_centroid_X, union_polygon_centroid_Y, union_polygon_centroid_Z = union_polygon_centroid.GetPoint()
					projection_center_X = round(union_polygon_centroid_X, 1)
					projection_center_Y = round(union_polygon_centroid_Y, 1)
					
					""" If a polygon intersects a dateline, it may happen that the centroid is wrongly identified (Polygon intersects dateline vs. 
					Polygon spans globe with a hole over dateline). If polygon cuts lines (which we assume to be datelines) at lon -179.8 and lon 179.8, we assume that a dateline 
					intersection is present. In this case, we set the projection center to (false projection center (if closer than 80 deg lon to central meridian) + 100 deg 
					lon). It is quite experimental to set a fixed lon value but it works in all cases we tested. """
					
					dateline_1 = ogr.Geometry(ogr.wkbLineString)
					dateline_1.AddPoint(-179, 90)
					dateline_1.AddPoint(-179, -90)
				
					dateline_2 = ogr.Geometry(ogr.wkbLineString)
					dateline_2.AddPoint(179, 90)
					dateline_2.AddPoint(179, -90)
					
					""" If both lines (lon +179.8 & -179.8) are intersecting the polygon, we assume a dateline intersection rather than a global polygon. """
					
					if union_polygon.Intersects(dateline_1) == True and union_polygon.Intersects(dateline_2) == True:
						if projection_center_X < 80 and projection_center_X > -80:
							projection_center_X = projection_center_X + 100
					
					""" Change reference meridian in geographic coordinate system to center of input polygon to avoid problems during dateline intersection. """
					
					geogr_sr_reprojection_text = re.sub('(PRIMEM)(.*)(,)', r'\1["Reference_Meridian",' + str(projection_center_X) + '],', str(geogr_sr))
					geogr_sr_reprojection = osr.SpatialReference(geogr_sr_reprojection_text)
					
					sr_reprojection_text = 'PROJCS["LAEA_REPROJECTION_EQ_AREA",'+str(geogr_sr_reprojection_text)+',PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["central_meridian",' + str(projection_center_X) + '],PARAMETER["latitude_of_origin",' + str(projection_center_Y) + '],UNIT["Meter",1.0]]'
					sr_reprojection = osr.SpatialReference(sr_reprojection_text)
					
					geogr_to_proj_reprojection = osr.CoordinateTransformation(geogr_sr, sr_reprojection)
					proj_reprojection_to_geogr_reprojection = osr.CoordinateTransformation(sr_reprojection, geogr_sr_reprojection)
					geogr_reprojection_to_proj_reprojection = osr.CoordinateTransformation(geogr_sr_reprojection, sr_reprojection) 
					proj_reprojection_to_geogr = osr.CoordinateTransformation(sr_reprojection, geogr_sr)
					
					""" Reproject input polygon to consider dateline intersections """
					
					union_polygon.Transform(geogr_to_proj_reprojection)
					union_polygon.Transform(proj_reprojection_to_geogr_reprojection)
					
					""" LINEARRING features must be used when MULTIPOLYGON (multiple areas) is present and LINEARRING must be taken when POLYGON (single area) 
					is present. """
					
					for area in union_polygon:
						number_of_inner_polygons = 0
						number_of_holes = 0
						
						#############################
						#	Step 1: CUT POLYGON		#
						#############################
						
						""" Determine number of polygons and number of holes to correctly split and buffer the data. """
						
						if area.GetGeometryName() == "LINEARRING":
							area_polygon = ogr.Geometry(ogr.wkbPolygon)
							area_polygon.AddGeometry(area)
						else:
							area_polygon = ogr.Geometry(ogr.wkbPolygon)
							area_polygon.AddGeometry(area.GetGeometryRef(0)) 
						
						if area.GetGeometryName() == "LINEARRING":
							number_of_holes = union_polygon.GetGeometryCount() - 1
							#print "There is one inner polygon and", number_of_holes, "holes for this polygon."
							
						if area.GetGeometryName() == "POLYGON":
							number_of_inner_polygons = union_polygon.GetGeometryCount()
							number_of_holes = area.GetGeometryCount() - 1 # When union_polygon becomes MULTIPOLYGON due to clip and buffer (formation of islands) 
							#print "There are", number_of_inner_polygons, "inner polygons and", area.GetGeometryCount() - 1, "holes for this polygon."
						
						""" Get inner rings: Get centroid of hole and split union_polygon or area into multiple parts according to holes. This way, only outlines need to be buffered (buffering 
						outlines is faster than buffering inner ring). """
						
						if number_of_holes > 0:
							if area.GetGeometryName() == "LINEARRING":
								union_polygon_2 = ogr.CreateGeometryFromWkt(union_polygon.ExportToWkt())
							
							if area.GetGeometryName() == "POLYGON":
								area_2 = ogr.CreateGeometryFromWkt(area.ExportToWkt())
							
							for hole_count in range(number_of_holes):
								if number_of_inner_polygons == 0:
									inner_ring_geometry = union_polygon.GetGeometryRef(hole_count + 1)
								if number_of_inner_polygons > 0:
									inner_ring_geometry = area.GetGeometryRef(hole_count + 1)
								
								""" Get center of inner ring. Inner ring geometry has to be transformed from LINEARRING to POLYGON geometry. """
								
								inner_ring_polygon_geometry = ogr.Geometry(ogr.wkbPolygon)
								inner_ring_polygon_geometry.AddGeometry(inner_ring_geometry)
								
								inner_ring_centroid = inner_ring_polygon_geometry.Centroid()
								inner_ring_centroid_X, inner_ring_centroid_Y, inner_ring_centroid_Z = inner_ring_centroid.GetPoint()
								
								""" Generate lines which intersect centroid of inner rings.  """
								
								cut_line = ogr.Geometry(ogr.wkbLineString)
								cut_line.AddPoint(0, 90)
								cut_line.AddPoint(inner_ring_centroid_X, inner_ring_centroid_Y)
								cut_line.AddPoint(0, -90)
								
								""" Generate splitted area from reference area and (buffered) split lines. OGR doesn't support polygon splitting 
								by lines. """
								
								buffered_cut_line = cut_line.Buffer(0.00000000001) 
								if area.GetGeometryName() == "LINEARRING":
									union_polygon_2 = union_polygon_2.Difference(buffered_cut_line)
									area_polygon = union_polygon_2
								if area.GetGeometryName() == "POLYGON":
									area_2 = area_2.Difference(buffered_cut_line)
									area_polygon = area_2
						
						#################################
						#	Step 2: BUFFER POLYGON		#
						#################################
						
						""" Get each linear ring in the area polygon (polygon outlines or splitted polygon parts when polygon has holes) and buffer outlines. 
						Holes are not present anymore. """
						
						for linear_ring in area_polygon:
							
							""" area_polygon is MULTIPOLYGON when polygon with holes (splitted polygon) is present. area_polygon is LINEARRING when no holes 
							are present (no splitted polygon). Using GetGeometryRef(0) we get the linear ring from the polygon when splitting was conducted. """
							
							if linear_ring.GetGeometryName() == "LINEARRING":
								linear_ring = linear_ring
							else:
								linear_ring = linear_ring.GetGeometryRef(0)
							
							splitted_buffered_polygon = ogr.Geometry(ogr.wkbPolygon)
							splitted_buffered_ring = ogr.Geometry(ogr.wkbLinearRing)
							
							no_of_polygon_vertices = linear_ring.GetPointCount()
							vertices_angle_list = []
							
							""" Get buffer points for outer polygon boundary. """
					
							for vertex in xrange(no_of_polygon_vertices):
								current_vertex_X, current_vertex_Y, z = linear_ring.GetPoint(vertex)
								
								""" Check if polygon is closed (first and last vertex share same coordinates) and modify neighboring vertices accordingly. """
								
								if linear_ring.GetPoint(0) == linear_ring.GetPoint(no_of_polygon_vertices - 1):
									if vertex == 0:
										vertex_position = "first"
										previous_vertex_X, previous_vertex_Y, z = linear_ring.GetPoint(no_of_polygon_vertices - 2)
										next_vertex_X, next_vertex_Y, z = linear_ring.GetPoint(vertex + 1)
									if vertex == no_of_polygon_vertices - 1:
										vertex_position = "last"
										previous_vertex_X, previous_vertex_Y, z = linear_ring.GetPoint(vertex - 1)
										next_vertex_X, next_vertex_Y, z = linear_ring.GetPoint(1)
									if vertex > 0 and vertex < no_of_polygon_vertices - 1:
										vertex_position = "middle" 
										previous_vertex_X, previous_vertex_Y, z = linear_ring.GetPoint(vertex - 1)
										next_vertex_X, next_vertex_Y, z = linear_ring.GetPoint(vertex + 1)
								else:
									if vertex == 0:
										vertex_position = "first"
										previous_vertex_X, previous_vertex_Y, z = linear_ring.GetPoint(no_of_polygon_vertices - 1)
										next_vertex_X, next_vertex_Y, z = linear_ring.GetPoint(vertex + 1)
									if vertex == no_of_polygon_vertices - 1:
										vertex_position = "last"
										previous_vertex_X, previous_vertex_Y, z = linear_ring.GetPoint(vertex - 1)
										next_vertex_X, next_vertex_Y, z = linear_ring.GetPoint(0)
									if vertex > 0 and vertex < no_of_polygon_vertices - 1:
										vertex_position = "middle" 
										previous_vertex_X, previous_vertex_Y, z = linear_ring.GetPoint(vertex - 1)
										next_vertex_X, next_vertex_Y, z = linear_ring.GetPoint(vertex + 1)
								
								""" Calculate angle between previous vertex - current vertex and current vertex - next vertex from vincenty's inverse formula 
								( calculation on a spheroid ). """
								
								inverse_vincenty(flattening, major_axis, previous_vertex_Y, previous_vertex_X, current_vertex_Y, current_vertex_X )
								angle_prev = direction12
								
								inverse_vincenty(flattening, major_axis, current_vertex_Y, current_vertex_X, next_vertex_Y, next_vertex_X )
								angle_next = direction12
								
								""" Ensure that angles remain within 0-360 deg range. """
								
								if angle_prev < 0:
									angle_prev = angle_prev + 360
								if angle_next < 0:
									angle_next = angle_next + 360
								
								angle_prev_BP = angle_prev - 90 
								angle_next_BP = angle_next - 90 
								
								""" Ensure that buffer points are perpendicular to reference area and remain within 0-360 deg range. """
								
								if angle_prev_BP < 0:
									angle_prev_BP = angle_prev_BP + 360  
								if angle_next_BP < 0:
									angle_next_BP = angle_next_BP + 360 
								
								""" Add to list in which coordinates and angles for buffer vertices calculation are stored """
								
								vertices_angle_list.append([current_vertex_X, current_vertex_Y, angle_prev_BP, angle_prev, angle_next])
								
								""" Remember coordinates of the buffer polygon's first vertex - used later to close polygon """
								
								if vertex == 0:
									X_0 = current_vertex_X
									Y_0 = current_vertex_Y
									angle_0 = angle_prev_BP
									angle_prev_0 = angle_prev
									angle_next_0 = angle_next
								
								""" Calculate buffer points between angle_prev-90 and angle_next-90 (used for round buffer edges). """
								
								if angle_next_BP > angle_prev_BP:
									if (angle_next_BP) - (angle_prev_BP) <= 180:
										
										""" Angles between previous and next buffer point - generate buffer points which are outside the original polygon """
										
										angles = numpy.arange(angle_prev_BP, angle_next_BP, 8)
										for angle in angles:
											vertices_angle_list.append([current_vertex_X, current_vertex_Y, angle, angle_prev, angle_next])
									if (angle_next_BP) - (angle_prev_BP) > 180:
										
										""" Scissor intersection: angles between next buffer point and 360 degrees and between 0 degrees and previous buffer point - generate buffer points which are outside the original polygon """
										
										angles = numpy.arange(angle_next_BP, 360, 8)
										for angle in angles:
											vertices_angle_list.append([current_vertex_X, current_vertex_Y, angle, angle_prev, angle_next])
										angles = numpy.arange(0, angle_prev_BP, 8)
										for angle in angles:
											vertices_angle_list.append([current_vertex_X, current_vertex_Y, angle, angle_prev, angle_next])
								if angle_next_BP < angle_prev_BP:
									if angle_next_BP - angle_prev_BP > -180:
										
										""" Angles between next and previous buffer point - generate buffer points which are outside the original polygon """
										
										angles = numpy.arange(angle_next_BP, angle_prev_BP, 8)
										for angle in angles:
											vertices_angle_list.append([current_vertex_X, current_vertex_Y, angle, angle_prev, angle_next])
									if angle_next_BP - angle_prev_BP <= -180:
										
										""" Scissor intersection - angles between previous buffer point and 360 degrees and between 0 degrees and next buffer point - generate buffer points which are outside the original polygon """
										
										angles = numpy.arange(angle_prev_BP, 360, 8)
										for angle in angles:
											vertices_angle_list.append([current_vertex_X, current_vertex_Y, angle, angle_prev, angle_next])
										angles = numpy.arange(0, angle_next_BP, 8)
										for angle in angles:
											vertices_angle_list.append([current_vertex_X, current_vertex_Y, angle, angle_prev, angle_next])
								
								vertices_angle_list.append([current_vertex_X, current_vertex_Y, angle_next_BP, angle_prev, angle_next])
							
							""" Close polygon using the first vertex. """
							
							vertices_angle_list.append([X_0, Y_0, angle_0, angle_prev_0, angle_next_0])
							
							""" Calculate coordinates of buffer points for outer polygon. """
							
							direct_vincenty(flattening, major_axis, vertices_angle_list, buffer_distance_polygon_BNSC, 1)
							
							for buffer_vertex in buffer_vertices_list:
								
								""" Add point to ring geometry. """
								
								splitted_buffered_ring.AddPoint(buffer_vertex[0], buffer_vertex[1])
							
							""" Add ring to polygon geometry. """
							
							splitted_buffered_polygon.AddGeometry(splitted_buffered_ring)
							
							""" Eliminate unwanted holes due to self-intersections on outer boundary using a planar buffer of zero distance. """
							
							splitted_buffered_polygon = splitted_buffered_polygon.Buffer(0)
							
							""" Errors may occur during Buffer(0) so that two polygons are formed from one polygon dur to severe self-intersections. 
							This would lead to an invalid geometry which could not be added to the BNSC_union_polygon. To prevent this, all geometries 
							in the splitted_buffered_polygon are investigated, Buffered (0) again and then added to the BNSC_union_polygon. """
							
							if splitted_buffered_polygon.IsValid() == False:
							
								for linear_ring_splitted_buffered_polygon in splitted_buffered_polygon:
									
									""" Add linear_ring_splitted_buffered_polygon to new polygon. """
									
									polygon_part_splitted_buffered_polygon = ogr.Geometry(ogr.wkbPolygon)
									polygon_part_splitted_buffered_polygon.AddGeometry(linear_ring_splitted_buffered_polygon)
									polygon_part_splitted_buffered_polygon = polygon_part_splitted_buffered_polygon.Buffer(0)
									
									if polygon_part_splitted_buffered_polygon.IsValid() == False:
										print "Error due to severe self-intersection during buffering. Please check the shapefile for errors."
										if write_logfile == True and multicore_operation == False:
											logfile.write("Error due to severe self-intersection during buffering. Please check the shapefile for errors.")
											logfile.flush()
										exit()
									
									BNSC_union_polygon = BNSC_union_polygon.Union(polygon_part_splitted_buffered_polygon)
							
							if splitted_buffered_polygon.IsValid() == True:
								
								BNSC_union_polygon = BNSC_union_polygon.Union(splitted_buffered_polygon)
						
						""" Special case: If only one research area with hole(s) is investigated, iteration in union_polygon would not consider polygon1, polygon2, 
						polygon3, etc., but ring1, ring2, ring3, etc. As the inner ring is already considered during iteration (because it is assumed that in 'for area in union_polygon' 
						area is a research area and not a linear ring), the function has to be stopped here. Otherwise, every further iteration in union polygon would consider 
						the linear rings which have already been considered in the function. This would lead to too many resulting polygons (rings*no_of_craters and not no_of_craters). """
						
						if len(Area_IDs) == 1 and number_of_inner_polygons == 0 and number_of_holes >= 1:
							break
					
					""" Project union_polygon & BNSC_union_polygon back to original spatial reference (to be used in NSC_BNSC_exclude_craters function). """
					
					BNSC_union_polygon.Transform(geogr_reprojection_to_proj_reprojection)
					BNSC_union_polygon.Transform(proj_reprojection_to_geogr)
					
					union_polygon.Transform(geogr_reprojection_to_proj_reprojection)
					union_polygon.Transform(proj_reprojection_to_geogr)
					
					if sr.IsProjected(): 
						union_polygon.Transform(geog_to_proj)
						BNSC_union_polygon.Transform(geog_to_proj)
		
		#########################
		#						#
		# END OF BNSC BUFFERING	#
		#						#
		#########################
		
		if approach == "NSC":
			union_polygon.Transform(sr_to_eq_area)
			geodesic_area = union_polygon.GetArea()/1000000
			union_polygon.Transform(eq_area_to_sr)
			
			""" Exclude craters outside reference area for CSFD measurement (NSC only). """
			
			if distance_crater_area > 0:
				crater_index += 1
				continue
				
			if generate_polygon_file == True:
				polygon_feature.SetGeometry(union_polygon) 
				polygon_feature.SetField('Size_sq_km', geodesic_area)
				polygon_feature.SetField('crater_ID', original_crater_id)
				layer_polygon.CreateFeature(polygon_feature)
			
			crater_diam = craters_for_counting_list[crater_index][1]
			crater_X = craters_for_counting_list[crater_index][2]
			crater_Y = craters_for_counting_list[crater_index][3]
			crater_area_list.append([crater_diam, crater_X, crater_Y, crater_diam * bufferfactor_crater, geodesic_area])
		
		if approach == "BNSC":
			
			""" 'If BNSC_union_polygon' condition is used to only write geometry to shapefile if buffering took place. Otherwise 
			 extra geometries would be added here during multicore processing.  """
			
			if 'BNSC_union_polygon' in locals() or 'BNSC_union_polygon' in globals():
				BNSC_union_polygon.Transform(sr_to_eq_area)
				geodesic_area = BNSC_union_polygon.GetArea()/1000000
				
				if geodesic_area == 0:
					
					""" This is a special case. If two large craters cut the reference area (the reference area is completely superposed 
					by crater + ejecta blanket), the smaller crater will not be used to buffer the area even though it cuts the reference 
					area and is located distant enough to the larger crater to be (theoretically) included in the counting. Since the area is 
					eliminated by the larger crater + ejecta blanket already, there is nothing to buffer (geodesic area size = 0). 
					This is a matter of dimensions. The suitability of the reference area for CSFD measurements should be investigated again. """
					
					continue
					
				BNSC_union_polygon.Transform(eq_area_to_sr)
				
				if generate_polygon_file == True:
					polygon_feature.SetGeometry(BNSC_union_polygon) 
					polygon_feature.SetField('Size_sq_km', geodesic_area)
					polygon_feature.SetField('crater_ID', original_crater_id)
					layer_polygon.CreateFeature(polygon_feature)
				
				crater_diam = craters_for_counting_list[crater_index][1]
				crater_X = craters_for_counting_list[crater_index][2]
				crater_Y = craters_for_counting_list[crater_index][3]
				crater_area_list.append([crater_diam, crater_X, crater_Y, crater_diam * bufferfactor_crater, geodesic_area])
				
				del BNSC_union_polygon
		
		crater_index += 1
	crater_area_out_q.put(crater_area_list)
	print process_count, "Done."
	
	""" Detailled logfile is only written in single-core mode """
	
	if write_logfile == True and multicore_operation == False:
		logfile.write("Done.\n")
		logfile.flush()

""" Calculate fractions from craters and individual buffer areas. Until here, if more than one (n) research area was selected, the crater_area_list 
 stores informartion about every crater for every individual research area (crater1 - geodesic_area_1, crater1 - geodesic_area_2). As fractions 
 are calculated using 'fraction_crater1 = total_area_size_non_buffered / (crater1_geodesic_area1 + crater1_geodesic_area2)', duplicate craters 
 must be identified in the list according to their centroid coordinates and their respective geodesic_area values must be added up. """

def create_crater_fractions_list(crater_area_list):
	from collections import defaultdict
	global crater_fractions_list
	
	crater_fractions_list = []
	crater_fractions_dict = defaultdict(int)
	
	for crater_diam, crater_x, crater_y, ejecta_diam, buffer_area in crater_area_list:
		crater_fractions_dict[(crater_diam, crater_x, crater_y, ejecta_diam)] += buffer_area
	
	crater_fractions_list_temp = [[crater_diam, crater_x, crater_y, ejecta_diam, buffer_area] for (crater_diam, crater_x, crater_y, ejecta_diam), buffer_area in crater_fractions_dict.items()]
	
	for crater_diam, crater_x, crater_y, ejecta_diam, buffer_area in crater_fractions_list_temp:
		fraction = total_area_size / buffer_area
		crater_fractions_list.append([crater_diam, fraction, crater_x, crater_y, ejecta_diam, buffer_area])


""" Write SCC/DIAM file header """

def write_crater_stats_file_header():
	global crater_stats_file
	
	now = datetime.datetime.now()
	date = str(now.day) + "." + str(now.month) + "." + str(now.year)
	
	crater_stats_file = open(path_to_craterstats_outfile, "w")
	
	if approach == "BCC":
		methodology_name = "Buffered crater count"
	
	if approach == "TRAD":
		methodology_name = "Traditional crater counting approach"
	
	if approach == "NSC":
		methodology_name = "Non-sparseness correction"
	
	if approach == "BNSC":
		methodology_name = "Buffered non-sparseness correction"
	
	if outfile_type == "SCC":
		crater_stats_file.write("# Spatial crater count for Craterstats - " + methodology_name + "\n"\
								"#\n"\
								"# Date of measurement = " + str(date) + "\n"\
								"#\n"\
								"# Ellipsoid axes: \n"\
								"a-axis radius = " + str(major_axis/1000) + " <km>\n"\
								"b-axis radius = " + str(minor_axis/1000) + " <km>\n"\
								"c-axis radius = " + str(minor_axis/1000) + " <km>\n"\
								"coordinate_system_name = " + str(sr_name) + "\n"\
								"#\n"\
								"# area_shapes:\n"\
								"unit_boundary = {vertex, sub_area, tag, lon, lat\n")
	if outfile_type == "DIAM":
		crater_stats_file.write("# Diam file for Craterstats - " + methodology_name + "\n"\
								"#\n"\
								"# Date of measurement export = " + str(date) + "\n"\
								"#\n")

""" Write area information to SCC/DIAM file. """

def write_crater_stats_file_area(n_area, area_feature, area_feature2, area_feature_SCC_File, union_polygon, Area_Names, Area_ID, an, write_logfile, logfile):
	global previous_vertex_count, vertices_list, inner_vertices_list, ring_index
	
	""" Get vertices from original polygon (area_feature_SCC_File) - to write original vertices to SCC File and segmentized polygon (area_feature) - to work with during processing. """
	
	vertices_list = []
	inner_vertices_list = []
	vertices_list_SCC_File = []
	inner_vertices_list_SCC_File = []
	
	area_shape = area_feature.GetGeometryRef()
	area_polygon = area_shape.GetGeometryRef(0)
	
	area_shape_SCC_File = area_feature_SCC_File.GetGeometryRef()
	area_polygon_SCC_File = area_shape_SCC_File.GetGeometryRef(0)

	""" Transform vertices to geographic coordinates if spatial reference is projected (used for vertices in SCC/DIAM file). """

	if sr.IsProjected:
		area_shape.Transform(proj_to_geog)
		area_shape_SCC_File.Transform(proj_to_geog)

	""" Get area size. """
	
	area_shape2 = area_feature2.GetGeometryRef()
	
	if area_shape2.GetGeometryCount() == 1: # Polygon without hole
		area_polygon2 = area_shape2.GetGeometryRef(0)
		get_area_size(area_polygon2)
	
	if area_shape2.GetGeometryCount() > 1: # Polygon with hole
		get_area_size(area_shape2)
	
	##################################
	# ORIGINAL VERTICES FOR SCC FILE #
	##################################
	
	no_of_polygon_vertices = area_polygon_SCC_File.GetPointCount()
	
	if no_of_polygon_vertices == 0:
		print "THIS IS STRANGE. I CAN'T READ ANY VERTICES. IS THIS A MULTIPART FEATURE? IF SO, PLEASE EXPLODE MULTIPART FEATURE AND TRY AGAIN."
		exit()
	
	""" Get outer ring vertices. """
	
	print "Area ", Area_ID, ":\n", "Number of Polygon vertices: ", no_of_polygon_vertices, "\n"
	if write_logfile == True:
		logfile.write("Area " + str(Area_ID) + ":\n" + "Number of Polygon vertices: " + str(no_of_polygon_vertices) + "\n")

	for vertex in xrange(no_of_polygon_vertices):
		current_vertex_X, current_vertex_Y, z = area_polygon_SCC_File.GetPoint(vertex)
		vertices_list_SCC_File.append([current_vertex_X, current_vertex_Y, n_area+ring_index])
	
	""" Get inner ring vertices. """
	
	if area_shape_SCC_File.GetGeometryCount() > 1:
		no_of_inner_rings = 0
		number_of_inner_polygons = area_shape_SCC_File.GetGeometryCount() - 1
		for inner_polygon_count in range(number_of_inner_polygons):
			ring_index += 1
			no_of_inner_rings += 1 
			inner_polygon = area_shape_SCC_File.GetGeometryRef(inner_polygon_count+1)
			no_of_inner_polygon_vertices = inner_polygon.GetPointCount()
			for vertex in xrange(no_of_inner_polygon_vertices):
				current_vertex_X, current_vertex_Y, z = inner_polygon.GetPoint(vertex)
				inner_vertices_list_SCC_File.append([current_vertex_X, current_vertex_Y, n_area + ring_index])
		
		""" Since n_area will be + 1 for the next area, number of ring_index will be subtracted by current number of rings. Otherwise rings are not correctly enumerated in SCC file. """
		
		ring_index -= no_of_inner_rings
	
	if outfile_type == "SCC":
		crater_stats_file.write("#\n"\
								"# Area_name " + str(n_area) + " = " + str(an) + "\n"\
								"# Area_size " + str(n_area) + " = " + str(Area_Size) + " <km^2>\n"\
								"#\n")
		if n_area == 1:
			vertex_count = 1
		else:
			vertex_count = previous_vertex_count
			
		for area_vertex in vertices_list_SCC_File:
			crater_stats_file.write(str(vertex_count) + "\t" + str(area_vertex[2]) + "\text\t" + str(area_vertex[0]) + "\t" + str(area_vertex[1]) + "\n")
			vertex_count += 1
		
		if len(inner_vertices_list_SCC_File) > 0:
			for area_vertex in inner_vertices_list_SCC_File:
				crater_stats_file.write(str(vertex_count) + "\t" + str(area_vertex[2]) + "\tint\t" + str(area_vertex[0]) + "\t" + str(area_vertex[1]) + "\n")
				vertex_count += 1
		
		previous_vertex_count = vertex_count
		
	if outfile_type == "DIAM":
		crater_stats_file.write("# " + str(an) + " = " + str(Area_Size) + " <km^2>\n")
	
	######################################
	# SEGMENTIZED VERTICES FOR BUFFERING #
	######################################
	
	""" Get outer ring vertices. """
	
	no_of_polygon_vertices = area_polygon.GetPointCount()
	for vertex in xrange(no_of_polygon_vertices):
		current_vertex_X, current_vertex_Y, z = area_polygon.GetPoint(vertex)
		vertices_list.append([current_vertex_X, current_vertex_Y, n_area+ring_index])

	if no_of_polygon_vertices == 0:
		print "THIS IS STRANGE. I CAN'T READ ANY VERTICES. IS THIS A MULTIPART FEATURE? IF SO, PLEASE EXPLODE MULTIPART FEATURE AND TRY AGAIN."
		exit()
	
	""" Get inner ring vertices. """
	
	if area_shape.GetGeometryCount() > 1:
		number_of_inner_polygons = area_shape.GetGeometryCount() - 1
		for inner_polygon_count in range(number_of_inner_polygons):
			ring_index += 1
			inner_polygon = area_shape.GetGeometryRef(inner_polygon_count+1)
			no_of_inner_polygon_vertices = inner_polygon.GetPointCount()
			for vertex in xrange(no_of_inner_polygon_vertices):
				current_vertex_X, current_vertex_Y, z = inner_polygon.GetPoint(vertex)
				inner_vertices_list.append([current_vertex_X, current_vertex_Y, n_area + ring_index])

""" Write results from CSFD measurements to SCC/DIAM file. """

def write_crater_stats_file_craters(crater_fractions_list):
	
	if outfile_type == "SCC":
		crater_stats_file.write("}\n"\
								"#\n"\
								"Total_area = " + str(total_area_size) + " <km^2>\n"\
								"#\n"\
								"# crater_diameters: \n"\
								"crater = {diam, fraction, lon, lat, ejecta_diam, buffer_area\n")
								
		for crater_diam, fraction, crater_x, crater_y, ejecta_diam, buffer_area in crater_fractions_list:
			crater_stats_file.write(str(crater_diam) + "\t" + str(fraction) + "\t" + str(crater_x) + "\t" + str(crater_y) + "\t" + str(ejecta_diam) + "\t" + str(buffer_area) + "\n")
		crater_stats_file.write("}")
		crater_stats_file.close()
	
	if outfile_type == "DIAM":
		crater_stats_file.write("#\n"\
								"# Area, <km^2>\n"\
								"area = " + str(total_area_size) + "\n"\
								"#\n"\
								"# crater_diameters: \n"\
								"crater = {diameter, fraction, lon, lat, ejecta_diam, buffer_area\n")
		for crater_diam, fraction, crater_x, crater_y, ejecta_diam, buffer_area in crater_fractions_list:
			crater_stats_file.write(str(crater_diam) + "\t" + str(fraction) + "\t" + str(crater_x) + "\t" + str(crater_y) + "\t" + str(ejecta_diam) + "\t" + str(buffer_area) + "\n")
		crater_stats_file.write("}")
		crater_stats_file.close()

""" Write results from CSFD measurements with traditional crater counting to SCC/DIAM file. """

def write_crater_stats_file_craters_TRAD(all_craters):
	
	if outfile_type == "SCC":
		crater_stats_file.write("}\n"\
								"#\n"\
								"Total_area = " + str(total_area_size) + " <km^2>\n"\
								"#\n"\
								"# crater_diameters: \n"\
								"crater = {diam, fraction, lon, lat, topo_scale_factor\n")
								
		for crater_id, crater_diam, crater_x, crater_y in all_craters:
			crater_stats_file.write(str(crater_diam) + "\t" + "1" + "\t" + str(crater_x) + "\t" + str(crater_y) + "\t" + "1" + "\n")
		crater_stats_file.write("}")
		crater_stats_file.close()
	
	if outfile_type == "DIAM":
		crater_stats_file.write("#\n"\
								"# Area, <km^2>\n"\
								"area = " + str(total_area_size) + "\n"\
								"#\n"\
								"# crater_diameters: \n"\
								"crater = {diameter, fraction, lon, lat, topo_scale_factor\n")
		for crater_id, crater_diam, crater_x, crater_y in all_craters:
			crater_stats_file.write(str(crater_diam) + "\t" + "1" + "\t" + str(crater_x) + "\t" + str(crater_y) + "\t" + "1" + "\n")
		crater_stats_file.write("}")
		crater_stats_file.close()

""" Show GUI """
	
def user_interface():
	app = QApplication(sys.argv)
	window = GUI_window()
	window.show()
	sys.exit(app.exec_())

if __name__ == '__main__':
	multiprocessing.freeze_support()
	
	""" Suppress error messages which may occur during multi-core operation """
	
	ctypes.windll.kernel32.SetErrorMode(3)
	
	""" Call GUI """
	
	user_interface()
	

