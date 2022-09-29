## Essential Climate Variable Analysis Tool V2
## Version 1.6.2 - 29/09/20 - 09:45
## Version 2.0   - 22/01/21 - 12:28
## Version 3     - 05/05/21 - 13:51
## Version 4     - 30/06/21
#	- Updated Calculations for efficiency and readability
#	- Added multiple types of graph solutions
#	- Upgraded to using np.arrays
#	- Added .pci as operand functionality (for batch submissions)
## Daniel Westwood - daniel.westwood@stfc.ac.uk

## NetCDF4
from netCDF4 import Dataset

## Matplotlib Packages
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Path, PathPatch
from matplotlib import colors
import matplotlib.widgets as wg
import matplotlib as m
import matplotlib.ticker as ticker
m.use('TkAgg') ## Faster rendering

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.path as Path

## System
import os
import sys
from datetime import datetime
from getopt import getopt

## Shape
import shapely.affinity
from shapely.geometry import Point, Polygon

## Numpy
import numpy as np
import numpy.ma as ma
import warnings
import math

# Standard dw package tools
STD_PY_LOC = '/home/users/dwest77/dwest77_RSG/std_py'
try:
	sys.path.append(STD_PY_LOC)
	import find_files as ff
	import pmath as pm
	import output_data as od
except: 
	print('ImportError: Std_py library missing: requires find_files.py, pmath.py and output_data.py as minimum')
	sys.exit()

## --- Global Variables --- ##
ECV_IDS = ['erb', 'srb', 'cld']

# Dictionary giving alignment to title
TITLE_DICT	= { 'fontsize':'large',
				  'fontweight':'normal',
				  'color':'black',
				  'verticalalignment':'baseline',
				  'horizontalalignment':'center'}
# Transform integer month value to different formats
MONTHS		= { 1:'01', 2:'02', 3:'03', 4:'04', 5:'05', 6:'06',
				  7:'07', 8:'08', 9:'09', 10:'10', 11:'11', 12:'12'}
CAL_MONTHS	= { 1:'Jan',2:'Feb',3:'Mar',4:'Apr',5:'May',6:'Jun',
				  7:'Jul',8:'Aug',9:'Sep',10:'Oct',11:'Nov',12:'Dec'}  
				  
MONTH_CONTENT = [31,28,31,30,31,30,31,31,30,31,30,31]

# Imported data from config giles
FILE_PATHS	= ff.get_from_file('config/','file_std_roots','.con')
ECV_SET_DICT  = ff.get_from_file('config/','ecv_sets','.con')

# General Variables
VERBOSE	   = False
LENFIG = 4
WIDFIG = 8

# Single set of figures with map and graph

# 1 - Single Fig (Map or Graph)
# 2 - Vertical Figs (Map and Graph)
# 3 - Horizontal Figs (Maps)

# 4 - Triple Fig (Maps and Graph)
# 5 - Quadruple Fig (Maps and Graphs)

MAP_FORMATS = [
[ [0.15, 0.2, 0.7, 0.7] ], # Single Fig
[ [0.05, 0.5, 0.9, 0.45] ], # Vertical Fig
[ [0.02,0.2,0.55,0.7], [0.57,0.2,0.45,0.7] ], # Horizontal Fig
[ [0.02,0.5,0.45,0.45], [0.52,0.5,0.45,0.45] ], # Triple Fig
[ [0.05,0.55,0.4,0.4], [0.55,0.55,0.4,0.4] ] ] # Quadruple Fig

CBAR_FORMATS = [
[ [0.25,0.12,0.5,0.02] ],
[ [0.8,0.4,0.15,0.15] ],
[ [0.15,0.15,0.3,0.02], [0.65,0.15,0.3,0.02] ], # Triple Fig - No HFig Yet
[ [0.12,0.47,0.25,0.02], [0.62,0.47,0.25,0.02] ], # Triple Fig
[ [0.12,0.5,0.25,0.02], [0.62,0.5,0.25,0.02] ] ] # Quadruple Fig


GRAPH_FORMATS = [
[ [0.15, 0.2, 0.7, 0.7] ], # Single Fig
[ [0.1,0.15,0.8,0.25] ], # Vertical Fig
[ [0.02,0.5,0.45,0.45], [0.52,0.5,0.45,0.45] ], # Triple Fig - No HFig Proper Yet
[ [0.1,0.15,0.8,0.25] ], # Triple Fig
[ [0.05,0.1,0.4,0.3], [0.55,0.1,0.4,0.3] ] ] # Quadruple Fig

RAL_LOGOS = [
[ [0.04,0.1,0.15,0.15] ], # Single
[ [0.05,0.4,0.15,0.15] ], # Vertical
[ [0.05,0.4,0.15,0.15] ], # No HFig Proper yet
[ [0.05,0.4,0.15,0.15] ], # No Triple Fig proper
[ [0.05,0.4,0.15,0.15] ] ] # No Quadruple Fig proper

CCI_LOGOS = [
[ [0.78,0.1,0.2,0.15] ], # Single
[ [0.8,0.4,0.15,0.15] ], # Vertical
[ [0.8,0.4,0.15,0.15] ], # No HFig Proper yet
[ [0.8,0.4,0.15,0.15] ], # No Triple Fig proper
[ [0.8,0.4,0.15,0.15] ] ] # No Quadruple Fig proper

## --- Standard Tools --- ##

def identify_c3s_instrument(year, month, instrument): # ECV Class
	# Takes year and month as input
	# Determines which instrument the data was recorded by due to the date range
	# Assumes c3s instrument - for other instruments, need more identify functions
	
	# Operating Ranges of three known c3s instruments
	atsr2 = ((year > 1995 and year < 2002) or (year == 1995 and month >= 6) or (year == 2002 and month <= 4))
	aatsr = ((year > 2002 and year < 2012) or (year == 2002 and month >= 5) or (year == 2012 and month <= 4))
	slstr = (year >= 2017)
	
	if slstr:
		if 'slstr' in instrument:
			return instrument
		else:
			return 'slstra'
	elif aatsr:
		return 'aatsr'
	elif atsr2:
		return 'atsr2'
	else:
		return None
		
def console_intro(): # keep
	# Simple console interface - welcome message
	
	print('''

 --- > Essential Climate Variable Analysis tool V4 < ---		 

 --- Dev: Daniel Westwood (29/09/2020)
 --- Updates: V2 (19/01/2021)
 ---		  V3 (05/05/2021)
 ---		  V4 (01/07/2021)
 
 Required: config files
		   pci file
  
	''')
	 
## --- End standard tools --- ##

## --- ECV Reading class --- ##

class ecv_reading(): # Single data scene
# Set general properties of ecv reading object
	def __init__(self, year, month, day, timeinc, isanom):
		# Defines properties of single ecv/scene
		# Runs import_data to obtain properties from corresponding netcdf file
		
		# Single value properties:
		# - year, month, day
		# - isanom
		# - date_index
		# - timeinc
		
		# Double value (multi-ecv) properties:
		# - family, instrument
		# - daynight, landsea
		# - ecv/variable
		# - label, spc,  pre_units
		# - ncf, basedata, baseunits
		# - is_valid
		# - lat/lon values (6)
		
		self.year	   = year
		self.month	  = month
		self.day		= day
		self.isanom	 = isanom
		self.date_index = None
		self.timeinc	= timeinc # Monthly or Daily

# Set other properties of ecv reading object
	def set_ecv_properties(self,ecvs, ecv_sets, families, instruments, labels, spcs, pre_units, vinfos):
		self.ecvs        = ecvs			   # actual variable to study
		self.ecv_sets    = ecv_sets	   # variable superset (srb, erb, cld)
		self.families    = families	   # family of instruments (c3s, cris, iasi)
		self.instruments = instruments # instrument subset of family (atsr2, aatsr, slstr etc.)
		self.labels      = labels		   # Additions to variable names (adjusted 'A/B')
		self.spcs        = spcs			   # Special case keys
		self.pre_units   = pre_units	 # Replacement/supplementary units
		self.vinfos      = vinfos
		
		self.is_valid  = [True for int in range(len(self.ecvs))]
		self.basedata  = [[] for int in range(len(self.ecvs))]
		self.baseunits = [None for int in range(len(self.ecvs))]
		self.latlist   = [None for int in range(len(self.ecvs))]
		self.lonlist   = [None for int in range(len(self.ecvs))]
		self.lat_min   = [None for int in range(len(self.ecvs))]
		self.lat_max   = [None for int in range(len(self.ecvs))]
		self.lon_min   = [None for int in range(len(self.ecvs))]
		self.lon_max   = [None for int in range(len(self.ecvs))]
		
		self.all_imports = []
	
	def set_map_properties(self, daynights, landseas):
		if daynights != None:
			self.daynights = daynights
		else:
			self.daynights = ['day' for int in range(len(self.ecvs))]
		self.landseas = landseas
	
# Sorting file paths and imports
	def sort_file_imports(self):
		base_path = ''
		
		# all imports [ [ 0, path/for/import, 'lat', 'lon],
		for index in range(len(self.ecvs)):
			daynight_original = self.daynights[index]
			if self.daynights[index] == 'day+night':
				dn = ['day','night']
			else:
				dn = [self.daynights[index]]
			for daynight in dn:
				
				self.daynights[index] = daynight
				variable = self.ecvs[index] + self.labels[index]
				if self.instruments[index] == 'saved':
					self.get_presaved_path(variable, index)
				elif self.families[index] == 'esa_cci':
					self.get_esa_cci_path(index)
				elif self.families[index] == 'iasi':
					self.get_iasi_path(index)
				elif 'l2grid' in self.families[index]:
					self.get_gridded_path(index)
				elif 'l2os' in self.families[index]:
					self.get_os_path(index)
				elif self.families[index] == 'cris':
					self.get_cris_path(index)
				elif self.families[index] == 'cs':
					self.get_custom_path(index)
				else:
					 print('Warning: Unrecognised instrument family: {} - add pathfinder function at ln 258'.format(self.families[index]))
			self.daynights[index] = daynight_original
		
		self.sort_paths()
		self.concat_all_basedata()
			
	def sort_paths(self):
		# Needs list of all indexes of variables
		# Corresponding list of paths, 1 path per index
		main_paths = {}
		main_paths_list = []
		main_lats = []
		main_lons = []
		latlonincs = []
		for entry in self.all_imports:
			index = entry[0]
			path = entry[1]
			try:
				main_paths[path].append(index)
			except:
				main_paths[path] = [index]
				main_paths_list.append(path)
				main_lats.append(entry[2])
				main_lons.append(entry[3])
				latlonincs.append(entry[4])
		for num, path in enumerate(main_paths_list):
			print('Importing from ',path)
			self.import_set_of_variables(path, main_paths[path], main_lats[num], main_lons[num], latlonincs[num])
	
# Extract and assemble data from different files
	def concat_all_basedata(self):
		for index in range(len(self.ecvs)):
			if len(self.basedata[index]) == 0:
				self.is_valid[index] = False
				self.basedata[index] = np.nan
			else:
				'''
				# Old method of averaging
				bd_sum = 0
				for repeat in range(len(self.basedata[index])):
					bd_sum += np.array(self.basedata[index][repeat])
				bd_sum = bd_sum/len(self.basedata[index])
				'''
				# New method of averaging
				self.basedata[index] = np.array(self.basedata[index])
				self.basedata[index][self.basedata[index] > 99999] = np.nan
				self.basedata[index][self.basedata[index] < -99999] = np.nan
				
				bd_sum = np.nanmean(self.basedata[index], axis=0)
				self.basedata[index] = bd_sum
	
	def import_set_of_variables(self, path, indexes, lat, lon, latloninc):
		try:
			ncf = Dataset(path, 'r', format = 'NETCDF4')
		except:
			print('Error: File not found {} - skipping'.format(path))
			indexes = []
		for index in indexes:
			variable = self.ecvs[index]
			self.latlist[index]   = ncf.variables[lat][:]
			self.lonlist[index]   = ncf.variables[lon][:]
			if self.spcs[index] == '20_nh3tpw': # Change later
				nh3 = ncf['nh3']
				tpw = ncf['tpw']
				new = np.array(nh3) - (0.0016*np.array(tpw) + 0.21)/1000 # UK region values
				self.basedata[index].append(new)
			elif self.spcs[index] == '21_nh3tpw': # Change later
				tpw = ncf['tpw']
				vza = ncf['satzen']
				new = np.array(tpw) / np.cos(np.array(vza) * (math.pi/180))
				self.basedata[index].append(new)
			elif self.spcs[index] == 'nh3tpwvza':
				print('Applied adj')
				nh3 = ncf['nh3']
				tpw = ncf['tpw']
				new = np.array(nh3) - (0.003*np.array(tpw) -0.06) # IASI gridded L2 value
				self.basedata[index].append(new)
				
			elif self.spcs[index] == '22_tpwvza':
				tpw = ncf['tpw']
				vza = ncf['satzen']
				new = np.array(tpw) / np.cos(np.array(vza) * (math.pi/180))
				self.basedata[index].append(new)	
				
			elif self.spcs[index] == '23_nh3tpwvza_nh3':
				nh3tpwvza = np.array(ncf['nh3tpwvza'])
				nh3 = np.array(ncf['nh3'])
				diff = nh3tpwvza-nh3
				self.basedata[index].append(diff)
			
			elif self.spcs[index] != 'None' and self.spcs[index] != 'XT' and self.spcs[index] != '':
				#variable = variable + self.labels[index]
				ncf_var = ncf.variables[variable][:]
				spc = int(self.spcs[index])
				self.basedata[index].append(np.array(ncf_var[spc]))
			else:
				nvar = np.array(ncf.variables[variable][:])
				self.basedata[index].append(np.array(ncf.variables[variable][:]))
			try:
				self.baseunits[index] = ncf.variables[variable].getncattr('units')
			except:
				self.baseunits[index] = ''
			if latloninc:
				self.lat_min[index] = ncf.getncattr('geospatial_lat_min')
				self.lat_max[index] = ncf.getncattr('geospatial_lat_max')
				self.lon_min[index] = ncf.getncattr('geospatial_lon_min')
				self.lon_max[index] = ncf.getncattr('geospatial_lon_max')
			else:
				self.lat_min[index] = self.latlist[index][0]
				self.lat_max[index] = self.latlist[index][len(self.latlist[index])-1]
				self.lon_min[index] = self.lonlist[index][0]
				self.lon_max[index] = self.lonlist[index][len(self.lonlist[index])-1] 
	
# Assembly of file paths by instrument			 
	def get_presaved_path(self, variable, index):
		base_path = '/home/users/dwest77/Documents/ECV_Images/CRIS/{}'.format(variable)
		main_path = base_path + '/{}_{}_{}_{}_{}{}.ncf'.format(self.families[index], variable, self.daynights[index], self.landseas[index], self.year, format(self.month,'02d'))
		new_entry = [index, main_path, 'lat','lon', False]
		self.all_imports.append(new_entry)
		
	def get_gridded_path(self, index):
		# Construct custom file path
		month = format(self.month, '02d')
		instrt = self.families[index].replace('_l2grid','')
		main_path = '/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/{}_l2/gridded_l2/{}/{}/{}/'.format(instrt, self.vinfos[index].split('_')[0], self.year, month)
		filename = '{}_gridded_l2_{}{}{}_{}.nc'.format(instrt, self.year, month, self.daynights[index], self.vinfos[index]) # Change sub version
		new_entry = [index, main_path+filename, 'lat','lon', False]
		self.all_imports.append(new_entry)
		
	def get_custom_path(self, index):
		# Construct custom file path
		daystring = format(self.day,'02d')
		month = format(self.month, '02d')
		day = format(self.day, '02d')
		main_path = '/home/users/dwest77/Documents/IASI_os_l3/NH3/{}/{}/IASI_nh3_Daily_cv9_v12_{}_{}_{}_day.nc'.format(self.year, month,
																													   self.year, month, day)
		new_entry = [index, main_path, 'lat','lon', False]
		self.all_imports.append(new_entry)		

	def get_iasi_path(self, index):   
		# Construct iasi file path
		iasi_root_file = FILE_PATHS['iasi_root_file']	 
		
		if int(self.year) > 2016:
			iasi_root_dir = FILE_PATHS['iasi_post16_dir']
		else:
			iasi_root_dir = FILE_PATHS['iasi_pre16_dir']
			
		ym = '{}{}'.format(self.year, MONTHS[self.month])
		main_path = iasi_root_dir + iasi_root_file.format(ym,self.daynights[index])
		new_entry = [index, main_path, 'latitude','longitude', False]
		self.all_imports.append(new_entry)
		
	def get_os_path(self, index):
		# Construct iasi os file path
		month = format(self.month, '02d')
		instrt = self.families[index].replace('_l2os','')
		main_path = '/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/{}_l2/oversampled_l2/{}/{}/{}/'.format(instrt, self.vinfos[index].split('_')[0], self.year, month)
		filename = '{}_oversampled_l2_{}{}_{}_{}.nc'.format(instrt, self.year, month, self.daynights[index], self.vinfos[index]) # Change sub version
		new_entry = [index, main_path+filename, 'lat','lon', False]
		self.all_imports.append(new_entry)
  
	def get_cris_path(self, index):
		# Construct cris file path
		cris_root_dir = FILE_PATHS['cris_root_dir']
		cris_root_file = FILE_PATHS['cris_root_file']
		ym = '{}{}'.format(self.year,MONTHS[self.month])
		main_path = cris_root_dir + cris_root_file.format(ym,self.daynights[index])
		new_entry = [index, main_path, 'latitude','longitude', False]
		self.all_imports.append(new_entry)
	   
	def get_esa_cci_path(self, index):
		# AATSR, ATSR2, SLSTRA, SLSTRB, SLSTRAB
		
		instrument = identify_c3s_instrument(self.year, self.month, self.instruments[index])
		doreplace = False
		if instrument == 'slstrab':
			instrument = 'slstra'
			doreplace = True
		elif instrument == None:
			self.is_valid[index] = False
			print('File outside of time range - {}/{} - skipping'.format(self.year, self.month))
			return None
		# time increment part
		tinc = self.timeinc.upper()
		time_inc_file = tinc + '-'
		time_inc_path = tinc.lower() + '/'
		file_date = '{}{}'.format(self.year,MONTHS[self.month])
		path_date = '{}/{}/'.format(self.year, MONTHS[self.month])
		
		if self.timeinc == 'Daily':
			day = format(self.day,'02d')
			file_date += '{}'.format(day)
			
		c3s_root_dir = FILE_PATHS['c3s_root_dir']
		c3s_root_file = FILE_PATHS['c3s_root_file']
		ext_path = FILE_PATHS['{}_path'.format(instrument)]
		
		ext_files = FILE_PATHS['{}_file'.format(instrument)] 
		# replace 3a with 3b
		if type(ext_files) != list:
			ext_files = [ext_files]
			
		# Fully assembled file dir
		main_dir = c3s_root_dir + self.ecv_sets[index] + ext_path + time_inc_path + path_date
		found_valid = False
		file_count = 0
		
		# Check all file names in list of possible names
		while not found_valid and file_count < len(ext_files):

			# Fully assembled file name
			main_file = c3s_root_file + time_inc_file + FILE_PATHS[self.ecv_sets[index]] + ext_files[file_count].format(file_date)
			main_path = main_dir + main_file
			if os.path.isfile(main_path):
				found_valid = True
			file_count += 1
		if not found_valid:
			print('File not found for {}/{} - {} - skipping'.format(self.year, self.month, self.index))
			self.is_valid[index] = False
			return None
		#
		new_entry = [index, main_path, 'lat','lon', True]
		self.all_imports.append(new_entry)
		if doreplace:
			new_entry = [index, main_path.replace('3a','3b'), 'lat','lon', True]
			all_imports.append(new_entry)
			
## --- End ECV reading class --- ##

## --- ECV Analysis class --- ##

class ecv_analysis():
# Init function for user specified inputs
	def __init__(self, plotconfig):
		# Create ecv reading list of ecv instances
		# Determine plot config info (default settings or from config file)
		   
		# -- Bools -- #
		self.use_ral = False
		self.use_cci = False
		
		self.output_figures = plotconfig['OutFigs'] # Temp
		
		# -- Strings -- # 
		self.time_format	= plotconfig['TimeFormat'] # Single, Range, Days, Months
		self.time_increment = plotconfig['TimeIncrement'] # Monthly/Daily
		self.calculate	  = plotconfig['Calculate']
		self.color		  = plotconfig['Color']
		
		# -- Nums -- #
		self.timegroupsize = int(plotconfig['TimeGroupSize']) # Number of readings to group for use in graph
		self.fignum		= int(plotconfig['FigureNumber'])
		self.alpha		 = float(plotconfig['Alpha'])
		
		# --- Time parameters --- #
		self.year_arr	   = [int(yr) for yr in plotconfig['Years']]
		self.month_arr	  = [int(mth) for mth in plotconfig['Months']]
		if self.time_increment == 'Daily':
			self.day_arr		= [int(day) for day in plotconfig['Days']]
		else:
			self.day_arr = ''
		if self.calculate == 'Anomaly':
			self.anomyear_arr   = [int(yr) for yr in plotconfig['AnomalyYears']]
			self.anommonth_arr  = [int(mth) for mth in plotconfig['AnomalyMonths']]
			if self.time_increment == 'Daily':
				self.anomdays_arr = [int(day) for day in plotconfig['AnomalyDays']]
			else:
				self.anomdays_arr = ''
		else:
			self.anomyear_arr   = ''
			self.anommonth_arr  = ''
			self.anomdays_arr   = ''
		
		if type(plotconfig['ecv']) != list:
			self.ecvs			= [plotconfig['ecv']]
		else:
			self.ecvs			= plotconfig['ecv']
			
		if type(plotconfig['ecv-set']) != list:
			self.ecv_sets		= [plotconfig['ecv-set']]
		else:
			self.ecv_sets		= plotconfig['ecv-set']
			
		if type(plotconfig['Family']) != list:
			self.families		= [plotconfig['Family']]
		else:
			self.families		= plotconfig['Family']
		
		# -- ECV Defined Arrays -- #
		self.ecv_reading_list	  = [] # Array of reading objects
		self.data_bounds		  = [[] for int in range(len(self.ecvs))] # Boundaries of data from lat/lon restrictions
		self.graph_data_bounds	  = [[] for int in range(len(self.ecvs))] # Lat/Lon Boundaries
		self.units				  = [None for int in range(len(self.ecvs))]
		
		self.graph_data			  = [[] for int in range(len(self.ecvs))] # Data for graph
		self.graph_data_count	  = [[] for int in range(len(self.ecvs))]
		self.graph_data_dates	  = [[] for int in range(len(self.ecvs))] # Dates applied to graph
		self.graph_data_errs	  = [[] for int in range(len(self.ecvs))]
		self.data_calculated	  = [[] for int in range(len(self.ecvs))] # Post-calculations data
		self.data_final			  = [[] for int in range(len(self.ecvs))] # Final Data set to be plotted
		
		self.running_max		  = [-99999 for int in range(len(self.ecvs))]
		self.running_min		  = [99999 for int in range(len(self.ecvs))]
		
		self.scale				  = [None for int in range(len(self.ecvs))]
		self.pre_units			  = [None for int in range(len(self.ecvs))]
		self.scale_off			  = [None for int in range(len(self.ecvs))]
		
		self.isvaluebound		  = [False for int in range(len(self.ecvs))]
		self.value_bounds		  = [[] for int in range(len(self.ecvs))]
			  
		try:
			self.extension = plotconfig['Extension']
		except:
			self.extension = ''
			
		try:
			self.fileexport = plotconfig['FileExport']
		except:
			self.fileexport = ''
			
		try:
			if type(plotconfig['Vinfo']) != list:
				self.vinfos = [plotconfig['Vinfo']]
			else:
				self.vinfos = plotconfig['Vinfo']
		except:
			self.vinfos = ['' for i in range(len(self.ecvs))]
			
		try: # Number of bins
			self.nbins = int(plotconfig['NBins'])
			if self.nbins == 0:
				self.togglebins = False
			else:
				self.togglebins = True
		except:
			self.nbins = 0
			self.togglebins = False
			
		try: # Graph Format
			self.graphformat = plotconfig['GraphFormat']
			if self.graphformat in ['','None','none']:
				self.graphformat = 'Straight'
		except:
			self.graphformat = 'Straight'
			
		# -- ECV User Arrays -- #
		# Special Case Keys
		try:
			if type(plotconfig['SpecialCaseKey']) != list:
				self.spcs = [plotconfig['SpecialCaseKey']]
			else:
				self.spcs = plotconfig['SpecialCaseKey']
			for spc in self.spcs:
				if spc in ['','None','none']:
					spc = None
		except:
			print('Warning: Special case keys missing')
			self.spcs = ['' for int in range(len(self.ecvs))]
			
		# Daynights   
		try:   
			if type(plotconfig['DayNight']) != list:
				self.daynights = [plotconfig['DayNight']]
			else:
				self.daynights = plotconfig['DayNight']
		except:
			print('Warning: Daynights missing')
			self.daynights = ['day' for int in range(len(self.ecvs))]
			
		# Landseas	
		try:   
			if type(plotconfig['LandSea']) != list:
				self.landseas = [plotconfig['LandSea']]
			else:
				self.landseas = plotconfig['LandSea']
		except:
			print('Warning: Landseas missing')
			self.landseas = ['land+sea' for int in range(len(self.ecvs))]
		
		# Instruments
		try:   
			if type(plotconfig['Instrument']) != list:
				self.instruments = [plotconfig['Instrument']]
			else:
				self.instruments = plotconfig['Instrument']
		except:
			print('Warning: Instruments missing')
			self.instruments = ['' for int in range(len(self.ecvs))]
			 
		try:   
			if type(plotconfig['Label']) != list:
				self.labels = [plotconfig['Label']]
			else:
				self.labels = plotconfig['Label']
		except:
			print('Warning: Labels missing')
			self.labels = ['' for int in range(len(self.ecvs))]	 
			 
			 
		try: # Graph Errors
			if plotconfig['GraphErrors'] == 'y':
				self.iserrorgraph	= True
			else:
				self.iserrorgraph	= False
		except:
			self.iserrorgraph = False
			
		try:
			if type(plotconfig['ScaleFactor']) == list:
				self.scale     = [float(sf) for sf in plotconfig['ScaleFactor']]
				self.pre_units = plotconfig['ScaleCoeff']
				self.scale_off = [float(so) for so in plotconfig['ScaleOffset']]
			else:
				self.scale     = [float(plotconfig['ScaleFactor'])]
				self.pre_units = [plotconfig['ScaleCoeff']]
				self.scale_off = [float(plotconfig['ScaleOffset'])]
		except:
			print('Warning: Scales, offsets and units missing')
			self.scale	   = [1 for int in range(len(self.ecvs))]
			self.pre_units = ['' for int in range(len(self.ecvs))]
			self.scale_off = [0 for int in range(len(self.ecvs))]

		# -- Bools -- #			
		self.mapgraph      = plotconfig['MapGraph']
		
		self.ismap		   = [False for int in range(len(self.ecvs))]
		self.ismapbound    = [False for int in range(len(self.ecvs))]
		self.isgraph	   = [False for int in range(len(self.ecvs))]
		self.isgraphbound  = [False for int in range(len(self.ecvs))]
		self.map_bounds    = [[] for index in range(len(self.ecvs))]
		self.user_map_bounds   = [[] for index in range(len(self.ecvs))]
		self.graph_bounds      = [[] for index in range(len(self.ecvs))]
		self.user_graph_bounds = [[] for index in range(len(self.ecvs))]
		
		## Set map and graph bounds ##
		
		# Assemble graph bounds
		if self.mapgraph in ['map+graph', 'graph']:
			self.isgraph = [True for int in range(len(self.ecvs))]
			bounds = plotconfig['GraphBounds']
			bounds_nest = []
			bds = []
			for idx, bound in enumerate(bounds):
				bds.append(float(bound))
				if (idx+1)%4 == 0:
					bounds_nest.append(bds)
					bds = []
			for index in range(len(self.ecvs)):
				if plotconfig['GraphCon'][index] == 'y':
					self.isgraphbound[index] = True
					self.user_graph_bounds[index] = bounds_nest[index]
				else:
					self.isgraphbound[index]   = False
					self.user_graph_bounds[index] = []
		
		# Assemble map bounds
		if self.mapgraph in ['map+graph','map']:
			self.ismap = [True for int in range(len(self.ecvs))]
			bounds = plotconfig['MapBounds']
			bounds_nest = []
			bds = []
			for idx, bound in enumerate(bounds):
				bds.append(float(bound))
				if (idx+1)%4 == 0:
					bounds_nest.append(bds)
					bds = []
			for index in range(len(self.ecvs)):
				if plotconfig['MapCon'][index] == 'y':
					self.ismapbound[index] = True
					self.user_map_bounds[index] = bounds_nest[index]
				else:
					self.ismapbound[index]   = False
					self.user_map_bounds[index] = []
		try:	
			if type(plotconfig['ValueCon']) != list:
				if plotconfig['ValueCon'] == 'y':
					self.isvaluebound = [True]
					self.value_bounds = [[float(bound) for bound in plotconfig['ValueBounds']]]
				else:
					self.isvaluebound = [False]
					self.value_bounds = [[0,0,0,0]]
			else:
				for ivbs in range(len(plotconfig['ValueCon'])):
					if plotconfig['ValueCon'][ivbs] == 'y':
						self.isvaluebound[ivbs] = True
						self.value_bounds[ivbs] = [float(plotconfig['ValueBounds'][ivbs*4 + i]) for i in range(4)]
					else:
						self.isvaluebound[ivbs] = False
						self.value_bounds[ivbs] = [0,0,0,0]
		except:
			self.isvaluebound = [False for int in range(len(self.ecvs))]
			self.value_bounds = [[0,0,0,0] for int in range(len(self.ecvs))]
			
		if plotconfig['ShowOutput'] == 'y':
			self.isshow	   = True
		else:
			self.isshow	   = False
			
		if plotconfig['PNGExport'] == 'y':
			self.ispngexport  = True
			try:
				self.saveto = plotconfig['SaveTo']
			except:
				self.saveto = ''
		else:
			self.ispngexport  = False
			self.saveto = ''
			
		self.doparts = plotconfig['DoParts']
		
		try:
			self.graphname = plotconfig['GraphName']
		except:
			self.graphname = ''
			
		self.graph_type = plotconfig['GraphType']
			
		#self.ecv_map_example = [None for int in range(len(self.ecvs))]

		if self.isshow or self.ispngexport:
			self.main()
		else:
			print('Output method is not defined, please specify png export or show output')
# Main function for analysis
	def main(self):
		# Defines output layout - config file later
		# Flow of processes:
		# - get ecvs (data)
		# - perform calculations
		# - output results
		
		# May later be deconstructed so not all files need to be loaded all together (rearrangement)
		
		print('Starting Main')
		
		if 'Calculate' in self.doparts:
			# Retrieve ecv readings
			is_abort = self.ecv_retrieval()
			if is_abort:
				print('Error: No valid ecv readings - no files selected could be read')
				return None
			# Organise readings into groups
			if self.calculate != '':
				self.create_groups()
			# Set map boundaries
			self.get_units()		
			
			self.ecv_examples = [self.find_valid_ecv_example(index) for index in range(len(self.ecvs))]
			
			# Perform calculations on ecv array - groups are irrelevant for calculations
			for index in range(len(self.ecvs)):
				self.calculations(index)
		if 'ImportGraph' in self.doparts:
			self.import_graph()
			self.get_units()
		if 'ExportGraph' in self.doparts:
			self.export_graph()
		if 'ExportMap' in self.doparts:
			self.export_map_as_netcdf()
		if 'ImportMap' in self.doparts:
			self.import_map_as_netcdf()
			self.get_units()
			for index in range(len(self.ecvs)):
				self.apply_value_bounds(index)
				self.map_bounds[index] = self.user_map_bounds[index]
				self.data_bounds[index] = [0, len(self.data_final[index]), 0, len(self.data_final[index][0])]
		if 'Plot' in self.doparts:
			self.assemble_outputs()
		return None
		
# Assembly of output figures				
	def assemble_outputs(self):
		# Single Fig (4x6)
		# Vertical Fig (8x6)
		# Horizontal Fig (4x12)
		# Triple Fig (8x12)
		# Quadruple Fig (8x12)		
		fig_options = ['Single','Vertical','Horizontal','Triple','Quadruple']
		double_fig_height = ['Vertical','Triple','Quadruple']
		double_fig_width = ['Horizontal','Triple','Quadruple']
		
		
		figure_index = None
		for index, option in enumerate(fig_options):
			if option == self.output_figures:
				figure_index = index
		if figure_index == None:
			print('Error: Invalid Figure Option in input -',self.output_figures)
			return None
		
		if self.output_figures in double_fig_height:
			figheight = LENFIG*2
		else:
			figheight = LENFIG
		if self.output_figures in double_fig_width:
			figwidth = WIDFIG*2
		else:
			figwidth = WIDFIG

		self.out_fig = plt.figure(figsize=(figwidth,figheight))
		self.out_fig.set_facecolor('white')
		self.out_fig.canvas.set_window_title('ECV L3 Map/Graph Tool')
		
		# Maps and Graphs
		if 'map' in self.mapgraph:
			self.map_control(figure_index)
		if 'graph' in self.mapgraph:
			self.graph_control(figure_index)
			
			
		# adding logos routine
		# self.show_logos(figure_index)
		
		if self.ispngexport:
			self.png_export()
		if self.isshow:
			plt.show()
		plt.close()
		return None
		
# Control of map plots
	def map_control(self, figure_index):
		two_maps = ['Horizontal','Triple','Quadruple']
		if self.output_figures in two_maps:
			number_of_maps = 2
		else:
			number_of_maps = 1
			
		for nm in range(number_of_maps):
			
			if self.calculate == 'Anomaly':
				unit_label = '{} Anomaly ({})'.format(self.ecvs[nm], '%')
			elif self.calculate == 'Trend':
				unit_label = '{} Trend ({})'.format(self.ecvs[nm], self.units[nm] + ' y-1')
			else:
				unit_label = '{} Mean ({})'.format(self.ecvs[nm], self.units[nm])
			
			res = (self.map_bounds[nm][1]-self.map_bounds[nm][0])/len(self.data_final[nm])
			
			self.out_fig = od.output_map(self.out_fig, self.data_final[nm],
										 self.map_bounds[nm], res,
										 self.data_bounds[nm], self.value_bounds[nm],
										 map_axes = MAP_FORMATS[figure_index][nm],
										 cbar_axes = CBAR_FORMATS[figure_index][nm],
										 map_title = self.map_title(nm),
										 cbar_title = unit_label,
										 color = self.color,
										 alpha = self.alpha,
										 landsea = self.landseas[nm])
		return None

# Control of graph plots
	def graph_control(self, figure_index):
		# Control Graph Inputs/Outputs
		# Can be expanded using graph_type variable later
		
		# Determines how many graphs are required
		two_graphs = ['Horizontal','Quadruple']
		if self.output_figures in two_graphs:
			number_of_graphs = 2
		else:
			number_of_graphs = 1
			
		# Plot all data on a single time series
		# Plot separate data on different time series
		# Plot two variable comparison
		# Plot three variable comparison
		
		# Graphtypes: SepTime, Time, Variable, 3D
		if 'Time' in self.graph_type:
			graph_titles = []
			if 'SepTime' in self.graph_type:
				for index in range(len(self.ecvs)):
					graph_titles.append( self.graph_title([index]))
			else:
				graph_titles.append(self.graph_title([i for i in range(len(self.ecvs))]))
			
			for index, graph_title in enumerate(graph_titles):
				time_array = self.graph_data_dates[index] # For Completeness
				ecv_label = self.ecvs[index] + self.labels[index] + ' ({})'.format(self.units[index])
				
				if len(graph_titles) == 1:
					graph_data = self.graph_data
					graph_errs = self.graph_data_errs
					axes = GRAPH_FORMATS[figure_index][0]
				else: # Separated graphs
					graph_data = [self.graph_data[index]]
					graph_errs = [self.graph_data_errs[index]]
					axes = GRAPH_FORMATS[figure_index][index]
				self.out_fig = od.output_time_series(self.out_fig,
				                         time_array, graph_data, graph_errs,
				                         graph_title=graph_title,
				                         graph_axes=axes, 
				                         graphformat=self.graphformat,
				                         do_errs=self.iserrorgraph, 
				                         ylabel=ecv_label)
				                         
		elif self.graph_type == 'Variable':
			graph_title = self.graph_title([0,1])
			xlabel = '{} ({})'.format(self.ecvs[0], self.units[0])
			ylabel = '{} ({})'.format(self.ecvs[1], self.units[1])
			self.out_fig = od.output_two_var_graph(self.out_fig,
			                         np.array(self.graph_data), np.array(self.graph_data_errs),
			                         graph_title=graph_title,
			                         graph_axes=GRAPH_FORMATS[figure_index][0],
			                         do_errs=self.iserrorgraph, 
			                         do_bins=self.togglebins,
			                         nbins=self.nbins,
			                         xlabel=xlabel,
			                         ylabel=ylabel)
		else:
			print('Error: Unrecognised format - {}'.format(self.graph_type))

	def export_graph(self):
		# Need to export:
		# self.graph_data[ecv_type]
		# self.graph_data_dates
		# self.graph_data_errs
		gdt = 'g'
		# Dates		nh3values-unit	tpwvalues-unit	nh3errs		tpwerrs
		for index in range(len(self.ecvs)):  
			ecv_type = self.ecvs[index]
			label	= self.labels[index]
			daynight = self.daynights[index]
			family   = self.families[index]
			landsea  = self.landseas[index]
			yrs_mths = str(self.year_arr[0]) + str(self.year_arr[1]) + '_' + str(self.month_arr[0]) + str(self.month_arr[1])
			
			base_path = '/home/users/dwest77/Documents/ECV_Images/output_files/graph_texts/{}'.format(ecv_type+label)
			if not os.path.isdir(base_path):
				os.makedirs(base_path)
			output_file = '/{}_{}_{}_{}_{}{}.txt'.format(family, ecv_type+label, daynight, landsea, yrs_mths, gdt)
			
			os.system('touch {}{}'.format(base_path, output_file))
			f = open(base_path+output_file, 'w')
			outstring = 'Dates, {}-{}, {}_errs\n'.format(ecv_type,self.units[index], ecv_type)
			for index2 in range(len(self.graph_data[index])):
				outstring += str(self.graph_data_dates[index][index2]) + ', '
				outstring += str(self.graph_data[index][index2]) + ', '
				outstring += str(self.graph_data_errs[index][index2]) + '\n'
			f.write(outstring)
			f.close()
		
	def import_graph(self):
		gdt = 'g'
		for index in range(len(self.ecvs)):  
			ecv_type = self.ecvs[index]
			label	= self.labels[index]
			daynight = self.daynights[index]
			family   = self.families[index]
			landsea  = self.landseas[index]
			yrs_mths = str(self.year_arr[0]) + str(self.year_arr[1]) + '_' + str(self.month_arr[0]) + str(self.month_arr[1])
		
			base_path = '/home/users/dwest77/Documents/ECV_Images/output_files/graph_texts/{}'.format(ecv_type+label)
			output_file = '/{}_{}_{}_{}_{}{}.txt'.format(family, ecv_type+label, daynight, landsea, yrs_mths, gdt)
		
			f = open(base_path+output_file, 'r')
			outstring = f.readlines()
			header = outstring[0].replace(' ','').split(',')
			self.units[index] = header[1].split('-')[1]
			for index2 in range(1,len(outstring)):
				line = outstring[index2].replace('\n','')
				line = line.replace(' ','').split(',')
			
				self.graph_data_dates[index].append(line[0])
				self.graph_data[index].append(np.float32(line[1]))
				self.graph_data_errs[index].append(np.float32(line[2]))
		
			f.close()	
	
# Titles and file names of outputs
	def output_filename(self):
		if self.time_format == 'Single':
			time_set = str(self.year_arr[0]) + format(self.month_arr[0],'02d')
			if self.time_increment == 'Daily':
				time_set += format(self.day_arr[0],'02d')
		else:
			time_set = str(self.year_arr[0]) +'-'+ str(self.year_arr[1])
			time_set += '_' + format(self.month_arr[0],'02d') +'-'+ format(self.month_arr[1],'02d')
			if self.time_increment == 'Daily':
				time_set += '_' + format(self.day_arr[0],'02d') +'-'+ format(self.day_arr[1],'02d')
			
		
		ecv_label = ''
		daynights = ''
		landseas = ''
		for index in range(len(self.ecvs)):
			ecv_label += self.families[index] + '-' + self.ecvs[index] + self.labels[index]
			if daynights != self.daynights[index]:
				daynights += self.daynights[index]
				
			if landseas != self.landseas[index]:
				landseas += self.landseas[index]
				
		mapgraph = self.mapgraph
		extension = self.extension
		filename = ecv_label + '_' + daynights + '_' + landseas + '_' + time_set + '_' + mapgraph + '_' + extension
		return filename
		
	def png_export(self):
		# Save figure as png in home directory
		CLOBBER = True
		family, ecv_label, ecv_desc = '', '', ''
		for index in range(len(self.ecvs)):
			family += self.families[index]
			ecv_label += self.ecvs[index]+self.labels[index]
			ecv_desc += '{}_{}{}{}_{}'.format(self.families[index], self.ecvs[index],self.labels[index], self.daynights[index], self.landseas[index])

		main_path = self.saveto + '/{}/{}/'.format(family.upper(), ecv_label.upper())
		if not(os.access(main_path, os.F_OK)):
			os.makedirs(main_path)
		main_path = main_path + self.output_filename() + '.png'
		self.out_fig.savefig(main_path)
			
	def graph_title(self, indexes):
		variable_title = ''
		for i, index in enumerate(indexes):
			var = self.ecvs[index]
			if self.labels[index] != '':
				var += '-{}'.format(self.labels[index])
			var += ' ({}) '.format(self.daynights[index])
			if i != len(indexes)-1:
				var += 'vs '
			variable_title += var
		
		minit = self.month_arr[0]
		mfin  = self.month_arr[1]	
		month_init = CAL_MONTHS[minit]
		month_fin  = CAL_MONTHS[mfin]
		
		year_init  = self.year_arr[0]
		year_fin   = self.year_arr[1]
			
		if self.time_format == 'Single':
			back = '{}-{}'.format( month_init, year_init )
		elif self.time_format == 'Range':
			back = '{}-{} to {}-{}'.format(month_init, year_init,
										   month_fin , year_fin)
		elif self.time_format == 'Months' or self.time_format == 'Days':
			back = ''
			if self.time_format == 'Days':
				if day_init != day_fin:
					back += '{}-{} '.format(day_init, day_fin)
				else:
					back += '{} '.format(day_init)

			if month_init != month_fin:
				back += '{}-{} '.format(month_init, month_fin)
			else:
				back += '{} '.format(month_init)
		
			if year_init != year_fin:
				back += '{}-{} '.format(year_init, year_fin)
			else:
				back += '{} '.format(year_init)
		variable_title += ' ' + back
		return variable_title

	def map_title(self, index):
		# Create map title given all parameters unique to this plot
		# Depends on specified times, ecv set etc.
		minit = self.month_arr[0]
		month_init = CAL_MONTHS[minit]
		year_init  = self.year_arr[0]
		day_init = ''
		if self.time_increment == 'Daily':
			day_init   = self.day_arr[0]
		
		if self.time_format != 'Single':
			
			mfin  = self.month_arr[1]
			
			if mfin > 12:
				mfin -= 12
				
			month_fin  = CAL_MONTHS[mfin]		
			year_fin   = self.year_arr[1]
			day_fin = ''
			if self.time_increment == 'Daily':
				day_fin   = self.day_arr[1]
		
		if self.calculate == 'Anomaly':
			anom_month_init = CAL_MONTHS[self.anommonth_arr[0]]
			anom_month_fin  = CAL_MONTHS[self.anommonth_arr[1]]
			
			anom_year_init  = self.anomyear_arr[0]
		start = ''
		start += self.ecvs[index]+self.labels[index]+' ({})'.format(self.daynights[index])
		start += ' '
		
		## Calculation type - Trend, Anomaly or Mean (Average and Single)
		if self.calculate == 'Trend':
			front = start + ' Trend for '
		elif self.calculate == 'Anomaly':
			front = start + ' Anom '
			if anom_month_init == anom_month_fin:
				front += anom_month_init
			else:
				front += anom_month_init+'-'+anom_month_fin
			front += ' '
			front += str(anom_year_init)
			front += ' from '
		else:
			front = start + ' Mean '
			
		## Dates used in calculations
		if self.time_format == 'Single':
			back = '{}-{}'.format( month_init, year_init )
		elif self.time_format == 'Range':
			back = '{}-{} to {}-{}'.format(month_init, year_init,
										   month_fin , year_fin)
		elif self.time_format == 'Months' or self.time_format == 'Days':
			back = ''
			if self.time_format == 'Days':
				if day_init != day_fin:
					back += '{}-{} '.format(day_init, day_fin)
				else:
					back += '{} '.format(day_init)

			if month_init != month_fin:
				back += '{}-{} '.format(month_init, month_fin)
			else:
				back += '{} '.format(month_init)
		
			if year_init != year_fin:
				back += '{}-{} '.format(year_init, year_fin)
			else:
				back += '{} '.format(year_init)
				
		elif self.time_format == 'Days':
			back = '{}-{}, {}-{} for {}'.format(year_init , year_fin,
													  month_init, month_fin, day_init)
		else:
			back = ' Other' # Unrecognised time format 
		plot_title = front + back
		print(plot_title)
		
		return plot_title
		
# Quick export and importing of saved data (Export/Import Map/Graph for DoParts variable)
	def export_map_as_netcdf(self):		
		for index in range(len(self.ecvs)):
			ecv_example = self.find_valid_ecv_example(index)
			latlist = ecv_example.latlist[index]
			lonlist = ecv_example.lonlist[index]
			
			new_latlist = []
			new_lonlist = []
			for lat in latlist:
				new_latlist.append(lat)
			for lon in lonlist:
				new_lonlist.append(lon)
			
			daynight = self.daynights[index]
			landsea  = self.landseas[index]
			family   = self.families[index]
			variable = self.ecvs[index] + self.labels[index]


			base_path = self.fileexport
			if not os.path.isdir(base_path):
				os.makedirs(base_path)

			output_file = self.output_filename() + '.nc'
			base_path = base_path + '/{}/{}/'.format(family.upper(), ecv_label.upper()) + output_file
			
			ncfile = Dataset(base_path, 'w', format = 'NETCDF4')
			lat_dim = ncfile.createDimension('lat',len(latlist))
			lon_dim = ncfile.createDimension('lon',len(lonlist))
			lat_var = ncfile.createVariable('lat', np.float32, ('lat',))
			lat_var.long_name = 'latitude'
			lat_var[:] = new_latlist
			lon_var = ncfile.createVariable('lon', np.float32, ('lon',))
			lon_var.long_name = 'longitude'
			lon_var[:] = new_lonlist
		
			data_var = ncfile.createVariable(variable, np.float32, ('lat','lon',))
			data_var.long_name = variable
			data_var.units = self.units[index]
			
			data_out = np.array(self.data_final[index])
			
			data_var[:,:] = data_out
		
			ncfile.close()
			
			print('Map Exported to', base_path, output_file)
			
	def import_map_as_netcdf(self):		
		for index in range(len(self.ecvs)):
			base_path = self.fileexport
			
			variable = self.ecvs[index] + self.labels[index]	
			input_file = self.output_filename() + '.nc'
			base_path = base_path + '/{}/{}/'.format(self.families[index].upper(), variable.upper()) + input_file
			ncfile = Dataset(base_path, 'r', format = 'NETCDF4')
			self.data_calculated[index] = np.array(ncfile[variable])
		
			ncfile.close()
			
# ECV retrieval functions - organising file retrievals					
	def get_readings(self, year, month, day):
		if VERBOSE:
			print('Retrieving ecv for {}/{}/{}'.format(day, month, year))
		# Create each ecv reading given a date (for all specified data sets)
		# Adds to graph dates so points can be added cumulatively later
		  
		isanom = self.check_anomaly(year, month, day)
		# Extend to use multiple variables at later date
		new_read = ecv_reading(year, month, day, self.time_increment, isanom)
		new_read.set_ecv_properties(self.ecvs, self.ecv_sets, self.families, self.instruments, self.labels, self.spcs, self.pre_units, self.vinfos)
		new_read.set_map_properties(self.daynights, self.landseas)
		new_read.sort_file_imports()
		addto = False
		for index in range(len(self.ecvs)):
			if new_read.is_valid[index] and len(new_read.basedata[index]) > 0:
				addto = True
		if addto:
			self.ecv_reading_list.append(new_read) # Add to list of ecv readings if its valid
			
	def ecv_retrieval(self):
		print('Starting ECV retrieval')
		# Function for managing the retrieval of all ecvs
		# Ecv retrieval is based on TimeFormat property (only gets the required files, not all files)
		# Ecv_reading_list contains all ecv instances
		# Grouping later - average groups of ecvs in ecv_reading list
		
		day_start, day_stop = 1, 1
		
		# Define year, month, day start and stop conditions
		year_current = self.year_arr[0]
		if len(self.year_arr) == 2:
			year_end = self.year_arr[1]
		
		month_start = self.month_arr[0]  
		if len(self.month_arr) == 2:
			month_stop = self.month_arr[1]
			
		if self.time_format == 'Days' or self.time_increment == 'Daily':
			day_start = self.day_arr[0]
			if len(self.day_arr) == 2:
				day_stop = self.day_arr[1]
				
		# Different time formats have different methods of retrieving all relevant files
		is_end = False
		if self.time_format == 'Single': # No complicated calculations
			print('Retrieving ECVS (1)')
			self.get_readings(year_current,month_start,day_start)

		elif self.time_format == 'Range': # Get readings for every date between two dates
			if self.time_increment == 'Monthly':
				num_ecvs = 12 - month_start + (year_end-year_current-1)*12 + month_end
				print('Retrieving ECVS ({})'.format(num_ecvs))
				
				self.find_all_months(month_start, month_stop, year_current, year_end, day_start)
			elif self.time_increment == 'Daily':
				num_ecvs = pm.date_to_value(day_stop, month_stop, year_end) - pm.date_to_value(day_start, month_start, year_current)
				print('Retrieving ECVS ({})'.format(num_ecvs))
				
				self.find_all_days(day_start, day_stop, month_start, month_stop, year_current, year_end)
			else:
				print('Error: {} not an accepted time format'.format(self.time_format))
				
		elif self.time_format == 'Months': # Get a period of months for each year
			num_ecvs = (month_stop - month_start)*(year_end - year_current)
			print('Retrieving ECVS ({})'.format(num_ecvs))
			
			while is_end == False:
				for month in range(month_start,month_stop+1): # Run through the period of months for each year
					if month > 12:
						self.get_readings(year_current+1, month-12, day_start)
					else:
						self.get_readings(year_current, month, day_start)
				# Stop Condition
				if year_current == year_end:
					is_end = True
				year_current += 1 # Run through the years with the same month each time
				
		elif self.time_format == 'Days': # Get a period of days in a period of months for each year
			num_ecvs = (day_stop-day_start)*(month_stop - month_start)*(year_end - year_current)
			print('Retrieving ECVS ({})'.format(num_ecvs))

			while is_end == False:
				for month in range(month_start,month_stop+1): # Run through the period of months and days for each year
					for day in range(day_start, day_stop+1):
						self.get_readings(year_current, month, day)
				# Stop Condition	
				if year_current == year_end:
				   is_end = True
				year_current += 1	  
		is_abort = True		
		for ecv in self.ecv_reading_list:			
			for index in range(len(self.ecvs)):
				if ecv.is_valid[index]:
					is_abort = False
		return is_abort
		
	def find_all_months(self, month_counter, month_stop, year_counter, year_stop, day):
		# Find all months in a range of years and retrieve ecvs for each date
		# Call get_readings with each specific month/year combination
		is_end = False
		while is_end == False:
			# End of Year Condition - Add to year_counter
			if month_counter == 13:
				month_counter = 1
				year_counter += 1
			# Retrieve single reading
			self.get_readings(year_counter, month_counter, day)
			# Stop Condition
			if year_counter == year_stop and month_counter == month_stop:
				is_end = True
			month_counter += 1 # Run through the whole range of dates between the specified values

	def find_all_days(self, day_counter, day_stop, month_counter, month_stop, year_counter, year_stop):
		# Find all days, months in a range of years and retrieve ecvs for each date
		# Call get_readings with each specific day/month/year combination
		is_end = False
		while is_end == False:
			# End of Month Condition - Add to month_counter
			if day_counter == MONTH_CONTENT[month_counter-1]+1:
				day_counter = 1
				month_counter += 1
			# End of Year Condition - Add to year_counter
			if month_counter == 13:
				month_counter = 1
				year_counter += 1
			# Retrieve single reading
			self.get_readings(year_counter, month_counter, day_counter)
			# Stop Condition
			if year_counter == year_stop and month_counter == month_stop and day_counter == day_stop:
				is_end = True
			day_counter += 1

	def create_groups(self):
		# New function for arranging ecvs into groups by user defined time range not simply time group size
		# Create arrays of grouped dates by time group size
		# Place ecvs into correct groups by matching grouped dates to the ecvs
			   
		# Define values for ease of use
		month_init = self.month_arr[0]
		month_fin  = self.month_arr[1]
		
		year_init  = self.year_arr[0]
		year_fin   = self.year_arr[1]
		
		if self.time_increment == 'Daily':
			day_init = self.day_arr[0]
			day_fin	 = self.day_arr[1]
		
		dates_array = []
		date_group = []
		counter = 0
		for yr in range(year_init, year_fin+1):
			for mth in range(month_init, month_fin+1):
				# Cycle months and years given (works for range or period)
				if mth > 12: # For month=14 i.e. Dec-feb = 12-14 months in same year
					month = mth - 12
					year = yr + 1
				else:
					month = mth
					year = yr
				
				# Separate daily loop	
				if self.time_increment == 'Daily':
					for day in range(day_init, day_fin+1):
						if counter >= int(self.timegroupsize):
							# Add assembled date to dates group
							dates_array.append(date_group)
							date_group = []
							counter = 0
						date_group.append([year, month, day])
						counter += 1
				else:
					if counter >= int(self.timegroupsize):
						# Add assembled date to dates group
						dates_array.append(date_group)
						date_group = []
						counter = 0
					date_group.append([year, month])
					counter += 1
		# Last catch
		dates_array.append(date_group)
		# Dates array contains all groups of dates together	   
		ecv_counter = 0
		self.dates_list = []
		for index, group in enumerate(dates_array):
			
			# Set ECV date indexes by matching with dates from dates_array
			for date_set in group:
				if ecv_counter < len(self.ecv_reading_list):
					
					# Match ecv date to the date in each group
					ecv = self.ecv_reading_list[ecv_counter]
					is_ym_match = (ecv.year == date_set[0] and ecv.month == date_set[1])
					if len(date_set) > 2:
						is_d_match = (ecv.day == date_set[2])
					else:
						is_d_match = True
						
					# Date of ecv matches date in group?
					if is_ym_match and is_d_match:
						# ECV matches date displayed
						ecv.date_index = index
						ecv_counter += 1
			
			# Assemble date string and add to dates_list			
			date_string = ''
			if self.time_increment == 'Daily':
				if self.timegroupsize > 1:
					date_string = '{}-'.format(group[0][2])
				date_string += '{}/{}/{}'.format(group[len(group)-1][2] , group[0][1], group[0][0])
												   
				value_in_days = pm.date_to_value(group[0][2], group[0][1], group[0][0])	
			else:
				if self.timegroupsize > 1:
					date_string += '{}-'.format(group[0][1])
				date_string += '{}/{}'.format(group[len(group)-1][1], group[0][0])
				value_in_days = pm.date_to_value(1, group[0][1], group[0][0])	
			self.dates_list.append([date_string, value_in_days])
			
			for index in range(len(self.ecvs)):
				if date_string not in self.graph_data_dates[index]: # New date to add
					self.graph_data_dates[index].append(date_string)
					self.graph_data[index].append(0) # Initial values
					self.graph_data_count[index].append(0)
					if self.iserrorgraph:
						self.graph_data_errs[index].append([])

# Main analysis and calculations - previously required 5 large functions, now condensed to 1 using numpy (v4)
	def calculations(self, index): 
		# Calculations hub to call all necessary functions for calculating specific parameters
		# Small function directs program to other calculating functions
		
		# self.map_bounds - User defined map boundaries
		# self.graph_bounds - User defined graph boundaries
		# map_bounds - actual boundaries to be applied
		# data_bounds - map_bounds acting on an array given that 0,0 in the array is lat_min,lon_min in coordinates
		# graph_data_bounds - actual graph data boundaries to be applied
		
		ecv_type = self.ecvs[index]
		ecv_example = self.ecv_examples[index]
		self.configure_bounds(index, ecv_example)
		print('** Starting Calculations for type: {} **'.format(ecv_type))
		
		# Concatenate all basedata to calculation array
		data_to_calculate = np.array([])
		for ecv in self.ecv_reading_list:
			if ecv.is_valid:
				# Apply sea/land mask here - extra step later
				data = ecv.basedata[index][self.data_bounds[index][0]:self.data_bounds[index][1], self.data_bounds[index][2]:self.data_bounds[index][3]]

				data_to_calculate = np.append(data_to_calculate, np.array(data)*self.scale[index] - self.scale_off[index])
		
		# Reshape array to 3d time/lat/lon array
		wid = self.data_bounds[index][1] - self.data_bounds[index][0]
		hig = self.data_bounds[index][3] - self.data_bounds[index][2]
		data_to_calculate = np.reshape(data_to_calculate, (len(self.ecv_reading_list),wid, hig))
		
		# Apply extremes filter
		data_to_calculate[data_to_calculate > 99999] = np.nan
		data_to_calculate[data_to_calculate < -99999] = np.nan
		
		#* Perform Calculations*#
		if self.calculate == '' or self.calculate == 'None' or self.calculate == 'Average':
			data_to_calculate[data_to_calculate==0] = np.nan
			self.data_calculated[index] = np.nanmean(data_to_calculate, axis=0)
		elif self.calculate == 'SDFM':
			# Standard Deviation Filtered Mean
			means = np.nanmean(data_to_calculate, axis=0)
			sds = np.nanstd(data_to_calculate, axis=0)
			#sds[:,:] = 0.09
			upperlims = means + 3*sds
			upperlims = np.reshape(np.repeat(upperlims, len(data_to_calculate)), np.shape(data_to_calculate))
			lowerlims = means - 3*sds
			lowerlims = np.reshape(np.repeat(lowerlims, len(data_to_calculate)), np.shape(data_to_calculate))
			data_to_calculate[data_to_calculate > upperlims] = np.nan
			data_to_calculate[data_to_calculate < lowerlims] = np.nan
			self.data_calculated[index] = np.nanmean(data_to_calculate, axis=0)	
		
		elif self.calculate == 'Total':
			self.data_calculated[index] = np.nansum(data_to_calculate, axis=0)
		elif self.calculate == 'Trend':
			dates_arr = np.arange(1,len(self.ecv_reading_list)+1)
			if self.time_increment == 'Daily':
				dates_arr = dates_arr/365
			else:
				dates_arr = dates_arr/12
			dates_arr = np.reshape(np.tile(dates_arr, len(data_to_calculate[0])*len(data_to_calculate[0][0])), 
			                       (len(data_to_calculate), len(data_to_calculate[0]), len(data_to_calculate[0][0])))
			self.data_calculated[index] = pm.np_linear_regressions(dates_arr, data_to_calculate)
		elif self.calculate == 'Anomaly':
			anom_data_to_calculate = np.array([])
			# Concatenate all anomaly data
			for ecv in self.ecv_reading_list:
				if ecv.is_anom and ecv.is_valid:
					# Apply sea/land mask here
					anom_data = ecv.basedata[index][self.data_bounds[index][0]:self.data_bounds[index][1], self.data_bounds[index][2]:self.data_bounds[index][3]]

					anom_data_to_calculate = np.append(anom_data_to_calculate, np.array(anom_data)*self.scale[index] - self.scale_off[index])
			
			# Reshape array to 3d time/lat/lon array
			anom_data_to_calculate = np.reshape(anom_data_to_calculate, (len(self.ecv_reading_list),wid, hig))
			avg_all_data = np.nanmean(data_to_calculate, axis=0)
			avg_anom_data = np.nanmean(anom_data_to_calculate, axis=0)
			
			self.data_calculated[index] = (avg_anom_data - avg_all_data)*(100/avg_all_data)
		
		if self.isgraph[index]:
			data_to_graph = np.array([])
			for ecv in self.ecv_reading_list:
				if ecv.is_valid:
					# Apply sea/land mask here
					data = ecv.basedata[index][self.graph_data_bounds[index][0]:self.graph_data_bounds[index][1], self.graph_data_bounds[index][2]:self.graph_data_bounds[index][3]]
				
					data_to_graph = np.append(data_to_graph, np.array(data)*self.scale[index] - self.scale_off[index])
			data_to_graph = np.reshape(data_to_graph, (len(self.ecv_reading_list),wid, hig))
			data_to_graph[data_to_graph > 99999] = np.nan
			data_to_graph[data_to_graph < -99999] = np.nan
			
			self.graph_data[index] = np.nanmean(data_to_graph, axis=(1,2))
			self.graph_data_errs[index] = np.nanstd(data_to_graph, axis=(1,2))/ (len(data_to_graph[0])*len(data_to_graph[0][0]))
		
		self.running_max[index] = np.nanmax(data_to_calculate)
		self.running_min[index] = np.nanmin(data_to_calculate)
		self.apply_value_bounds(index)
	
# General utility functions
	def find_valid_ecv_example(self, index):
		is_found = False
		counter = 0
		while is_found == False and counter < len(self.ecv_reading_list):
			ecv_example = self.ecv_reading_list[counter]
			if ecv_example.is_valid[index] == True:
				is_found = True
			counter += 1
		if not is_found:
			return None
		else:
			return ecv_example

	def get_units(self):
		for index in range(len(self.ecvs)):
			try:
				ecv_example = self.find_valid_ecv_example(index)
				units = ecv_example.baseunits[index]
			except:
				units = ''
			
			if self.pre_units[index] != None:
				if '_' in self.pre_units[index]:
					# Replace units (after conversion)
					self.units[index] = self.pre_units[index].replace('_','')
				else:
					# Addition to existing units (e.g yr-1)
					self.units[index] = units + self.pre_units[index]
			else:
				self.units[index] = units

	def show_logos(self, figure_index):
		## RAL and Cloud CCI Logos
		ral_ax = RAL_LOGOS[figure_index][0]
		cci_ax = CCI_LOGOS[figure_index][0]
		if self.use_ral:
			ral_ax = self.out_fig.add_axes(ral_ax) # RAL Logo
			# Remove axes (no x-y for images)
			ral_ax.set_axis_off()
			try:
				img = mpimg.imread('/gws/nopw/j04/cds_c3s_cloud/dwest77/ECV/archive/RAL_Space.jpg')
				imgplot = plt.imshow(img)
				imgplot.set_clim(0.0,7.0)
			except:
				print('ImageError: RAL Space Logo not found')
		if self.use_cci:
			cci_ax = self.out_fig.add_axes(cci_ax) # Cloud Logo
			cci_ax.set_axis_off()
			try:
				img = mpimg.imread('/gws/nopw/j04/cds_c3s_cloud/dwest77/ECV/archive/cci_cloud.png')
				imgplot = plt.imshow(img)
				imgplot.set_clim(0.0,7.0)
			except:
				print('ImageError: Cloud CCI Logo not found')

	def check_anomaly(self, year, month, day):
		# Check if this date is inside anomaly range given as input
		# Return true or false for setting isanom value for ecv
		
		if self.calculate != 'Anomaly':
			return False
		else:
			# Anomaly conditions
			try:
				is_day = ((self.anomday_arr[0] <= day and self.anomday_arr[1] <= day) or (self.time_increment == 'Monthly'))
			except:
				is_day = True
			is_month = (self.anommonth_arr[0] <= month and self.anommonth_arr[1] <= month)
			is_year = (self.anomyear_arr[0] <= year and self.anomyear_arr[1] <= year)
			return (is_day and is_month and is_year) # Boolean AND sum of all three

	def get_data_constraints(self, reset_latlons, bounds, ecv_index): 
		# Conversion between latitude/longitude and data coordinates
		# Takes the ecv reading list to convert the internal lat/lonlists of each ecv object
		ecv_example = self.find_valid_ecv_example(ecv_index)
		# Upgrade to find first valid example with correct type
	
		row_bound = len(ecv_example.basedata[ecv_index])	   # Rows of data points - size
		column_bound = len(ecv_example.basedata[ecv_index][0]) # Columns of data points - size

		offset_lat = ecv_example.lat_min[ecv_index]	 # for transforming origin to bottom corner
		offset_lon = ecv_example.lon_min[ecv_index]

		location_constraints = np.float32(bounds) # Convert to np format
		# Adjust constraints set by user to origin at bottom corner
		adjusted_constraints = [0,0,0,0] 
		for coord in range(2):
			adjusted_constraints[coord] = location_constraints[coord] - offset_lat
			adjusted_constraints[coord + 2] = location_constraints[coord + 2] - offset_lon

		# Adjust map boundaries by origin as well
		abs_lat = np.float32(ecv_example.lat_max[ecv_index] - ecv_example.lat_min[ecv_index])
		abs_lon = np.float32(ecv_example.lon_max[ecv_index] - ecv_example.lon_min[ecv_index])


		# Constrain Lat/Lon values for each ecv by user defined location
		# Only do when considering map constraints, not graph constraints (so not done twice)
		if reset_latlons:
			for index, ecv_read in enumerate(self.ecv_reading_list):
				new_lats = []
				new_lons = []
				for lat in ecv_read.latlist[ecv_index]: # Transform lat and lon list for the new constraints
					if not (lat < location_constraints[0] or lat > location_constraints[1]):
						new_lats.append(lat)
				for lon in ecv_read.lonlist[ecv_index]:
					if not (lon < location_constraints[2] or lon > location_constraints[3]):
						new_lons.append(lon)
				self.ecv_reading_list[index].latlist[ecv_index] = new_lats
				self.ecv_reading_list[index].lonlist[ecv_index] = new_lons # Reset lat and lon list for new constraints
				
		data_bounds = [0,0,0,0]
		for coord in range(2):
			data_bounds[coord] = int(round((adjusted_constraints[coord]/abs_lat)*row_bound))
			data_bounds[coord + 2] = int(round((adjusted_constraints[coord + 2]/abs_lon)*column_bound))
		# Get data boundaries corresponding to the established lat lon constraints set by the user
		return data_bounds		   

	def apply_value_bounds(self, index):
		# Constrain values inside/outside a range specified by the user
		# Values above or below are set to those values.
		
		if self.isvaluebound[index]:
			print('bounding values')	
			calculate = np.array(self.data_calculated[index])
			max_outer = np.float32(self.value_bounds[index][1])
			min_outer = np.float32(self.value_bounds[index][0])
			max_inner = np.float32(self.value_bounds[index][2])
			min_inner = np.float32(self.value_bounds[index][3])
			
			calculate[calculate > max_outer] = max_outer
			calculate[calculate < min_outer] = min_outer
			
			min_inner_f = calculate > min_inner
			max_inner_f = calculate < max_inner
			calculate[min_inner_f & max_inner_f] = np.float32(0)
			
			self.data_final[index] = calculate
			self.running_min[index] = self.value_bounds[index][0]
			self.running_max[index] = self.value_bounds[index][1]
		else:
			self.data_final[index] = self.data_calculated[index]
			self.value_bounds[index] = [self.running_min[index], self.running_max[index],0,0]

		
	def configure_bounds(self, index, ecv_example):
		# Determine map/graph bounds and data bounds to use for map and graph plotting
		for index in range(len(self.ecvs)):
			if self.ismapbound[index]:
				self.map_bounds[index] = self.user_map_bounds[index]
				self.data_bounds[index] = self.get_data_constraints(True, self.map_bounds[index], index)
			else:
				self.map_bounds[index] = [ecv_example.lat_min[index], ecv_example.lat_max[index], ecv_example.lon_min[index], ecv_example.lon_max[index]]
				self.data_bounds[index] = [0,len(ecv_example.basedata[index]),0,len(ecv_example.basedata[index][0])]
			#
			if self.isgraph[index]:
				if not self.isgraphbound[index]:
					self.graph_data_bounds[index] = self.data_bounds[index] # Graph matches the map boundaries
				else: # Graph uses a specified range smaller than the map boundaries
					self.graph_data_bounds[index] = self.get_data_constraints(not self.ismapbound[index], self.user_graph_bounds[index], index)
					
			# May be required for reinstating is_near_land function from v3 - removes non land/sea data	
			#self.ecv_map_example[index] = Basemap(projection='cyl',llcrnrlat=self.map_bounds[index][0], urcrnrlat=self.map_bounds[index][1], 
									   #llcrnrlon=self.map_bounds[index][2], urcrnrlon=self.map_bounds[index][3], lat_ts=20, resolution='l')   
	  
## --- End ECV analysis class --- ##

## Global function calls
options, operands = getopt(sys.argv[1:], "", ["pcifile="])
# Determine options and processing method
if operands != []:
	pcifile = operands[0].replace('pcifile=','')
	plotconfig = ff.get_from_file('pcis/',pcifile,'')
	status = ecv_analysis(plotconfig)
	if status == 0:
		print('Program exited with code 0')
elif __name__ == "__main__":
	console_intro()

	usedefault = ff.accept_input('Use Default Settings? (y/n): ',accepted=['y','yes','n','no'])
	if usedefault == 'y' or usedefault == 'yes':
		plotconfig = ff.get_from_file('pcis/','default.pci','')
	else:
		configfile = ff.accept_input('Plot Config File: ', atype='file',fileroot='pcis/')
		plotconfig = ff.get_from_file('pcis/',configfile,'')
	
	status = ecv_analysis(plotconfig)
	if status == 0:
		print('Program exited with code 0')
		
else:
	print('Warning: Program not on main and pci file not given')
	
			
			
