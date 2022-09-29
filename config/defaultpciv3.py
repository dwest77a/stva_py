## Python default.py
# Daniel Westwood - daniel.westwood@stfc.ac.uk
#
# Usage:
#   - Default parameters for pci content

    # ----- Time Settings - Specify groupings of ecvs ----- #
    # # TimeFormat     : Single, Range, Months, days      # #
    # # TimeIncrement  : Monthly, Daily                   # #
    # # TimeGroupSize  : Number of points in group        # #

    # Years, AnomalyYear, Months, AnomalyMonths - ints                - Dropdowns
    
    # ecvs - which ECVs are being used                                 - Dropdown
    # ecv-set - list of all sets being used                            
    # erb, srb, cld - groupings of sets into the ecv products          - Dropdowns
    
    # Calculate, TimeFormat - Specify calculations                     - Dropdowns
    # Orders - equation format for a composite ecv (eg. sis - srs)     - Textbox
    # Color - viridis or others from the color dict                    - Dropdowns
    # MapCon, Graph, ValueCon, PNGExport, Composite                    - Tickbox
    # MapBounds, GraphBounds, ValueBounds                              - Tickboxes

    ## ATSR2 and AATSR Radiation Datasets as ecv-set(Jasmin: /gws/nopw/j04/cds_c3s_cloud...)

    ## Surface Radiation Budget Datasets
    
    # - Surface Incoming Shortwave (sis)
    # - Surface Outgoing Shortwave (srs) 'Reflected'
    # - Surface Net Shortwave (sns)
    # - Surface Downwelling Longwave (sdl)
    # - Surface Outgoing Longwave (sol)
    # - Surface Net Longwave (snl)

    ## Earth Radiation Budget Datasets

    # - Cloud Count (cloud_count)
    # - Reflected Solar Flux (rsf)
    # - Outgoing Longwave Radiation (olr)

    ## Cloud Properties

    # - Liquid count month or day (liquid_count_month)
    # - Ice count month or day (ice_count_month)
    # - Cloud Fractional Cover (cfc)
    # - Cloud Top Pressure mean (ctp)
    # - Cloud Top Height mean (cth)
    # - Cloud Top Temperature (ctt)
    # - Cloud Effective Radius mean (cer)
    # - Cloud Optical Thickness (cot)
    # - Liquid Water Path mean (lwp)
    # - Ice Water Path mean (iwp)


default_plotconfig = {
    'TimeFormat': 'Single',
    'TimeIncrement':'Monthly', # Mainly for when using range (monthly or daily)
    'TimeGroupSize':1,
    'Calculate':'None',
    'MapGraph':'map',
    'FigureNumber':1,
    'DoParts':['Calculate','ExportMap'],

    'DayNight':['day','day'],
    'Years':[2016, 2018],
    'Months':[1,2],                          
    'Days':[1,1],        
    'AnomalyYears':[],
    'AnomalyMonths':[],
    'AnomalyDays':[],

    'ecv':['nh3'],
    'ecv-set':['nh3'],
    'Family':['cris'],
    'Instruments':['',''],
    'SpecialCaseKey':['None','None'],

    'MapCon':'n',
    'GraphCon':'n',
    'MapBounds':[49.5,60.5,-10,2],
    'GraphBounds':[],
    'LandSea':['land+sea','land+sea'],
    'ShowOutput':'n',
    'PNGExport':'y',
    'GraphType':'Variable',
    'GraphFormat':'Straight',
    'GraphErrors':'y',
    'NBins':20,
    
    'Color':'viridis',
    'Alpha':0.7,
    'ScaleFactor':[1000,1],
    'ScaleOffset':[0,0],
    'ScaleCoeff':['_ppbv',''],
    'ValueCon':'n',
    'ValueBounds':[],
    'SaveTo':'/home/users/dwest77/Documents/ECV_Images/Output'}
    
