#!/usr/bin/env python3

import os
import sys

from collections import OrderedDict
from datetime    import timedelta


def generate_echo(inpfile, ladas_cpl = 0):
   """
    Echo generator of inpfile, ignore line starts with "## "
   """
   if ladas_cpl == 0 :
      use_rc_defaults = 'GEOSldas=>'    # use defaults for LDAS
   else :
      use_rc_defaults = 'GEOSagcm=>'    # use defaults for AGCM

   with open (inpfile) as fin :
      for line in fin:
         if use_rc_defaults in line:  
            line = line.split(use_rc_defaults)[1]
         yield line

def parseInputFile(inpfile, ladas_cpl=0):
   """ 
     Parse the input file and return a dict of options
     Inpfile : input file
     Output  : dict, if inpfile does not exist, return empty {}
   """ 
   if not os.path.exists(inpfile):
     assert ladas_cpl > 0, " No exeinput file only if ladas_cpl > 0"
     return {}
        
   inpdict = OrderedDict()
   errstr  = "line [%d] of [%s] is not in the form 'key: value'"

   echoLines = generate_echo(inpfile, ladas_cpl=ladas_cpl)
         
   linenum = 0
   for line in echoLines:
     line = line.strip()
     linenum += 1
     if not line:
        continue
     # handle comments
     position = line.find('#')
     if position==0: # comment line
        continue
     if position>0:  # strip out comment
        line = line[:position]
     # we expect a line to be of the form
     # key : value
     assert ':' in line, errstr % (linenum, inpfile)
                                                             
     key, val = line.split(':',1)
     key = key.strip()
     val = val.strip()
     if not key or not val:
       print ("WARNING: " + errstr % (linenum, inpfile))
       continue
     #raise Exception(errstr % (linenum, inpfile))
     if key in inpdict:
       raise Exception('Duplicate key [%s] in [%s]' % (key, inpfile))
     inpdict[key] = val.strip()
   return inpdict

def echoInputFile(inpfile, ladas_cpl = 0):
   """
    Echo inpfile, ignore line starts with "## "
   """
   echoLines = generate_echo(inpfile, ladas_cpl=ladas_cpl)

   for line in echoLines:
      sys.stdout.write(line) 
      sys.stdout.flush() 

def printExeInputSampleFile():
   """
    Print sample exeinp file to screen
   """
   print ('####################################################################################')
   print ('#                                                                                  #')
   print ('#                             REQUIRED INPUTS                                      #')
   print ('#                                                                                  #')
   print ('#   The inputs to ldas_setup can be provided in an "exeinp" file such as this      #')
   print ('#     sample file.                                                                 #')
   print ('#   When ldas_setup is called by fvsetup to set up an experiment with the coupled  #')
   print ('#     land-atm data assimilation system, the inputs to ldas_setup are provided as  #')
   print ('#     optional command line arguments.                                             #')
   print ('#                                                                                  #')
   print ('####################################################################################')
   print ()
   print ('############################################################')
   print ('#                                                          #')
   print ('#                 EXPERIMENT INFO                          #')
   print ('#                                                          #')
   print ('# Format for start/end times is yyyymmdd hhmmss.           #')
   print ('#                                                          #')
   print ('############################################################')
   print ()
   print ('EXP_ID:')
   print ('EXP_DOMAIN:')
   print ('NUM_LDAS_ENSEMBLE:')
   print ('BEG_DATE:')
   print ('END_DATE:')
   print ()
   print ('############################################################')
   print ('#                                                          #')
   print ('#                 RESTART INFO                             #')
   print ('#                                                          #')
   print ('# (i) Select "RESTART" option from [1,2,3,M]:              #')
   print ('#                                                          #')
   print ('# Use one of the following options to restart from a       #')
   print ('#   GEOSldas experiment:                                   #')
   print ('#                                                          #')
   print ('# RESTART: 1 = No remapping; zoom optional.                #')
   print ('#   ---------------------------------------                #')
   print ('#   New experiment has SAME tile space (grid), boundary    #')
   print ('#   conditions, and snow model parameter (WEMIN) as        #')
   print ('#   experiment that provides restart files.                #')
   print ('#   Domain of new experiment can be same as or smaller     #')
   print ('#   than domain of restart experiment.                     #')
   print ('#                                                          #')
   print ('# RESTART: 2 = Remapping of bcs version and/or resolution. #')
   print ('#   ------------------------------------------------------ #')
   print ('#   New experiment has DIFFERENT tile space (grid) or      #')
   print ('#   DIFFERENT boundary conditions or DIFFERENT snow model  #')
   print ('#   parameter (WEMIN).                                     #')
   print ('#   New and restart experiments *must* be GLOBAL.          #')
   print ('#                                                          #')
   print ('# -------------------------------------------------------- #')
   print ('#                                                          #')
   print ('# Use the following option if your restart files are from  #')
   print ('#   a GEOSgcm experiment:                                  #')
   print ('#                                                          #')
   print ('# RESTART: 3 = No remapping; global only.                  #')
   print ('#   -------------------------------------                  #')
   print ('#   New experiment has SAME tile space (grid), boundary    #')
   print ('#   conditions, and snow model parameter (WEMIN) as the    #')
   print ('#   GCM experiment that provides restart files.            #')
   print ('#   New experiment *must* be GLOBAL.                       #')
   print ('#                                                          #')
   print ('# -------------------------------------------------------- #')
   print ('#                                                          #')
   print ('# Use the following option if you do not have restart      #')
   print ('#   files.                                                 #')
   print ('#                                                          #')
   print ('# RESTART: M = Remap from archived MERRA-2 restart files.  #')
   print ('#   ----------------------------------------------------   #')
   print ('#   Must run on machine with archived MERRA-2 restarts.    #')
   print ('#   New experiment *must* be GLOBAL.                       #')
   print ('#   Restart date/time must match that of an archived       #')
   print ('#   MERRA-2 restart.                                       #')
   print ('#                                                          #')
   print ('# ======================================================== #')
   print ('# IMPORTANT:                                               #')
   print ('#   Except for RESTART=1, SPIN-UP is REQUIRED in almost    #')
   print ('#   all cases.                                             #')
   print ('# ======================================================== #')
   print ('#                                                          #')
   print ('#                                                          #')
   print ('# (ii) Specify experiment ID/location of restart files     #')
   print ('#                                                          #')
   print ('# For RESTART = 1, 2, or 3:                                #')
   print ('#   -----------------------                                #')
   print ('#   Specify RESTART_ID, RESTART_PATH, RESTART_DOMAIN.      #')
   print ('#                                                          #')
   print ('# For RESTART = 3:                                         #')
   print ('#   --------------                                         #')
   print ('#   Place restarts into the following directory tree:      #')
   print ('#     RESTART_PATH/RESTART_ID/output/RESTART_DOMAIN/rs/    #')
   print ('#                                                          #')
   print ('# For RESTART = M:                                         #')
   print ('#   --------------                                         #')
   print ('#   No need to specify RESTART_ID, RESTART_PATH,           #')
   print ('#   and RESTART_DOMAIN.                                    #')
   print ('#                                                          #')
   print ('############################################################')
   print ()
   print ('RESTART:')
   print ('#RESTART_ID:')
   print ('#RESTART_PATH:')
   print ('#RESTART_DOMAIN:')
   print ()
   print ('############################################################')
   print ('#                                                          #')
   print ('#         SURFACE METEOROLOGICAL FORCING                   #')
   print ('#                                                          #')
   print ('#  Surface meteorological forcing time step is in seconds. #')
   print ('#                                                          #')
   print ('#  NOTE:                                                   #')
   print ('#    When forcing is on cube-sphere (CS) grid, must use:   #')
   print ('#    - Model tile space (BCS) derived from same CS grid.   #')
   print ('#    - Nearest-neighbor interpolation (MET_HINTERP: 0).    #')
   print ('#                                                          #')
   print ('#  For more information, see:                              #')
   print ('#    GEOSldas/doc/README.MetForcing_and_BCS.md             #')
   print ('#                                                          #')
   print ('############################################################')
   print ()
   print ('MET_TAG:')
   print ('MET_PATH:')
   print ('FORCE_DTSTEP:')
   print ()
   print ('############################################################')
   print ('#                                                          #')
   print ('#          LAND BOUNDARY CONDITIONS (BCS)                  #')
   print ('#                                                          #')
   print ('#  Path to and (atmospheric) resolution of BCS.            #')
   print ('#    Path includes BCS_VERSION.                            #')
   print ('#                                                          #')
   print ('#  For more information, see:                              #')
   print ('#    GEOSldas/doc/README.MetForcing_and_BCS.md             #')
   print ('#    [..]/GEOSsurface_GridComp/Utils/Raster/make_bcs       #')
   print ('#                                                          #')
   print ('############################################################')
   print ()
   print ('BCS_PATH:')
   print ('BCS_RESOLUTION:')
   print ()
   print ('############################################################')
   
   current_directory = os.path.dirname(__file__)    # path where ldas_setup is
   # rc template files are in [current_directory]/../etc/
   # add defaults from GEOSldas_LDAS.rc
   _fn = current_directory+'/../etc/GEOSldas_LDAS.rc'
   echoInputFile(_fn)
   _fn = current_directory+'/../etc/GEOS_SurfaceGridComp.rc'
   echoInputFile(_fn)


def getExeKeys(option):
   # required keys for exe input file depending on restart option
   rqdExeInpKeys = {'1' : ['EXP_ID',        'EXP_DOMAIN',      'NUM_LDAS_ENSEMBLE',  'BEG_DATE',
                           'END_DATE',      'RESTART_PATH',    'RESTART_DOMAIN',     'RESTART_ID',
                           'MET_TAG',       'MET_PATH',        'FORCE_DTSTEP',       'BCS_PATH',
                           'BCS_RESOLUTION'],
                    '0' : ['EXP_ID',        'EXP_DOMAIN',      'NUM_LDAS_ENSEMBLE',  'BEG_DATE',
                           'END_DATE',
                           'MET_TAG',       'MET_PATH',        'FORCE_DTSTEP',       'BCS_PATH',
                           'BCS_RESOLUTION']
                }

   assert option == '0' or option == '1', '"%s" option is not recognized ' % option
   return rqdExeInpKeys[option]

def getResourceKeys(option):
   # ------
   # Required resource manager input fields
   # ------
   RmInpKeys ={'required' : ['account', 'walltime', 'ntasks_model'],
               'optional' : ['job_name', 'qos', 'oserver_nodes', 'writers-per-node', 'ntasks-per-node', 'constraint']
              }
   assert option in ['required', 'optional'],' "%s" option is not supported' % option

   return RmInpKeys[option]
 

def printResourceInputSampleFile():
   """
    print sample resource manager input file
   """
   requiredKeys = getResourceKeys('required')
   optionalKeys = getResourceKeys('optional')

   print ('#')
   print ('# REQUIRED inputs')
   print ('#')
   print ('# NOTE:')
   print ('# - account          = computational project number')
   print ('#                      [At NCCS: Use command "getsponsor" to see available account number(s).]' )
   print ('# - walltime         = walltime requested; format is HH:MM:SS (hours/minutes/seconds)')
   print ('# - ntasks_model     = number of processors requested for the model (typically 126; output server is not included)')
   print ('#')
   for key in requiredKeys:
       print (key + ':')
   print ()
   print ('#')
   print ('# OPTIONAL inputs')
   print ('#')
   print ('# NOTE:')
   print ('# - job_name         = name of experiment; default is "exp_id"')
   print ('# - qos              = quality-of-service; do not specify by default; specify "debug" for faster but limited service')
   print ('# - oserver_nodes    = number of nodes for oserver ( default is 0, for future use )')
   print ('# - writers-per-node = tasks per oserver_node for writing ( default is 5, for future use );')
   print ('#                        IMPORTANT REQUIREMENT: total #writers = writers-per-node * oserver_nodes >= 2;')
   print ('#                        jobs will hang when oserver_nodes = writers-per-node = 1.')
   print ('# - ntasks-per-node  = requesting fewer ntasks-per-node than total number of cores per node increases allocated memory;')
   print ('#                        ntasks_model should be a multiple of ntasks-per-node') 
   print ('# - constraint       = name of chip set(s) (NCCS default is "[mil|cas]", NAS default is "cas_ait")')
   print ('#')
   for key in optionalKeys:
       print ('#'+key + ':')

def printInputSampleFile(cmdLineArgs):
   '''
     ./ldas_setup sample ...
    
     "sample" sub-command:
       '--exeinp' and '--batinp' are mutually exclusive command line arguments.
       Specifying one will set it to True and set the other one to False.
       That is, we can have either: {'exeinp': False, 'batinp': True }
                                or: {'exeinp': True,  'batinp': False} 
   '''
   if cmdLineArgs['exeinp']:
      printExeInputSampleFile()
   elif cmdLineArgs['batinp']:
      printResourceInputSampleFile()
   else:
      raise Exception('unrecognized sample option')

def printDictionary(d):
    """     
    Private method: print a 'flat' dictionary
    """     
            
    for key, val in d.items():
        print (key.ljust(24), ':', val)

def hours_to_hhmmss(hours):

    # Convert hours to timedelta
    td = timedelta(hours=hours)

    # Extract hours, minutes, seconds
    total_seconds = int(td.total_seconds())
    hours, remainder = divmod(total_seconds, 3600)
    minutes, seconds = divmod(remainder, 60)

    # Format as HHMMSS
    return f"{hours:02d}{minutes:02d}{seconds:02d}"
