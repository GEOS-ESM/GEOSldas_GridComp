#!/usr/bin/env python3

import os
import sys
from collections    import OrderedDict

def parseInputFile(inpfile, ladas_cpl=0):
   """ 
     Private method: parse input file and return a dict of options
     Input: input file
     Output: dict
   """ 
   if not os.path.exists(inpfile):
     assert ladas_cpl > 0, " No exeinput file only if ladas_cpl > 0"
     return {}
        
   inpdict = OrderedDict()
   errstr  = "line [%d] of [%s] is not in the form 'key: value'"
   # determine which default values to pick from GEOS_SurfaceGridComp.rc
   # must be in the format '# GEOSldas=>' or # GEOSagcm=>
   if ladas_cpl == 0 :
      use_rc_defaults = '# GEOSldas=>'    # use defaults for LDAS
   else :
      use_rc_defaults = '# GEOSagcm=>'    # use defaults for AGCM
            
   linenum = 0
   with open(inpfile) as fin:
     for line in fin:
       line = line.strip()
       linenum += 1
       if not line:
          continue
       if line.startswith(use_rc_defaults):
          line = line.split(use_rc_defaults)[1]
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
   if ladas_cpl == 0 :
      use_rc_defaults = '# GEOSldas=>'    # use defaults for LDAS
   else :
      use_rc_defaults = '# GEOSagcm=>'    # use defaults for AGCM

   with open (inpfile) as fin :
     for line in fin:
        if line.startswith("## "):
           continue
        if line.startswith(use_rc_defaults):
           line = line.split(use_rc_defaults)[1]
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
   print ('#   To run ldas_setup, the paremeters can be provided in an exeinput file as       #')
   print ('#   this sample file. However, if no such file is provided in case GEOSldas is     #')
   print ('#   coupling with GEOSadas, the inputs should be provided as optional arguments    #')
   print ('#   in the command line:                                                           #')
   print ('#   ./ldas_setup setup /run/folder no_exeinp batinp --                             #')
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
   print ('# (i) Select "RESTART" option:                             #')
   print ('#                                                          #')
   print ('# Use one of the following options if you *have* a         #')
   print ('#   GEOSldas restart file:                                 #')
   print ('#                                                          #')
   print ('# RESTART: 1                                               #')
   print ('#   YES, have restart file from GEOSldas                   #')
   print ('#        in SAME tile space (grid) with SAME boundary      #')
   print ('#        conditions and SAME snow model parameter (WEMIN). #')
   print ('#        The restart domain can be for the same or         #')
   print ('#        a larger one.                                     #')
   print ('#                                                          #')
   print ('# RESTART: 2                                               #')
   print ('#   YES, have restart file from GEOSldas but               #')
   print ('#        in a DIFFERENT tile space (grid) or with          #')
   print ('#        DIFFERENT boundary conditions or DIFFERENT snow   #')
   print ('#        model parameter (WEMIN).                          #')
   print ('#        Restart *must* be for the GLOBAL domain.          #')
   print ('#                                                          #')
   print ('# Use one of the following options if you DO NOT have a    #')
   print ('#   GEOSldas restart file                                  #')
   print ('#   (works for global domain ONLY!):                       #')
   print ('#                                                          #')
   print ('# RESTART: 0                                               #')
   print ('#   Cold start from some old restart for Jan 1, 0z.        #')
   print ('#                                                          #')
   print ('# RESTART: M                                               #')
   print ('#   Re-tile from archived MERRA-2 restart file.            #')
   print ('#                                                          #')
   print ('# -------------------------------------------------------- #')
   print ('# IMPORTANT:                                               #')
   print ('#   Except for RESTART=1, SPIN-UP is REQUIRED in almost    #')
   print ('#   all cases.                                             #')
   print ('# -------------------------------------------------------- #')
   print ('#                                                          #')
   print ('#                                                          #')
   print ('# (ii) Specify experiment ID/location of restart file:     #')
   print ('#                                                          #')
   print ('# For RESTART=1 or RESTART=2:                              #')
   print ('#   Specify RESTART_ID, RESTART_PATH, RESTART_DOMAIN with  #')
   print ('#   restarts stored as follows:                            #')
   print ('#     RESTART_PATH/RESTART_ID/output/RESTART_DOMAIN/rs/    #')
   print ('#                                                          #')
   print ('# For RESTART=0 or RESTART=M:                              #')
   print ('#   There is no need to specify RESTART_ID, RESTART_PATH,  #')
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
   # required keys for exe input file depending on restart options
   rqdExeInpKeys = {'1' : ['EXP_ID',        'EXP_DOMAIN',      'NUM_LDAS_ENSEMBLE',  'BEG_DATE',
                         'END_DATE',      'RESTART_PATH',    'RESTART_DOMAIN',     'RESTART_ID',
                         'MET_TAG',       'MET_PATH',        'FORCE_DTSTEP',       'BCS_PATH',
                         'BCS_RESOLUTION'],
                    '0' : ['EXP_ID',    'EXP_DOMAIN',       'NUM_LDAS_ENSEMBLE',  'BEG_DATE',
                         'END_DATE',       'MET_TAG',         'MET_PATH',           'FORCE_DTSTEP',
                         'BCS_PATH', 'BCS_RESOLUTION']
                }

   assert option == '0' or option == '1', '"%s" option is not recognized ' % option
   return rqdExeInpKeys[option]

def verifyExeInpKeys(ExeInputs):
   option = '1'
   if (ExeInputs['RESTART'] == 'G' or  ExeInputs['RESTART'] == '0'):
      option = '0'
 
   rqdExeInpKeys = getExeKeys(option)
   for key in rqdExeInpKeys:
      assert key in ExeInputs,' "%s" is required in the inputs ( from exeinpfile or command line) ' % (key)   

def getResourceKeys(option):
   # ------
   # Required resource manager input fields
   # ------
   RmInpKeys ={'required' : ['account', 'walltime', 'ntasks_model'],
               'optional' : ['job_name', 'qos', 'oserver_nodes', 'writers-per-node', 'ntasks-per-node', 'constraint']
              }
   assert option in ['required', 'optional'],' "%s" option is not supported' % option

   return RmInpKeys[option]
 
def verifyResourceInputs(ResourceInputs):
   #-----
   # verify resource input keys are correct
   #-----
   rqdRmInpKeys    = getResourceKeys('required')
   optSlurmInpKeys = getResourceKeys('optional')
   allKeys = rqdRmInpKeys + optSlurmInpKeys
   for key in rqdRmInpKeys:
     assert key in ResourceInputs,' "%s" is required in the inputs ( from batinpfile or command line) ' % (key)   

   for key in ResourceInputs:
     assert key in allKeys, ' "%s" is not recognized ' % key 
   
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
   print ('# - qos              = quality-of-service; do not specify by default; specify "debug" for faster but limited service.')
   print ('# - oserver_nodes    = number of nodes for oserver ( default is 0, for future use )')
   print ('# - writers-per-node = tasks per oserver_node for writing ( default is 5, for future use ),')
   print ('#      IMPORTANT REQUIREMENT: total #writers = writers-per-node * oserver_nodes >= 2')
   print ('#      Jobs will hang when oserver_nodes = writers-per-node = 1.')
   print ('# - ntasks-per-node  , allocate ntasks per nodes. Usually it is less then total cores per nodes to get more memory.')
   print ('#                      make ntasks_model be a multiple of ntasks-per-node') 
   print ('# - constraint       = name of chip set(s) (NCCS default is "[mil|cas]", NAS default is "cas_ait").')
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

if __name__=='__main__':
    inpfile = '/gpfsm/dnb34/wjiang/develop_ldas/GEOSldas_nc4/install-SLES15/etc/GEOS_SurfaceGridComp.rc'
    inpdict = parseInputFile(inpfile)
    print(inpdict)

    printExeInputSampleFile()
     
