#!/usr/bin/env python

import sys
import os
import shutil
import csv
import imp
import utils
import numpy as np
from collections import OrderedDict


#######################################
### --- Import MadGraph modules --- ###
#######################################

# - We will load the following two modules:
#   madevent_interface: launches event generation (among other things)
#   madgraph.various.banner: editing run_card (among other things)

# - Paths to MadGraph - #
# Note: You need to edit this to correspond to the paths in your system!
MG5_rootDir    = "/home/de3u14/lib/build/hep/MadGraph/MG5_aMC_v2_4_2"
MG5_procDir    = "/scratch/de3u14/MWTC/samples/MWTC_pp_Zh_LO_MG5_242"
MWT_Calculator = "/home/de3u14/lib/projects/MWTC/MWT_Calculator/MWT_Calculator"

sys.path.append(os.path.join(MG5_procDir,'bin','internal'))
sys.path.append( MG5_rootDir )

run_card_path   = os.path.join(MG5_procDir,  "./Cards/run_card.dat")
param_card_path = os.path.join(MG5_procDir,  "./Cards/param_card.dat")

# - Import MadGraph modules - #
import madevent_interface  as ME
import madgraph.various.banner as banner_mod      # run_card
import check_param_card as param_card_mod         # param_card
import logging, logging.config, coloring_logging

# - Logging - #
logging.config.fileConfig(os.path.join(MG5_procDir, 'bin', 'internal', 'me5_logging.conf'))
logging.root.setLevel(logging.ERROR)
logging.getLogger('madevent').setLevel(logging.ERROR)
#logging.getLogger('madgraph').setLevel(logging.ERROR)

#######################################################

# - Output of this script - #
output = "xsec_scan_{}.dat".format( sys.argv[1] )

###################################################################################

# - Access to MadEvent prompt - #
utils.silentrm( os.path.join(MG5_procDir,'RunWeb') )
launch = ME.MadEventCmd(me_dir=MG5_procDir)

# - Cards - #
run_card    = banner_mod.RunCard(       run_card_path   ) 

# - Define empty lists for holding the results - #
param_pts_list = []
xsec_list      = []
xsec_unc_list  = []
result_list    = []

##################################
### --- Cross section scan --- ###
##################################

### --- Parameters --- ###

PS = 0.3
rs = 0.0
mh = 125.0

gt_min   = 1.0
gt_max   = 9.0
gt_nBins = 9

mA_min   =  300.0
mA_max   = 2000.0
mA_nBins = 28

#S_min  = 0.0
#S_max  = 0.5


### ------------------ ###

# - Start of the for loop for the scan - #
#for cba in numpy.linspace(cba_min, cba_max, cba_bins):
#for  tb in numpy.linspace(tb_min,   tb_max, tb_bins):

# - Turn pythia=OFF
# - Turn delphes=OFF
utils.silentrm( os.path.join(MG5_procDir, 'Cards/pythia_card.dat')  )
utils.silentrm( os.path.join(MG5_procDir, 'Cards/delphes_card.dat') )

# - Write data to file - #
run_card['nevents'] = 1000
run_card.write(run_card_path)

data = OrderedDict(
        [
        ('mA'       , None),
        ('gt'       , None),
        ('PS'       , None),
        ('rs'       , None),
        ('mh'       , None),
        ('M1N'      , None),
        ('M2N'      , None),
        ('M1C'      , None),
        ('M2C'      , None),
        ('wh'       , None),
        ('w1N'      , None),
        ('w1C'      , None),
        ('w2N'      , None),
        ('w2C'      , None),
        ('xsec'     , None),
        ('xsec_unc' , None)
        ])

result_str = ""
for i in range(0,len(data)):
    result_str += "{:.3e} "
result_str += "\n"

with open(output, 'w') as f_out:

    for key in data.keys():
        f_out.write( "{} ".format( key ) )

    f_out.write( "\n" )
    f_out.flush()

    for mA in np.linspace( mA_min, mA_max, mA_nBins ):
        for gt in np.linspace( gt_min, gt_max, gt_nBins ):

            # - Start calculation
            utils.silentrm( param_card_path )
            os.system("{} {} {} {} {} {} > {}".format( MWT_Calculator, gt, mA, PS, rs,
                mh, param_card_path ) )

            launch.run_cmd('generate_events -f')
            print( "Event generation now should be finished." )
            print( "Run name: %s" % launch.run_name )

            # - Get results
            data['xsec']     = launch.results.current['cross']
            data['xsec_unc'] = launch.results.current['error']
            
            param_card  = param_card_mod.ParamCard( param_card_path ) 

            data['gt']  = param_card['tcinput'].param_dict[ (1,)  ].value
            data['mA']  = param_card['tcinput'].param_dict[ (2,)  ].value
            data['PS']  = param_card['tcinput'].param_dict[ (3,)  ].value
            data['rs']  = param_card['tcinput'].param_dict[ (4,)  ].value
            data['mh']  = param_card['tcinput'].param_dict[ (5,)  ].value

            data['M1N'] = param_card['mass'].param_dict[ (50,)  ].value
            data['M2N'] = param_card['mass'].param_dict[ (51,)  ].value
            data['M1C'] = param_card['mass'].param_dict[ (52,)  ].value
            data['M2C'] = param_card['mass'].param_dict[ (53,)  ].value

            data['wh']  = param_card['decay'].param_dict[ (25,)  ].value
            data['w1N'] = param_card['decay'].param_dict[ (50,)  ].value
            data['w2N'] = param_card['decay'].param_dict[ (51,)  ].value
            data['w1C'] = param_card['decay'].param_dict[ (52,)  ].value
            data['w2C'] = param_card['decay'].param_dict[ (53,)  ].value

            result = result_str.format( *data.values() )

            # - Store results in the list
            f_out.write( result )
            f_out.flush()

            print( data )

            param_pts_list.append( ( data['gt'], data['mA']) )
            xsec_list.     append( data['xsec']     )
            xsec_unc_list. append( data['xsec_unc'] )

            ### --- WARNING! --- ###
            # - This command deletes and entire directory tree
            # - Be very cautious with this command
            shutil.rmtree(os.path.join(MG5_procDir, 'Events', 'run_01'), ignore_errors=True )
            # - End of for loop - #
 

## - Exit from MadEventCmd - #
launch.run_cmd('quit')
