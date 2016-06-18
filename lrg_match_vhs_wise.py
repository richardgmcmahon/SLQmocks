from __future__ import print_function, division

"""

WARNING: logging sandpit testing in progress

TODO: Outer join maybe to determine why some LRGs do not match with 
WISE and/or VHS.


needs /home/sr525/ etc

python2.7 wise_vhs_match.py filename ra_col_name dec_col_name 

Based on:
python2.7 /home/sr525/Python_Code/wise_vhs_match.py RA DEC \
/data/des/RedMagic_LRGs/sva1_gold_1.0.2_run_redmapper_v6.3.3_redmagic_0.5-10.fit



HISTORY:

201603xx


"""

import os
import sys
import socket
import math
import time
t0 = time.clock()

import inspect
import imp
import traceback

import glob
import logging
import logging.handlers

now = now = time.localtime(time.time())
datestamp = time.strftime("%Y%m%d",now)

module_name =  os.path.splitext(__file__)[0]
print('Module name: ', module_name)
prefix = module_name 

from argparse import ArgumentParser
parser = ArgumentParser()

parser.set_defaults(nomatch=False)
parser.add_argument("--nomatch", action='store_true',
    dest = 'nomatch', help = "skip matching")

parser.set_defaults(debug=False)
parser.add_argument("--debug", action='store_true',
    dest='debug', help="debug option")

parser.set_defaults(showplots=False)
parser.add_argument("--showplots", action='store_true',
    dest='showplots', help="showplots option")

parser.set_defaults(verbose=False)
parser.add_argument("--verbose", action='store_true',
    dest='verbose', help="verbose option")

parser.set_defaults(xkcd=False)
parser.add_argument("--xkcd", action='store_true',
    dest='xkcd', help="xkcd cartoon plot style")

args = parser.parse_args()

debug = args.debug
print('debug: ', debug)

showplots = args.showplots
print('showplots: ', showplots)

verbose = args.verbose
print('verbose: ', verbose)

LOG_FILENAME = prefix + '_' + datestamp + '.log'

# Set up a specific logger with our desired output level
my_logger = logging.getLogger('MyLogger')
my_logger.setLevel(logging.DEBUG)

# Add the log message handler to the logger
handler = logging.handlers.RotatingFileHandler(LOG_FILENAME,
                                               maxBytes=20,
                                               backupCount=5)
my_logger.addHandler(handler)

# Log some messages
for i in range(20):
    my_logger.debug('i = %d' % i)

 
# add filemode="w" to overwrite
# add a date or time stamp
logging.basicConfig(filename="sample.log", filemode="w", level=logging.INFO)

logging.debug("This is a debug message")
logging.info("Informational message")
logging.error("An error has happened!")


import inspect

import matplotlib as mpl
print('matplotlib: ', mpl.__version__)
print('matplotlib.matplotlib_fname(): ',mpl.matplotlib_fname())

# Force matplotlib to not use any Xwindows backend
# useful when you want to turn off plt.show()
try: mpl.use('Agg')
except: pass

# Fonts, latex:
mpl.rc('font',**{'family':'serif', 'serif':['TimesNewRoman']})
mpl.rc('text', usetex=True)
#mpl.rcParams['text.usetex']=True


import matplotlib.pyplot as plt
import numpy as np
print('numpy: ', np.__version__)

import scipy
print('scipy: ', scipy.__version__)
from scipy.ndimage import interpolation 

from scipy import signal, ndimage

import astropy
print('astropy: ', astropy.__version__)
print(inspect.getfile(astropy))

from astropy.table import Table
from astropy.io import fits

sys.path.append('/home/rgm/soft/OM10/OM10/')
import om10

sys.path.append('/home/sr525/Python_Code/')
import srpylib as srl

# modulised version of Sophie's wise_vhs_match
from wise_vhs_match import wise_vhs_match
print(inspect.getfile(wise_vhs_match))
# help(wise_vhs_match)

# turn tex off since turned on by om10
mpl.rc('text', usetex=False)

try:
  import psfex
except ImportError:
  print("ImportError: psfex not available")
  print()

from astroML.plotting import hist

sys.path.append('/home/rgm/soft/python/lib/')
import librgm as rgm
from librgm.table_index_column import table_index_column
from librgm.plotid import plotid
from librgm.table_stats import table_stats
#help(rgm.table_stats)


logging.basicConfig(level=logging.DEBUG, 
    format='%(asctime)s - %(levelname)s - % (message)s',
                    datefmt='%m/%d/%Y %I:%M:%S %p')
pid = os.getpid()
print('Current working directory: %s' % (os.getcwd()))
print('User: ', os.getenv('USER'))
print('Hostname: ', os.getenv('HOSTNAME'))
print('Host:     ', os.getenv('HOST'))
print('Hostname: ', socket.gethostname())
print()
print('__file__: ', __file__)
print('__name__: ', __name__)
#print('__module__', __module__)
#print('__init__', __init__)

print(inspect.getsource(inspect.getsource))
help(inspect.getargspec)
arg_spec = inspect.getargspec(__name__)
print('NAMES   :', arg_spec[0])
print('*       :', arg_spec[1])
print('**      :', arg_spec[2])
print('defaults:', arg_spec[3])


print('inspect.stack depth: ', len(inspect.stack()))
print('Running: ',inspect.stack()[0][3])
if len(inspect.stack()) > 1: print('Caller: ',inspect.stack()[1][3])
print()

for level in inspect.stack():
    frame, filename, line_num, func, src_code, src_index = level
    print(level)
    print(inspect.getargvalues(frame))
    print()

frame = inspect.currentframe()
print('inspect.currentframe')
args, _, _, values = inspect.getargvalues(frame)
print('function name "%s"' % inspect.getframeinfo(frame)[2])
for i in args:
    print("    %s = %s" % (i, values[i]))


(frame, filename, line_number, function_name, lines, index) = \
    inspect.getouterframes(inspect.currentframe())[0]
print(frame, filename, line_number, function_name, lines, index)


trace = traceback.extract_stack()
progname_str=os.path.basename(trace[0])
print('len(trace): ', len(trace))
for each in trace:
    print('trace: ', each)


infile = '/data/des/RedMagic_LRGs/sva1_gold_1.0.2_run_redmapper_v6.3.3_redmagic_0.5-10.fit'



t1 = time.clock()
print(time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime()))
print()
print('Preamble complete')
print("Elapsed Time:", t1 - t0)
print('Starting wise_vhs_match')
#nrows=1000
nrows=None
result = wise_vhs_match(infile, 'RA','DEC', 
    prefix='lrg_match_vhs_wise', nrows=nrows)

result.meta['infile']= infile

result.info()
result.info('stats')


# See what files are created
logfiles = glob.glob('%s*' % LOG_FILENAME)
for filename in logfiles:
    print(filename)

print(time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime()))
t1 = time.clock()
print("Elapsed Time:", t1 - t0)


