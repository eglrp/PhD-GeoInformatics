## investigate if modis brdf has lin rel to dmc brdf for sliding win
## basic approach is to assume both brdf fns are kernel fns and then eval kernel fns for different possible angles
## one thing to note is that the modis angles do not vary, so esentially the assumption must be that the brdf adjustment
## due to differing viewing angles remains constant inside the sliding window

import numpy as np
import pylab
import os
import scipy.stats
import matplotlib.pyplot as pyplot
import matplotlib as mpl

