"""
run_experiments.py
"""
from utils import *
from triple_collocation import *
import numpy as np
import sys
"""
Exp 1.

Run the whole time-series thing for both FireCCI products

This should reveal whether estimates remain stable between
the two datasets. Also will show improvement of new fireCCI?

For this just use 2005-2012 as before..
"""




if __name__ =="__main__":



	tile = sys.argv[1]
	outputs, tile_ds  = do_TriCol(tile, N=120, start_year=2005, end_year=2011, fireCCI5=True)
	# save it
	save_experiment("fire_cci50_exp1", outputs, tile, tile_ds)

	# and 41
	outputs, tile_ds = do_TriCol(tile, N=120, start_year=2005, end_year=2011, fireCCI5=False)
	# save it
	save_experiment("fire_cci41_exp1", outputs, tile, tile_ds)

	# run the full time-series for FIRECCI50 too...

	"""
	Exp 1.

	Run the whole time-series thing for both FireCCI products

	This should reveal whether estimates remain stable between
	the two datasets. Also will show improvement of new fireCCI?

	For this just use 2005-2012 as before..
	"""
	outputs, tile_ds = do_TriCol(tile, N=120, start_year=2002, end_year=2014, fireCCI5=True)
	# save it
	save_experiment("FULL_RECORD", outputs, tile, tile_ds)
