"""
run_experiments.py
"""



from triple_collocation import *




"""
Exp 1.

Run the whole time-series thing for both FireCCI products

This should reveal whether estimates remain stable between
the two datasets. Also will show improvement of new fireCCI?

For this just use 2005-2012 as before..
"""
tile = "h30v10"
outputs = do_TriCol(tile, N=120, start_year=2005, end_year=2012, fireCCI5=True)
# save it
save_experiment("fire_cci50_exp1", outputs)

# and 41
outputs = do_TriCol(tile, N=120, start_year=2005, end_year=2012, fireCCI5=False)
# save it
save_experiment("fire_cci41_exp1", outputs)
