#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 21:56:23 2021

@author: fabian
"""

from linopy import Model
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import pypsa
import subprocess
import sys

import logging

logging.basicConfig(level=logging.INFO)

n = pypsa.examples.scigrid_de()
m = n.optimize.create_model()

N = 24
blocks = np.repeat(np.arange(1, N + 1), len(n.snapshots) / N)
m.blocks = xr.DataArray(blocks, [n.snapshots])

Filepath = '/home/ken/Desktop/pypsa_with_PIPS/pypsa-model'
m.to_block_files(Filepath)

m.solve()
original_stdout = sys.stdout
#with open('/home/ken/Desktop/pypsa_with_PIPS/pypsa-model/Python_Solution.txt', 'w') as f:
#    sys.stdout = f # Change the standard output to the file we created.
#    print(m.solution.variables)
#    sys.stdout = original_stdout # Reset the standard output to its original value
result = subprocess.run(["/home/ken/Desktop/Linopy_PIPS/PIPS-IPMpp/build/pipsipmLinopyCallback", "Lückenfüller", str(N), Filepath])
