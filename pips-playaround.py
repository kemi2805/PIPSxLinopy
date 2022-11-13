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

result = subprocess.run(["mpirun","/home/ken/Desktop/PIPS-IPMpp/build/pipsipmLinopyCallbackExample", "-5", str(N), Filepath])
