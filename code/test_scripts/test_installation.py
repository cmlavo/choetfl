"""
This script verifies the installation of the required libraries for the project.
"""

import platform
import numpy
import pandas
import pytfa
import etfl
import gurobipy

import importlib.metadata

print("Testing library versions:\n" \
f"Python - {platform.python_version()}\n" \
f"NumPy - {numpy.__version__}\n" \
f"Pandas - {pandas.__version__}\n" \
f"PyTFA - {importlib.metadata.version('pytfa')}\n" \
f"ETFL - {importlib.metadata.version('etfl')}\n" \
f"GurobiPy - {importlib.metadata.version('gurobipy')} (This will vary depending on the installed solver)\n" \
)