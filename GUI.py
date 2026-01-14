import functions
import tkinter
from tkinter import filedialog
import os
import requests
import sys
import numpy
import pandas
import pathlib
import pubchempy as pcp
import pandas as pds
import PySide6.QtGui
import matplotlib.pyplot as plt
import functions
if __name__ == '__main__':
    compound = pcp.get_compounds('2-Octene,(Z)-','name')[0]
    print(compound.molecular_formula)
    print(compound.canonical_smiles)
    print(compound.iupac_name)


