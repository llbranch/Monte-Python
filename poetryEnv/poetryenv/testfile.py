""" Imports """
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import find_peaks, peak_prominences, peak_widths
import sys


""" Pass Data in an argument when the file is run"""
ltoutputch1 = str(sys.argv[1])
ltoutputch4 = str(sys.argv[2])

print(f"arg 1:{ltoutputch1}")
print(f"arg 2:{ltoutputch4}")

ch1data = pd.read_csv(ltoutputch1, sep="\t")
ch4data = pd.read_csv(ltoutputch4, sep="\t")

ch1data.columns = ["Time", "V(compout)-CH1"]
ch4data.columns = ["Time", "V(compout)-CH4"]

print(ch1data)
