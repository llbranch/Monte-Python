import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
df = pd.read_csv('monte_carlo_input4ch4_07_10_2023.txt', names=['time','current'], sep=' ')
extra = pd.read_csv('monte_carlo_extradata4chT4_07_10_2023.txt')
fig0, ax0 = plt.subplots()
ax0.plot(df['time']*1e9,df['current'], marker='.')
ax0.set_xlabel('time [ns]')
ax0.set_ylabel('current [A]')
ax0.grid()
plt.show()
fig1, ax1 = plt.subplots()
ax1.scatter(extra['T4_prop_times']/1e3, extra['T4_endpoint_dist'], marker='.')
ax1.set_xlabel('time [ns]')
ax1.set_ylabel('distance [cm]')
ax1.grid()
plt.show()