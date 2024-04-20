import numpy as np
import matplotlib.pyplot as plt
import os
from helper import *

def microstrip_branch_line_coupler(data_dir):
  freq = np.load(os.path.join(data_dir, 'frequencies.npy'))
  s11 = np.load(os.path.join(data_dir, 's11.npy'))
  s21 = np.load(os.path.join(data_dir, 's21.npy'))
  s31 = np.load(os.path.join(data_dir, 's31.npy'))
  s41 = np.load(os.path.join(data_dir, 's41.npy'))

  plt.figure(figsize=(10, 5))
  plt.plot(freq/1e9, 20*np.log10(np.abs(s11)), label='S11')
  plt.plot(freq/1e9, 20*np.log10(np.abs(s21)), label='S21')
  plt.plot(freq/1e9, 20*np.log10(np.abs(s31)), label='S31')
  plt.plot(freq/1e9, 20*np.log10(np.abs(s41)), label='S41')
  plt.xlabel('Frequency (GHz)')
  plt.ylabel('Magnitude (dB)')
  plt.legend()
  plt.xlim([freq[0]/1e9, freq[-1]/1e9])
  plt.ylim([-40, 0])
  plt.grid()
  plt.figure(figsize=(10, 5))
  plt.plot(freq/1e9, np.angle(s11), label='S11')
  plt.plot(freq/1e9, np.angle(s21), label='S21')
  plt.plot(freq/1e9, np.angle(s31), label='S31')
  plt.plot(freq/1e9, np.angle(s41), label='S41')
  plt.xlabel('Frequency (GHz)')
  plt.ylabel('Phase (rad)')
  plt.legend()
  plt.xlim([freq[0]/1e9, freq[-1]/1e9])
  plt.ylim([-np.pi, np.pi])
  plt.grid()
  
  plt.show()

if __name__ == '__main__':
  microstrip_branch_line_coupler(data_dir='../data/microstrip_branch_line_coupler')