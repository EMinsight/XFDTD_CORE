import numpy as np
import matplotlib.pyplot as plt
import os

def microstrip_low_pass(data_dir:str):
  freq_path = os.path.join(data_dir, 'frequencies.npy')
  s11_path = os.path.join(data_dir, 's11.npy')
  s21_path = os.path.join(data_dir, 's21.npy')

  freq = np.load(freq_path)
  s11 = np.load(s11_path)
  s21 = np.load(s21_path)

  fig, ax = plt.subplots(2, 1, figsize=(8, 6))
  ax[0].plot(freq/1e9, 20*np.log10(np.abs(s11)), label='$S_{11}$')
  ax[0].set_xlabel('Frequency (GHz)')
  ax[0].set_ylabel('Magnitude (dB)')
  ax[0].set_ylim(-40, 5)
  ax[0].legend()
  ax[0].grid(True)
  ax[1].plot(freq/1e9, 20*np.log10(np.abs(s21)), label='$S_{21}$')
  ax[1].set_xlabel('Frequency (GHz)')
  ax[1].set_xlabel('Frequency (GHz)')
  ax[1].set_ylabel('Magnitude (dB)')
  ax[1].set_ylim(-50, 5)
  ax[1].legend()
  ax[1].grid(True)

  plt.show()




if __name__ == "__main__":
  microstrip_low_pass("../data/microstrip_low_pass")
