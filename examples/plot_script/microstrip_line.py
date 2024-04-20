import numpy as np
import matplotlib.pyplot as plt
import os

def microstrip_line(data_dir:str):
  freq = np.load(os.path.join(data_dir, 'frequencies.npy'))
  s11 = np.load(os.path.join(data_dir, 's11.npy'))
  v1 = np.load(os.path.join(data_dir, 'v1.npy'))
  c1 = np.load(os.path.join(data_dir, 'c1.npy'))

  fig,ax = plt.subplots(2,1,figsize=(10, 5))
  ax[0].plot(c1[0,:]*1e9, c1[1,:], label='Current')
  ax[0].set_xlabel('time (ns)')
  ax[0].set_ylabel('Current (A)')
  ax[0].legend()
  ax[0].grid()
  ax[1].plot(v1[0,:]*1e9, v1[1,:], label='Voltage')
  ax[1].set_xlabel('time (ns)')
  ax[1].set_ylabel('Voltage (V)')
  ax[1].legend()
  ax[1].grid()

  plt.figure(figsize=(10, 5))
  plt.plot(freq/1e9, 20*np.log10(np.abs(s11)), label='S11')
  plt.xlabel('Frequency (GHz)')
  plt.ylabel('Magnitude (dB)')
  plt.legend()
  plt.xlim([freq[0]/1e9, freq[-1]/1e9])
  plt.ylim([-40, 0])
  plt.grid()

  # calculate Z0
  def dft(time, data, freq):
    res = np.zeros(freq.shape, dtype=np.complex128)
    for i in range(len(freq)):
      res[i] = np.sum(data*np.exp(-1j*2*np.pi*freq[i]*time))
    return res

  i_fft = dft(c1[0,:], c1[1,:], freq)
  v_fft = dft(v1[0,:], v1[1,:], freq)
  z0 = v_fft/i_fft
  plt.figure(figsize=(10, 5))
  plt.plot(freq/1e9, np.real(z0), label='Z0_real')
  plt.plot(freq/1e9, np.imag(z0), label='Z0_imag')
  plt.xlabel('Frequency (GHz)')
  plt.ylabel('Impedance (Ohm)')
  plt.legend()
  plt.xlim([freq[0]/1e9, freq[-1]/1e9])
  plt.grid()

  plt.show()

if __name__ == '__main__':
  microstrip_line(data_dir='../data/microstrip_line')