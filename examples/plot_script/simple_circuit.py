import numpy as np
import matplotlib.pyplot as plt
import os

def simple_circuit(data_dir:str):
    i_data_path = os.path.join(data_dir, 'i_monitor.npy')
    v_data_path = os.path.join(data_dir, 'v_monitor.npy')
    i_data = np.load(i_data_path)
    v_data = np.load(v_data_path)

    fig, ax = plt.subplots(2, 2, figsize=(9,9))
    ax[0][0].plot(i_data[0,:]*1e9, i_data[1,:])
    ax[0][0].set_title('Current')
    ax[0][0].set_xlabel('Time(ns)')
    ax[0][0].set_ylabel('Current')
    ax[0][0].set_xlim(0, i_data[0,-1]*1e9)
    ax[0][0].set_ylim(np.min(i_data[1,:]), np.max(i_data[1,:]))
    ax[0][0].grid(True)
    ax[0][1].plot(v_data[0,:]*1e9, v_data[1,:])
    ax[0][1].set_title('Voltage')
    ax[0][1].set_xlabel('Time(ns)')
    ax[0][1].set_ylabel('Voltage')
    ax[0][1].set_xlim(0, v_data[0,-1]*1e9)
    ax[0][1].set_ylim(np.min(v_data[1,:]), np.max(v_data[1,:]))
    ax[0][1].grid(True)
    # plot fft
    freq = np.arange(0,1000e6,10e6)
    def dft(time, data, freq):
        result = np.zeros_like(freq, dtype = 'complex_')
        for i in range(len(freq)):
            for j in range(len(time)):
                result[i] += data[j]*np.exp(-1j*2*np.pi*freq[i]*time[j])
        return result
    i_fft = dft(i_data[0,:], i_data[1,:], freq)
    v_fft = dft(v_data[0,:], v_data[1,:], freq)
    # i_fft = np.fft.fft(i_data[1,:])
    # v_fft = np.fft.fft(v_data[1,:])
    # freq = np.fft.fftfreq(len(i_data[1,:]), i_data[0,1] - i_data[0,0])
    ax[1][0].plot(freq/1e6, np.abs(i_fft))
    ax[1][0].set_title('Current')
    ax[1][0].set_xlabel('Frequency (MHz)')
    ax[1][0].set_ylabel('Current')
    ax[1][0].set_xlim(0, 1000)
    # ax[1][0].set_ylim(0, np.max(np.abs(i_fft)))
    ax[1][0].grid(True)
    ax[1][1].plot(freq/1e6, np.abs(v_fft))
    ax[1][1].set_title('Voltage')
    ax[1][1].set_xlabel('Frequency (MHz)')
    ax[1][1].set_ylabel('Voltage')
    ax[1][1].set_xlim(0, 1000)
    # ax[1][1].set_ylim(0, np.max(np.abs(v_fft)))
    ax[1][1].grid(True)

    i_max_idx = np.argmax(np.abs(i_fft))
    v_max_idx = np.argmax(np.abs(v_fft))
    ax[1][0].plot(abs(freq[i_max_idx]/1e6), np.abs(i_fft[i_max_idx]), 'ro')
    ax[1][1].plot(abs(freq[v_max_idx]/1e6), np.abs(v_fft[v_max_idx]), 'ro')
    ax[1][0].text(0.1, 0.8, f"Max freq: {abs(freq[i_max_idx]/1e6)}", transform=ax[1][0].transAxes)
    ax[1][1].text(0.1, 0.8, f"Max freq: {abs(freq[v_max_idx]/1e6)}", transform=ax[1][1].transAxes)
    r = np.abs(v_fft[v_max_idx]/i_fft[i_max_idx])
    print(f"Resistance: {r}")
    plt.show()
    
    
if __name__ == "__main__":
    simple_circuit('../data/simple_circuit')
