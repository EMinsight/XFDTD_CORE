import numpy as np
import matplotlib.pyplot as plt
import os
from helper import *


def half_wave_dipole(data_dir: str, save_dir: str, save: bool = False):
    fig_arr = []
    title_arr = []
    freq_path = os.path.join(data_dir, 'frequencies.npy')
    s11_path = os.path.join(data_dir, 's11.npy')

    freq = np.load(freq_path)
    s11 = np.load(s11_path)

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.plot(freq/1e9, 20*np.log10(np.abs(s11)), label='$S_{11}$')
    ax.set_xlabel('Frequency (GHz)')
    ax.set_ylabel('Magnitude (dB)')
    ax.set_xlim(5, 10)
    ax.set_ylim(-30, 5)
    ax.legend()
    ax.set_title('S11')
    ax.grid(True)
    title_arr.append('S11')
    fig_arr.append(fig)
    

    c1 = np.load(os.path.join(data_dir, 'c1.npy'))
    v1 = np.load(os.path.join(data_dir, 'v1.npy'))

    i_fft = list(dft(c1[0, :], c1[1, :], f) for f in freq)
    v_fft = list(dft(v1[0, :], v1[1, :], f) for f in freq)
    z0 = np.array([v_fft[i] / i_fft[i] for i in range(len(freq))])
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.plot(freq/1e9, np.real(z0), label='Z0_real')
    ax.plot(freq/1e9, np.imag(z0), label='Z0_imag')
    ax.set_xlabel('Frequency (GHz)')
    ax.set_ylabel('Impedance (Ohm)')
    ax.legend()
    ax.set_xlim([freq[0]/1e9, freq[-1]/1e9])
    ax.set_ylim([-1000, 2000])
    ax.grid(True)
    fig_arr.append(fig)
    title_arr.append('Z0')
    
    far_field_freq = np.load(os.path.join(data_dir, 'power.npy'))[0, :]
    radiation_power = np.load(os.path.join(data_dir, 'power.npy'))[1, :]
    
    for i in range(len(far_field_freq)):
        xy_data = FarFieldData(
            freq=far_field_freq[i],
            source_power=radiation_power[i],
            data_dir=os.path.join(data_dir, 'xy'),
            angle=np.linspace(-np.pi, np.pi, 360)
        )

        fig, ax = plot_directivity(xy_data, const_theta=True)
        ax.set_title('XY Plane')
        fig_arr.append(fig)
        title_arr.append(xy_data.format_freq_str() + ' XY Plane')
        
        xz_data = FarFieldData(
            freq=far_field_freq[i],
            source_power=radiation_power[i],
            data_dir=os.path.join(data_dir, 'xz'),
            angle=np.linspace(-np.pi, np.pi, 360)
        )
        
        fig, ax = plot_directivity(xz_data, const_theta=False)
        ax.set_title('XZ Plane')
        fig_arr.append(fig)
        title_arr.append(xz_data.format_freq_str() + ' XZ Plane')
        
        yz_data = FarFieldData(
            freq=far_field_freq[i],
            source_power=radiation_power[i],
            data_dir=os.path.join(data_dir, 'yz'),
            angle=np.linspace(-np.pi, np.pi, 360)
        )
        
        fig, ax = plot_directivity(yz_data, const_theta=False)
        ax.set_title('YZ Plane')
        fig_arr.append(fig)
        title_arr.append(yz_data.format_freq_str() + ' YZ Plane')
    
    return title_arr, fig_arr


if __name__ == '__main__':
    data_dir = '../data/half_wave_dipole'
    args = my_plot_arg_parser(
        default_data_dir=data_dir, default_save=False, default_save_dir='half_wave_dipole')

    save = args.save
    save_dir = args.save_dir

    if not os.path.exists(data_dir):
        raise FileNotFoundError(f'{data_dir} does not exist')
    if save and not os.path.exists(save_dir):
        print(f'{save_dir} does not exist, creating...')
        os.makedirs(save_dir)

    title_arr, fig_arr = half_wave_dipole(data_dir, save_dir, save)
    if save:
        print(f'Saved figures to {save_dir}')
        for t,f in zip(title_arr, fig_arr):
            f.savefig(os.path.join(save_dir, t + '.png'))
    else:
        plt.show()
