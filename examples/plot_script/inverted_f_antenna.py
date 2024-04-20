import numpy as np
import matplotlib.pyplot as plt
import os
from helper import *

def inverted_f_antenna(data_dir: str, save_dir: str, save: bool = False):
    # v1_path = os.path.join(data_dir, 'v1.npy')
    # c1_path = os.path.join(data_dir, 'c1.npy')
    freq_path = os.path.join(data_dir, 'frequencies.npy')
    s11_path = os.path.join(data_dir, 's11.npy')
    power_path = os.path.join(data_dir, 'power.npy')

    # v1 = np.load(v1_path)
    # c1 = np.load(c1_path)
    freq = np.load(freq_path)
    s11 = np.load(s11_path)

    fig_arr = []
    title_arr = []

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.plot(freq/1e9, 20*np.log10(np.abs(s11)), label='$S_{11}$')
    ax.set_xlabel('Frequency (GHz)')
    ax.set_ylabel('Magnitude (dB)')
    ax.set_xlim(0, np.max(freq)/1e9)
    ax.set_ylim(-30, 5)
    ax.legend()
    ax.set_title('S11')
    ax.grid(True)
    fig_arr.append(fig)
    title_arr.append('S11')

    far_field_freq = np.load(power_path)[0, :]
    radiation_power = np.load(power_path)[1, :]

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
        
    if save:
        for i in range(len(fig_arr)):
            fig_arr[i].savefig(os.path.join(save_dir, title_arr[i] + '.png'))

if __name__ == '__main__':
    data_dir = '../data/inverted_f_antenna'
    args = my_plot_arg_parser(default_data_dir=data_dir,default_save_dir='inverted_f_antenna', default_save=False)
    data_dir = args.data_dir
    save_dir = args.save_dir
    save = args.save
    if save and not os.path.exists(save_dir):
        os.makedirs(save_dir)
    inverted_f_antenna(data_dir, save_dir, save)
    if not save:
        plt.show()
    else:
        full_save_path = os.path.realpath(save_dir)
        print(f"Images saved in {full_save_path}")
        
