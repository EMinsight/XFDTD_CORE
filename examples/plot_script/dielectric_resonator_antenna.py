import numpy as np
import matplotlib.pyplot as plt
import os
from helper import *


def plot_each_freq(data_dir: str):
    title_list = []
    fig_list = []

    power_path = os.path.join(data_dir, 'power.npy')
    freq = np.load(power_path)[0, :]
    radiation_power = np.load(power_path)[1, :]
    for i in range(len(freq)):
        xy_far_field = FarFieldData(
            freq=freq[i],
            source_power=radiation_power[i],
            data_dir=os.path.join(data_dir, 'xy_plane'),
            angle=np.linspace(-np.pi, np.pi, 360)
        )

        xz_far_field = FarFieldData(
            freq=freq[i],
            source_power=radiation_power[i],
            data_dir=os.path.join(data_dir, 'xz_plane'),
            angle=np.linspace(-np.pi, np.pi, 360)
        )

        yz_far_field = FarFieldData(
            freq=freq[i],
            source_power=radiation_power[i],
            data_dir=os.path.join(data_dir, 'yz_plane'),
            angle=np.linspace(-np.pi, np.pi, 360)
        )

        xy_fig, xy_ax = plot_directivity(xy_far_field, const_theta=True)
        xy_ax.set_title('XY Plane')
        title_list.append(f'XY Plane {freq[i]/1e9:.2f} GHz')
        fig_list.append(xy_fig)

        xz_fig, xz_ax = plot_directivity(xz_far_field, const_theta=False)
        xz_ax.set_title('XZ Plane')
        title_list.append(f'XZ Plane {freq[i]/1e9:.2f} GHz')
        fig_list.append(xz_fig)

        yz_fig, yz_ax = plot_directivity(yz_far_field, const_theta=False)
        yz_ax.set_title('YZ Plane')
        title_list.append(f'YZ Plane {freq[i]/1e9:.2f} GHz')
        fig_list.append(yz_fig)
    return title_list, fig_list


def dielectric_resonator_antenna(data_dir: str):
    title_list = []
    fig_list = []

    freq = np.load(os.path.join(data_dir, 'frequencies.npy'))
    s11 = np.load(os.path.join(data_dir, 's11.npy'))

    fig, ax = plt.subplots()
    ax.plot(freq/1e9, 20*np.log10(np.abs(s11)))
    ax.set_xlabel('Frequency (GHz)')
    ax.set_ylabel('S11 (dB)')
    ax.set_title('S11')
    ax.grid()
    title_list.append('S11')
    fig_list.append(fig)

    v1 = np.load(os.path.join(data_dir, 'v1.npy'))
    c1 = np.load(os.path.join(data_dir, 'c1.npy'))

    fig, ax = plt.subplots(2, 1)
    ax[0].plot(v1[0, :]*1e9, v1[1, :], label='v1')
    ax[0].set_xlabel('Time (ns)')
    ax[0].set_ylabel('Voltage (V)')
    ax[0].grid()
    ax[0].legend()
    ax[0].set_xlim(0, 5)
    ax[1].plot(c1[0, :]*1e9, c1[1, :], label='c1')
    ax[1].set_xlabel('Time (ns)')
    ax[1].set_ylabel('Current (A)')
    ax[1].grid()
    ax[1].legend()
    ax[1].set_xlim(0, 2)
    fig.suptitle('Voltage and Current')
    title_list.append('Voltage and Current')
    fig_list.append(fig)

    t_l, f_l = plot_each_freq(data_dir)

    title_list.extend(t_l)
    fig_list.extend(f_l)

    return title_list, fig_list


if __name__ == '__main__':

    data_dir = '../data/dielectric_resonator_antenna'
    args = my_plot_arg_parser(default_data_dir=data_dir, default_save=False,
                              default_save_dir='./dielectric_resonator_antenna')
    if not os.path.exists(data_dir):
        raise FileNotFoundError(f'{data_dir} does not exist')

    if args.save and not os.path.exists(args.save_dir):
        print(f'Creating {args.save_dir}')
        os.makedirs(args.save_dir, exist_ok=True)

    title_list, fig_list = dielectric_resonator_antenna(
        data_dir=args.data_dir)
    
    if args.save:
        print(f'Saving to {args.save_dir}')
        for t,f in zip(title_list, fig_list):
            f.savefig(os.path.join(args.save_dir, f'{t}.png'))
    else :
        plt.show()
