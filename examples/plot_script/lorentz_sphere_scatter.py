from time import sleep
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
from helper import *


def plot_relative_permittivity(data_dir: str):
    data = np.load(os.path.join(data_dir, "relative_permittivity.npy"))
    freq = data[0, :]
    eps_r = data[1, :]
    real_data = np.real(eps_r)
    imag_data = np.imag(eps_r)
    sigma = np.abs(imag_data) * 2 * np.pi * freq * 8.85e-12

    print(data.shape)
    f, ax = plt.subplots(2, 1)
    ax[0].plot(freq/1e9, real_data, label='Real')
    ax[0].plot(freq/1e9, imag_data, label='Imaginary')
    ax[0].set_ylim([-10, 10])
    ax[0].set_xlabel('Frequency (GHz)')
    ax[0].set_ylabel('Relative permittivity')
    ax[0].grid()
    ax[0].legend()
    ax[0].set_title(" Relative Permittivity")
    ax[1].plot(freq/1e9, sigma, label='Conductivity')
    ax[1].set_xlabel('Frequency (GHz)')
    ax[1].set_ylabel('Conductivity (S/m)')
    ax[1].grid()
    ax[1].legend()


def dielectric_sphere_scatter(data_dir: str, save_dir: str, save: bool = False):
    x = input("plot move? (y/n)")
    if x == 'y':
        def get_sorted_files(dir: str):
            files = os.listdir(dir)
            files = [f for f in files if f.endswith('.npy')]
            files.sort()
            return [os.path.join(dir, f) for f in files]

        def read_data(file: str): return np.load(file)

        def process_data(data):
            d = np.squeeze(data)
            d[np.abs(d < 1e-4)] = 1e-4
            d = 10 * np.log10(np.abs(d))
            return d

        files = get_sorted_files(os.path.join(data_dir, "movie_ex_xz"))
        files = files[:len(files)]
        ani_ex_xz = plot2d_gif(files, read_data, process_data)
        # files = get_sorted_files(os.path.join(data_dir, "movie_ex_yz"))
        # files = files[:len(files)//2]
        # ani_ex_yz = plot2d_gif(files, read_data, process_data)
        if save:
            ani_ex_xz.save(os.path.join(save_dir, "ex_xz.gif"))
            # ani_ex_yz.save(os.path.join(save_dir, "ex_yz.gif"))
        plt.show()
        exit(0)

    time = np.load(os.path.join(data_dir, "time.npy"))
    incident_wave = np.load(os.path.join(data_dir, "incident_wave.npy"))
    freq = 1e9

    a = dft(time, incident_wave, np.array([freq]))
    power_i = 0.5 * (np.abs(a) ** 2) / (120 * np.pi)
    # xy_data = FarFieldData(
    #     freq=freq,
    #     source_power=power_i,
    #     data_dir=os.path.join(data_dir, "xy"),
    #     angle=np.linspace(-np.pi, np.pi, 360)
    # )

    xz_data = FarFieldData(
        freq=freq,
        source_power=power_i,
        data_dir=os.path.join(data_dir, "xz"),
        angle=np.linspace(-np.pi, np.pi, 360)
    )

    f, ax = plot_rcs_dBsm(xz_data, const_angle_degree=0,
                          const_angle_symbol=r'\varphi', plot_phi=False)
    ax.set_title('XZ Plane')
    ax.set_xlim([0, np.pi])

    yz_data = FarFieldData(
        freq=freq,
        source_power=power_i,
        data_dir=os.path.join(data_dir, "yz"),
        angle=np.linspace(-np.pi, np.pi, 360)
    )

    f, ax = plot_rcs_dBsm(yz_data, const_angle_degree=90,
                          const_angle_symbol=r'\varphi', plot_theta=False)
    ax.set_title('YZ Plane')
    ax.set_xlim([0, np.pi])

    try:
        mie = 10 * \
            np.log10(np.loadtxt(os.path.join(data_dir, "bistatic_rcs.dat")))

        plt.figure()
        plt.plot(np.linspace(0, 180, 180), mie[1, :], label="Mie")
        plt.plot(np.linspace(-180, 180, 360),
                 xz_data.data_theta_dB,  label="FDTD")
        plt.xlabel('Angle (degree)')
        plt.ylabel('RCS (dBsm)')
        plt.xlim([0, 180])
        plt.ylim([-45, 0])
        plt.grid(True)
        plt.legend()
        plt.title('RCS in XZ Plane')
        err = np.abs(mie[1, :] - xz_data.data_theta_dB[-180:])
        plt.figure()
        plt.plot(np.linspace(0, 180, 180), err)
        plt.xlabel('Angle (degree)')
        plt.ylabel('Error (dB)')
        plt.grid(True)
        plt.xlim([0, 180])
        print('mean error = ', np.mean(err))
    except:
        print('Mie data not found')


if __name__ == "__main__":
    data_dir = "../data/dispersive_material_scatter/lorentz_medium"

    args = my_plot_arg_parser(
        default_data_dir=data_dir, default_save_dir='lorentz_medium', default_save=False)

    # plot_relative_permittivity(args.data_dir)
    dielectric_sphere_scatter(args.data_dir, args.save_dir, args.save)
    plt.show()

    # dielectric_sphere_scatter(args.data_dir, args.type)
