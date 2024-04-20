import numpy as np
import matplotlib.pyplot as plt
import os
from helper import *


def plot_relative_permittivity(data_dir: str):
    data = np.load(os.path.join(data_dir, "relative_permittivity.npy"))
    freq = data[0, :]
    eps_r = data[1, :]
    real_data = np.real(eps_r)
    imag_data = np.imag(eps_r)
    sigma = -imag_data * 2 * np.pi * freq * 8.85e-12

    f, ax = plt.subplots(2, 1, figsize=(6, 10))
    ax[0].plot(freq/1e9, real_data, label='Real')
    ax[0].plot(freq/1e9, imag_data, label='Imaginary')
    ax[0].set_ylim([-5, 5])
    ax[0].set_title('Relative Permittivity')
    ax[0].set_xlabel('Frequency (GHz)')
    ax[0].set_ylabel('Relative permittivity')
    ax[0].grid()
    ax[0].legend()

    ax[1].plot(freq/1e9, sigma, label='Conductivity')
    ax[1].set_xlabel('Frequency (GHz)')
    ax[1].set_ylabel('Conductivity (S/m)')
    ax[1].set_title('Conductivity')
    ax[1].grid()

    return f, ax


def dispersive_sphere_scatter(d_data_dir: str, n_data_dir: str, save_dir: str, m_data_dir: str = None, save: bool = False):
    t_list = []
    f_list = []
    f, _ = plot_relative_permittivity(d_data_dir)
    t_list.append("Relative Permittivity")
    f_list.append(f)

    time = np.load(os.path.join(d_data_dir, "time.npy"))
    incident_wave = np.load(os.path.join(d_data_dir, "incident_wave.npy"))
    freq = 1e9

    a = dft(time, incident_wave, np.array([freq]))
    power_i = 0.5 * (np.abs(a) ** 2) / (120 * np.pi)

    d_xz_data = FarFieldData(
        freq=freq,
        source_power=power_i,
        data_dir=os.path.join(d_data_dir, "xz"),
        angle=np.linspace(-np.pi, np.pi, 360)
    )
    n_xz_data = FarFieldData(
        freq=freq,
        source_power=power_i,
        data_dir=os.path.join(n_data_dir, "xz"),
        angle=np.linspace(-np.pi, np.pi, 360)
    )

    f, ax = plt.subplots(1, 1)
    ax.plot(d_xz_data.angle, d_xz_data.data_theta_dB, label='Dispersive FDTD')
    ax.plot(n_xz_data.angle, n_xz_data.data_theta_dB,
            label='Non-dispersive FDTD')
    if m_data_dir != None:
        mie = np.load(os.path.join(m_data_dir, "mie_xz.npy"))
        mie_rcs = 10 * np.log10(mie[1, :])
        ax.plot(mie[0, :], mie_rcs, label='Mie')
    ax.set_xlabel('Angle (radian)')
    ax.set_ylabel('RCS (dBsm)')
    ax.set_xlim([0, np.pi])
    ax.set_title('RCS in XZ Plane')
    ax.legend()
    ax.grid()

    t_list.append("RCS in XZ Plane")
    f_list.append(f)

    d_yz_data = FarFieldData(
        freq=freq,
        source_power=power_i,
        data_dir=os.path.join(d_data_dir, "yz"),
        angle=np.linspace(-np.pi, np.pi, 360)
    )

    n_yz_data = FarFieldData(
        freq=freq,
        source_power=power_i,
        data_dir=os.path.join(n_data_dir, "yz"),
        angle=np.linspace(-np.pi, np.pi, 360)
    )

    f, ax = plt.subplots(1, 1)
    ax.plot(d_yz_data.angle, d_yz_data.data_phi_dB, label='Dispersive FDTD')
    ax.plot(n_yz_data.angle, n_yz_data.data_phi_dB,
            label='Non-dispersive FDTD')
    if m_data_dir != None:
        mie = np.load(os.path.join(m_data_dir, "mie_yz.npy"))
        mie_rcs = 10 * np.log10(mie[1, :])
        ax.plot(mie[0, :], mie_rcs, label='Mie')
    ax.set_xlabel('Angle (radian)')
    ax.set_ylabel('RCS (dBsm)')
    ax.set_xlim([0, np.pi])
    ax.set_title('RCS in YZ Plane')
    ax.legend()
    ax.grid()

    t_list.append("RCS in YZ Plane")
    f_list.append(f)

    return t_list, f_list


if __name__ == '__main__':

    parser = arg.ArgumentParser()
    parser.add_argument('--d_data_dir', type=str)
    parser.add_argument('--n_data_dir', type=str)
    parser.add_argument('--m_data_dir', type=str, default=None)
    parser.add_argument('--save', type=bool, default=False)
    parser.add_argument('--save_dir', type=str, default='./save')

    args = parser.parse_args()
    save = args.save
    save_dir = args.save_dir
    d_data_dir = args.d_data_dir
    n_data_dir = args.n_data_dir
    m_data_dir = args.m_data_dir

    if not os.path.exists(d_data_dir):
        raise FileNotFoundError(f"Data directory {d_data_dir} not found")
    if not os.path.exists(n_data_dir):
        raise FileNotFoundError(f"Data directory {n_data_dir} not found")
    if m_data_dir != None and not os.path.exists(m_data_dir):
        raise FileNotFoundError(f"Data directory {m_data_dir} not found")
    if not os.path.exists(save_dir) and save:
        print(f"Creating save directory {save_dir}")
        os.makedirs(save_dir)

    t_arr, f_arr = dispersive_sphere_scatter(
        d_data_dir=d_data_dir, n_data_dir=n_data_dir, m_data_dir=m_data_dir, save_dir=save_dir, save=save)

    if save:
        for t, f in zip(t_arr, f_arr):
            f.savefig(os.path.join(save_dir, f"{t}.png"))
        print(f"Saved figures to {save_dir}")
    else:
        plt.show()
