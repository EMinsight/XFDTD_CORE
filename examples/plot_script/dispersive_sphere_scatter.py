import numpy as np
import matplotlib.pyplot as plt
import os
from helper import *


def plot_relative_permittivity(data_dir: str):
    data = np.load(os.path.join(data_dir, "relative_permittivity.npy"))
    freq = np.real(data[0, :])
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


def plot_mono_rcs(data_dir: str, incident_wave: np.ndarray, time: np.ndarray):
    freq_unit = 1e9
    dt = time[1] - time[0]

    w_theta = np.load(os.path.join(data_dir, "w_theta.npy"))
    u_phi = np.load(os.path.join(data_dir, "u_phi.npy"))
    w_phi = np.load(os.path.join(data_dir, "w_phi.npy"))
    u_theta = np.load(os.path.join(data_dir, "u_theta.npy"))

    def e_t(w_theta, u_phi, nta: float = 120*np.pi):
        return -nta * w_theta - u_phi

    def e_p(w_phi, u_theta, nta: float = 120*np.pi):
        return -nta * w_phi + u_theta

    e_theta = e_t(w_theta, u_phi)
    e_phi = e_p(w_phi, u_theta)

    e_theta_fft = np.fft.fftshift(np.fft.fft(e_theta))
    e_phi_fft = np.fft.fftshift(np.fft.fft(e_phi))

    # expand incident wave to the same length as the simulation
    incident_wave = np.pad(incident_wave, (0, len(
        w_theta) - len(incident_wave)), 'constant')
    incident_wave_fft = np.fft.fftshift(np.fft.fft(incident_wave))

    freq = np.fft.fftshift(np.fft.fftfreq(len(w_theta), dt))
    print(f"Mono RCS Frequency Resolution: {(freq[1]-freq[0])/freq_unit} GHz")

    nta = 120 * np.pi
    p_i = 0.5 * np.abs(incident_wave_fft) ** 2 / nta
    p_s = 0.5 * (np.abs(e_theta_fft) ** 2 + np.abs(e_phi_fft) ** 2) / nta
    rcs = 10 * np.log10(4 * np.pi * p_s / p_i)

    f, ax = plt.subplots()
    ax.plot(freq/freq_unit, rcs, label='FDTD', linewidth=1.5, marker='v')
    ax.set_xlabel('Frequency (GHz)')
    ax.set_ylabel('RCS (dBsm)')
    ax.set_title('Mono RCS')
    ax.set_xlim([0, 4])
    ax.set_ylim([-50, 0])
    ax.grid()

    if not os.path.exists(os.path.join(data_dir, "mie_mono_rcs.npy")):
        return ['Mono RCS'], [f]

    mie_mono_rcs = np.load(os.path.join(data_dir, "mie_mono_rcs.npy"))

    mie_freq = np.real(mie_mono_rcs[0, :])
    mie_rcs = 10 * np.log10(np.abs(mie_mono_rcs[1, :]))

    ax.plot(mie_freq/freq_unit, mie_rcs,
            label='Mie', linewidth=1.5, marker='o')
    ax.legend()

    return ['Mono RCS'], [f]


def plot_bistatic_rcs(d_data_dir: str, n_data_dir: str, incident_wave: np.ndarray, time: np.ndarray):
    t_list = []
    f_list = []

    freq = 1e9

    a = dft(time, incident_wave, np.array([freq]))
    power_i = 0.5 * (np.abs(a) ** 2) / (120 * np.pi)

    f, ax = plt.subplots(1, 1)
    if os.path.exists(os.path.join(d_data_dir)):
        d_xz_data = FarFieldData(
            freq=freq,
            source_power=power_i,
            data_dir=os.path.join(d_data_dir, 'xz'),
            angle=np.linspace(-np.pi, np.pi, 360)
        )

        ax.plot(d_xz_data.angle, d_xz_data.data_theta_dB,
                label='Dispersive FDTD')

    if os.path.exists(os.path.join(n_data_dir)):
        n_xz_data = FarFieldData(
            freq=freq,
            source_power=power_i,
            data_dir=os.path.join(n_data_dir, 'xz'),
            angle=np.linspace(-np.pi, np.pi, 360)
        )

        ax.plot(n_xz_data.angle, n_xz_data.data_theta_dB,
                label='Non-dispersive FDTD')

    if os.path.exists(os.path.join(d_data_dir, 'xz', 'mie_bistatic_rcs.npy')):
        mie_bistatic_rcs = np.load(os.path.join(
            d_data_dir, 'xz', 'mie_bistatic_rcs.npy'))
        mie_angle = np.real(mie_bistatic_rcs[0, :])
        mie_rcs = 10 * np.log10(np.abs(mie_bistatic_rcs[1, :]))
        ax.plot(mie_angle, mie_rcs, label='Mie')

    ax.set_xlabel('Angle (radian)')
    ax.set_ylabel('RCS (dBsm)')
    ax.set_xlim([0, np.pi])
    ax.set_title('Bistatic RCS in XZ Plane')
    ax.legend()
    ax.grid()

    t_list.append("Bistatic RCS in XZ Plane")
    f_list.append(f)

    return t_list, f_list


def dispersive_sphere_scatter(d_data_dir: str, n_data_dir: str):
    t_list = []
    f_list = []

    f, _ = plot_relative_permittivity(d_data_dir)
    t_list.append("Relative Permittivity")
    f_list.append(f)

    time = np.load(os.path.join(d_data_dir, "time.npy"))
    incident_wave = np.load(os.path.join(d_data_dir, "incident_wave.npy"))

    t, f = plot_bistatic_rcs(os.path.join(
        d_data_dir, 'fd'), os.path.join(n_data_dir, 'fd'), incident_wave, time)

    t_list.extend(t)
    f_list.extend(f)

    t, f = plot_mono_rcs(os.path.join(d_data_dir, "td"), incident_wave, time)

    t_list.extend(t)
    f_list.extend(f)

    return t_list, f_list


if __name__ == '__main__':

    parser = arg.ArgumentParser()
    parser.add_argument('--d_data_dir', type=str)
    parser.add_argument('--n_data_dir', type=str)
    parser.add_argument('--save', type=bool, default=False)
    parser.add_argument('--save_dir', type=str, default='./save')

    args = parser.parse_args()
    save = args.save
    save_dir = args.save_dir
    d_data_dir = args.d_data_dir
    n_data_dir = args.n_data_dir

    if not os.path.exists(d_data_dir):
        raise FileNotFoundError(f"Data directory {d_data_dir} not found")

    if not os.path.exists(save_dir) and save:
        print(f"Creating save directory {save_dir}")
        os.makedirs(save_dir)

    t_arr, f_arr = dispersive_sphere_scatter(
        d_data_dir=d_data_dir, n_data_dir=n_data_dir)

    if save:
        for t, f in zip(t_arr, f_arr):
            t = t.replace(" ", "_")
            f.savefig(os.path.join(save_dir, f"{t}.png"))
        print(f"Saved figures to {save_dir}")
    else:
        plt.show()
