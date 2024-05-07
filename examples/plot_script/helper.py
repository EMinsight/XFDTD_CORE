import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.animation as animation
import argparse as arg


def plot_pattern_with_constant_theta(phi, pat_theta, pat_phi, max_val, step_size, n_ring, info_pat_theta, info_pat_phi):
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    plot_range = step_size * n_ring
    min_val = max_val - plot_range
    theta = np.linspace(0, 2 * np.pi, 50)
    circle_x = np.cos(theta)
    circle_y = np.sin(theta)

    # plot ring
    for i in range(n_ring):
        r = (1 / n_ring) * (i + 1)
        ax.plot(r * circle_x, r * circle_y, 'k--', linewidth=0.5)
        ax.text(0.09, r, str(min_val + step_size*(i+1)), fontsize=12,
                verticalalignment='center', horizontalalignment='center')

    # plot line
    r = np.linspace(0, 1, 20)
    for i in range(12):
        theta = i * np.pi / 6
        ax.plot(r * np.cos(theta), r * np.sin(theta), 'k--', linewidth=0.5)
        ax.text(r[-1]*1.08 * np.cos(theta), r[-1]*1.08 * np.sin(theta), str(i * 30) +
                '$^\circ$', fontsize=12, verticalalignment='center', horizontalalignment='center')

    pat_theta[pat_theta < min_val] = min_val
    pat_theta = (pat_theta - min_val) / plot_range
    pat_phi[pat_phi < min_val] = min_val
    pat_phi = (pat_phi - min_val) / plot_range

    x1 = pat_theta * np.cos(phi)
    y1 = pat_theta * np.sin(phi)

    x2 = pat_phi * np.cos(phi)
    y2 = pat_phi * np.sin(phi)

    ax.plot(x1, y1, label=info_pat_theta)
    ax.plot(x2, y2, label=info_pat_phi)

    # hide axis
    ax.axis([-1.3, 1.3, -1.3, 1.3])
    ax.axis('off')
    ax.axis_equal = True
    ax.legend(loc='upper right', fontsize=8)

    return fig, ax


def plot_pattern_with_constant_phi(theta, pat_theta, pat_phi, max_val, step_size, n_ring, info_pat_theta, info_pat_phi):
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    plot_range = step_size * n_ring
    min_val = max_val - plot_range
    circle_x = np.cos(np.linspace(0, 2 * np.pi, 50))
    circle_y = np.sin(np.linspace(0, 2 * np.pi, 50))

    # plot ring
    for i in range(n_ring):
        r = (1 / n_ring) * (i + 1)
        ax.plot(r * circle_x, r * circle_y, 'k--', linewidth=0.5)
        ax.text(0.09, r, str(min_val + step_size*(i+1)), fontsize=12,
                verticalalignment='center', horizontalalignment='center')

    r = np.arange(-1, 1.1, 0.1)
    for i in range(3, 9):
        phi = i * np.pi / 6
        ax.plot(r * np.cos(phi), r * np.sin(phi), 'k--', linewidth=0.5)
        ax.text(r[-1]*1.08 * np.cos(phi), r[-1]*1.08 * np.sin(phi), str((i - 3) * 30) +
                '$^\circ$',
                horizontalalignment='center', fontsize=12, verticalalignment='center')
        ax.text(r[0]*1.08 * np.cos(phi), r[-1]*1.08 * np.sin(phi), str((i - 3) * 30) +
                '$^\circ$',
                horizontalalignment='center', fontsize=12, verticalalignment='center')
    ax.text(0, r[1]*1.18, str(180) +
            '$^\circ$', horizontalalignment='center',
            fontsize=12, verticalalignment='center')

    pat_theta[pat_theta < min_val] = min_val
    pat_theta = (pat_theta - min_val) / plot_range
    pat_phi[pat_phi < min_val] = min_val
    pat_phi = (pat_phi - min_val) / plot_range

    x1 = -pat_theta * np.cos(theta + np.pi / 2)
    y1 = pat_theta * np.sin(theta + np.pi / 2)

    x2 = -pat_phi * np.cos(theta + np.pi / 2)
    y2 = pat_phi * np.sin(theta + np.pi / 2)

    ax.plot(x1, y1, label=info_pat_theta)
    ax.plot(x2, y2, label=info_pat_phi)

    ax.axis([-1.3, 1.3, -1.3, 1.3])
    ax.axis('off')
    ax.axis_equal = True
    ax.legend(loc='upper right', fontsize=8)

    return fig, ax


def plot2d_gif(sorted_data_files, read_func, process_func, vmin=-25, vmax=5, cmap='jet'):
    f, ax = plt.subplots()
    data = process_func(read_func(sorted_data_files[0]))

    im = ax.imshow(data, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.grid(True)
    f.colorbar(im)

    def update(i):
        im.set_data(process_func(read_func(sorted_data_files[i])))
        # show frame number
        ax.set_title(sorted_data_files[i].split('/')[-1])
        return im

    ani = animation.FuncAnimation(f, update, frames=len(
        sorted_data_files), repeat=True, interval=0.5*2000/10, repeat_delay=1000)
    return ani


class FarFieldData:
    def __init__(self,
                 freq,
                 source_power,
                 data_dir,
                 angle: np.ndarray = None) -> None:
        self.freq = freq
        self.source_power = source_power
        self.angle = angle
        self.data_dir = data_dir
        self.load_data()

    def load_data(self):
        prefix = self.format_freq_str()+"_"

        self.a_theta = np.load(os.path.join(
            self.data_dir, prefix + 'a_theta.npy'))
        self.a_phi = np.load(os.path.join(self.data_dir, prefix + 'a_phi.npy'))
        self.f_theta = np.load(os.path.join(
            self.data_dir, prefix + 'f_theta.npy'))
        self.f_phi = np.load(os.path.join(self.data_dir, prefix + 'f_phi.npy'))

    def format_freq_str(self):
        return f'{self.freq/1e9:.2f}GHz'

    @property
    def frequency(self):
        return self.freq

    @property
    def k(self):
        return 2 * np.pi * self.freq / 3e8

    @property
    def radiation_power(self):
        return (8 * np.pi * 120 * np.pi * self.source_power) / self.k ** 2

    @property
    def data_tehta(self):
        return np.abs(120 * np.pi * self.a_theta + self.f_phi) ** 2 / self.radiation_power

    @property
    def data_phi(self):
        return np.abs(-120 * np.pi * self.a_phi + self.f_theta) ** 2 / self.radiation_power

    @property
    def data_theta_dB(self):
        return 10 * np.log10(self.data_tehta)

    @property
    def data_phi_dB(self):
        return 10 * np.log10(self.data_phi)


def plot_rcs_dBsm(data: FarFieldData, const_angle_degree: float, const_angle_symbol: str, plot_phi: bool = True, plot_theta: bool = True, freq_unit: str = 'GHz'):

    fig, ax = plt.subplots()
    rcs_dBsm_theta = data.data_theta_dB
    rcs_dBsm_phi = data.data_phi_dB

    label_prefix = f'{data.freq/1e9:.2f}'
    label_prefix = str(label_prefix) + freq_unit + " "

    if plot_phi:
        ax.plot(data.angle, rcs_dBsm_phi, label=label_prefix +
                r'$RCS_\varphi$' + rf' at ${const_angle_symbol}$=' + str(const_angle_degree) + ' degree')
    if plot_theta:
        ax.plot(data.angle, rcs_dBsm_theta, label=label_prefix +
                r'$RCS_\theta$' + rf' at ${const_angle_symbol}$=' + str(const_angle_degree) + ' degree')
    ax.set_xlabel('Angle (degree)')
    ax.set_ylabel('RCS (dBsm)')
    ax.legend()
    ax.grid()
    return fig, ax


def plot_directivity(data: FarFieldData, const_theta: bool, freq_unit: str = 'GHz', step_size=10, n_ring=4):
    pat1 = data.data_theta_dB
    pat2 = data.data_phi_dB
    max_val = np.max(np.array([pat1, pat2]))
    max_val = step_size * np.ceil(max_val / step_size)
    label_prefix = f'{data.freq/1e9:.2f}'
    label_prefix = str(label_prefix) + freq_unit + " "
    if const_theta:
        fig, ax = plot_pattern_with_constant_theta(
            data.angle, pat1, pat2, max_val, step_size, n_ring, label_prefix + r"$D_{\theta}$", label_prefix + r"$D_{\varphi}$")
    else:
        fig, ax = plot_pattern_with_constant_phi(
            data.angle, pat1, pat2, max_val, step_size, n_ring, label_prefix + r"$D_{\theta}$", label_prefix + r"$D_{\varphi}$")
    return fig, ax


def dft(time, data, freq):
    res = 0+0j
    dt = time[1] - time[0]
    res = np.sum(data*np.exp(-1j*2*np.pi*freq*time))
    return res * dt


def my_plot_arg_parser(default_data_dir: str = 'data', default_save: bool = False, default_save_dir: str = 'save'):
    parser = arg.ArgumentParser()
    parser.add_argument('--data_dir', type=str, default=default_data_dir)
    parser.add_argument('--save', type=bool, default=default_save)
    parser.add_argument('--save_dir', type=str, default=default_save_dir)
    return parser.parse_args()
