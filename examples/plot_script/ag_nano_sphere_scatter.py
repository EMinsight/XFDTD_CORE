from helper import my_plot_arg_parser, plot2d_gif, FarFieldData, dft
import os
import numpy as np
import matplotlib.pyplot as plt


class MyTHz:

    @staticmethod
    def __repr__():
        return 'THz'

    @staticmethod
    def unit():
        return 1e12


class MyGHz:

    @staticmethod
    def __repr__():
        return 'GHz'

    @staticmethod
    def unit():
        return 1e9


def plot_relative_permittivity(data_dir: str, freq_unit=MyGHz()):

    data = np.load(os.path.join(data_dir, "relative_permittivity.npy"))
    freq = np.real(data[0, :])
    eps_r = data[1, :]
    real_data = np.real(eps_r)
    imag_data = np.imag(eps_r)
    sigma = -imag_data * 2 * np.pi * freq * 8.85e-12

    f, ax = plt.subplots(2, 1, figsize=(6, 10))
    ax[0].plot(freq/freq_unit.unit(), real_data, label='Real')
    ax[0].plot(freq/freq_unit.unit(), imag_data, label='Imaginary')
    ax[0].set_ylim([-50, 5])
    ax[0].set_xlim([300, 1200])

    ax[0].set_title('Relative Permittivity')
    ax[0].set_xlabel(f'Frequency ({freq_unit})')
    ax[0].set_ylabel('Relative permittivity')
    ax[0].grid()
    ax[0].legend()

    ax[1].plot(freq/freq_unit.unit(), sigma, label='Conductivity')
    ax[1].set_xlabel(f'Frequency ({freq_unit})')
    ax[1].set_ylabel('Conductivity (S/m)')
    ax[1].set_title('Conductivity')
    ax[1].set_xlim([300, 1200])
    ax[1].set_ylim([0, 1e4])
    ax[1].grid()

    return f, ax


def plot_mono_rcs(drude_data_dir: str, m_lorentz_data_dir: str, incident_wave: np.ndarray, time: np.ndarray, freq_unit=MyGHz()):

    dt = time[1] - time[0]

    nta = 120 * np.pi

    def e_field(data_dir: str, nta: float = 120*np.pi):
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

        return e_theta, e_phi

    if os.path.exists(os.path.join(drude_data_dir)):
        e_theta, e_phi = e_field(drude_data_dir)
        e_theta_fft_drude = np.fft.fftshift(np.fft.fft(e_theta))
        e_phi_fft_drude = np.fft.fftshift(np.fft.fft(e_phi))
        p_s_drude = 0.5 * (np.abs(e_theta_fft_drude) ** 2 +
                           np.abs(e_phi_fft_drude) ** 2) / nta

    if os.path.exists(os.path.join(m_lorentz_data_dir)):
        e_theta, e_phi = e_field(m_lorentz_data_dir)
        e_theta_fft_lor = np.fft.fftshift(np.fft.fft(e_theta))
        e_phi_fft_lor = np.fft.fftshift(np.fft.fft(e_phi))
        p_s_lor = 0.5 * (np.abs(e_theta_fft_lor) ** 2 +
                         np.abs(e_phi_fft_lor) ** 2) / nta

    incident_wave = np.pad(incident_wave, (0, len(
        e_theta) - len(incident_wave)), 'constant')
    incident_wave_fft = np.fft.fftshift(np.fft.fft(incident_wave))

    freq = np.fft.fftshift(np.fft.fftfreq(len(e_theta), dt))
    print(
        f"Mono RCS Frequency Resolution: {(freq[1]-freq[0])/freq_unit.unit()} {freq_unit}")

    nta = 120 * np.pi
    p_i = 0.5 * np.abs(incident_wave_fft) ** 2 / nta

    f, ax = plt.subplots()
    ax.set_xlabel(f'Frequency ({freq_unit})')
    ax.set_ylabel('RCS (dBsm)')
    ax.set_title('Mono RCS')
    ax.set_xlim([300, 1200])
    ax.set_ylim([-165, -130])
    ax.grid()
    
    if os.path.exists(os.path.join(drude_data_dir)):
        rcs = 10 * np.log10(4 * np.pi * p_s_drude / p_i)
        ax.plot(freq/freq_unit.unit(), rcs, label='Drude Model',
                linewidth=1.5, marker='v')
        

    if os.path.exists(os.path.join(m_lorentz_data_dir)):
        rcs = 10 * np.log10(4 * np.pi * p_s_lor / p_i)
        ax.plot(freq/freq_unit.unit(), rcs, label='MLorentz Model',
                linewidth=1.5, marker='.')
    

    if not os.path.exists(os.path.join(drude_data_dir, "ag_drude_mie_mono_rcs.npy")):
        return ['mono_rcs'], [f]

    mie_mono_rcs = np.load(os.path.join(
        drude_data_dir, "ag_drude_mie_mono_rcs.npy"))

    mie_freq = np.real(mie_mono_rcs[0, :])
    mie_rcs = 10 * np.log10(np.abs(mie_mono_rcs[1, :]))

    ax.plot(mie_freq/freq_unit.unit(), mie_rcs,
            label='Mie', linewidth=1.5)
    ax.legend()

    return ['mono_rcs'], [f]


def au_nano_sphere_scatter(drude_data_dir, m_lor_data_dir):
    t_list = []
    f_list = []

    freq_unit = MyTHz()

    f, ax = plot_relative_permittivity(drude_data_dir, freq_unit=freq_unit)
    t_list.append('relative_permittivity')
    f_list.append(f)

    incident_wave = np.load(os.path.join(drude_data_dir, 'incident_wave.npy'))
    i_time = incident_wave[0, :]
    i_value = incident_wave[1, :]

    # plot sampled electric field
    drude_ex_point = np.load(os.path.join(drude_data_dir, 'ex_point.npy'))
    m_lor_ex_point = np.load(os.path.join(m_lor_data_dir, 'ex_point.npy'))

    ex_time = drude_ex_point[0, :]
    drude_ex_value = drude_ex_point[1, :]
    m_lor_ex_value = m_lor_ex_point[1, :]

    f, ax = plt.subplots(2, 1, figsize=(8, 10))
    ax[0].plot(i_time*1e12, i_value, label='Incident Wave', linewidth=2)
    ax[0].plot(ex_time*1e12, drude_ex_value, label='Drude Model',
               linestyle='dashdot', linewidth=2)
    ax[0].plot(ex_time*1e12, m_lor_ex_value,
               label='Lorentz Model', linestyle=':', linewidth=2)

    ax[0].set_title('Ag Nano Sphere Scatter Time Domain')
    ax[0].set_xlabel('Time (fs)')
    ax[0].set_ylabel('Amplitude')
    ax[0].set_xlim([i_time[0]*1e12, i_time[-1]*1e12])
    ax[0].set_ylim([-1, 1])
    ax[0].grid()
    ax[0].legend()

    i_fft = np.fft.fftshift(np.fft.fft(i_value))
    drude_ex_fft = np.fft.fftshift(np.fft.fft(drude_ex_value))
    m_lor_ex_fft = np.fft.fftshift(np.fft.fft(m_lor_ex_value))
    freq = np.fft.fftshift(np.fft.fftfreq(len(i_value), i_time[1] - i_time[0]))
    print(f"Resolution: {(freq[1] - freq[0])/freq_unit.unit()} {freq_unit}")
    freq_min = 0
    freq_max = freq[-1]
    # find index
    freq_min_idx = np.argmin(np.abs(freq - freq_min))
    freq_max_idx = np.argmin(np.abs(freq - freq_max))
    freq = freq[freq_min_idx:freq_max_idx]

    i_fft = i_fft[freq_min_idx:freq_max_idx]
    drude_ex_fft = drude_ex_fft[freq_min_idx:freq_max_idx]
    m_lor_ex_fft = m_lor_ex_fft[freq_min_idx:freq_max_idx]

    max_ex = np.max(np.abs(i_fft))

    ax[1].plot(freq/freq_unit.unit(), np.abs(i_fft)/max_ex,
               label='Incident Wave', linewidth=2)
    ax[1].plot(freq/freq_unit.unit(), np.abs(drude_ex_fft)/max_ex,
               label='Drude Model', linestyle='dashdot', linewidth=2)
    ax[1].plot(freq/freq_unit.unit(), np.abs(m_lor_ex_fft)/max_ex,
               label='Lorentz Model', linestyle=':', linewidth=2)
    ax[1].set_title('Normalized Amplitude Spectrum')
    ax[1].set_xlabel(f'Frequency ({freq_unit})')
    ax[1].set_ylabel('Amplitude')
    ax[1].set_xlim([300, 1200])
    ax[1].set_ylim([0, 1.25])
    ax[1].set_yticks(np.arange(0, 1.25, 0.25))
    ax[1].grid()
    ax[1].legend()
    t_list.append('sampled_electric_field')
    f_list.append(f)

    # plot mono rcs
    t, f = plot_mono_rcs(os.path.join(
        drude_data_dir, 'td'), os.path.join(
        m_lor_data_dir, 'td'), incident_wave[1, :], i_time, freq_unit=MyTHz())
    t_list.extend(t)
    f_list.extend(f)

    # plot bistatic rcs
    f, ax = plt.subplots()
    freq = 800 * freq_unit.unit()
    i_dft = dft(i_time, i_value, freq)
    
    if os.path.exists(os.path.join(drude_data_dir, 'fd', 'xz')):
        xz_data_drude = FarFieldData(
            freq=freq,
            source_power=0.5 * np.abs(i_dft) ** 2 / 120 / np.pi,
            data_dir=os.path.join(drude_data_dir, 'fd', 'xz'),
            angle=np.linspace(-np.pi, np.pi, 360))

        ax.plot(xz_data_drude.angle, xz_data_drude.data_theta_dB, label='FDTD', linewidth=1.5, marker='v')
    
    if os.path.exists(os.path.join(m_lor_data_dir, 'fd', 'xz')):
        xz_data_lor = FarFieldData(
            freq=freq,
            source_power=0.5 * np.abs(i_dft) ** 2 / 120 / np.pi,
            data_dir=os.path.join(m_lor_data_dir, 'fd', 'xz'),
            angle=np.linspace(-np.pi, np.pi, 360))
        
        ax.plot(xz_data_lor.angle, xz_data_lor.data_theta_dB, label='MLorentz', marker='.' ,linewidth=1.5)

    if os.path.exists(os.path.join(drude_data_dir, 'fd', 'xz', 'ag_drude_mie_bistatic_rcs.npy')):
        mie_bistatic_rcs = np.load(os.path.join(
            drude_data_dir, 'fd', 'xz', 'ag_drude_mie_bistatic_rcs.npy'))
        mie_angle = np.real(mie_bistatic_rcs[0, :])
        mie_rcs = 10 * np.log10(np.abs(mie_bistatic_rcs[1, :]))
        ax.plot(mie_angle, mie_rcs, label='Mie', linestyle='--', linewidth=1.5)

    ax.set_title('Bistatic RCS')
    ax.set_xlabel('Angle (rad)')
    ax.set_ylabel('RCS (dBsm)')
    ax.set_xlim([0, np.pi])
    ax.set_xticks(np.linspace(0, np.pi, 5))
    ax.set_ylim([-165, -130])
    ax.legend()
    ax.grid()

    t_list.append('bistatic_rcs')
    f_list.append(f)

    # plot 2d gif
    x = input('Input y to plot the gif: ')
    if x != 'y':
        return t_list, f_list

    def sorted_files(data_dir):
        files = os.listdir(data_dir)
        data = sorted([os.path.join(data_dir, file)
                      for file in files if file.endswith('.npy')])
        min_i = 0
        # max_i = 7*len(data)//10
        data = data[min_i:]
        data = [data[i] for i in range(0, len(data), 1)]
        return data

    def process_func(data):
        data = np.squeeze(data)
        data[np.abs(data) < 1e-4] = 1e-4
        return 20*np.log10(np.abs(data))

    def read_func(file):
        return np.load(file)

    ani = plot2d_gif(sorted_files(os.path.join(
        drude_data_dir, 'movie_xz')), read_func, process_func, cmap='inferno')

    plt.show()

    return t_list, f_list


if __name__ == '__main__':
    # args = my_plot_arg_parser(default_data_dir='../data/ag_nano_sphere_scatter',
    #                           default_save_dir='./ag_nano_sphere_scatter')

    drude_data_dir = '../data/ag_nano_sphere_scatter_ag_drude'
    m_lor_data_dir = '../data/ag_nano_sphere_scatter_ag_m_lor'
    save_dir = './ag_nano_sphere_scatter'

    if not os.path.exists(drude_data_dir):
        raise FileNotFoundError(f'{drude_data_dir} does not exist')
    if not os.path.exists(m_lor_data_dir):
        raise FileNotFoundError(f'{m_lor_data_dir} does not exist')

    save = False
    if save and not os.path.exists(save_dir):
        print(f'Creating {save_dir}')
        os.makedirs(save_dir)

    t_list, f_list = au_nano_sphere_scatter(
        drude_data_dir=drude_data_dir, m_lor_data_dir=m_lor_data_dir)

    if save:
        print(f'Saving to {save_dir}')
        for t, f in zip(t_list, f_list):
            f.savefig(os.path.join(save_dir, f'{t}.png'))
    else:
        plt.show()
