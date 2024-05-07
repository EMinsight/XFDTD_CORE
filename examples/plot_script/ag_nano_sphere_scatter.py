from helper import my_plot_arg_parser, plot2d_gif
import os
import numpy as np
import matplotlib.pyplot as plt


def au_nano_sphere_scatter(drude_data_dir, m_lor_data_dir):
    t_list = []
    f_list = []
    incident_wave = np.load(os.path.join(drude_data_dir, 'incident_wave.npy'))
    i_time = incident_wave[0, :]
    i_value = incident_wave[1, :]
    
    drude_ex_point = np.load(os.path.join(drude_data_dir, 'ex_point.npy'))
    m_lor_ex_point = np.load(os.path.join(m_lor_data_dir, 'ex_point.npy'))
    
    ex_time = drude_ex_point[0, :]
    drude_ex_value = drude_ex_point[1, :]
    m_lor_ex_value = m_lor_ex_point[1, :]
    
    f, ax = plt.subplots(2, 1, figsize=(8, 10))
    ax[0].plot(i_time*1e12, i_value, label='Incident Wave', linewidth=2)
    ax[0].plot(ex_time*1e12, drude_ex_value, label='Drude Model', linestyle='dashdot', linewidth=2)
    ax[0].plot(ex_time*1e12, m_lor_ex_value, label='Lorentz Model', linestyle=':', linewidth=2)
    
    ax[0].set_title('Ag Nano Sphere Scatter Time Domain')
    ax[0].set_xlabel('Time (fs)')
    ax[0].set_ylabel('Amplitude')
    ax[0].set_xlim([i_time[0]*1e12, i_time[-1]*1e12])
    ax[0].grid()
    ax[0].legend()
    
    # fft
    # extend i_value to the length is 20000
    # i_value = np.pad(i_value, (0, 20000 - len(i_value)), constant_values=0)
    # drude_ex_value = np.pad(drude_ex_value, (0, 20000 - len(drude_ex_value)), constant_values=0)
    # m_lor_ex_value = np.pad(m_lor_ex_value, (0, 20000 - len(m_lor_ex_value)), constant_values=0)
    
    i_fft = np.fft.fftshift(np.fft.fft(i_value))
    drude_ex_fft = np.fft.fftshift(np.fft.fft(drude_ex_value))
    m_lor_ex_fft = np.fft.fftshift(np.fft.fft(m_lor_ex_value))
    freq = np.fft.fftshift(np.fft.fftfreq(len(i_value), i_time[1] - i_time[0]))
    print(f"Resolution: {(freq[1] - freq[0])/1e12} THz")
    freq_min = 0
    freq_max = freq[-1]
    # find index
    freq_min_idx = np.argmin(np.abs(freq - freq_min))
    freq_max_idx = np.argmin(np.abs(freq - freq_max))
    freq = freq[freq_min_idx:freq_max_idx]
    
    i_fft = i_fft[freq_min_idx:freq_max_idx]
    drude_ex_fft = drude_ex_fft[freq_min_idx:freq_max_idx]
    m_lor_ex_fft = m_lor_ex_fft[freq_min_idx:freq_max_idx]
    
    freq_unit = 1e12  # THz
    max_ex = np.max(np.abs(i_fft))
    
    
    ax[1].plot(freq/freq_unit, np.abs(i_fft)/max_ex, label='Incident Wave', linewidth=2)
    ax[1].plot(freq/freq_unit, np.abs(drude_ex_fft)/max_ex, label='Drude Model', linestyle='dashdot', linewidth=2)
    ax[1].plot(freq/freq_unit, np.abs(m_lor_ex_fft)/max_ex, label='Lorentz Model', linestyle=':', linewidth=2)
    ax[1].set_title('Normalized Amplitude Spectrum')
    ax[1].set_xlabel('Frequency (THz)')
    ax[1].set_ylabel('Amplitude')
    ax[1].set_xlim([300, 1200])
    ax[1].set_ylim([0, 1.25])
    # ax[1].set_xticks(np.arange(300, 1200, 50))
    ax[1].set_yticks(np.arange(0, 1.25, 0.25))
    ax[1].grid()
    ax[1].legend()
    t_list.append('incident_wave')
    f_list.append(f)

    # # # plot relative permittivity
    # # relative_permittivity = np.load(os.path.join(
    # #     data_dir, 'relative_permittivity.npy'))
    # # freq = np.real(relative_permittivity[0, :])
    # # relative_permittivity_data = np.real(relative_permittivity[1, :])
    # # f, ax = plt.subplots(2, 1, figsize=(8, 6))
    # # ax[0].plot(freq[10:]/1e12,
    # #            np.real(relative_permittivity_data[10:]), label='Real')
    # # ax[0].plot(freq[10:]/1e12, np.imag(relative_permittivity_data[10:]),
    # #            label='Imaginary')
    # # ax[0].set_xlim([freq[10]/1e12, freq[-1]/1e12])
    # # ax[0].set_title('Relative Permittivity')
    # # ax[0].grid()
    # # ax[0].set_xlabel('Frequency (THz)')
    # # ax[0].legend()
    # # # sigma
    # # sigma = -np.imag(relative_permittivity_data) * 2 * np.pi * freq * 8.85e-12
    # # ax[1].plot(freq[10:]/1e12, sigma[10:], label='Conductivity')
    # # ax[1].set_xlabel('Frequency (THz)')
    # # ax[1].set_ylabel('Conductivity (S/m)')
    # # ax[1].set_xlim([freq[10]/1e12, freq[-1]/1e12])
    # # ax[1].grid()
    # # t_list.append('relative_permittivity')
    # # ax[0].set_title('Relative Permittivity')


    def sorted_files(data_dir):
        files = os.listdir(data_dir)
        data = sorted([os.path.join(data_dir, file) for file in files if file.endswith('.npy')])
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
    # ani.save(os.path.join(data_dir, 'movie_yz.gif'), writer='ffmpeg', fps=5)
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

    t_list, f_list = au_nano_sphere_scatter(drude_data_dir=drude_data_dir, m_lor_data_dir=m_lor_data_dir)

    if save:
        print(f'Saving to {save_dir}')
        for t, f in zip(t_list, f_list):
            f.savefig(os.path.join(save_dir, f'{t}.png'))
    else:
        plt.show()
