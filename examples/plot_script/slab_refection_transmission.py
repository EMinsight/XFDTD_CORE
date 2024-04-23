import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.animation import FuncAnimation
from helper import my_plot_arg_parser, dft


def refection_coefficient(epsilon_1, epsilon_2):
    return np.abs((np.sqrt(epsilon_1) - np.sqrt(epsilon_2))/(np.sqrt(epsilon_1) + np.sqrt(epsilon_2)))


def transmission_coefficient(epsilon_1, epsilon_2):
    return np.abs(2*np.sqrt(epsilon_1)/(np.sqrt(epsilon_1) + np.sqrt(epsilon_2)))


if __name__ == "__main__":
    args = my_plot_arg_parser(default_data_dir='../data/slab_reflection_transmission',
                              default_save=False, default_save_dir='./slab_reflection_transmission')

    if not os.path.exists(args.data_dir):
        raise FileNotFoundError(f'{args.data_dir} does not exist')
    if args.save and not os.path.exists(args.save_dir):
        print(f'Creating {args.save_dir}')
        os.makedirs(args.save_dir)

    epsilon_data = np.load(os.path.join(args.data_dir, 'epsilon.npy'))
    freq = np.abs(epsilon_data[0, :])
    epsilon = epsilon_data[1, :]
    real_data = np.real(epsilon)
    imag_data = np.imag(epsilon)
    sigma = -imag_data * 2 * np.pi * freq * 8.85e-12

    f, ax = plt.subplots(2, 1, figsize=(6, 10))
    ax[0].plot(freq/1e9, real_data, label='Real')
    ax[0].plot(freq/1e9, imag_data, label='Imaginary')
    ax[0].set_xlim([freq[0]/1e9, freq[-1]/1e9])
    ax[0].set_title('Relative Permittivity')
    ax[0].set_xlabel('Frequency (GHz)')
    ax[0].set_ylabel('Relative permittivity')
    ax[0].grid()
    ax[0].legend()

    ax[1].plot(freq/1e9, sigma, label='Conductivity')
    ax[1].set_xlim([freq[0]/1e9, freq[-1]/1e9])
    ax[1].set_xlabel('Frequency (GHz)')
    ax[1].set_ylabel('Conductivity (S/m)')
    ax[1].set_title('Conductivity')
    ax[1].grid()

    reflect = np.load(os.path.join(args.data_dir, 'reflect_ex.npy'))
    transmit = np.load(os.path.join(args.data_dir, 'transmit_ex.npy'))
    incident = np.load(os.path.join(args.data_dir, 'incident.npy'))

    f, ax = plt.subplots(2, 1, figsize=(10, 10))
    ax[0].plot(reflect[0, :], reflect[1, :], label='Reflect')
    ax[0].plot(transmit[0, :], transmit[1, :], label='Transmit')
    ax[0].plot(incident[0, :], incident[1, :], label='Incident')
    ax[0].set_xlim([0, reflect[0, -1]])
    ax[0].set_title('Reflection and Transmission')
    ax[0].legend()
    ax[0].grid()
    ax[0].set_xlabel('Time (s)')

    r_fft = np.fft.fftshift(np.fft.fft(reflect[1, :]))
    t_fft = np.fft.fftshift(np.fft.fft(transmit[1, :]))
    i_fft = np.fft.fftshift(np.fft.fft(incident[1, :]))
    time = reflect[0, :]
    freq_fft = np.fft.fftshift(np.fft.fftfreq(len(time), time[1] - time[0]))

    freq_min = 0
    freq_max = freq_fft[-1]
    # find index
    freq_min_idx = np.argmin(np.abs(freq_fft - freq_min))
    freq_max_idx = np.argmin(np.abs(freq_fft - freq_max))
    freq_fft = freq_fft[freq_min_idx:freq_max_idx]
    r_fft = r_fft[freq_min_idx:freq_max_idx]
    t_fft = t_fft[freq_min_idx:freq_max_idx]
    i_fft = i_fft[freq_min_idx:freq_max_idx]

    ax[1].plot(freq_fft/1e9, np.abs(r_fft), label='Reflect')
    ax[1].plot(freq_fft/1e9, np.abs(t_fft), label='Transmit')
    ax[1].plot(freq_fft/1e9, np.abs(i_fft), label='Incident')
    ax[1].set_xlim([freq[0]/1e9, freq[-1]/1e9])
    ax[1].set_title('Reflection and Transmission')
    ax[1].legend()
    ax[1].grid()
    ax[1].set_xlabel('Frequency (GHz)')

    f, ax = plt.subplots(2, 1, figsize=(10, 10))
    theory_r = refection_coefficient(epsilon_1=1, epsilon_2=epsilon)
    theory_t = transmission_coefficient(epsilon_1=1, epsilon_2=epsilon)
    ax[0].plot(freq_fft/1e9, np.abs(r_fft/i_fft),
               label='Reflect', linestyle='-.')  # dotted line
    ax[0].plot(freq_fft/1e9, np.abs(t_fft/i_fft),
               label='Transmit', linestyle=':')  # circle line
    ax[0].plot(freq/1e9, np.abs(theory_r), label='Theory Reflect')
    ax[0].plot(freq/1e9, np.abs(theory_t), label='Theory Transmit')
    ax[0].set_xlim([freq[0]/1e9, freq[-1]/1e9/5])
    ax[0].set_ylim([0, 1])
    ax[0].grid()
    ax[0].set_xlabel('Frequency (GHz)')
    ax[0].legend()
    ax[0].set_title('Reflection and Transmission Coefficient')

    interpolate_freq = np.linspace(freq[0], freq[-1], 1000)
    theory_r_abs = np.abs(theory_r)
    theory_t_abs = np.abs(theory_t)
    interp_theory_r = np.interp(interpolate_freq, freq, theory_r_abs)
    interp_theory_t = np.interp(interpolate_freq, freq, theory_t_abs)
    interp_r = np.interp(interpolate_freq, freq_fft, np.abs(r_fft/i_fft))
    interp_t = np.interp(interpolate_freq, freq_fft, np.abs(t_fft/i_fft))
    err_rel_t = np.abs((np.abs(interp_t) - np.abs(interp_theory_t)) /
                       np.abs(interp_theory_t))
    err_rel_r = np.abs((np.abs(interp_r) - np.abs(interp_theory_r)) /
                       np.abs(interp_theory_r))

    ax[1].plot(interpolate_freq/1e9, err_rel_r, label='Reflect')
    ax[1].plot(interpolate_freq/1e9, err_rel_t, label='Transmit')
    ax[1].set_xlim([freq[0]/1e9, freq[-1]/1e9/5])
    ax[1].grid()
    ax[1].set_xlabel('Frequency (GHz)')
    ax[1].legend()
    ax[1].set_title('Relative Error')

    def get_sorted_files(dir: str):
        files = os.listdir(dir)
        files = [f for f in files if f.endswith('.npy')]
        files.sort()
        return [os.path.join(dir, f) for f in files]

    files = get_sorted_files(os.path.join(args.data_dir, 'ex_movie'))

    def read_data(file_name: str):
        data = np.load(file_name)
        data = np.squeeze(data)
        # print(data.shape)
        return data

    f, ax = plt.subplots()
    ax.plot(read_data(files[0]))
    ax.set_ylim(-2, 2)
    ax.add_patch(plt.Rectangle((50, -2), 150, 4, fill=None, edgecolor='red'))

    def update_plot(i):
        data = read_data(files[i])
        ax.clear()
        ax.plot(data)
        ax.set_ylim(-2, 2)
        ax.add_patch(plt.Rectangle(
            (50, -2), 150, 4, fill=None, edgecolor='red'))
        ax.set_title(files[i].split('/')[-1])
        return ax

    ani = FuncAnimation(f, update_plot, frames=len(files),
                        interval=100, repeat=True)

    plt.show()
