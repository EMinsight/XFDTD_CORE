import numpy as np
import matplotlib.pyplot as plt
import os
from helper import my_plot_arg_parser

def simple_capacitor(data_dir:str):
    v_data_path = os.path.join(data_dir, 'v_monitor.npy')
    v_data = np.load(v_data_path)
    t_data = v_data[0,:]
    v_data = v_data[1,:]
    exact_v_data_path = os.path.join(data_dir, 'exact.npy')
    f, [ax1,ax2] = plt.subplots(2, 1)
    f.set_size_inches(9, 9)
    ax1.plot(t_data, v_data, label='Simulation')
    exact_data = np.load(exact_v_data_path)
    ax1.plot(t_data, exact_data, label='Exact')
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Voltage (V)')
    ax1.set_title('Voltage')
    ax1.grid()
    ax1.legend()
    ax2.plot(t_data, v_data-exact_data)
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('Error (V)')
    ax2.set_title('Error')
    ax2.grid()
    return ['Voltage'], [f]
    


if __name__ == "__main__":
    args = my_plot_arg_parser(default_data_dir='../data/simple_capacitor', default_save=False, default_save_dir='./simple_capacitor')
    title_list, fig_list = simple_capacitor(args.data_dir)
    if not os.path.exists(args.data_dir):
        raise FileNotFoundError(f'{args.data_dir} does not exist')
    if args.save and not os.path.exists(args.save_dir):
        print(f'Creating {args.save_dir}')
        os.makedirs(args.save_dir)
    
    if args.save:
        for i, fig in enumerate(fig_list):
            fig.savefig(os.path.join(args.save_dir, f'{title_list[i]}.png'))
    else:
        plt.show()
