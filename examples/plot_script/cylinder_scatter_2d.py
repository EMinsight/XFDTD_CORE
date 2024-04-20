from helper import *

if __name__ == "__main__":
    data_dir = "../data/cylinder_scatter_2d"
    # data_dir = "../data/microstrip_low_pass/movie"
    save_dir = "./cylinder_scatter_2d"
    args = my_plot_arg_parser(
        default_data_dir=data_dir, default_save_dir='./', default_save=False)
    data_dir = args.data_dir
    save_dir = args.save_dir

    if not os.path.exists(data_dir):
        raise FileNotFoundError(f'{data_dir} does not exist')

    if args.save and not os.path.exists(save_dir):
        print(f'Creating {save_dir}')
        os.makedirs(save_dir)

    def get_sorted_files(dir: str):
        files = os.listdir(dir)
        files = [f for f in files if f.endswith('.npy')]
        files.sort()
        return [os.path.join(dir, f) for f in files]

    def read_data(file: str): return np.load(file)

    def process_data(data):
        data = np.squeeze(data)
        data[np.abs(data) < 1e-4] = 1e-4
        data = 10 * np.log10(np.abs(data))
        return data

    files = get_sorted_files(data_dir)
    ani = plot2d_gif(files, read_data, process_data, -40, 10)
    if args.save:
        print(f"Saving to {os.path.join(args.save_dir, 'ex_xz.gif')}")
        ani.save(os.path.join(args.save_dir, "ex_xz.gif"), writer='ffmpeg')
    else:
        plt.show()
