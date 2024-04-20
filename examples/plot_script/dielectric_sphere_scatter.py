import numpy as np
import matplotlib.pyplot as plt
import os
from helper import *


Mie = [
    0.46682631277049674,
    0.46659445551718964,
    0.46589940308337674,
    0.46474271192906497,
    0.46312697110159523,
    0.46105579468596547,
    0.4585338112935272, 
    0.45556665064639923, 0.45216092733002616, 0.44832422180045917, 0.4440650587460003, 0.43939288291459794, 0.4343180325287672, 0.42885171041865777, 0.42300595301112154, 0.41679359731822346, 0.4102282460725267, 0.4033242311586347, 0.39609657549101196, 0.3885609534869538, 0.3807336502809274, 0.37263151982238635, 0.3642719419937568, 0.3556727788787115, 0.3468523303032827, 0.33782928876401386, 0.3286226938483601, 0.31925188624321993, 0.3097364614179217, 0.3000962230585164, 0.290351136320995, 0.2805212809622884, 0.27062680439984477, 0.26068787474335225, 0.2507246338360325, 0.24075715033795445, 0.23080537288020841, 0.2208890833166056, 0.2110278500989372, 0.20124098180279457, 0.19154748083352566, 0.1819659973461243, 0.17251478341863757, 0.16321164752599773, 0.15407390936995174, 0.14511835513082108, 0.1363611932180693, 0.1278180106088668, 0.11950372987686878, 0.1114325670269996, 0.10361799026597232, 0.09607267985227774, 0.088808489183221, 0.08183640728998196, 0.07516652292436067, 0.06880799043257196, 0.06276899762191346, 0.05705673583507474, 0.0516773724540636, 0.046636026060945926, 0.04193674448563349, 0.037582485971623905, 0.033575103688729506, 0.029915333817310156, 0.0266027874212347, 0.023635946316684735, 0.02101216313093404, 0.018727665729406064, 0.01677756617065708, 0.015155874327533842, 0.013855516288726277, 0.012868357628414753, 0.012185231602890166, 0.011795972302103526, 0.011689452751334319, 0.011853627923815502, 0.01227558258951908, 0.01294158388870697, 0.013837138481622333, 0.014947054088191922, 0.016255505194192204, 0.017746102663372655, 0.01940196695890443, 0.021205804642604315, 0.02313998778704308, 0.02518663590424562, 0.027327699965573906, 0.029545048060886123, 0.03182055222149742, 0.03413617591111838, 0.036474061672075384, 0.038816618400957306, 0.04114660771858095, 0.04344722889399109, 0.045702201781234596, 0.04789584723096428, 0.05001316444658468, 0.05203990476666059, 0.05396264137164211, 0.05576883443353906, 0.057446891251902006, 0.05898622094818214, 0.060377283323057834, 0.06161163151741549, 0.06268194815709172, 0.06358207470393622, 0.06430703378091732, 0.06485304428652561, 0.06521752916325202, 0.06539911573605425, 0.06539762858905078, 0.0652140750017941, 0.06485062301992937, 0.06431057228842062, 0.06359831782836771, 0.06271930699033258, 0.06167998986759157, 0.06048776350142561, 0.059150910257043034, 0.05767853079260778, 0.05608047208475665, 0.05436725101157895, 0.0525499740279915, 0.05064025349847188, 0.048650121277961624, 0.04659194015319121, 0.0444783137735217, 0.04232199571249317, 0.04013579830850357, 0.03793250193533908, 0.0357247653506094, 0.03352503776250695, 0.03134547324276424, 0.029197848096308345, 0.02709348177603389, 0.025043161904500837, 0.02305707393340884, 0.021144735936643212, 0.019314938993795723, 0.017575693578635546, 0.015934182321369537, 0.01439671946504122, 0.012968717285444797, 0.011654659690877965, 0.010458083163328967, 0.00938156514672192, 0.008426719931059992, 0.007594202024151614, 0.006883716945528182, 0.006294039320601216, 0.00582303809750669, 0.005467708654873951, 0.00522421151635962, 0.005087917337606842, 0.005053457783720463, 0.0051147818707537, 0.005265217303431957, 0.0054975363037093215, 0.005804025391057881, 0.006176558545887372, 0.006606673162410721, 0.007085648176804148, 0.007604583740814004, 0.00815448180015853, 0.008726326931241363, 0.009311166788880961, 0.009900191521968419, 0.010484811523162942, 0.011056732892846127, 0.011608030016473031, 0.012131214678030832, 0.01261930116036241, 0.013065866815412702, 0.013465107623761305, 0.013811888302830744, 0.014101786566597886, 0.014331131186148765, 0.014497033549644359, 0.014597412471822958, 0.01463101205666085
]


def dielectric_sphere_scatter(data_dir: str, save_dir: str, save: bool = False):
    x = input("plot move? (y/n)")
    if x == 'y':
        print(f"Direcotry: {data_dir}")

        def get_sorted_files(dir: str):
            files = os.listdir(dir)
            files = [f for f in files if f.endswith('.npy')]
            files.sort()
            return [os.path.join(dir, f) for f in files]

        def read_data(file: str): return np.load(file)

        def process_data(data):
            d = np.squeeze(data)
            d[(np.abs(d) < 1e-4)] = 1e-4
            d = 10 * np.log10(np.abs(d))
            return d

        ani_ex_xz = plot2d_gif(get_sorted_files(os.path.join(
            data_dir, "movie_ex_xz")), read_data, process_data)
        if save:
            ani_ex_xz.save(os.path.join(save_dir, "ex_xz.gif"))
        plt.show()

    time = np.load(os.path.join(data_dir, "time.npy"))
    incident_wave = np.load(os.path.join(data_dir, "incident_wave.npy"))
    freq = 1e9

    a = dft(time, incident_wave, np.array([freq]))
    power_i = 0.5 * (np.abs(a) ** 2) / (120 * np.pi)
    xy_data = FarFieldData(
        freq=freq,
        source_power=power_i,
        data_dir=os.path.join(data_dir, "xy"),
        angle=np.linspace(-np.pi, np.pi, 360)
    )

    f, ax = plot_directivity(xy_data, const_theta=True)
    ax.set_title('XY Plane')
    if save:
        f.savefig(os.path.join(save_dir, "xy_directivity.png"))

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
    if save:
        f.savefig(os.path.join(save_dir, "xz_rcs.png"))

    yz_data = FarFieldData(
        freq=freq,
        source_power=power_i,
        data_dir=os.path.join(data_dir, "yz"),
        angle=np.linspace(-np.pi, np.pi, 360)
    )

    _, ax = plot_rcs_dBsm(yz_data, const_angle_degree=90,
                          const_angle_symbol=r'\varphi', plot_theta=False)
    ax.set_title('YZ Plane')
    ax.set_xlim([0, np.pi])
    if save:
        f.savefig(os.path.join(save_dir, "yz_rcs.png"))

    # compare with mie theory
    mie =np.array(Mie)
    print(f"mie.shape: ${mie.shape}")
    mie = 10 * np.log10(mie)

    f = plt.figure()
    plt.plot(np.linspace(0, 180, 180), mie, label="Mie")
    plt.plot(np.linspace(-180, 180, 360),
             xz_data.data_theta_dB,  label="FDTD")
    plt.title('RCS in XZ Plane (Mie vs FDTD)')
    plt.xlabel('Angle (degree)')
    plt.ylabel('RCS (dBsm)')
    plt.xlim([0, 180])
    plt.ylim([-28, 0])
    plt.grid(True)
    plt.legend()
    plt.title('RCS in XZ Plane')
    if save:
        f.savefig(os.path.join(save_dir, "xz_rcs_compare.png"))

    err = np.abs(mie - xz_data.data_theta_dB[-180:])
    f = plt.figure()
    plt.plot(np.linspace(0, 180, 180), err)
    plt.title('Error between Mie and FDTD')
    plt.xlabel('Angle (degree)')
    plt.ylabel('Error (dB)')
    plt.grid(True)
    print('mean error = ', np.mean(err))
    if save:
        f.savefig(os.path.join(save_dir, "xz_rcs_error.png"))


if __name__ == "__main__":
    data_dir = r"../data/dielectric_sphere_scatter/"
    args = my_plot_arg_parser(
        default_data_dir=data_dir, default_save_dir='dielectric_sphere_scatter', default_save=False)

    save = args.save
    save_dir = args.save_dir
    data_dir = args.data_dir

    if not os.path.exists(data_dir):
        raise FileNotFoundError(f"Data directory {data_dir} not found")
    if not os.path.exists(save_dir) and save:
        print(f"Creating save directory {save_dir}")
        os.makedirs(save_dir)

    dielectric_sphere_scatter(data_dir=data_dir, save_dir=save_dir, save=save)
    if save:
        full_path = os.path.realpath(save_dir)
        print(f"Results are saved in {full_path}")
    else:
        plt.show()
