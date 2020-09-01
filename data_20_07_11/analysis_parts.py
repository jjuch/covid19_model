import numpy as np
import datetime
import matplotlib.pyplot as plt

from scipy import optimize
from pyswarms.single.global_best import GlobalBestPSO
from pyswarms.utils.plotters import plot_cost_history

from analysis import load_and_process_data
from fracDiffDelay import FOPFDD

def calculate_error_PSO(x, ref_data, Ts, Tmax):
    K_list = x[:, 0]
    tau_list = x[:, 1]
    alpha_list = x[:, 2]
    L_list = x[:, 3]
    res = []
    for i in range(len(L_list)):
        K, tau, alpha, L = K_list[i], tau_list[i], alpha_list[i], L_list[i]
        
        model = FOPFDD(K, tau, alpha, L)
        t, y = model.step_response(Ts, Tmax, verbose=False)
        weight = np.ones(int(Tmax))
        weight = [1/t_el for t_el in t]
        kwad_difference_per_sample = [(r - y_el)**2 * w for r, y_el, w in zip(ref_data, y, weight)]
        # print(K, ', ', tau, ', ', alpha, ', ', L, ' : ', sum(kwad_difference_per_sample))
        res.append(sum(kwad_difference_per_sample))
    return np.array(res)


if __name__ == "__main__":
    file_name = 'cum_cases_flanders.csv'
    t_fl, data_fl = load_and_process_data(file_name, plot=False)

    # Important dates National Security Board
    start_date = datetime.date(2020, 1, 24) # https://ec.europa.eu/info/live-work-travel-eu/health/coronavirus-response/timeline-eu-action_en
    dates = [datetime.date(2020, 3, 18), datetime.date(2020, 5, 18)]
    dates.insert(0, start_date)
    dates_converted = [(d - start_date).days for d in dates]
    print(dict(zip(dates, dates_converted)))
    dates_converted.append(len(t_fl))

    t_cut = [t_fl[dates_converted[i]:dates_converted[i + 1]] for i in range(len(dates_converted[1:]))]
    data_cut = [data_fl[dates_converted[i]:dates_converted[i + 1]] for i in range(len(dates_converted[1:]))]

    # plt.figure()
    # plt.plot(t_fl, data_fl, linewidth=4, label='original')
    # for part in range(len(t_cut)):
    #     plt.plot(t_cut[part], data_cut[part], label='part {}'.format(part + 1))
    # plt.legend()
    # plt.show()

    ####### Part 1
    K = 1.2
    tau = 26
    alpha = 0.8
    L = 95
    fopfdd1 = FOPFDD(K, tau, alpha, L)
    t_1, y_1 = fopfdd1.step_response(1, t_cut[0][-1], verbose=True)

    if False:
        plt.figure()
        plt.plot(t_cut[0], data_cut[0], label='data')
        plt.plot(t_1, y_1, label='model')
        plt.legend()
        plt.xlabel('Time [days]')
        plt.ylabel('Cumulative cases')
        plt.title('Flanders')
        plt.show()


    ################ FOPFDD model - PSO
    # # Create bounds
    # K_min, K_max = 1, 1.5
    # tau_min, tau_max = 1, 100
    # alpha_min, alpha_max = 0.75, 0.85
    # L_min, L_max = 50, 150
    # bounds = (np.array([K_min, tau_min, alpha_min, L_min]), np.array([K_max, tau_max, alpha_max, L_max]))

    # # Initialize swarm
    # options = {'c1': 0.5, 'c2': 0.3, 'w': 0.9}
    # kwargs = {"ref_data": data_cut[0] , "Ts": 1 , "Tmax": t_cut[0][-1]}
    # optimizer = GlobalBestPSO(n_particles=10, dimensions=4, options=options, bounds=bounds)
    # cost, pos = optimizer.optimize(calculate_error_PSO, iters=50, **kwargs)
    # plot_cost_history(cost_history=optimizer.cost_history)
    # plt.show()

    # pos = np.array([1.2, 38,  0.81572044, 90.25755211])
    # fopfdd1_opt = FOPFDD(*pos.tolist())
    # t1_opt, data1_opt = fopfdd1_opt.step_response(1, t_cut[0][-1], verbose=True)
    # plt.figure()
    # plt.plot(t_cut[0], data_cut[0], label='data')
    # plt.plot(t1_opt, data1_opt, label='model')
    # plt.legend()
    # plt.xlabel('Time [days]')
    # plt.ylabel('Cumulative cases')
    # plt.title('Flanders')
    # plt.show()


################ FOPFDD model - PSO
    # # Create bounds
    # K_min, K_max = 1, 1.5
    # tau_min, tau_max = 1, 100
    # alpha_min, alpha_max = 0.75, 0.85
    # L_min, L_max = 50, 150
    # bounds = (np.array([K_min, tau_min, alpha_min, L_min]), np.array([K_max, tau_max, alpha_max, L_max]))

    # # Initialize swarm
    # options = {'c1': 0.5, 'c2': 0.3, 'w': 0.9}
    # kwargs = {"ref_data": data_cut[0] , "Ts": 1 , "Tmax": t_cut[0][-1]}
    # optimizer = GlobalBestPSO(n_particles=10, dimensions=4, options=options, bounds=bounds)
    # cost, pos = optimizer.optimize(calculate_error_PSO, iters=50, **kwargs)
    # plot_cost_history(cost_history=optimizer.cost_history)
    # plt.show()

    # pos = np.array([1.2, 38,  0.81572044, 90.25755211])
    # fopfdd1_opt = FOPFDD(*pos.tolist())
    # t1_opt, data1_opt = fopfdd1_opt.step_response(1, t_cut[0][-1], verbose=True)
    # plt.figure()
    # plt.plot(t_cut[0], data_cut[0], label='data')
    # plt.plot(t1_opt, data1_opt, label='model')
    # plt.legend()
    # plt.xlabel('Time [days]')
    # plt.ylabel('Cumulative cases')
    # plt.title('Flanders')
    # plt.show()



####### Part 2
    K = 6
    tau = 38
    alpha = 0.974
    L = 45
    fopfdd2 = FOPFDD(K, tau, alpha, L)
    t_2, y_2 = fopfdd2.step_response(1, 200, verbose=True)

    if True:
        plt.figure()
        plt.plot(t_cut[1], data_cut[1], label='data')
        plt.plot(t_2, y_2, label='model')
        plt.legend()
        plt.xlabel('Time [days]')
        plt.ylabel('Cumulative cases')
        plt.title('Flanders')
        plt.show()

