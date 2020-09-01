import matplotlib.pyplot as plt
import numpy as np
import math
import datetime

from scipy import optimize
from pyswarms.single.global_best import GlobalBestPSO
from pyswarms.utils.plotters import plot_cost_history

from fracDiffDelay import FDD, FOPFDD, FDDVarA
from data_vis import read_data



def convert_to_int_list(data):
    return [float(el) for el in data.tolist()]

def convert_csv_raw_data(data):
    data_new = np.array(data)

def load_and_process_data(fileName, plot=False):
    data_list = read_data(fileName)
    data = np.array(data_list)
    time = data[:, 0]
    cum_data = data[:, 1]
    t = convert_to_int_list(time)
    cum_da = convert_to_int_list(cum_data)
    if plot:
        plt.figure()
        plt.plot(t, cum_da)
        plt.show()
    return t, cum_da

def calulate_error(z, *params):
    K, tau, alpha, L = z
    ref_data, Ts, Tmax = params
    model = FOPFDD(K, tau, alpha, L)
    t, y = model.step_response(Ts, Tmax, verbose=False)
    kwad_difference_per_sample = [(r - y_el)**2 for r, y_el in zip(ref_data, y)]
    print('r: ', len(ref_data), ' - y: ', len(y), ' - test: ', kwad_difference_per_sample[-3:])
    print(K, ', ', tau, ', ', alpha, ', ', L, ' : ', sum(kwad_difference_per_sample))
    return sum(kwad_difference_per_sample)


def calculate_error_PSO(x, ref_data, Ts, Tmax):
    K_list = x[:, 0]
    tau_list = x[:, 1]
    alpha_list = x[:, 2]
    L_list = x[:, 3]
    res = []
    for i in range(len(L_list)):
        K, tau, alpha, L = K_list[i], tau_list[i], alpha_list[i], L_list[i]
        if 61**(1/alpha)*0.8 < L < 61**(1/alpha)*1.2:
            print('It is in the range')
            model = FOPFDD(K, tau, alpha, L)
            t, y = model.step_response(Ts, Tmax, verbose=False)
            weight = np.ones(int(Tmax))
            # weight = [1/t_el for t_el in t]
            kwad_difference_per_sample = [(r - y_el)**2 * w for r, y_el, w in zip(ref_data, y, weight)]
            # print(K, ', ', tau, ', ', alpha, ', ', L, ' : ', sum(kwad_difference_per_sample))
            res.append(sum(kwad_difference_per_sample))
        else:
            res.append(1000)
    return np.array(res)

def calculate_error_PSO_VarA(x, ref_data, Ts, Tmax):
    K_list = x[:, 0]
    tau_list = x[:, 1]
    alpha1_list = x[:, 2]
    alpha2_list = x[:, 3]
    alpha3_list = x[:, 4]
    L_list = x[:, 5]
    res = []
    for i in range(len(L_list)):
        K, tau, alpha1, alpha2, alpha3, L = K_list[i], tau_list[i], alpha1_list[i], alpha2_list[i], alpha3_list[i], L_list[i]
        # if (61**(1/alpha1)*0.8 < L < 61**(1/alpha1)*1.2) and (61**(1/alpha2)*0.8 < L < 61**(1/alpha2)*1.2) and (61**(1/alpha3)*0.8 < L < 61**(1/alpha3)*1.2):
        if True:
            print('It is in the range')
            alpha = create_list_alpha(alpha1, alpha2, alpha3 , dates_converted)
            model = FOPFDD(K, tau, alpha, L)
            t, y = model.step_response(Ts, Tmax, verbose=False)
            weight = np.ones(int(Tmax))
            # weight = [1/t_el for t_el in t]
            kwad_difference_per_sample = [(r - y_el)**2 * w for r, y_el, w in zip(ref_data, y, weight)]
            # print(K, ', ', tau, ', ', alpha, ', ', L, ' : ', sum(kwad_difference_per_sample))
            res.append(sum(kwad_difference_per_sample))
        else:
            res.append(1000)
    return np.array(res)

def create_list_alpha(a1, a2, a3, dates):
    alpha = dates[0] * [a1] + list(np.linspace(a1, a2, 14)) + (dates[1] - dates[0] - 14) * [a2] + list(np.linspace(a2, a3, 14)) + (dates[2] - dates[1] - 14) * [a3]
    return alpha




if __name__ == '__main__':
    file_name = 'cum_cases_flanders.csv'
    t_fl, data_fl = load_and_process_data(file_name, plot=False)

    ################## FDD model - manual
    # fdd = FDDVarA(0.40, 3000)
    # t_mod, data_mod = fdd.step_response(1, 1000, K=1.2, plot=False, verbose=True)

    # plt.figure()
    # plt.plot(t_fl, data_fl, label='data')
    # plt.plot(t_mod, data_mod, label='model')
    # plt.legend()
    # plt.xlabel('Time [days]')
    # plt.ylabel('Cumulative cases')
    # plt.title('Flanders')
    # plt.show()
    # exit()


    ################## FOPFDD model - manual
    # fopfdd = FOPFDD(1.2, 40, 0.77, 120)
    # z = ( 1.28991665e+00, 6.53533072e+01, 4.75404563e-01, 1.00729913e+04)
    # params = (data_fl, 1, t_fl[-1]) 
    # fopfdd = FOPFDD(*z)
    # t_mod2, data_mod2 = fopfdd.step_response(1, t_fl[-1], verbose=True)
    # res = calulate_error(z, *params)
    
    # plt.figure()
    # plt.plot(t_fl, data_fl, label='data')
    # plt.plot(t_mod2, data_mod2, label='model')
    # plt.legend()
    # plt.xlabel('Time [days]')
    # plt.ylabel('Cumulative cases')
    # plt.title('Flanders')
    # plt.show()
    # exit()

    ################# FOPFDD model - differential evolution
    # bounds = [(1.1, 1.3), (30, 50), (0.7, 0.8), (200, 300)]
    # params = (data_fl, 1, t_fl[-1])
    # res = optimize.differential_evolution(calulate_error, bounds, args=params, maxiter=150)
    # print(res.x)
    # print(res.message)
    # print(res.jac)

    ################ FOPFDD model - PSO
    # # Create bounds
    # K_min, K_max = 1, 1.5
    # tau_min, tau_max = 1, 100
    # alpha_min, alpha_max = 0.40, 0.95
    # L_min, L_max = 61**(1/0.95)*0.8, 61**(1/0.40)*1.2
    # bounds = (np.array([K_min, tau_min, alpha_min, L_min]), np.array([K_max, tau_max, alpha_max, L_max]))

    # # Initialize swarm
    # options = {'c1': 0.5, 'c2': 0.3, 'w': 0.9}
    # kwargs = {"ref_data": data_fl , "Ts": 1 , "Tmax": t_fl[-1]}
    # optimizer = GlobalBestPSO(n_particles=10, dimensions=4, options=options, bounds=bounds)
    # cost, pos = optimizer.optimize(calculate_error_PSO, iters=50, **kwargs)
    # plot_cost_history(cost_history=optimizer.cost_history)
    # plt.show()


    ############# Time-Varying FOPFDD model - PSO
    # Important dates National Security Board
    start_date = datetime.date(2020, 1, 24) # https://ec.europa.eu/info/live-work-travel-eu/health/coronavirus-response/timeline-eu-action_en
    dates = [datetime.date(2020, 3, 18), datetime.date(2020, 5, 18)]
    dates_converted = [(d - start_date).days for d in dates]
    print(dict(zip(dates, dates_converted)))
    dates_converted.append(len(t_fl))

    # Create bounds
    K_min, K_max = 1, 1.5
    tau_min, tau_max = 1, 100
    alpha1_min, alpha1_max = 0.50, 0.75
    alpha2_min, alpha2_max = 0.65, 0.95
    alpha3_min, alpha3_max = 0.50, 0.75
    # L_min, L_max = 61**(1/0.95)*0.8, 61**(1/0.40)*1.2
    L_min, L_max = 61**(1/0.95)*0.8, 750
    bounds = (np.array([K_min, tau_min, alpha1_min, alpha2_min, alpha3_min, L_min]), np.array([K_max, tau_max, alpha1_max, alpha2_max, alpha3_max, L_max]))
    print(L_min, L_max)

    # Initialize swarm
    options = {'c1': 0.5, 'c2': 0.3, 'w': 0.9}
    kwargs = {"ref_data": data_fl , "Ts": 1 , "Tmax": t_fl[-1]}
    optimizer = GlobalBestPSO(n_particles=5, dimensions=6, options=options, bounds=bounds)
    cost, pos = optimizer.optimize(calculate_error_PSO_VarA, iters=50, **kwargs)
    plot_cost_history(cost_history=optimizer.cost_history)
    plt.show()

    alpha = create_list_alpha(pos[2], pos[3], pos[4], dates_converted)
    params = (data_fl, 1, t_fl[-1]) 
    fopfdd = FOPFDD(pos[0], pos[1], alpha, pos[5])
    t_mod2, data_mod2 = fopfdd.step_response(1, t_fl[-1], verbose=True)
    res = calculate_error_PSO_VarA(np.array([pos]), *params)
    
    plt.figure()
    plt.plot(t_fl, data_fl, label='data')
    plt.plot(t_mod2, data_mod2, label='model')
    plt.legend()
    plt.xlabel('Time [days]')
    plt.ylabel('Cumulative cases')
    plt.title('Flanders')
    plt.show()
    exit()
