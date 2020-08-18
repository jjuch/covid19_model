############################################################################   Fractional Diffusive Delay:                                           #
#                                                                         #
#       Calculate the impulse response of the Fractional Diffusive Delay  #
###########################################################################

import numpy as np
import scipy.special as sc
import math
from mpmath import mp, mpf, isinf, fsum, fmul, fdiv, power, gamma, factorial, rgamma
import mpmath
import matplotlib.pyplot as plt
import warnings
import sys
import control


class FDD:
    def __init__(self, alpha, L):
        if isinstance(L, float) or isinstance(L, int):
            if L > 0:
                self.L = L
            else:
                raise AssertionError('Parameter L should be larger than 0.')
        else:
            raise TypeError('The parameter L should be an int or a float.')
        if isinstance(alpha, float):
            if 0 < alpha < 1:
                self.alpha = alpha
            else:
                raise AssertionError('The parameter alpha should be between 0 and 1.')
        elif self.__iterable__(alpha):
            type_not_respected = [not isinstance(a, float) for a in alpha]
            print(type_not_respected)
            if any(type_not_respected):
                raise TypeError('The parameter alpha should be a list of floats')
            else:
                bounds_respected = [0 < a < 1 for a in alpha]
                if all(bounds_respected):
                    self.alpha = alpha
                else:
                    raise AssertionError('The parameters alpha should be between 0 and 1.')
        else:
            raise TypeError('The parameter alpha should be a float or a list of floats')

        
    def __iterable__(self, obj):
        try:
            iter(obj)
        except Exception:
            return False
        else:
            return True


    def impulse_response(self, Ts:float, Tmax:float, K=1.0, N:int=200, P=20, plot=False, verbose=False):
        '''
        This function returns FDL's exp(-(L*s)^alpha) time response for the time vector t=[Ts:Ts:Tmax].
        
        Parameters:
        -----------
            : L: the dead time [s]
            : alpha: the fractional power ]0,1[
            : Ts: sampling time [s]
            : Tmax: the end value of the time vector [s]
            : N: Quality of the output (The number of summation terms)
            : P: Estimated peak size

        Returns:
        --------
            : t: time vector [s]
            : I: impulse response vector
        '''
        # Produce time axis: the time power will converge as t > 1
        t = np.arange(Ts, Tmax, Ts)

        if isinstance(self.alpha, float):
            alpha = self.alpha
        else:
            alpha = self.alpha[0]

        # Init parameters of loop
        i = 0
        summed_terms = len(t) * [0]

        while i < N:
            # Gamma-function of integer argument is infinite
            arg = i * alpha
            if arg.is_integer():
                single_term = len(t) * [0]
            else:
                ti = time.time()
                try:
                    # Calculate the different terms of the single_term
                    factorial_i = math.factorial(i)
                    gamma_ia = math.gamma(-i * alpha)
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        time_power = [t_el**(-i * alpha - 1) for t_el in t]
                    L_ai = (-1)**i * self.L**(i*alpha)
                    const = L_ai/(factorial_i * gamma_ia)
                    if math.isinf(const):
                        raise OverflowError
                    single_term = [tp * const for tp in time_power]

                except OverflowError:
                    # Check if there are any overflows
                    mp.dps = 65
                    mp.prec = 100
                    factorial_i_hp = factorial(i)
                    gamma_ia_hp = gamma(-i * alpha)
                    with warnings.catch_warnings():
                        warnings.filterwarnings('error')
                        try:
                            L_ai_hp = mpf((-1)**i * self.L**(i * alpha))
                        except:
                            if verbose:
                                print('Overflow stop at {}%, due to large L'.format(math.ceil(i/N*100)))
                            break
                    const_hp = mpf(L_ai_hp/(factorial_i_hp * gamma_ia_hp))
                    single_term = [tp * float(const_hp) for tp in time_power]


                    # Check if time power is not infinite at t < 1
                    q = [i for i, val in enumerate(time_power) if math.isinf(val)]
                    
                    for j in range(len(q)):
                        time_power_temp = power(t[q[j]], (-alpha * i - 1))
                        quotient_hp = time_power_temp * const_hp
                        if math.isinf(float(quotient_hp)):
                            single_term[q[j]] = sys.float_info.max
                        else:
                            single_term[q[j]] = float(quotient_hp)

            # Add the current iteration to the previous iterations
            with warnings.catch_warnings():
                warnings.filterwarnings('error')
                try:
                    summed_terms = [sum_el + single_el if not math.isinf(sum_el + single_el) else sys.float_info.min for sum_el, single_el in zip(summed_terms, single_term)]                   
                except:
                    summed_terms = [fsum([sum_el, single_el]) if not isinf(fsum([sum_el, single_el])) else sys.float_info.min for sum_el, single_el in zip(summed_terms, single_term)]

            # Update iteration counter
            i += 1

            if verbose and (i % 5 == 0):
                # Print progress
                print('Progress: {:d}%'.format(math.floor(i/N * 100)), end='\r', flush=True)
        if verbose:
            print('', end='\n', flush=True)

        # Polish off the error due to cutting of the summation
        q1 = [i for i, val in enumerate(summed_terms) if abs(val) > P]
        q2 = [i for i, val in enumerate(summed_terms) if val < 0]
        
        I_norm = summed_terms
        if len(q1) is not 0:
            #TODO: improve cleaning, no magic number 10, look at derivative?
            if q1[-1] + 10 < len(I_norm):
                I_norm[0:(q1[-1] + 10)] = (q1[-1] + 10) * [0]
        if len(q2) is not 0:
            I_norm[0:q2[-1] + 1] = (q2[-1] + 1) * [0]
        
        I = [K * I_norm_el for I_norm_el in I_norm]

        if plot:
            plt.figure()
            plt.plot(t, I)
            plt.show()
        return t, I

    def step_response(self, Ts:float, Tmax:float, K=1.0, N:int=200, P=10**3, plot=False, verbose=False):
        t, I = self.impulse_response(Ts, Tmax, K=K, N=N, P=P, verbose=verbose)
        I_sum = np.cumsum(I)
        if plot:
            plt.figure()
            plt.plot(t, I_sum)
            plt.show()
        return t, I_sum



class FDDVarA:
    def __init__(self, alpha, L):
        if isinstance(L, float) or isinstance(L, int):
            if L > 0:
                self.L = L
            else:
                raise AssertionError('Parameter L should be larger than 0.')
        else:
            raise TypeError('The parameter L should be an int or a float.')
        if isinstance(alpha, float):
            if 0 < alpha < 1:
                self.alpha = alpha
            else:
                raise AssertionError('The parameter alpha should be between 0 and 1.')
        elif self.__iterable__(alpha):
            type_not_respected = [not isinstance(a, float) for a in alpha]
            if any(type_not_respected):
                raise TypeError('The parameter alpha should be a list of floats')
            else:
                bounds_respected = [0 < a < 1 for a in alpha]
                if all(bounds_respected):
                    self.alpha = alpha
                else:
                    raise AssertionError('The parameters alpha should be between 0 and 1.')
        else:
            raise TypeError('The parameter alpha should be a float or a list of floats')

        
    def __iterable__(self, obj):
        try:
            iter(obj)
        except Exception:
            return False
        else:
            return True
    

    def __calculate_single_term_var_alpha__(self, t, i, alpha):
        i_alpha = [i * a for a in alpha]
        sign = (-1) ** i
        i_fact = math.factorial(i)

        def catch(func, func_hp, *args, exceptions=Exception, filter_name='error',  **kwargs):
                with warnings.catch_warnings():
                    warnings.filterwarnings(filter_name)
                    try:
                        res = func(*args, **kwargs)
                        if res is None:
                            raise exceptions
                        else:
                            return res
                    except exceptions:
                        return func_hp(*args, **kwargs)

        def L_ai_func(ia):
            try:
                i_fact_recip = i_fact**(-1)
                L_i_alpha = self.L**ia
                return sign * L_i_alpha * i_fact_recip
            except:
                return 

        def L_ai_hp_func(ia):
            return fdiv(sign * power(self.L, ia), i_fact)

        L_ai = [0 if ia.is_integer() else catch(L_ai_func, L_ai_hp_func, ia,exceptions=(OverflowError, Exception)) for ia in i_alpha]
        rec_gamma = [1 if ia.is_integer() else rgamma(-ia) if math.isinf(sc.rgamma(-ia)) else sc.rgamma(-ia) for ia in i_alpha]
        const = [Lai * rec_g for Lai, rec_g in zip(L_ai, rec_gamma)]

        def single_term_func(ti, ia, c):
            if c < sys.float_info.max:
                return ti**(-ia - 1) * c
            else:
                return fmul(ti**(-ia - 1), c)
                
        
        def single_term_hp_func(ti, ia, c):
            return fmul(power(ti, (-ia - 1)), c)

        single_term = [0 if ia.is_integer() else catch(single_term_func, single_term_hp_func, time, ia, c) for time, ia, c in zip(t, i_alpha, const)]
        return single_term


    def __calculate_single_term_cst_alpha__(self, t, i, alpha):
        i_alpha = i * alpha
        if i_alpha.is_integer():
            return len(t) * [0]
        else:
            sign = (-1) ** i
            i_fact = math.factorial(i)

            def catch(func, func_hp, *args, exceptions=Exception, filter_name='error',  **kwargs):
                with warnings.catch_warnings():
                    warnings.filterwarnings(filter_name)
                    try:
                        res = func(*args, **kwargs)
                        if res is None:
                            raise exceptions
                        else:
                            return res
                    except exceptions:
                        return func_hp(*args, **kwargs)
            
            def L_ai_func():
                try:
                    i_fact_recip = i_fact**(-1)
                    L_i_alpha = self.L**i_alpha
                    return sign * L_i_alpha * i_fact_recip
                except:
                    return 

            def L_ai_hp_func():
                return fdiv(sign * power(self.L, i_alpha), i_fact)
            
            L_ai = catch(L_ai_func, L_ai_hp_func, exceptions=OverflowError)
            rec_gamma = rgamma(-i_alpha) if math.isinf(sc.rgamma(-i_alpha)) else sc.rgamma(-i_alpha)
            const = L_ai * rec_gamma
            

            def single_term_func(ti):
                if const < sys.float_info.max:
                    return ti**(-i_alpha - 1) * const
                else:
                    return fmul(ti**(-i_alpha - 1), const)
                    
            
            def single_term_hp_func(ti):
                return fmul(power(ti, (-i_alpha - 1)), const)

            single_term = [catch(single_term_func, single_term_hp_func, time) for time in t]
            
            return single_term



    def impulse_response(self, Ts:float, Tmax:float, K=1.0, N:int=200, P=20, plot=False, verbose=False):
        '''
        This function returns FDL's exp(-(L*s)^alpha) time response for the time vector t=[Ts:Ts:Tmax].
        
        Parameters:
        -----------
            : L: the dead time [s]
            : alpha: the fractional power ]0,1[
            : Ts: sampling time [s]
            : Tmax: the end value of the time vector [s]
            : N: Quality of the output (The number of summation terms)
            : P: Estimated peak size

        Returns:
        --------
            : t: time vector [s]
            : I: impulse response vector
        '''
        # Produce time axis: the time power will converge as t > 1
        t = np.arange(Ts, Tmax, Ts)

        # Prepare alpha vector
        diff_start_idx = [0]
        diff_stop_idx = []
        if isinstance(self.alpha, float):
            alpha = len(t) * [self.alpha]
        else:
            alpha = self.alpha
            diff_alpha = [1 if not math.isclose(a_diff, 0.0) else 0 for a_diff in np.diff(alpha, prepend=alpha[0])]
            
            add_idx_bool = True
            for idx, el in enumerate(diff_alpha):
                if not math.isclose(el, 0.0):
                    if (not add_idx_bool) and (idx + 1 is not len(alpha)) and (math.isclose(diff_alpha[idx + 1], 0)):
                        add_idx_bool = True 
                    if add_idx_bool:
                        diff_stop_idx.append(idx)
                        diff_start_idx.append(idx)
                        add_idx_bool = False
                elif not add_idx_bool:
                    add_idx_bool = True
        diff_stop_idx.append(len(alpha))

        type_alpha = [] # 0: cst| 1: var
        for k in range(len(diff_start_idx)):
            if sum(diff_alpha[diff_start_idx[k]: diff_stop_idx[k]]) == diff_stop_idx[k] - diff_start_idx[k]:
                type_alpha.append(1)
            else:
                type_alpha.append(0)

        # Init parameters of loop
        i = 0
        summed_terms = len(t) * [0]

        # Set precisions
        mp.dps = 65
        mp.prec = 100

        while i < N:
            # Init one term vector
            single_term = len(t) * [0]

            for k in range(len(diff_start_idx)):
                if type_alpha[k]:
                    single_term[diff_start_idx[k]:diff_stop_idx[k]] = self.__calculate_single_term_var_alpha__(t[diff_start_idx[k]:diff_stop_idx[k]], i, alpha[diff_start_idx[k]:diff_stop_idx[k]])
                    
                else:
                    single_term[diff_start_idx[k]:diff_stop_idx[k]] = self.__calculate_single_term_cst_alpha__(t[diff_start_idx[k]:diff_stop_idx[k]], i, alpha[diff_start_idx[k]])
                
            # Add the current iteration to the previous iterations
            with warnings.catch_warnings():
                warnings.filterwarnings('error')
                try:
                    summed_terms = [sum_el + single_el if not math.isinf(sum_el + single_el) else sys.float_info.min for sum_el, single_el in zip(summed_terms, single_term)]                   
                except:
                    summed_terms = [fsum([sum_el, single_el]) if not isinf(fsum([sum_el, single_el])) else sys.float_info.min for sum_el, single_el in zip(summed_terms, single_term)]


            # Update iteration counter
            i += 1

            if verbose and (i % 5 == 0):
                # Print progress
                print('Progress: {:d}%'.format(math.floor(i/N * 100)), end='\r', flush=True)
        if verbose:
            print('', end='\n', flush=True)

        # Polish off the error due to cutting of the summation
        q1 = [i for i, val in enumerate(summed_terms) if abs(val) > P]
        q2 = [i for i, val in enumerate(summed_terms) if val < 0]
        
        I_norm = summed_terms
        if len(q1) is not 0:
            #TODO: improve cleaning, no magic number 10, look at derivative?
            if q1[-1] + 10 < len(I_norm):
                I_norm[0:(q1[-1] + 10)] = (q1[-1] + 10) * [0]
        if len(q2) is not 0:
            I_norm[0:q2[-1] + 1] = (q2[-1] + 1) * [0]
        
        I = [K * I_norm_el for I_norm_el in I_norm]

        if plot:
            plt.figure()
            plt.plot(t, I)
            plt.show()
        return t, I

    def step_response(self, Ts:float, Tmax:float, K=1.0, N:int=200, P=10**3, plot=False, verbose=False):
        t, I = self.impulse_response(Ts, Tmax, K=K, N=N, P=P, verbose=verbose)
        I_sum = np.cumsum(I)
        if plot:
            plt.figure()
            plt.plot(t, I_sum)
            plt.show()
        return t, I_sum



class FOPFDD(FDD):
    def __init__(self, K, tau, alpha, L):
        super().__init__(alpha, L)
        self.K = K
        self.tau = tau

    def step_response(self, Ts:float, Tmax:float, N:int=200, P=20, plot=False, verbose=False):
        t_fdd, I_fdd = super().impulse_response(Ts, Tmax, N=N, P=P, verbose=verbose)
        sys_fo = control.tf(self.K, [self.tau, 1])
        t, y_fo = control.step_response(sys_fo, T=t_fdd)
        y_fofdd_full = np.convolve(y_fo, I_fdd)
        y_fofdd = y_fofdd_full[0:int(np.floor((len(y_fofdd_full) + 1)/2))]
        y = [Ts * y_el for y_el in y_fofdd]
        if plot:
            plt.figure()
            plt.plot(t, y, label='fopfdd')
            plt.show()
        return t, y


if __name__ == '__main__':
    import time

    # t = time.time()
    # a_var = 2000 * [0.5] + list(np.linspace(0.5, 0.6, 24)) + 1975 * [0.6]
    # f = FDDVarA(a_var, 1)
    # f.impulse_response(0.0005, 2, plot=False, verbose=True)
    # # f.step_response(0.0005, 2, plot=True, verbose=True)
    # elapsed = time.time() - t
    # print('toc: ', elapsed)

    t = time.time()
    a_var = 238 * [0.85] + list(np.linspace(0.85, 0.95, 24)) + 238 * [0.95]
    # a_var = 500 * [0.95]
    f = FDDVarA(a_var, 100)
    # f.impulse_response(0.0005, 2, plot=True, verbose=True)
    # f.impulse_response(0.5, 250, plot=False, verbose=True)
    f.step_response(0.5, 250, plot=True, verbose=True)
    elapsed = time.time() - t
    # print('toc: ', elapsed)

    t = time.time()
    f2 = FDD(0.75, 100)
    # f2.impulse_response(0.0005, 2, plot=True, verbose=True)
    # f2.impulse_response(0.05, 250, plot=False, verbose=True)
    elapsed = time.time() - t
    # print('toc: ', elapsed)

    # fdd1 = FDD(0.5, 1)
    # fdd1.impulse_response(0.0005, 2, plot=True, verbose=True)
    # fdd2 = FDD(0.85, 1)
    # fdd2.cumulative_impulse_response(0.005, 2, plot=True, verbose=True)

    # fdd1 = FDD(0.5, 1)
    # t1, I1 = fdd1.impulse_response(0.0005, 2, N=200, verbose=True)
    # fdd2 = FDD(0.85, 1)
    # t2, I2 = fdd2.impulse_response(0.0005, 2, N=200, verbose=True)
    # fdd3 = FDD(0.95, 1)
    # t3, I3 = fdd3.impulse_response(0.0005, 2, N=600, P=4.4, verbose=True)
    # plt.figure()
    # plt.plot(t1, I1, t2, I2, t3, I3)
    # plt.show()

    # fopfdd = FOPFDD(1, 1, 0.5, 1)
    # fopfdd.step_response(0.0005, 2, verbose=True, plot=True)

    # fdd4 = FOPFDD(1.373242801772693, 1.6472210441695725, 0.7997688363521038, 850.3613823603612)
    # t, y = fdd4.step_response(1, 695.0)
    # print(y[220:700])