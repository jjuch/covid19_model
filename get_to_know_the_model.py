from data_20_07_11.fracDiffDelay import FDD, FDDVarA
import matplotlib.pyplot as plt

L = 10
alpha = [0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.96, 0.98, 0.99]
# alpha = [0.96, 0.98]
# alpha = [0.99, 0.995]

Ts = 0.005
Tmax = 50

plt.figure()
for a in alpha:
    print(a)
    # fdd = FDD(a, L)
    fdd = FDDVarA(a, L)
    t, y = fdd.impulse_response(Ts, Tmax, verbose=True)
    plt.plot(t, y, label=str(a))
    plt.plot(L**a, y[round(L**a / Ts)], 'r+', markersize=4)

plt.legend()
plt.show()
