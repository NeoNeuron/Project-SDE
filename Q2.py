#!/usr/bin/python 
# ================
# Author: Kai Chen
# Date: 2021-01-08
# ================

from Q1 import *
import time

if __name__ == '__main__':
    mu,sigma = -0.5, 1
    T = 1
    dt0 = 0.5
    dts = dt0*0.5**np.arange(10)
    n_trials = 10000
    np.random.seed(0)
    traces = []
    t0 = time.time()
    T_max = np.arange(0, T, dts[-1])
    dW = np.random.randn(n_trials, T_max.shape[0])*np.sqrt(dts[-1])
    for i, dt in enumerate(dts):
        t = np.arange(0, T+dt, dt)
        traces_buffer = np.ones((3,n_trials, len(t)))
        dWi = dW.reshape((n_trials,-1, 2**(len(dts)-1-i))).sum(2)
        for idx, f in enumerate([EM, Milstein, RK]):
            for ti in range(len(t)-1):
                traces_buffer[idx, :, ti+1] = f(mu, sigma, B, Sigma, dSigma, traces_buffer[idx,:,ti], dt, dWi[:, ti])
        traces.append(traces_buffer.copy())
    
    print(f'Evolve SDE took {time.time()-t0:.3f} s')
    
    convergence = np.zeros((3, 2, len(dts)-1))
    for i in range(convergence.shape[2]):
        convergence[:,0,i] = np.mean(np.abs(traces[i]-traces[-1][:,:,::2**(len(dts)-1-i)]), axis=1).max(1)
        convergence[:,1,i] = np.abs(traces[i].mean(1)-traces[-1][:,:,::2**(len(dts)-1-i)].mean(1)).max(1)
    np.save('convergence_data.npy', convergence)

    labels = ['Strong', 'Weak']
    titles = ['Euler-Maruyama Scheme', 'Milstein Scheme', 'Runge-Kutta Scheme']
    colors = ['navy', 'orange', 'springgreen']
    fig, axes = plt.subplots(1,3,figsize=(14,4))
    for idx, ax in enumerate(axes):
        lines = ax.loglog(dts[:-1], convergence[idx,:,:].T, alpha=1.0)
        [line.set_color(colors[i]) for i, line in enumerate(lines)]
        # plot first order reference line
        line_half = ax.loglog(dts[:-1], 5e-1*dts[:-1]**0.5, colors[0], ls='--')[0]
        line_one = ax.loglog(dts[:-1], 1e-1*dts[:-1]**1, colors[1], ls='--')[0]
        # ax.legend(lines, labels)
        ax.legend([*lines, line_half, line_one], [*labels, r'O($\Delta t^{1/2}$)', r'O($\Delta t)$'])
        ax.set_xlabel('Time-Step Size', fontsize=16)
        ax.set_ylabel(f'Error at $X_t$={T:.1f}', fontsize=16)
        ax.set_title(titles[idx], fontsize=18)
        ax.grid(ls='--')
    plt.tight_layout()
    plt.savefig('numerical_convergence.png')