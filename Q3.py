#!/usr/bin/python 
# ================
# Author: Kai Chen
# Date: 2021-01-08
# ================

from Q1 import *
import time

if __name__ == '__main__':
    # Mean Square Stability
    mu,sigma = -3, 1
    T = 20
    dt0 = 1.0
    dts = dt0*0.5**np.arange(4)
    n_trials = 50000
    np.random.seed(0)
    traces = []
    t0 = time.time()
    T_max = np.arange(0, T, dts[-1])
    dW = np.random.randn(n_trials, T_max.shape[0])*np.sqrt(dts[-1])
    for i, dt in enumerate(dts):
        t = np.arange(0, T+dt, dt)
        traces_buffer = np.ones((n_trials, len(t)))
        dWi = dW.reshape((n_trials,-1, 2**(len(dts)-1-i))).sum(2)
        for ti in range(len(t)-1):
            traces_buffer[:, ti+1] = EM(mu, sigma, B, Sigma, dSigma, traces_buffer[:,ti], dt, dWi[:, ti])
        traces.append(traces_buffer.copy())
    
    print(f'Evolve SDE took {time.time()-t0:.3f} s')
    
    mean_stability = np.zeros((len(dts)))
    for i in range(mean_stability.shape[0]):
        mean_stability[i] = np.mean(traces[i][:,-1]**2)
    np.save('mean_stability.npy', mean_stability)
    fig, (ax1, ax2) = plt.subplots(1,2,figsize=(10,4))
    ax1.loglog(dts, mean_stability, 'navy', label='mean stability')
    ax1.set_xlabel('Time-Step Size',fontsize=14)
    ax1.set_ylabel(r'$E\left[X_n^2\right]$', fontsize=14)
    ax1.set_title('Mean Stability')

    # Asymptotic Stability
    np.random.seed(0)
    t0 = time.time()
    T = 20
    dts = dt0*0.5**np.arange(4)
    T_max = np.arange(0, T, dts[-1])
    dW = np.random.randn(T_max.shape[0])*np.sqrt(dts[-1])
    for i, dt in enumerate(dts):
        t = np.arange(0, T+dt, dt)
        traces_buffer = np.ones(len(t))
        dWi = dW.reshape((-1, 2**(len(dts)-1-i))).sum(1)
        for ti in range(len(t)-1):
            traces_buffer[ti+1] = EM(mu, sigma, B, Sigma, dSigma, traces_buffer[ti], dt, dWi[ti])
        ax2.semilogy(t, np.abs(traces_buffer), color='navy', alpha=dt/dts.max(), label=f'$\Delta$t = {dt:.3f} s',)
    print(f'Evolve SDE took {time.time()-t0:.3f} s')
    ax2.set_xlabel('Time', fontsize=14)
    ax2.set_ylabel(r'$\left|X_n\right|$', fontsize=14)
    ax2.set_title('Asymptotically Stablility')
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig('stability.png')