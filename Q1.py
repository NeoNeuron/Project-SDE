#!/usr/bin/python 
# ================
# Author: Kai Chen
# Date: 2021-01-08
# ================

import numpy as np
import matplotlib.pyplot as plt

def B(mu, X):
    return mu*X

def Sigma(sigma, X):
    return sigma*X

def dSigma(sigma, X):
    return sigma

def EM(mu, sigma, B, Sigma, dSigma, X, dt, dW):
        return X + B(mu, X)*dt + Sigma(sigma, X)*dW

def Milstein(mu, sigma, B, Sigma, dSigma, X, dt, dW):
        return X + B(mu, X)*dt + Sigma(sigma, X)*dW + 0.5*Sigma(sigma, X)*dSigma(sigma, X)*(dW**2-dt)


def RK(mu, sigma, B, Sigma, dSigma, X, dt, dW):
    X1 = X + Sigma(sigma, X)*np.sqrt(dt)
    return X + B(mu, X)*dt + Sigma(sigma, X)*dW + 0.5*(Sigma(sigma, X1)-Sigma(sigma, X))*(dW**2-dt)/np.sqrt(dt)

if __name__ == '__main__':
    mu,sigma = -0.5, 1
    T = 10
    dt = 0.1
    t = np.arange(0, T, dt)
    traces = np.ones((3,len(t)))
    labels = ['Euler-Maruyama Scheme', 'Milstein Scheme', 'Runge-Kutta Scheme']
    colors = ['navy', 'orange', 'springgreen']
    np.random.seed(0)
    for idx, f in enumerate([EM, Milstein, RK]):
        for ti in range(len(t)-1):
            dW = np.random.randn()*np.sqrt(dt)
            traces[idx, ti+1] = f(mu, sigma, B, Sigma, dSigma, traces[idx,ti], dt, dW)
        
    lines = plt.plot(t, traces.T)
    [line.set_color(colors[i]) for i, line in enumerate(lines)]
    plt.legend(lines, labels)
    plt.xlabel('Time', fontsize=16)
    plt.ylabel('Value', fontsize=16)
    plt.tight_layout()
    plt.savefig('numerical_trace.png')