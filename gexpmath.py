import math
from decimal import *
import random
from datetime import datetime
import numpy as np
import pandas as pd
from scipy.stats import chi2
from scipy.optimize import minimize

setcontext(Context(prec=100, rounding=ROUND_HALF_EVEN, Emin=-999999, Emax=999999, capitals=1, clamp=0, flags=[Inexact, FloatOperation, Rounded], traps=[InvalidOperation, DivisionByZero, Overflow]))
D = Decimal
PI = D('3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172')


class GNBDistribution(object):
    """docstring for GNBDistribution"""

    def __init__(self, v, q, theta, max_n=10):
        self.alpha = theta
        self.p = q
        self.u = v
        self.max_n = max_n
        self.probabilities = []
        self.cum_probabilities = []
        for i in range(self.max_n):
            prob = GNB_Probability(i, self.u, self.p, self.alpha)
            self.probabilities.append(prob)
            if i == 0:
                self.cum_probabilities.append(prob)
            else:
                self.cum_probabilities.append(self.cum_probabilities[i - 1] + prob)
            if self.cum_probabilities[i] > D("0.999"):
                self.max_n = i + 1
                break

        self.sample = None
        self.num_of_outliers = None

    def generate_sample(self, sample_size=1000):
        self.sample = [None] * sample_size
        random.seed(datetime.now())
        self.sample_size = sample_size
        self.num_of_outliers = 0
        for i in range(sample_size):
            pr = random.random()
            j = 0
            while j < len(self.cum_probabilities) and self.cum_probabilities[j] < pr:
                j += 1

            self.sample[i] = j
            if j == self.max_n:
                self.num_of_outliers += 1

        print(self.probabilities)
        print("-" * 30)
        print(sum(self.probabilities))
        return self.sample

    def make_hist_dict(self, no_outliers=False, density=False):
        if self.sample is None:
            raise(Exception("No sample in object. Call object.generate_sample(sample_size) before calling this function"))

        hist = np.histogram(self.sample, self.max_n + 1, range=(0, self.max_n + 1))

        x = hist[1][:-1]
        y = hist[0]
        if no_outliers and self.max_n in x:
            x = x[:-1]
            y = y[:-1]

        if density:
            y = np.array(y) / sum(y)

        return {'x': x, 'y': y}

    def chi_square(self, q=0.05, hist=None):
        if hist is None:
            hist = self.make_hist_dict(no_outliers=True)

        sample_bins = hist['y']
        print(sample_bins)

        chi2_stat = 0
        n = sum(sample_bins)
        print(np.array(self.probabilities) * n)
        for i in range(len(sample_bins)):
            p = self.probabilities[i]
            m = sample_bins[i]
            num = (m - p * n) * (m - p * n)
            den = n * p
            chi2_stat += num / den

        df = len(sample_bins)
        crit_val = chi2.ppf(1 - q, df)
        pvalue = 1 - chi2.cdf(np.float64(chi2_stat), df)
        # result = pvalue > q
        result = chi2_stat < crit_val
        return result, chi2_stat, crit_val, pvalue

    def __str__(self):
        return f"GNB(v={self.u}, q={self.p}, \u03B8={self.alpha})"


def poisson(k, lamb):
    if lamb == 0:
        return D(1 / math.factorial(k))
    else:
        return D(math.exp(-lamb)) * (lamb ** k) / math.factorial(k)


def generalized_gamma(x, params):
    u, p, theta = params
    numerator = abs(u) * D(math.exp(-((x / theta) ** u))) * (x ** (u * p - 1))
    denumenator = theta ** (u * p) * D(gamma(p))

    return numerator / denumenator


def GexpElement_k(x, alpha, beta, k):

    return ((x ** k) / D(math.factorial(k))) * D(gamma(alpha * k + beta))


def Gexp(x, alpha, beta):
    cumsum = D(0)
    prevCumsum = D(0)
    i = 0
    current_element = D(0)
    prev_element = D(0)
    while (i < 10000):
        prevCumsum += current_element
        current_element = (GexpElement_k(D(x), D(alpha), D(beta), D(i)))
        cumsum += current_element
        if i > 0 and i % 10 == 0:
            delta = D((i + 1) / i) * (abs(current_element) / abs(prev_element))
            eps = D('0.001')
            if abs(current_element) < eps and delta < 1:
                break
            # if (i % 500 == 0):
            #     print('howdy', i, delta, current_element)
        prev_element = current_element
        i += 1

    return cumsum



def GNB_Probability(n, u, p, alpha):
    n = D(n)
    u = D(u)
    p = D(p)
    alpha = D(alpha)
    if u > 1:
        multiplier = (alpha ** n) / (D(gamma(p)) * D(math.factorial(n)))

        gexp = Gexp(-alpha, 1 / u, n / u + p)
        # print('gexp params: x = ', -alpha, ', alpha = ', 1/u, ', beta = ', n/u + p, sep='')
        # print('multiplier =', multiplier, 'gexp =', gexp)
        return multiplier * gexp
    elif u < 1:
        multiplier = u / (D(alpha ** (u * p)) * D(gamma(p)) * D(math.factorial(n)))
        gexp = Gexp(-(alpha ** (-u)), u, n + u * p)

        return multiplier * gexp
    elif u == 1:
        dr1 = gamma(n + p) / (gamma(p) * math.factorial(n))
        dr2 = (alpha / (alpha + 1)) ** n
        dr3 = 1 / (alpha + 1) ** p

        return dr1 * dr2 * dr3
 

# from cmath import sin, sqrt, pi, exp

_p = [D('676.5203681218851')
    ,D('-1259.1392167224028')
    ,D('771.32342877765313')
    ,D('-176.61502916214059')
    ,D('12.507343278686905')
    ,D('-0.13857109526572012')
    ,D('9.9843695780195716e-6')
    ,D('1.5056327351493116e-7')
    ]

EPSILON = 1e-07


def gamma(z):
    z = D(z)
    if z < D('0.5'):
        y = PI / (decimal_sin(PI*z) * gamma(1-z))  # Reflection formula
    else:
        z -= 1
        x = D('0.99999999999980993')
        for (i, pval) in enumerate(_p):
            x += pval / (z+i+1)
        t = z + len(_p) - D('0.5')
        y = D(2*PI).sqrt() * t**(z+D('0.5')) * D(-t).exp() * x
    return y
# """
# The above use of the reflection (thus the if-else structure) is necessary, even though 
# it may look strange, as it allows to extend the approximation to values of z where 
# Re(z) < 0.5, where the Lanczos method is not valid.
# """


def digamma(x):
    if x > 2:
        return digamma(x - 1) + 1 / x
    elif x < 1:
        return digamma(x + 1) - 1 / x
    else:
        eul_mac = D('0.577215664901532860606512090082')
        res = 0
        for n in range(1, 200):
            res += (x - 1) / (n * n + x * n - n)
        res -= eul_mac
        return res


def decimal_gamma(x):
    """
    Lanczos approximation for gamma function. All 
    """
    if x < 0.5:
        # Reflection formula
        gamma = PI / (decimal_sin(PI * x) * decimal_gamma(1 - x))
    else:
        x -= 1
        g = D('1')
        gamma = D(D('2') * PI).sqrt() * ((x + g + D('0.5')) ** (x + D('0.5'))) * D(-(x + g + D('0.5'))).exp() * A(x, g)
    return gamma

def A(x, g):
    res = 0
    steps = 15
    numerator = 0
    denumenator = 0
    for i in range(steps):
        if i == 0:
            numerator = 1
            denumenator = 2
        elif i == 1: 
            numerator = x
            denumenator = x + 1
        else:
            numerator *= (x - i + 1)
            denumenator *= (x + i)
        
        res += p(g, i) * numerator / denumenator
    return res



def p(g, k):
    coef = D('2').sqrt() / PI
    _sum = 0
    # Gamma(-1/2) = -2sqrt(PI)
    gamma_val = -2 * PI.sqrt()
    for l in range(k + 1):
        _sum += cheb(2 * k + 1, 2 * l + 1) * gamma_val * (l + g + D('0.5')) ** (-l - D('0.5')) * D(l + g + D('0.5')).exp()
        # Gamma(x + 1) = x * Gamma(x)
        gamma_val *= (l - D('0.5'))

    return coef * _sum

def cheb(n, m):
    """
    Recursive formula for Chebyshev coefficients. n >= m always
    """
    
    if n == 1 and m == 1:
        return 1
    elif n == 2 and m == 2:
        return 1
    elif m == 1:
        return -cheb(n - 2, 1)
    elif m == n:
        return 2 * cheb(n - 1, m - 1)
    else:
        return 2 * cheb(n - 1, m - 1) - cheb(n - 2, m)


def decimal_sin(x):
    """
    Maclaurin series for sin(x)
    """
    sin = D('0')
    factorial = D('1')
    steps = 100
    for i in range(steps):
        sign = -1 if i % 2 else 1
        sin += sign * x ** (2 * i + 1) / factorial
        # factorial = (2(i + 1) + 1)
        factorial = factorial * (2 * (i + 1)) * (2 * (i + 1) + 1)
    return sin


def l2optimization(h):
    x0 = init_approx(h)
    # 
    res = minimize(l2metric, x0, args=(len(h), h), method='Nelder-Mead', tol=5 * 1e-6, callback=callbackF)
    # v, q, theta
    return res['x']

Nfeval = 0
def callbackF(Xi):
    global Nfeval
    print('{0:4d}   {1: 3.6f}   {2: 3.6f}   {3: 3.6f}'.format(Nfeval, Xi[0], Xi[1], Xi[2]))
    Nfeval += 1

def l2metric(params, max_n, h):
    l2 = D('0')
    v = params[0]
    q = params[1]
    theta = params[2]
    if len(h) > max_n:
        h = h[:max_n]
    if len(h) < max_n:
        h = h + [0] * (max_n - len(h))

    for i in range(max_n):
        l2 += (D(GNB_Probability(i, v, q, theta)) - D(h[i])) * (D(GNB_Probability(i, v, q, theta)) - D(h[i]))
    l2 = l2.sqrt()
    return l2


def init_approx(x):
    Ex = 0
    for i in range(len(x)):
        Ex += i * x[i]
    Ex /= sum(x)
    Dx = 0
    for i in range(len(x)):
        Dx += i * i * x[i]
    Dx /= sum(x)
    Dx -= Ex * Ex
    if (Dx > Ex):
        theta = 1 / (Dx / Ex - 1)
        q = theta * theta / (1 + theta) * Ex
        init_v = [0.5, 0.9, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
        minv = 0
        minl2 = -1
        probs = np.array(x) / sum(x)
        for v in init_v:
            l2 = l2metric((v, q, theta), len(x), probs)
            if minl2 == -1 or minl2 > l2:
                minl2 = l2
                minv = v
        v = minv
    else:
        minl2 = -1
        minv = 0.4
        minq = 0.4
        mintheta = 0.4
        probs = np.array(x) / sum(x)
        for v in np.arange(0.0, 4, 0.4):
            for p in np.arange(0.0, 4, 0.4):
                for theta in np.arange(0.0, 4, 0.4):
                    l2 = l2metric((v, q, theta), len(x), probs)
                    if minl2 == -1 or minl2 > l2:
                        minl2 = l2
                        minv = v
                        minq = q
                        mintheta = theta
        v = minv
        q = minq
        theta = mintheta

    return (v, q, theta)

def l1v(v, q, theta, N, m):
    res = 0
    for i in range(N):
        num = m[i] * theta ** i
        denum = gamma(q) * GNB_Probability(i, v, q, theta) * math.factorial(i)
        cumsum = 0
        for k in range(10000):
            current_element = (-theta) ** k / math.factorial(k) * (-k - i) / (v * v) * gamma((k + i) / v + q) * digamma((k + i) / v + q)
            cumsum += current_element
            if k > 0 and k % 10 == 0:
                delta = D((k + 1) / k) * (abs(current_element) / abs(prev_element))
                eps = D('0.01')
                if abs(current_element) < eps and delta < 1:
                    break
                
            prev_element = current_element
        res += num / denum * cumsum
        print(i)
    print('l1v')
    return res


def l1q(v, q, theta, N, m):
    res = 0

    for i in range(N):
        num = m[i] * theta ** i
        denum = gamma(q) * GNB_Probability(i, v, q, theta) * math.factorial(i)
        cumsum = 0
        current_element = 0
        for k in range(10000):
            current_element = (GexpElement_k(D(-theta), D(1 / v), D(i / v + q), D(k))) * digamma(k/v + i/v + q)
            cumsum += current_element
            if k > 0 and k % 10 == 0:
                delta = D((k + 1) / k) * (abs(current_element) / abs(prev_element))
                eps = D('0.01')
                if abs(current_element) < eps and delta < 1:
                    break
                
            prev_element = current_element
        res += num / denum * cumsum
        print(i)
    print('l1q')
    return res - digamma(q) * sum(m)


def l1t(v, q, theta, N, m):
    r0 = m[0] / (GNB_Probability(0, v, q, theta) * gamma(q))
    cumsum0 = 0
    for k in range(1, 10000):
        current_element = (-1) ** k * k * theta ** (k - 1) / math.factorial(k) * gamma(k / v + q)
        cumsum0 += current_element

        if k > 0 and k % 10 == 0:
            delta = D((k + 1) / k) * (abs(current_element) / abs(prev_element))
            eps = D('0.01')
            if abs(current_element) < eps and delta < 1:
                break
            
        prev_element = current_element
    r0 = r0 * cumsum0


    r = 0    
    for i in range(1, N):
        mult = m[i] / (GNB_Probability(i, v, q, theta) * gamma(q) * math.factorial(i))
        cumsum = 0 
        for k in range(10000):
            current_element = (-1) ** k * (k + i) * theta ** (k + i - 1) / math.factorial(k) * gamma(k / v + i / v + q)
            cumsum += current_element

            if k > 0 and k % 10 == 0:
                delta = D((k + 1) / k) * (abs(current_element) / abs(prev_element))
                eps = D('0.01')
                if abs(current_element) < eps and delta < 1:
                    break
                
            prev_element = current_element
        r += mult * cumsum
        print(i)
    print('l1t')
    return r0 + r

was_there = False
def l2v(v, q, theta, N, m):
    res = 0
    for i in range(N):
        mult = m[i] / (GNB_Probability(i, v, q, theta) * gamma(q) * math.factorial(i))
        ge = Gexp(-(theta ** (-v)), v, vq + i)
        r1 = 1 - q * v * theta.ln()

        cumsum = 0
        for k in range(10000):
            current_element = (-theta ** (-v)) ** k / math.factorial(k) * (k + q) * gamma(v * k + v * q + i) * digamma(v * k + v * q + i)
            cumsum += current_element

            if k > 0 and k % 10 == 0:
                delta = D((k + 1) / k) * (abs(current_element) / abs(prev_element))
                eps = D('0.01')
                if abs(current_element) < eps and delta < 1:
                    break
                
            prev_element = current_element

        r2 = v * cumsum

        cumsum = 0
        for k in range(1, 10000):
            current_element = (-1) ** (k + 1) * theta ** (-v * k) * k * theta.ln() / math.factorial(k) * gamma(v * k + v * q + i)
            cumsum += current_element

            if k > 0 and k % 10 == 0:
                delta = D((k + 1) / k) * (abs(current_element) / abs(prev_element))
                eps = D('0.01')
                if abs(current_element) < eps and delta < 1:
                    break
                
            prev_element = current_element
        r3 = v * cumsum
        res += mult * (r1 + r2 + r3)
        print(i)
    print('l2v')
    return res

def l2q(v, q, theta, N, m):
    res = 0
    for i in range(N):
        mult = m[i] / (GNB_Probability(i, v, q, theta) * gamma(q) * math.factorial(i))
        cumsum = 0
        for k in range(0, 10000):
            current_element = (-(theta ** -v)) ** k / math.factorial(k) * v * gamma(v * k + v * q + i) * digamma(v * k + v * q + i)
            cumsum += current_element

            if k > 0 and k % 10 == 0:
                delta = D((k + 1) / k) * (abs(current_element) / abs(prev_element))
                eps = D('0.01')
                if abs(current_element) < eps and delta < 1:
                    break

            prev_element = current_element
        res += mult * ((cumsum - (v * theta.ln() + digamma(q)) * Gexp(-theta ** (-v), v, v * q + i)) / theta ** (v * q))
        print(i)
    print('l2q')
    return res

def l2t(v, q, theta, N, m):
    res = 0
    for i in range(N):
        mult = m[i] / (GNB_Probability(i, v, q, theta) * gamma(q) * math.factorial(i))
        cumsum = 0
        for k in range(1, 10000):
            current_element = (-1) ** (k + 1) * theta ** (-v * k - 1) * k * v / math.factorial(k) * gamma(v * k + v * q + i)
            cumsum += current_element

            if k > 0 and k % 10 == 0:
                delta = D((k + 1) / k) * (abs(current_element) / abs(prev_element))
                eps = D('0.01')
                if abs(current_element) < eps and delta < 1:
                    break

            prev_element = current_element
        res += mult * ((cumsum - v * q * theta ** (v * q - 1) * Gexp(-theta ** (-v), v, v * q + i)) / theta ** (v * q))
        print(i)
    print('l2t')
    return res

def Fun(x, args):
    v = x[0]
    q = x[1]
    theta = x[2]
    N = args[0]
    m = args[1]
    if v < D('1.1') and v > D('0.99'):
        v = 0.98
        global was_there
        was_there = True
    if v > 1:
        return (l1v((v, q, theta, N, m)), l1q((v, q, theta, N, m)), l1t((v, q, theta, N, m)))
    else:
        return (l2v((v, q, theta, N, m)), l2q((v, q, theta, N, m)), l2t((v, q, theta, N, m)))

def Jac(x, args):
    v = x[0]
    q = x[1]
    theta = x[2]
    N = args[0]
    m = args[1]
    h = D('0.001')
    if v < D('1.1') and v > D('0.99'):
        v = 0.98
        global was_there
        was_there = True
    if v < 1:
        return [[(l1v((v + h, q, theta, N, m)) - l1v((v - h, q, theta, N, m))) /  2 * h, (l1v((v, q + h, theta, N, m)) - l1v((v, q - h, theta, N, m))) /  2 * h, (l1v((v, q, theta + h, N, m)) - l1v((v, q, theta - h, N, m))) /  2 * h],
                [(l1q((v + h, q, theta, N, m)) - l1q((v - h, q, theta, N, m))) /  2 * h, (l1q((v, q + h, theta, N, m)) - l1q((v, q - h, theta, N, m))) /  2 * h, (l1q((v, q, theta + h, N, m)) - l1q((v, q, theta - h, N, m))) /  2 * h],
                [(l1t((v + h, q, theta, N, m)) - l1t((v - h, q, theta, N, m))) /  2 * h, (l1t((v, q + h, theta, N, m)) - l1t((v, q - h, theta, N, m))) /  2 * h, (l1t((v, q, theta + h, N, m)) - l1t((v, q, theta - h, N, m))) /  2 * h]]
    else:
        return [[(l2v((v + h, q, theta, N, m)) - l2v((v - h, q, theta, N, m))) /  2 * h, (l2v((v, q + h, theta, N, m)) - l2v((v, q - h, theta, N, m))) /  2 * h, (l2v((v, q, theta + h, N, m)) - l2v((v, q, theta - h, N, m))) /  2 * h],
                [(l2q((v + h, q, theta, N, m)) - l2q((v - h, q, theta, N, m))) /  2 * h, (l2q((v, q + h, theta, N, m)) - l2q((v, q - h, theta, N, m))) /  2 * h, (l2q((v, q, theta + h, N, m)) - l2q((v, q, theta - h, N, m))) /  2 * h],
                [(l2t((v + h, q, theta, N, m)) - l2t((v - h, q, theta, N, m))) /  2 * h, (l2t((v, q + h, theta, N, m)) - l2t((v, q - h, theta, N, m))) /  2 * h, (l2t((v, q, theta + h, N, m)) - l2t((v, q, theta - h, N, m))) /  2 * h]]

def ml_gnb(hist):
    x0 = init_approx(hist)
    N = len(hist)
    m = hist

    result = newton_system(Fun, J, x0, eps=D('0.01'), args=(N, m))
    if result[0] < D('1.1') and result[0] > D('0.99') or result[0] == 0.98:
        result[0] = 1
    was_there = False
    return result

def newton_system(F, J, x, eps, args):
    """
    Solve nonlinear system F=0 by Newton's method.
    J is the Jacobian of F. Both F and J must be functions of x.
    At input, x holds the start value. The iteration continues
    until ||F|| < eps.
    """
    F_value = F(x, args)
    F_norm = np.linalg.norm(F_value, ord=2)  # l2 norm of vector
    iteration_counter = 0
    while abs(F_norm) > eps and iteration_counter < 200:
        print()
        print('-' * 30)
        print(x)
        print('-' * 30)
        print()
        print(datetime.now())
        delta = np.linalg.solve(J(x, args), -F_value)
        print(datetime.now())
        x = x + delta
        F_value = F(x, args)
        F_norm = np.linalg.norm(F_value, ord=2)
        iteration_counter += 1
    if iteration_counter == 200:
        return None

    return x
