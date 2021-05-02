## **Maksimchik Daniil** 

- daniil.mm@icloud.com
- +375295517635
- vk: https://vk.com/makcimchik_daniil_com


|   Skills      |    Level       |
|:-------------:|:---------------|
| C, C++           | two semester of study at the university           |
| Git           | semester of study at the university             |
| Java, python, scala | semester of study at the university|
| Assembler   | semester of study at the university|

### Education
The third year of study at bsuir at the faculty of computer systems and networks.

| Language      | Level           | 
|:------------- |:---------------:| 
| English       |  B1 Streamline  | 
| Russian       |      Native     |

### Code examples

from math import cos as cos
from math import sin as sin
from math import exp
from cmath import pi as pi
import numpy as np

precision = 2


def function(x):
    y = sin(x) + cos(4 * x)
    return y


def discrete_x_scale(N):
    x = []
    for i in range(N):
        x.append(2 * pi / N * i)
    return x


def discrete_y_values(function_, x_scale):
    y = []
    N = len(x_scale)
    for n in range(N):
        y.append(round(function_(x_scale[n]), 5))
    return y


def rounding(value):
    if abs(value) < pow(10, -1 * precision):
        return 0
    else:
        return round(value, precision)


def fir_filter(fft_complex_values, N):
    M = 800
    H = np.zeros(M)

    A = np.zeros(N)
    B = np.zeros(N)

    A = calculation_rationing(A, M, 0.196)
    B = calculation_rationing(B, M, 0.204)
    B = conversion(B, M)

    for i in range(M):
        H[i] = A[i] + B[i]

    H = conversion(H, M)

    fir = []
    for i in range(N):
        fir.append(complex(0.0, 0.0j))

    j = M
    while j < N - 1:
        fir[j] = complex(0.0, 0.0j)
        for i in range(M):
            fir[j] = fir[j] + fft_complex_values[j - i] * H[i]
        j = j + 1
    return fir


def calculation_rationing(mus, M, FC):
    for i in range(M):
        if (i == M / 2):
            mus[i] = 2 * pi * FC
        else:
            mus[i] = sin(2 * pi * FC * (i - M / 2)) / (i - M / 2)
        mus[i] = mus[i] * (0.42 - 0.5 * cos(2 * pi * i / M) + 0.08 * cos(4 * pi * i / M))

    sum = 0
    for i in range(M):
        sum = sum + mus[i]

    for i in range(M):
        mus[i] = mus[i] / sum

    return mus


def conversion(mus, M):
    for i in range(M):
        mus[i] = -1 * mus[i]
    mus[400] = mus[400] + 1
    return mus


def iir_filter(discrete_values, N, fc):
    noise_discrete_values = np.zeros(N)
    for i in range(N):
        noise_discrete_values[i] = discrete_values[i] + np.random.uniform(-0.5, 0.5)

    x = exp(-14.445 * fc)
    a0 = (1 - x) ** 4
    b1 = 4 * x
    b2 = -6 * x ** 2
    b3 = 4 * x ** 3
    b4 = -1 * x ** 4

    iir = np.zeros(N)

    j = 4
    while j < N:
        iir[j] = a0 * noise_discrete_values[j] + b1 * iir[j - 1] + b2 * iir[j - 2] + b3 * iir[j - 3] + b4 * iir[j - 4]
        j = j + 1

    return noise_discrete_values, iir


def fastFourierTransform(a, N, direction):
    if N == 1:
        return a

    wN = complex(cos(2 * pi / N), direction * sin(2 * pi / N))
    w = complex(1.0, 0.0j)

    firstPart = []
    secoundPart = []
    for n in range(N // 2):
        firstPart.append(a[n])
    for n in range(N // 2, N, 1):
        secoundPart.append(a[n])

    y = []
    transformResult1 = []
    transformResult2 = []

    for j in range(N // 2):
        firstResult = complex(rounding((a[j] + a[j + N // 2]).real), rounding((a[j] + a[j + N // 2]).imag))
        secoundResult = complex(rounding((a[j] - a[j + N // 2]).real), rounding((a[j] - a[j + N // 2]).imag)) * w
        transformResult1.append(firstResult)
        transformResult2.append(secoundResult)
        w = w * wN

    transformResult1 = fastFourierTransform(transformResult1, N // 2, direction)
    transformResult2 = fastFourierTransform(transformResult2, N // 2, direction)

    for n in range(N // 2):
        y.insert(2 * n, transformResult1[n])
        y.insert(2 * n + 1, transformResult2[n])

    return y
