import numpy as np

def gerk(a, b, c, initial, terminal, timesteps, func, enforce_rules=False):

    try:
        class GerkError(Exception):
            __module__ = "builtins"
        assert type(timesteps) is int and timesteps > 0
        assert len(b) == len(c) == len(a) + 1 
        if enforce_rules:
            assert sum(b) == 1
            assert round(sum([bi * ci for bi, ci in zip(b, c)]), 1) == 0.5
            for i in range(len(c)):
                assert sum(a[i]) == c[i]
    except AssertionError:
        raise GerkError("A condition has failed. Please check your Runge-Kutta coefficients.")
    
    k = np.zeros(len(b))
    h = (terminal - initial[0]) / timesteps
    x_n, y_n = initial
    x, y = [x_n], [y_n]
    x_k, y_k = x_n, y_n

    for _ in range(timesteps):
        k[0] = func(x_k, y_k)
        for i in range(1, len(a) + 1):
            x_k += h * c[i]
            y_k += h * sum([a[i - 1][j] * k[j] for j in range(len(a[i - 1]))])
            k[i] = func(x_k, y_k)

        x_n += h
        y_n += h * sum(b * k)
        x.append(x_n)
        y.append(y_n)
        x_k, y_k = x_n, y_n

    return x, y

def adaptive_gerk(a, b_1, b_2, c, initial, terminal, timesteps, func, enforce_rules=False, tolerance=1e-4):

    try:
        class GerkError(Exception):
            __module__ = "builtins"
        assert type(timesteps) is int and timesteps > 0
        assert len(b_1) == len(b_2) == len(c) == len(a) + 1
        if enforce_rules:
            assert sum(b_1) == 1
            assert round(sum([b1i * ci for b1i, ci in zip(b_1, c)]), 1) == 0.5
            assert round(sum([b2i * ci for b2i, ci in zip(b_2, c)]), 1) == 0.5
            for i in range(len(c)):
                assert sum(a[i]) == c[i]
    except AssertionError:
        raise GerkError("A condition has failed. Please check your Runge-Kutta coefficients.")
    
    k = np.zeros(len(b_1))
    arr_b_1 = np.array(b_1)
    arr_b_2 = np.array(b_2)
    h = (terminal - initial[0]) / timesteps
    order = min([x for x in b_1 if x], [x for x in b_2 if x])
    x_n, y_n = initial
    x, y = [x_n], [y_n]
    x_k, y_k = x_n, y_n

    for _ in range(timesteps):
        k[0] = func(x_k, y_k)
        for i in range(1, len(a) + 1):
            x_k += h * c[i]
            y_k += h * sum([a[i - 1][j] * k[j] for j in range(len(a[i - 1]))])
            k[i] = func(x_k, y_k)

        if h * np.dot((arr_b_1 - arr_b_2), k) > tolerance:
            h *= 0.9 ** (1 / (order + 1))
            continue

        x_n += h
        y_n += h * sum(b_1 * k)
        x.append(x_n)
        y.append(y_n)
        x_k, y_k = x_n, y_n

    return x, y
