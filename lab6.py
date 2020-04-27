import matplotlib.pyplot as plt
import numpy as np

def func(x):
	return x[0] ** 2 / 3 + x[0] * x[1] / 6 - 127 * x[0] / 6 + x[1] ** 2 / 3 - 91 * x[1] / 6 + 239 / 2

def argmin3(x1, x2, x3):
	if (func(x1) < func(x2)):
		if (func(x1) < func(x3)):
			return x1
		else:
			return x3
	else:
		if (func(x2) < func(x3)):
			return x2
		else:
			return x3

def argmin2(x1, x2):
	if (func(x1) < func(x2)):
		return x1
	else:
		return x2
    
def show(path):
    plt.title("Searching trajectory")
    plt.plot(path[:, 0], path[:, 1], 'r')
    levels = [-225, -100, 20]
    for i in range(len(path)):
        levels.append(func(path[i]))
    levels.sort()
    
    x = np.arange(5, 40, 0.1)
    y = np.arange(-10, 40, 0.1)
    X, Y = np.meshgrid(x, y)
    Z = func(np.array([X, Y], float))
    plt.contour(X, Y, Z, levels)
    
    plt.show()

def slidingWindow(x, p, h):
    x1 = x.copy()
    x2 = x.copy()
    interval = np.zeros(2, dtype = float)

    if (p[0] != 0):
        while (func(x1) < func(x) or func(x) > func(x2)):
            x1 = x.copy()
            x2 = x.copy()
            x1[0] -= h
            x2[0] += h
            x1[1] -= p[1] / p[0] * h
            x2[1] += p[1] / p[0] * h
            interval[0] = x1[0]
            interval[1] = x2[0]
            
            if (func(x1) > func(x2)):
                x[0] += h / 2
                x[1] += p[1] / p[0] * h / 2
            elif (func(x1) < func(x2)):
                x[0] -= h / 2.
                x[1] -= p[1] / p[0] * h / 2
    else:
        while (func(x1) < func(x) or func(x) > func(x2)):
            x1 = x.copy()
            x2 = x.copy()
            x1[1] -= h
            x2[1] += h
            interval[0] = x1[1]
            interval[1] = x2[1]
            
            if (func(x1) > func(x2)):
                x[1] += h / 2
            elif (func(x1) < func(x2)):
                x[1] -= h / 2
                
    return interval

def squareInterp(x, p, h0, eps):
    h = 2
    x1 = x.copy()
    x2 = x.copy()
    x3 = x.copy()
    X = x.copy()
    interval = slidingWindow(x, p, h0);
    xmin = np.zeros(2, dtype = float)

    if (p[0] != 0):
        x1[0] = (interval[0] + interval[1]) / 2
        x1[1] += p[1] / p[0] * (x1[0] - x[0])
        while (np.fabs(func(xmin) - func(X)) > eps and np.fabs(xmin[0] - X[0]) > eps):
            x2[0] = x1[0] + h
            x2[1] += p[1] / p[0] * (x2[0] - x[0])
            if (func(x1) > func(x2)):
                x3[0] = x1[0] + 2 * h
            else:
                x3[0] = x1[0] - h
            x3[1] += p[1] / p[0] * (x3[0] - x[0])
            X[0] = (x2[0] * x2[0] - x3[0] * x3[0]) * func(x1) + (x3[0] * x3[0] - x1[0] * x1[0]) * func(x2) + (x1[0] * x1[0] - x2[0] * x2[0]) * func(x3)
            X[0] /= 2 * ((x2[0] - x3[0]) * func(x1) + (x3[0] - x1[0]) * func(x2) + (x1[0] - x2[0]) * func(x3))
            X[1] += p[1] / p[0] * (X[0] - x[0])
            xmin = argmin3(x1, x2, x3)
            if (X[0] >= x1[0] and X[0] <= x3[0]):
                x1 = argmin2(xmin, X)
            else:
                x1[0] = X[0]
            x1[1] += p[1] / p[0] * (x1[0] - x[0])
    else:
        x1[1] = (interval[0] + interval[1]) / 2
        while (np.fabs(func(xmin) - func(X)) > eps and np.fabs(xmin[1] - X[1]) > eps):
            x2[1] = x1[1] + h
            if (func(x1) > func(x2)):
                x3[1] = x1[1] + 2 * h
            else:
                x3[1] = x1[1] - h
            X[1] = (x2[1] * x2[1] - x3[1] * x3[1]) * func(x1) + (x3[1] * x3[1] - x1[1] * x1[1]) * func(x2) + (x1[1] * x1[1] - x2[1] * x2[1]) * func(x3)
            X[1] /= 2 * ((x2[1] - x3[1]) * func(x1) + (x3[1] - x1[1]) * func(x2) + (x1[1] - x2[1]) * func(x3))
            xmin = argmin3(x1, x2, x3)
            if (X[1] >= x1[1] and X[1] <= x3[1]):
                x1 = argmin2(xmin, X)
            else:
                x1[1] = X[1]
                
    return X

def coordinateDescent(x, eps):
    iter = 0
    coef = 0
    p = np.zeros(2, dtype = float)
    x0 = np.zeros(2, dtype = float)
    path = np.array(x, float)
    
    while (np.fabs(func(x0) - func(x)) > eps and np.fabs(x0[coef] - x[coef]) > eps):
        x0 = x.copy()
        if (iter % 2 == 0):
            coef = 0
        else:
            coef = 1
        p[coef] = 1
        x = squareInterp(x, p, 10, eps)
        path = np.append(path, x)
        iter += 1
        p[coef] = 0
    
    print(iter, "iterations")
    path = path.reshape((iter + 1, 2))
    show(path)
    return x

eps = 1e-7
x = np.array([9, 25], float)

print("Coordinate descent method")
extr = coordinateDescent(x, eps);
print("Maximum is on point (", extr[0], ",", extr[1], ")\nfmax =", -func(extr))