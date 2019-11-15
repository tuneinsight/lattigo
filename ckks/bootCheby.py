from math import cos, pi, sin
from functools import reduce
from operator import mul






class Chebyshev:
    """
    Chebyshev(a, b, n, func)
    Given a function func, lower and upper limits of the interval [a,b],
    and maximum degree n, this class computes a Chebyshev approximation
    of the function.
    Method eval(x) yields the approximated function value.
    """
    def __init__(self, a, b, n, func):
        self.a = a
        self.b = b
        self.func = func

        bma = 0.5 * (b - a)
        bpa = 0.5 * (b + a)

        K = 9
        deg = [3, 3, 3, 4, 5, 4, 3, 3, 3]
        nodes = []
        e = 0.01
        for i in range(K):
            for j in range(1, deg[i]+1):
                nodes += [i-(K>>1) + e * cos(((j-0.5)/deg[i])*pi)]
        #nodes = [cos(pi * (k + 0.5) / n) * bma + bpa for k in range(n)]

        
        f = [func(nodes[k]) for k in range(n)]
        fac = 2.0 / n
        self.c = [fac * sum([f[k] * cos(pi * j * (k + 0.5) / n) for k in range(n)]) for j in range(n)]

        #self.c = self.lagrange(nodes, func)

        self.c = [self.Lagrange_polynomial(f[i], i, nodes) for i in range(len(nodes))]

        #print(self.c(2.001))


    def Lagrange_polynomial(self, x, i, points):
        p = 1
        for k in range(len(points)):
            if k != i:
                p *= (x - points[k])/(points[i] - points[k])
        return p

    def lagrange(self, pts, f):
        """ Return a Lagrange interpolating polynomial for f
        at the points pts. """
     
        def basis(others):
            # w_{i}(x) = (x - x_{1}) * ... * (x - x_{i-1}) * (x - x_{i+1}) * ...
            return lambda x : reduce(mul, [x - o for o in others])
     
        # Get each basis polynomial and coefficients
        bases = [basis(pts[:i] + pts[i + 1:]) for i in range(len(pts))]
        # a_{i} = f(x_{i}) / w_{i}(x_{i})
        weights = [f(pt) / bases[i](pt) for i, pt in enumerate(pts)]

        return lambda x : sum([weights[i] * bases[i](x) for i in range(len(pts))])

    def eval(self, x):
        a,b = self.a, self.b
        assert(a <= x <= b)
        y = (2.0 * x - a - b) * (1.0 / (b - a))
        y2 = 2.0 * y
        (d, dd) = (self.c[-1], 0)             # Special case first step for efficiency
        for cj in self.c[-2:0:-1]:            # Clenshaw's recurrence
            (d, dd) = (y2 * d - dd + cj, d)
        return y * d - dd + 0.5 * self.c[0]   # Last step is different

def f(x):
    return 1.0/(2*pi) * sin(2*pi*x)
ch = Chebyshev(-4, 4, 31, f)
print(ch.eval(0.001))
print(f(2.001))
