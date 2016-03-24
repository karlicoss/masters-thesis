reset()

from itertools import product
from pprint import pprint
import string

from sage.misc.viewer import viewer
viewer.pdf_viewer("atril")


def rvar(name):
    return var(name, domain='real')

_view_later = []
def view_later(formula):
    _view_later.append(formula)

def view_all():
    view(_view_later, tightpage=True, debug=False)


# B is R
# C is T
A, B, C, D = var('A B C D', domain='complex')
k = var('k', domain='complex')

f_inc(x) = A * exp(i * k * x) + B * exp(-i * k * x)
f_out(x) = C * exp(i * k * x) + D * exp(-i * k * x)

# f_inc(x) = i * (A - B) * sin(k * x) + (A + B) * cos(k * x)
# f_out(x) = i * (C - D) * sin(k * x) + (C + D) * cos(k * x)
# f_inc(x) = A * sin(k * x) + B * cos(k * x)
# f_out(x) = C * sin(k * x) + D * cos(k * x)

f_inc_d = f_inc.derivative(x)
f_out_d = f_out.derivative(x)


def symmetric_Smatrix(R, T):
    return matrix([
        [R, T], 
        [T, R]
    ])


def asymmetric_Smatrix(Bs, Cs):
    return matrix([
        [Bs.coefficient(A), Bs.coefficient(D)],
        [Cs.coefficient(A), Cs.coefficient(D)],
    ])


# Checks the matrix that it is a valid S-matrix
def sanity_checks(S): # variable: k
    for v in [-10, -1, 1, 10]:
        show("Unitary on {}:".format(v))
        show(n(S(k=v) * S(k=v).H))
        
    for v in [1 + i, 1 - i, -1 - i, -1 + i]:
        show("S-matrix property on {}:".format(v))
        show(n(S(k=v) * S(k=v.conjugate()).H))


class PopovSolver(object):
    # wires: integer
    # a: symbolic/float
    def __init__(self, wires, a=rvar('a')):
        super(PopovSolver, self).__init__()
        self.wires = wires
        self.a = a

    def spectrum(self, L_val=1):
        return [pi * k / L_val for k in range(1, 10)]

    def solve_analytic(self):
        a = self.a
        W = self.wires
        L = self.L
        # to make formulas look pretty

        num = W * k * cos(k) + (a + 2 * i * k) * sin(k)
        den = W * k * cos(k) + (a - 2 * i * k) * sin(k)
        return num / den


    def solve_symbolic(self, L_val=1):
        L = rvar('L')
        letters = string.uppercase[7:][:self.wires]
        ones = [var(p + '1') for p in letters]
        twos = [var(p + '2') for p in letters]

        wavefunctions  = [o * sin(k * x) + t * cos(k * x) for o, t in zip(ones, twos)]
        wavefunctionsd = [wf.derivative(x) for wf in wavefunctions]

        equations = []

        for wf in wavefunctions:
            equations.append(wf(x=L) == 0)

        for wf in wavefunctions + [f_out]:
            equations.append(f_inc(x=0) == wf(x=0))

        equations.append(-f_inc_d(x=0) + sum(wfd(x=0) for wfd in wavefunctionsd) + f_out_d(x=0) == a * f_out(x=0))

        solutions = solve(
            equations,
            B, C, *(ones + twos),
            solution_dict=True
        )

        Bs = solutions[0][B].full_simplify()
        Cs = solutions[0][C].full_simplify()
        SM = asymmetric_Smatrix(Bs, Cs)(L=L_val)
        return SM.det()


class LoopSolver(object):
    # wires: integer
    # a: symbolic/float
    def __init__(self, wires, a=rvar('a'), L=rvar('L')):
        super(LoopSolver, self).__init__()
        self.wires = wires
        self.a = a
        self.L = L
    
    def spectrum(self, L_val=1):
        return [k / (L_val / (2 * pi)) for k in range(1, 10)]

    def solve_analytic(self):
        a = self.a
        W = self.wires
        L = self.L
        # # to make formulas look pretty

        num = 2 * W * k * cos(L * k) + (a + 2 * i * k) * sin(L * k) - 2 * W * k
        den = 2 * W * k * cos(L * k) + (a - 2 * i * k) * sin(L * k) - 2 * W * k
        return num / den


    def solve_symbolic(self):
        a = self.a
        L = self.L


        letters = string.uppercase[7:][:self.wires]
        ones = [var(p + '1') for p in letters]
        twos = [var(p + '2') for p in letters]

        wavefunctions  = [o * sin(k * x) + t * cos(k * x) for o, t in zip(ones, twos)]
        wavefunctionsd = [wf.derivative(x) for wf in wavefunctions]

        equations = []

        equations.append(f_inc(x=0) == f_out(x=0))

        for wf in wavefunctions:
            equations.append(f_inc(x=0) == wf(x=0))
            equations.append(wf(x=L) == f_out(x=0))

        equations.append(-f_inc_d(x=0) + sum(wfd(x=0) - wfd(x=L) for wfd in wavefunctionsd) + f_out_d(x=0) == a * f_inc(x=0))

        solutions = solve(
            equations,
            B, C, *(ones + twos),
            solution_dict=True
        )

        Bs = solutions[0][B].full_simplify()
        Cs = solutions[0][C].full_simplify()
        SM = asymmetric_Smatrix(Bs, Cs)
        return SM.det()


# TODO look at eigenvalues
class IntervalSolver(object):
    # wires: integer
    # a: symbolic/float
    # b: symbolic/float
    def __init__(self, wires, a=rvar('a'), b=rvar('b')):
        super(IntervalSolver, self).__init__()
        self.wires = wires
        self.a = a
        self.b = b


    def solve_analytic(self):
        a = self.a
        b = self.b
        W = self.wires
        # to make formulas look pretty

        L = 1
        coeff = 0 if W == 0 else W**2 + 1
        # TODO 2 * k might be becasue the actual length is 2 * L instead of L :)
        num = W * (a + b) * k * cos(2 * k) + (a * b - coeff * k**2) * sin(2 * k) + 2 * W * i * k**2 * cos(2 * k) + i * (a + b) * k * sin(2 * k)
        den = W * (a + b) * k * cos(2 * k) + (a * b - coeff * k**2) * sin(2 * k) - 2 * W * i * k**2 * cos(2 * k) - i * (a + b) * k * sin(2 * k)
        return num / den

    def solve_symbolic(self, L_val=1):
        a = self.a
        b = self.b


        L = rvar('L')
        letters = string.uppercase[7:][:self.wires]
        ones = [var(p + '1') for p in letters]
        twos = [var(p + '2') for p in letters]

        wavefunctions  = [o * sin(k * x) + t * cos(k * x) for o, t in zip(ones, twos)]
        wavefunctionsd = [wf.derivative(x) for wf in wavefunctions]

        equations = []

        if self.wires == 0:
            equations.append(f_inc(x=0) == f_out(x=0))
            equations.append(-f_inc_d(x=0) + f_out_d(x=0) == a * f_inc(x=0))
        else:
            last = f_inc(x=0)
            for wf in wavefunctions:
                cur = wf(x=-L)
                equations.append(last == cur)
                last = cur

            last = f_out(x=0)
            for wf in wavefunctions:
                cur = wf(x=L)
                equations.append(last == cur)
                last = cur

            derl = -f_inc_d(x=0) + sum(wfd(x=-L) for wfd in wavefunctionsd) == a * f_inc(x=0)
            derr =  f_out_d(x=0) - sum(wfd(x= L) for wfd in wavefunctionsd) == b * f_out(x=0)

            equations.append(derl)
            equations.append(derr)

        solutions = solve(
            equations,
            B, C, *(ones + twos),
            solution_dict=True
        )

        Bs = solutions[0][B].full_simplify()
        Cs = solutions[0][C].full_simplify()
        SM = asymmetric_Smatrix(Bs, Cs)(L=L_val)
        return SM.det()

# TODO look at eigenvalues
class FractalSolver(object):
    # wires: integer
    # a: symbolic/float
    # b: symbolic/float
    def __init__(self, wires, a=rvar('a')):
        super(FractalSolver, self).__init__()
        self.wires = wires
        self.a = a


    def solve_analytic(self):
        a = self.a
        W = self.wires
        # to make formulas look pretty

        L = 1
        if W == 0:
            num = a + 2 * i * k
            den = a - 2 * i * k
            return num / den
        elif W == 1: 
            # coeff = 0 if W == 0 else W**2 + 1
            # num = W * (a + b) * k * cos(2 * k) + (a * b - coeff * k**2) * sin(2 * k) + 2 * W * i * k**2 * cos(2 * k) + i * (a + b) * k * sin(2 * k)
            # den = W * (a + b) * k * cos(2 * k) + (a * b - coeff * k**2) * sin(2 * k) - 2 * W * i * k**2 * cos(2 * k) - i * (a + b) * k * sin(2 * k)
            num = 2 * a * k * cos(k) + (a**2 - 2 * k**2) * sin(k) + 2 * i * k**2 * cos(k) + 2 * i * a * k * sin(k)
            den = 2 * a * k * cos(k) + (a**2 - 2 * k**2) * sin(k) - 2 * i * k**2 * cos(k) - 2 * i * a * k * sin(k)
            return num / den
        elif W == 2:
            num = -(4 * k ** 3 - 4 * a**2 * k) * cos(k) * sin(k) + (-6 * a * k**2 + a**3) * sin(k)**2 + 3 * a * k**2 + 2 * i * k**3 + 6 * i * a * k**2 * cos(k) * sin(k) - (4 * i * k**3 - 2 * i * a**2 * k) * sin(k)**2
            den = -(4 * k ** 3 - 4 * a**2 * k) * cos(k) * sin(k) + (-6 * a * k**2 + a**3) * sin(k)**2 + 3 * a * k**2 - 2 * i * k**3 - 6 * i * a * k**2 * cos(k) * sin(k) + (4 * i * k**3 - 2 * i * a**2 * k) * sin(k)**2
            return num / den
        elif W == 3:
            return None
        else:
            return None

    def solve_symbolic(self, L_val=1):
        L = rvar('L')
        letters = string.uppercase[7:][:self.wires]
        ones = [var(p + '1') for p in letters]
        twos = [var(p + '2') for p in letters]

        wavefunctions  = [o * sin(k * x) + t * cos(k * x) for o, t in zip(ones, twos)]
        wavefunctionsd = [wf.derivative(x) for wf in wavefunctions]

        equations = []

        all_wfs  = wavefunctions + [f_out]
        all_wfds = wavefunctionsd + [f_out_d]
        last = f_inc(x=0)
        lastd = f_inc_d(x=0)
        curL = L
        for cur, curd in zip(all_wfs, all_wfds):
            equations.append(last == cur(x=0))
            equations.append(-lastd + curd(x=0) == a * cur(x=0))
            last = cur(x=curL) # TODO adjust length
            lastd = curd(x=curL)
            # curL = curL / 2

        solutions = solve(
            equations,
            B, C, *(ones + twos),
            solution_dict=True
        )

        Bs = solutions[0][B].full_simplify()
        Cs = solutions[0][C].full_simplify()
        SM = asymmetric_Smatrix(Bs, Cs)
        return SM(L=L_val).det()

def icayley(x):
    return i * (1 + x) / (1 - x)

DPI = 200

def plot_all(Sdet, suffix="", rrange=(-15, 15), irange=(-5, 5), points=500):
    complex_plot(Sdet, rrange, irange, plot_points=points).save('plot{}.png'.format(suffix), figsize=[12, 2])
    complex_plot(abs(Sdet), rrange, irange, plot_points=points).save('plot_abs{}.png'.format(suffix), figsize=[12, 2])
    complex_plot(ln(abs(Sdet)), rrange, irange, plot_points=points).save('plot_ln{}.png'.format(suffix), figsize=[12, 2])
    # unit_circle = circle((0, 0), 1)
    # (complex_plot(ln(abs(Sdet(k=cayley(k)))), (-1, 1), (-1, 1), plot_points=points) + unit_circle).save('plot_circle.png', dpi=DPI)

# S = solve_popov_analytic(a_val=10)
# S = solve_double_loop(a_val=5, b_val=5)

def contour_integral_analysis(expr):
    ig(k) = ln(abs(expr(k=k))) / (k - 1) ** 2
    for R in [1, 5, 10, 20, 50, 100, 200, 400, 500, 1000, 10000]:
        pig(t) = ig(k = R * exp(i * t) + R * i) * R * i * exp(i * t)
        show(numerical_integral(lambda q: pig(t=q).real(), 0, pi))
        show(numerical_integral(lambda q: pig(t=q).imag(), 0, pi))    


def region_analysis(expr):
    x, y = var('x y')
    region_plot(lambda x, y: 0 <= abs(expr(k=x + i * y)) <= 0.5, (x, -10, 10), (y, -10, 10)).save('region.png')


def test_matrices(S, Sa):
    from numpy import arange
    for rp in arange(0, 3, 0.4):
        for ip in arange(0, 3, 0.4):
            k = rp + i * ip
            try:
                diff = n(S(k=k)) - n(Sa(k=k))
                if abs(diff) < 0.0001:
                    print("OK")
                else:
                    print("ERROR")
            except ValueError as err:
                if "division by zero" in err.message:
                    print("computation failed on {}".format(k))
                else:
                    raise err



# contour_integral_analysis(S)
# S = solve_double_loop(a_val=5)

# S = solve_interval()(a=a)
# view_later(S)


a = var('a', domain='real')
b = var('b', domain='real')
# a = 2
# b = -2
# !!!! case 0 only works for a = b


# for wires in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
#     S = IntervalSolver(wires, a=a, b=b).solve_symbolic()
#     Sa = IntervalSolver(wires, a=a, b=b).solve_analytic()
#     # view_later(S)
#     test_matrices(S, Sa)

# view_all()
    

# view_all()



# region_analysis(S)




from numpy import arange
from sage.plot.colors import rainbow
# for q in arange(0, 3, 1):
    # S = solve_double_loop(a_val=q, b_val=-q)
    # plot_all(S, str(q))

def plot_bound_energies():
    p = var('p')
    plots = []
    values = arange(-3, 3, 0.2)
    for q, col in zip(values, rainbow(len(values))):
        print("Plotting for " + str(q))
        Sdet = solve_double_loop_analytic(a_val=q)
        plots.append(plot(abs(Sdet(k = p * i)), (p, 0, 4), ymin=0, ymax=10, plot_points=30, color=col, legend_label="q={}".format(q)))

    sum(plots).save('plot_fewfef.png', dpi=DPI)


def check_zero(Sdet):
    t = var('t')
    for R in [1, 2, 4, 16, 32, 64, 128, 256]:
        ig = ln(abs(Sdet)) * 1 / (k) ** 2
        ff = ig(k = R * exp(i * t) + i * R) * R * i * exp(i * t)
        print("R = {}".format(R))
        show(ff)
        show(integral(ff, t, 0, pi))
        show(numerical_integral(lambda q: ff(t=q).real(), 0, pi))
        show(numerical_integral(lambda q: ff(t=q).imag(), 0, pi))    
    
# a = 2
# b = 1
# solver = IntervalSolver(1, a=a, b=b)
# Sd = solver.solve_analytic()
# plot_all(Sd)
# check_zero(Sd)

# a = 1
a = rvar('a')
L = rvar('L')

a = 2
L = 3


for w in [1, 2, 3]:
    solver = LoopSolver(w, a=a, L=L)
    S = solver.solve_symbolic()
    Sa = solver.solve_analytic()
    test_matrices(S, Sa)
    # view_later(S)
    
    # view_later(Sa)

# view_all()

# for w in [1, 2, 3, 4]:
#     solver = PopovSolver(w, a=0)
#     # S = solver.solve_symbolic()
#     Sa = solver.solve_analytic()
#     # view_later(S)
#     view_later(Sa)
# view_all()

# for w in [0, 1, 2, 3, 4]:
#     solver = PopovSolver(w, a=a)
#     # S = solver.solve_symbolic()
#     Sa = solver.solve_symbolic()
#     # view_later(S)
#     view_later(Sa)
# view_all()

# loop = LoopSolver(1, a=a)
# print(loop.spectrum(L_val=2))
# plot_all(loop.solve_symbolic(L_val=2), suffix='loop')

# popov = PopovSolver(1, a=a)
# print(popov.spectrum())
# plot_all(popov.solve_symbolic(L_val=1), suffix='popov')


    # test_matrices(S, Sa)
# view_all()

# solver = FractalSolver(0, a=a)
# S = solver.solve_symbolic()
# view_later(S)
# Sa = solver.solve_analytic()
# view_later(Sa)

# solver = FractalSolver(1, a=a)
# S = solver.solve_analytic()
# Sa = solver.solve_symbolic()
# # test_matrices(S, Sa)
# # S = solver.solve_symbolic(L_val=1)
# # test_matrices(S, Sa)
# view_later(S)
# view_later(Sa)

# solver = FractalSolver(2, a=a)
# Sa = solver.solve_analytic()
# Sa = solver.solve_symbolic()
# # S = solver.solve_symbolic(L_val=1)
# # test_matrices(S, Sa)
# view_later(Sa)


# solver = FractalSolver(3, a=a)
# Sa = solver.solve_symbolic()
# view_later(Sa)

# solver = FractalSolver(4, a=a)
# Sa = solver.solve_symbolic()
# view_later(Sa)

# solver = FractalSolver(2, a=a)
# S = solver.solve_symbolic(L_val=1)
# view_later(S)

# test_matrices(S, Sa)

# for wires in range(0, 4):
#     solver = FractalSolver(wires, a=a)
#     Sd = solver.solve_symbolic(L_val=1)
#     view_later(Sd)
#     # plot_all(Sd, suffix=str(wires))
# view_all()
# v = SM.eigenvalues()[0]
# t = var('t')
# plot(ln(v(k=i*t).real()), (t, 0, 10)).save('plot.png')