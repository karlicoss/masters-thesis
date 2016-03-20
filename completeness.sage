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
A, B, C, D = var('A B C D')
k = var('k')
assume(k, 'imaginary')

f_inc(x) = A * exp(i * k * x) + B * exp(-i * k * x)
f_out(x) = C * exp(i * k * x) + D * exp(-i * k * x)

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


def solve_popov_analytic(a_val=0):
    a = a_val
    rp = k * cos(k) + a * sin(k)
    ip = 2 * k * sin(k)
    Sd = (rp + i * ip) / (rp - i * ip)
    return Sd

def solve_popov(a_val=0, L_val=1):
    P = var('P')
    a, L = var('a L')
    assume(a, 'real')
    assume(L, 'real')

    f2(x) = P * sin(k * x)
    f2d = f2.derivative(x)

    solutions = solve([f_inc(0) == f2(L), f2(L) == f_out(0), -f_inc_d(0) - f2d(L) + f_out_d(0) == a * f_inc(0)], B, C, P, solution_dict=True)

    Bs = solutions[0][B].full_simplify()
    Cs = solutions[0][C].full_simplify()
    SM = asymmetric_Smatrix(Bs, Cs)
    # view_later(SM(L=1).eigenvalues())
    return SM(a=a_val, L=L_val)


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
        L = 1
        sprimes = sum(primes_first_n(self.wires))
        num = self.wires * (self.a + self.b) * k * cos(2 * k) + 2 * self.wires * i * k**2 * cos(2 * k) + i * (self.a + self.b) * k * sin(2 * k) + (self.a * self.b - sprimes * k**2) * sin(2 * k)
        den = self.wires * (self.a + self.b) * k * cos(2 * k) - 2 * self.wires * i * k**2 * cos(2 * k) - i * (self.a + self.b) * k * sin(2 * k) + (self.a * self.b - sprimes * k**2) * sin(2 * k)
        return num / den

    def solve_symbolic(self, L_val=1):
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
        SM = asymmetric_Smatrix(Bs, Cs)
        return SM(L=L_val).det()


def solve_loop():
    A, B = var('A B') # TODO FIXME
    a, L = var('a L')
    assume(a, 'real')
    assume(L, 'real')

    f2(x) = A * sin(k * x) + B * cos(k * x) # exp(i * k * x) + B * exp(- i * k * x)
    f2d = f2.derivative(x)

    solutions = solve([f_inc(0) == f2(0), f2(0) == f2(L), f2(0) == f_out(0), -f_inc_d(0) + f2d(0) - f2d(L) + f_out_d(0) == a * f_inc(0)], R, T, A, B, solution_dict=True)
    
    Rs = solutions[0][R].full_simplify()
    Ts = solutions[0][T].full_simplify()

    show("R = ", Rs)
    show("T = ", Ts)
    return Rs(a = 1, L = 1), Ts(a = 1, L = 1)


rrange = (0, 1000)
irange = (0, 10)
points = 500


def icayley(x):
    return i * (1 + x) / (1 - x)

DPI = 200

def plot_all(Sdet, suffix=""):
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

a = 2
b = 2
# !!!! case 0 only works for a = b
# a = var('a', domain='real')
# b = var('b', domain='real')


S = IntervalSolver(0, a=a, b=b).solve_symbolic()
Sa = IntervalSolver(0, a=a, b=b).solve_analytic()
test_matrices(S, Sa)

# view_all()

S = IntervalSolver(1, a=a, b=b).solve_symbolic()
Sa = IntervalSolver(1, a=a, b=b).solve_analytic()
test_matrices(S, Sa)


# S = solve_double_loop().rational_simplify(algorithm='noexpand')(a=a, b=b)
S = IntervalSolver(2, a=a, b=b).solve_symbolic()
Sa = IntervalSolver(2, a=a, b=b).solve_analytic()
test_matrices(S, Sa)

S = IntervalSolver(3, a=a, b=b).solve_symbolic()
Sa = IntervalSolver(3, a=a, b=b).solve_analytic()
test_matrices(S, Sa)

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
        ig = ln(abs(Sdet)) * 1 / (k - 1) ** 2
        ff = ig(k = R * exp(i * t) + i * R) * R * i * exp(i * t)
        print("R = {}".format(R))
        show(numerical_integral(lambda q: ff(t=q).real(), 0, pi))
        show(numerical_integral(lambda q: ff(t=q).imag(), 0, pi))    
    

# S = solve_popov(a_val=0)

# S = solve_double_loop(a_val=-2, b_val=-2)
# Sd = S.det()
# y = var('y')
# solve_double_loop(a_val=1)
# Sd = solve_double_loop_analytic(a_val=y)
# check_zero(Sd)
# plot_bound_energies()

# plot_all(Sd)

# view_later(Sd)
# solve_popov()
# solve_interval()
# view_all()
# draw_double_loop()
# SM = solve_popov(a_val=0, L_val=1)
# v = SM.eigenvalues()[0]
# t = var('t')
# plot(ln(v(k=i*t).real()), (t, 0, 10)).save('plot.png')