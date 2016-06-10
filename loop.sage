# Imports for IDE
from sage.all_cmdline import reset, exp, view, matrix, i, sin, cos, solve, norm, n
from sage.plot import complex_plot
from sage.symbolic.assumptions import assume
from sage.symbolic.expression_conversions import fast_callable
from sympy import var, latex, pretty_print, diff, pi, CC, RR, ln, oo

# apparently, resetting is necessary for script mode =/
reset()

from itertools import product
from numpy import arange
### stuff for latex viewing

from sage.misc.viewer import viewer

viewer.pdf_viewer("atril")

_view_later = []


def view_later(formula):
    _view_later.append(formula)


def view_all():
    view(_view_later, tightpage=True, debug=False)


###

def safe_subst(expr, v, value):
    try:
        res = n(expr.subs({v: value}))
        return res
    except ZeroDivisionError:
        return "division by zero"
    except ValueError as err:
        if "division by zero" in err.message:
            return "division by zero"
        else:
            raise err


def cvar(name):
    return var(name, domain=CC)


def rvar(name):
    return var(name, domain=RR)


def pvar(name):
    return var(name, domain='positive')


def symmetric_Smatrix(Rs, Ts):
    return matrix([
        [Rs, Ts],
        [Ts, Rs]
    ])


class BaseSolver(object):  # symmetric?
    def __init__(self):
        super(BaseSolver, self).__init__()
        self.R, self.T = var('R'), cvar('T')
        self.k = cvar('k')
        self.x = rvar('x')
        self.f_inc = exp(i * self.k * self.x) + self.R * exp(-i * self.k * self.x)
        self.f_out = self.T * exp(i * self.k * self.x)
        self.f_inc_d = self.f_inc.derivative(self.x)
        self.f_out_d = self.f_out.derivative(self.x)


# Checks the matrix that it is a valid S-matrix
def check_Smatrix_properties(S):  # variable: k
    for v in [-10, -1, 1, 10]:
        print("Unitary on real axis: {}".format(v))
        m = n(S(k=v) * S(k=v).H)
        # show(m)
        if norm(m - 1) > 0.01:
            print "!!!!!!!!!!ERROR!!!!!!!!!"
        else:
            print "OK"

    for q in [1, 2, 4, 8]:
        for v in [1 + i, 1 - i, -1 - i, -1 + i]:
            w = q * v
            print("S-matrix property on {}:".format(w))
            m = n(S(k=w) * S(k=w.conjugate()).H)
            # show(m)
            if norm(m - 1) > 0.01:
                print "!!!!!!!!!!ERROR!!!!!!!!!"
            else:
                print "OK"


def check_determinants_same(Sdet_s, Sdet_a):
    for rp, ip in product(arange(-3, 3, 0.5), arange(-3, 3, 0.5)):
        k = rp + i * ip
        try:
            diff = n(Sdet_s(k=k) - Sdet_a(k=k))
            if abs(diff) < 0.0001:
                print("OK")
            else:
                print("ERROR")
        except ValueError as err:
            if "division by zero" in err.message:
                print("computation failed on {}".format(k))
            else:
                raise err
                # except RuntimeError as err:
                #     print(err)


class LoopSolver(BaseSolver):
    # a: symbolic/float
    def __init__(self, a=rvar('a'), L=rvar('L')):
        super(LoopSolver, self).__init__()
        self.a = a
        self.L = L

    def analytic_Sdet(self):
        a = self.a
        L = self.L
        k = self.k
        num = 2 * k * cos(L * k) + (a + 2 * i * k) * sin(L * k) - 2 * k
        den = 2 * k * cos(L * k) + (a - 2 * i * k) * sin(L * k) - 2 * k
        return num / den

    def symbolic_S(self):
        R, T = self.R, self.T
        f_inc, f_out = self.f_inc, self.f_out
        f_inc_d, f_out_d = self.f_inc_d, self.f_out_d
        x = self.x
        k = self.k
        a = self.a
        L = self.L

        P, Q = cvar('P'), cvar('Q')

        wf = P * sin(k * x) + Q * cos(k * x)
        wfd = wf.derivative(x)

        equations = [
            f_inc(x=0) == f_out(x=0),
            f_inc(x=0) == wf(x=0),
            wf(x=L) == f_out(x=0),
            -f_inc_d(x=0) + wfd(x=0) - wfd(x=L) + f_out_d(x=0) == a * f_out(x=0)  # TODO f_inc?
        ]

        solution = solve(
            equations,
            R, T, P, Q,
            solution_dict=True,
        )[0]

        Rs = solution[R].full_simplify()
        Ts = solution[T].full_simplify()
        return symmetric_Smatrix(Rs, Ts)

    def symbolic_Sdet(self):
        return self.symbolic_S().det()

    @staticmethod
    def test():
        for a in [-1, 0, 1]:
            print("a = " + str(a))
            solver = LoopSolver(a=a, L=1)
            symbolic = solver.symbolic_S()
            check_Smatrix_properties(symbolic)
            analytic = solver.analytic_Sdet()
            check_determinants_same(symbolic.det(), analytic)

    @staticmethod
    def test_abs_convergence():
        for a in [-1, 0, 1]:
            print("a = " + str(a))
            solver = LoopSolver(a=a, L=1)
            analytic = solver.analytic_Sdet()
            for rp in [-5, 0, 5]:
                print("RP = " + str(rp))
                for ip in [0, 10, 25, 50, 100]:
                    k = rp + i * ip
                    print(safe_subst(analytic, solver.k, k))  # ugh
            print("===================")


class IntervalLoopSolver(BaseSolver):
    # wires: integer
    # a: symbolic/float
    def __init__(self, L=rvar('L')):
        super(IntervalLoopSolver, self).__init__()
        self.L = L

    def analytic_Sdet(self):
        L = self.L
        k = self.k
        num = (2 * cos(k) + 2 * i * sin(k)) * cos(L * k) + (+2 * i * cos(k) - 3 * sin(k)) * sin(L * k) - 2
        den = (2 * cos(k) - 2 * i * sin(k)) * cos(L * k) + (-2 * i * cos(k) - 3 * sin(k)) * sin(L * k) - 2
        return num / den

    def symbolic_S(self):
        R, T = self.R, self.T
        f_inc, f_out = self.f_inc, self.f_out
        f_inc_d, f_out_d = self.f_inc_d, self.f_out_d
        x = self.x
        k = self.k
        L = self.L

        P, Q = cvar('P'), cvar('Q')
        wf = P * sin(k * x) + Q * cos(k * x)
        wfd = wf.derivative(x)

        PL, QL = cvar('PL'), cvar('QL')
        wfl = PL * sin(k * x) + QL * cos(k * x)
        wfld = wfl.derivative(x)

        equations = [
            f_inc(x=0) == wf(x=0),
            wf(x=1) == f_out(x=0),
            f_inc(x=0) == wfl(x=0),
            wfl(x=L) == f_out(x=0),

            -f_inc_d(x=0) + wfd(x=0) + wfld(x=0) == 0,
            f_out_d(x=0) - wfd(x=1) - wfld(x=L) == 0,
        ]

        solution = solve(
            equations,
            R, T, P, Q, PL, QL,
            solution_dict=True,
        )[0]

        Rs = solution[R].full_simplify()
        Ts = solution[T].full_simplify()
        return symmetric_Smatrix(Rs, Ts)

    def symbolic_Sdet(self):
        return self.symbolic_S().det()

    @staticmethod
    def test():
        for L in [1, 0.5, 2]:
            print("L = " + str(L))
            solver = IntervalLoopSolver(L=L)
            symbolic = solver.symbolic_S()
            check_Smatrix_properties(symbolic)
            analytic = solver.analytic_Sdet()
            check_determinants_same(symbolic.det(), analytic)


class IntervalsSolver(BaseSolver):
    # wires: integer
    # a: symbolic/float
    def __init__(self, W=rvar('W'), a=rvar('a'), b=rvar('b'), L=rvar('L')):
        super(IntervalsSolver, self).__init__()
        self.W = W
        self.a = a
        self.b = b
        self.L = L

    def analytic_Sdet(self):
        W = self.W
        k = self.k
        a = self.a
        b = self.b
        L = self.L

        coeff = W ** 2 + 1
        num = W * (a + b) * k * cos(L * k) + (a * b - coeff * k ** 2) * sin(L * k) + 2 * W * i * k ** 2 * cos(L * k) + i * (a + b) * k * sin(L * k)
        den = W * (a + b) * k * cos(L * k) + (a * b - coeff * k ** 2) * sin(L * k) - 2 * W * i * k ** 2 * cos(L * k) - i * (a + b) * k * sin(L * k)
        return num / den



def plot_all(det, suffix=None, rrange=(-20, 20), irange=(-2, 2), points=1500):
    plot_abs = complex_plot(abs(det), rrange, irange, plot_points=points)  # , aspect_ratio=)
    plot_abs.save("abs" + ("" if suffix is None else "_" + suffix) + ".png", figsize=[8, 4])


# LoopSolver.test_abs_convergence()

def divide_both(expr, d):
    num, den = [(e / d).expand() for e in expr.numerator_denominator()]
    print(num)
    print(den)
    print(num.maxima_methods().divide(den))
    # print((den / d).expand())
    # term = Term(1, numer=(num / d).expand(), denom=(den / d).expand())
    print(term)


def icayley(z):
    return i * (1 + z) / (1 - z)


def cayley(z):
    return (z - i) / (z + i)

CCC(r) = ((icayley(r) + icayley(-r)) / 2).imag()
RRR(r) = ((icayley(r) - icayley(-r)) / 2).imag()


def changevar(f, eqn, newvar):
    dx = diff(eqn.rhs(), newvar)
    # print(f)
    # print(f.subs(eqn))
    # print(dx)
    return f.subs(eqn) * dx


def complex_integral(expr, x, ff, tt):
    rp = fast_callable(expr.real(), vars=[x], domain=CC)
    ip = fast_callable(expr.imag(), vars=[x], domain=CC)
    eps = 1e-02
    rpart, _ = numerical_integral(rp, ff, tt, eps_rel=eps, rule=1)
    impart, _ = numerical_integral(ip, ff, tt, eps_rel=eps, rule=1)
    return rpart + i * impart


# In Cayley space
def loop_integral_symbolic():
    k = cvar('k')
    f = k.imag()
    t = rvar('t')
    R = pvar('R')
    assume(abs(R) - 1 < 0)
    pig = changevar(f(k=icayley(k)), k == R * exp(i * t), t).full_simplify()
    # print(complex_integral(pig(R=0.99), 0, 2 * pi))
    print(pig.integrate(t, 0, 2 * pi))


# TODO SIGN DISCREPANCY :(
# in ordinary space
def loop_integral_symbolic_2():
    k = cvar('k')
    f = k.imag()
    t = rvar('t')
    R = pvar('R')

    ig = f * diff(cayley(k), k)
    pig = changevar(ig, k == R * exp(i * t) + R * i, t).full_simplify()
    print(complex_integral(pig(R=5), 0, 2 * pi))
    # print(pig.integrate(t, 0, 2 * pi))
    # pig = f(k=)


def stuff_for_paper():
    rrange = (-50, 50)
    irange = (0, 50)

    # solver = LoopSolver(L=1)
    # Sdet = solver.analytic_Sdet()(a=1)
    # Sdet = IntervalsSolver(L=1, a=1, b=-2, W=3).analytic_Sdet()
    Sdet = IntervalLoopSolver(L=0.0001).analytic_Sdet()
    # print(Sdet)
    r = 0.9
    (complex_plot(abs(Sdet), (-15, 15), (0, 20)) + circle((0, CCC(r)), RRR(r))).save('plot_a.png', figsize=[12, 6], aspect_ratio='equal')
    # view_later("S-matrix for a")
    # view_later(Sdet)
    # view_later(divide_both(Sdet, 2 * solver.k))
    # solver0 = LoopSolver(a=0, L=1)
    # solver1 = LoopSolver(a=10, L=1)
    #
    # for solver in [solver0, solver1]:
    #     Sdet = solver.analytic_Sdet()
    #     # TODO black labels on black plot!!
    #     plot_all(Sdet, suffix="loop_%s" % solver.a, rrange=rrange, irange=irange)
    #     view_later("S-matrix for a = %s" % solver.a)
    #     view_later(Sdet)
    # view_all()


# In Cayley space
def calculate_integral_symbolic_cayley():
    t = rvar('t')
    R = pvar('R')
    solver = LoopSolver(L=1, a=1)
    k = solver.k
    f = ln(abs(solver.analytic_Sdet()))
    # k = cvar('k')
    # f = ln(abs((cos(k) + (1/k + i) * sin(k) - 1) / (cos(k) + (1/k - i) * sin(k) - 1)))
    assume(abs(R) - 1 < 0)
    pig = changevar(f(k=icayley(k)), k == R * exp(i * t), t)  # .full_simplify()
    # print(pig.limit(R=1))
    # print(complex_integral(pig(R=0.99), 0, 2 * pi))
    # print(pig.integrate(t, 0, 2 * pi))
    print(pig)
    for tt in arange(0.0001, 0.0001 * 2 * pi, 0.0001 * pi / 20):
        print(pig(R=1)(t=tt).n())
        # print(complex_integral(pig(R=1), t, 0, 2 * pi))


def calculate_integral_symbolic():
    t = rvar('t')
    R = pvar('R')
    solver = LoopSolver(L=1, a=1)
    k = solver.k
    f = ln(abs(solver.analytic_Sdet()))
    # k = cvar('k')
    # f = k.imag()
    print(f)
    ig = f * diff(cayley(k), k)
    pig = changevar(ig, k == R * exp(i * t) + R * i, t).full_simplify()
    print(pig(t=0.4999 * pi).simplify()(R=1000).n())  # .limit(R=oo))
    # for rr in [1, 2, 5, 10, 15, 20, 25]:
    #     print(complex_integral(pig(R=rr), t, 0, 2 * pi))
    # print(pig.integrate(t, 0, 2 * pi))
    # pig = f(k=)


stuff_for_paper()
