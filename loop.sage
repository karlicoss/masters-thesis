# Imports for IDE
from sage.all_cmdline import reset, exp, view, matrix, i, sin, cos, solve, norm, n
from sage.plot import complex_plot
from sympy import var
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
    return var(name, domain='complex')


def rvar(name):
    return var(name, domain='real')


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

    # def solve_analytic(self):
    #     a = self.a
    #     W = self.wires
    #     L = self.L
    #     # # to make formulas look pretty
    #
    #     num = 2 * W * k * cos(L * k) + (a + 2 * i * k) * sin(L * k) - 2 * W * k
    #     den = 2 * W * k * cos(L * k) + (a - 2 * i * k) * sin(L * k) - 2 * W * k
    #     return num / den
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
                    print(safe_subst(analytic, solver.k, k)) # ugh
            print("===================")


# LoopSolver.test()

def plot_all(det, suffix=None, rrange=(-2, 20), irange=(-2, 2), points=1500):
    plot_abs = complex_plot(abs(det), rrange, irange, plot_points=points)
    plot_abs.save("abs" + ("" if suffix is None else "_" + suffix) + ".png", figsize=[12, 2])


LoopSolver.test_abs_convergence()

# solver = LoopSolver(a=10, L=1)
# plot_all(solver.analytic_Sdet(), suffix="loop_new")
# view_later(solver.symbolic_Sdet())
# view_later(solver.analytic_Sdet())
# view_all()
