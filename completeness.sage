from sage.all_cmdline import *
from sympy import var

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


# det = -1 means 100% transmission?
def symmetric_Smatrix(R, T):
    return matrix([
        [R, T], 
        [T, R]
    ])


def asymmetric_Smatrix(Bs, Cs):
    R = Bs.coefficient(A)
    T = Bs.coefficient(D)
    view_later(R)
    view_later(T)
    return matrix([
        # [T, R],
        # [-conjugate(R), conjugate(T)],
        [Bs.coefficient(A), Bs.coefficient(D)],
        [Cs.coefficient(A), Cs.coefficient(D)],
    ])


# Checks the matrix that it is a valid S-matrix
def sanity_checks(S): # variable: k
    for v in [-10, -1, 1, 10]:
        show("Unitary on real axis: {}".format(v))
        m = n(S(k=v) * S(k=v).H)
        # show(m)
        if norm(m - 1) > 0.1:
            print "ERROR!!!"
        else:
            print "OK"
        
    for q in [1, 2, 4, 8]:
        for v in [1 + i, 1 - i, -1 - i, -1 + i]:
            w = q * v
            show("S-matrix property on {}:".format(w))
            m = n(S(k=w) * S(k=w.conjugate()).H)
            # show(m)
            if norm(m - 1) > 0.1:
                print "ERROR!!!"
            else:
                print "OK"


class PopovSolver(object):
    # wires: integer
    # a: symbolic/float
    def __init__(self, wires, a=rvar('a'), L=rvar('L')):
        super(PopovSolver, self).__init__()
        self.wires = wires
        self.a = a
        self.L = L

    def spectrum(self, L_val=1):
        return [pi * k / L_val for k in range(1, 10)]

    # CORRECT SOLUTION (FROM THE PAPER)
    def solve_analytic_2(self):
        a = self.a
        W = self.wires
        L = self.L
        # to make formulas look pretty

        T = 2 * i * k / (2 * i * k - a - k * cot(k))
        R = (a + k * cot(k)) / (2 * i * k - a - k * cot(k))
        # num = W * k * cos(L * k) + (a + 2 * i * k) * sin(L * k)
        # den = W * k * cos(L * k) + (a - 2 * i * k) * sin(L * k)
        return matrix([[R, T],[T, R]]).det()
        # return num / den

    def solve_analytic(self):
        a = self.a
        W = self.wires
        L = self.L
        # to make formulas look pretty

        num = W * k * cos(L * k) + (a + 2 * i * k) * sin(L * k)
        den = W * k * cos(L * k) + (a - 2 * i * k) * sin(L * k)
        return num / den

    def solve_symbolic_S(self):
        W = self.wires
        a = self.a
        L = self.L

        letters = string.uppercase[7:][:W]
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
        SM = asymmetric_Smatrix(Bs, Cs)
        return SM

    def solve_symbolic(self):
        return self.solve_symbolic_S().det()


# TODO look at eigenvalues
class PopovSolver2(object):
    # wires: integer
    # a: symbolic/float
    def __init__(self, wires, a=rvar('a'), L=rvar('L')):
        super(PopovSolver2, self).__init__()
        self.wires = wires
        self.a = a
        self.L = L


    def solve_symbolic_S(self):
        a = self.a
        L = self.L
        W = self.wires

        letters = string.uppercase[7:][:3] # TODO
        ones = [var(p + '1') for p in letters]
        twos = [var(p + '2') for p in letters]

        wavefunctions  = [o * sin(k * x) + t * cos(k * x) for o, t in zip(ones, twos)]
        wavefunctionsd = [wf.derivative(x) for wf in wavefunctions]

        equations = []

        wire1, wire1_d = wavefunctions[0], wavefunctionsd[0]
        wire2, wire2_d = wavefunctions[2], wavefunctionsd[2]
        conn, conn_d = wavefunctions[1], wavefunctionsd[1]

        left = 0
        right = L

        equations.append(wire1(x=L) == 0)
        equations.append(wire2(x=L) == 0)

        equations.append(f_inc(x=left) == wire1(x=0))
        equations.append(wire1(x=0) == conn(x=0))

        equations.append(conn(x=L) == wire2(x=0))
        equations.append(wire2(x=0) == f_out(x=right))

        equations.append(-f_inc_d(x=left) + wire1_d(x=0) + conn_d(x=0) == 0) # TODO a
        equations.append(-conn_d(x=L) - wire2_d(x=0) + f_out_d(x=right) == 0) # TODO a

        solutions = solve(
            equations,
            B, C, *(ones + twos),
            solution_dict=True
        )

        Bs = solutions[0][B].full_simplify()
        Cs = solutions[0][C].full_simplify()
        SM = asymmetric_Smatrix(Bs, Cs)
        return SM

    def solve_symbolic(self):
        return self.solve_symbolic_S().det()


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


    def solve_symbolic_S(self):
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
        return SM

    def solve_symbolic(self):
        return self.solve_symbolic_S().det()


class IntervalLoopSolver(object):
    # wires: integer
    # a: symbolic/float
    def __init__(self, L=rvar('L')):
        super(IntervalLoopSolver, self).__init__()
        self.L = L


    def solve_symbolic_S(self):
        L = self.L


        letters = string.uppercase[7:][:2]
        ones = [var(p + '1') for p in letters]
        twos = [var(p + '2') for p in letters]

        wavefunctions  = [o * sin(k * x) + t * cos(k * x) for o, t in zip(ones, twos)]
        wavefunctionsd = [wf.derivative(x) for wf in wavefunctions]

        vf, vfd = wavefunctions[0], wavefunctionsd[0]
        wf, wfd = wavefunctions[1], wavefunctionsd[1]
        
        equations = [
            f_inc(x=0) == vf(x=0),
            vf(x=0)    == f_out(x=0),
            -f_inc_d(x=0) + vfd(x=0) + f_out_d(x=0) == 0,

            vf(x=L) == wf(x=0),
            vf(x=L) == wf(x=1),
            -vfd(x=L) + wfd(x=0) - wfd(x=1) == 0
        ]
        # for e in equations:
            # view_later(e)

        solutions = solve(
            equations,
            B, C, *(ones + twos),
            solution_dict=True
        )

        Bs = solutions[0][B].full_simplify()
        Cs = solutions[0][C].full_simplify()
        SM = asymmetric_Smatrix(Bs, Cs)
        return SM

    def solve_symbolic(self):
        return self.solve_symbolic_S().det()


# TODO look at eigenvalues
class IntervalSolver(object):
    # wires: integer
    # a: symbolic/float
    # b: symbolic/float
    def __init__(self, wires, a=rvar('a'), b=rvar('b'), L=rvar('L')):
        super(IntervalSolver, self).__init__()
        self.wires = wires
        self.a = a
        self.b = b
        self.L = L


    def solve_analytic(self, W=None):
        W = len(self.wires) if W is None else W # TODO ???
        a = self.a
        b = self.b
        L = self.L
        # to make formulas look pretty

        coeff = W**2 + 1
        # TODO 2 * k might be because the actual length is 2 * L instead of L :)
        ### CORRECT VERSION!!!!
        num = W * (a + b) * k * cos(L * k) + (a * b - coeff * k**2) * sin(L * k) + 2 * W * i * k**2 * cos(L * k) + i * (a + b) * k * sin(L * k)
        den = W * (a + b) * k * cos(L * k) + (a * b - coeff * k**2) * sin(L * k) - 2 * W * i * k**2 * cos(L * k) - i * (a + b) * k * sin(L * k)
        ###

        ### CORRECT VERSION WITH TAN !!!!
        # num = W * (a + b) * k  + (a * b - coeff * k**2) * tan(L * k) + 2 * W * i * k**2 + i * (a + b) * k * tan(L * k)
        # den = W * (a + b) * k  + (a * b - coeff * k**2) * tan(L * k) - 2 * W * i * k**2 - i * (a + b) * k * tan(L * k)
        ###

        # den = ((W + 1) * k + i)**2 - ((W - 1) * k - i)**2 * exp(2 * i * k)
        # num = - (W * 4)**2 * k**4 * exp(2 * i * k) \
        #     + (a * b - i * ((W - 1) * a + (W + 1) * b) * k - (W - 1) * (W + 1) * k**2 - (a * b + i * ((W + 1) * a + (W - 1) * b) * k - (W - 1) * (W + 1) * k**2) * exp(2 * i * k)) \
        #     * (a * b - i * ((W + 1) * a + (W - 1) * b) * k - (W + 1) * (W - 1) * k**2 - (a * b + i * ((W - 1) * a + (W + 1) * b) * k - (W + 1) * (W - 1) * k**2) * exp(2 * i * k))
        # num =  ((a * b - i * ((W - 1) * a + (W + 1) * b) * k - (W - 1) * (W + 1) * k**2 - (a * b + i * ((W + 1) * a + (W - 1) * b) * k - (W - 1) * (W + 1) * k**2) * exp(2 * i * k)) - 4 * W * k**2 * exp(i * k)) \
        #      * ((a * b - i * ((W + 1) * a + (W - 1) * b) * k - (W + 1) * (W - 1) * k**2 - (a * b + i * ((W - 1) * a + (W + 1) * b) * k - (W + 1) * (W - 1) * k**2) * exp(2 * i * k)) + 4 * W * k**2 * exp(i * k))
        return num / den

    def solve_symbolic_S(self):
        wires = self.wires
        a = self.a
        b = self.b
        L = self.L

        letters = string.uppercase[7:][:len(wires)]
        ones = [var(p + '1') for p in letters]
        twos = [var(p + '2') for p in letters]

        wavefunctions  = [o * sin(k * x) + t * cos(k * x) for o, t in zip(ones, twos)]
        # wavefunctions  = [o * exp(i * k * x) + t * exp(-i * k * x) for o, t in zip(ones, twos)]
        wavefunctionsd = [wf.derivative(x) for wf in wavefunctions]

        equations = []

        left = 0
        right = 0

        print(wires)
        for wlen, wf in zip(wires, wavefunctions):
            equations.append(f_inc(x=left) == wf(x=0))
            equations.append(wf(x=wlen)    == f_out(x=right))

        derl = -f_inc_d(x=left)  + sum(wfd(x=0)    for wfd       in wavefunctionsd)             == a * f_inc(x=left)
        derr =  f_out_d(x=right) - sum(wfd(x=wlen) for wfd, wlen in zip(wavefunctionsd, wires)) == b * f_out(x=right)

        equations.append(derl)
        equations.append(derr)

        pprint(equations)

        solutions = solve(
            equations,
            B, C, *(ones + twos),
            solution_dict=True
        )

        pprint(solutions[0])

        Bs = solutions[0][B].full_simplify()
        Cs = solutions[0][C].full_simplify()

        SM = asymmetric_Smatrix(Bs, Cs)
        return SM

    def solve_symbolic(self):
        return self.solve_symbolic_S().det()    


def interval_solver_nonuniform(wires, a, b):
    return IntervalSolver(wires, a=a, b=b, L=1)


def interval_solver_uniform(count, a, b):
    return interval_solver_nonuniform([1 for _ in range(count)], a, b)


# TODO look at eigenvalues
class FractalSolver(object):
    # wires: integer
    # a: symbolic/float
    # b: symbolic/float
    def __init__(self, wires, a=rvar('a'), L=rvar('L')):
        super(FractalSolver, self).__init__()
        self.wires = wires
        self.a = a
        self.L = L


    def solve_analytic(self):
        W = self.wires
        a = self.a
        L = self.L
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

    def solve_symbolic(self):
        a = self.a
        L = self.L
        W = self.wires

        letters = string.uppercase[7:][:W]
        ones = [var(p + '1') for p in letters]
        twos = [var(p + '2') for p in letters]

        wavefunctions  = [o * sin(k * x) + t * cos(k * x) for o, t in zip(ones, twos)]
        wavefunctionsd = [wf.derivative(x) for wf in wavefunctions]

        equations = []

        left = 0
        right = L * W

        all_wfs  = wavefunctions # + [f_out]
        all_wfds = wavefunctionsd # + [f_out_d]
        last = f_inc(x=left)
        lastd = f_inc_d(x=left)
        for i, cur, curd in zip(range(1, W + 1), all_wfs, all_wfds): #TODO i is overridden !!!!!!!!
            equations.append(last == cur(x=0))
            equations.append(-lastd + curd(x=0) == i * a * cur(x=0))
            last = cur(x=L)
            lastd = curd(x=L)
    
        equations.append(last == f_out(x=right))
        equations.append(-lastd + f_out_d(x=right) == (W + 1) * a * f_out(x=right))


        pprint(equations)

        solutions = solve(
            equations,
            B, C, *(ones + twos),
            solution_dict=True
        )

        Bs = solutions[0][B].full_simplify()
        Cs = solutions[0][C].full_simplify()
        SM = asymmetric_Smatrix(Bs, Cs)
        return SM.det()




class GridSolver(object):
    # wires: integer
    # a: symbolic/float
    # b: symbolic/float
    def __init__(self, wires, a=rvar('a'), L=rvar('L')):
        super(GridSolver, self).__init__()
        self.wires = wires
        self.a = a
        self.L = L

    def solve_analytic(self):
        a = self.a
        L = self.L
        W = self.wires
        # TODO a = 0
        # TODO L = 1

        if W == 1:
            num = exp(4 * i * k) + 9
            den = exp(4 * i * k) - 9
            return num / den


    def solve_symbolic_S(self):
        a = self.a
        L = self.L
        W = self.wires

        irange = range(W + 1)
        jrange = range(W + 1)

        from itertools import product

        vertices = set(product(irange, jrange))

        edges = set()

        for q, w in vertices:
            neighbours = [p for p in [(q - 1, w), (q, w - 1), (q + 1, w), (q, w + 1)] if p in vertices]
            for nb in neighbours:
                edge = frozenset([nb, (q, w)])
                edges.add(edge)

        def encode(edge):
            [(q, w), (e, r)] = list(edge)
            return "e" + ''.join([str(y) for y in [q, w, e, r]])

        def var1(edge):
            return var(encode(edge) + "_1")

        def var2(edge):
            return var(encode(edge) + "_2")

        print(edges)

        variables = [var1(e) for e in edges] + [var2(e) for e in edges]

        functions   = {e: var1(e) * exp(-i * k * x) + var2(e) * exp(i * k * x) for e in edges}
        functions_d = {e: functions[e].derivative(x) for e in edges} 

        equations = []
        for q, w in vertices:

            left  = (q - 1, w)
            right = (q + 1, w)
            up    = (q, w + 1)
            down  = (q, w - 1)

            def get_function(v):
                edge = frozenset([v, (q, w)])
                return functions[edge]

            def get_function_d(v):
                edge = frozenset([v, (q, w)])
                return functions_d[edge]

            values_at_vertex = []
            sum_at_vertex = []
            if left in vertices:
                print("Left: " + str(left))
                values_at_vertex.append(get_function(left)(x=L))
                sum_at_vertex.append(-get_function_d(left)(x=L))
            if right in vertices:
                values_at_vertex.append(get_function(right)(x=0))
                sum_at_vertex.append(get_function_d(right)(x=0))
            if down in vertices:
                print("Down: " + str(down))
                values_at_vertex.append(get_function(down)(x=L))
                sum_at_vertex.append(-get_function_d(down)(x=L))
            if up in vertices:
                values_at_vertex.append(get_function(up)(x=0))
                sum_at_vertex.append(get_function_d(up)(x=0))

            if (q, w) == (0, 0):
                values_at_vertex.append(f_inc(x=0))
                sum_at_vertex.append(-f_inc_d(x=0))

            if (q, w) == (W, W):
                values_at_vertex.append(f_out(x=0))
                sum_at_vertex.append(f_out_d(x=0))

            for e1, e2 in zip(values_at_vertex, values_at_vertex[1:]):
                equations.append(e1 == e2)

            equations.append(sum(sum_at_vertex) == a * values_at_vertex[0])

        # pprint(equations)
        # pprint(len(equations))
        # pprint(len(variables))

        solutions = solve(
            equations,
            B, C, *variables,
            solution_dict=True
        )

        pprint(solutions)

        Bs = solutions[0][B].full_simplify()
        Cs = solutions[0][C].full_simplify()
        SM = asymmetric_Smatrix(Bs, Cs)
        return SM

    def solve_symbolic(self):
        return self.solve_symbolic_S().det()


class TriangleSolver(object):
    # wires: integer
    # a: symbolic/float
    # b: symbolic/float
    def __init__(self, a=rvar('a'), L=rvar('L')):
        super(TriangleSolver, self).__init__()
        self.a = a
        self.L = L

    def solve_symbolic(self):
        a = self.a
        L = self.L

        letters = string.uppercase[7:10]
        ones = [var(p + '1') for p in letters]
        twos = [var(p + '2') for p in letters]


        wavefunctions  = [o * sin(k * x) + t * cos(k * x) for o, t in zip(ones, twos)]
        # wavefunctions  = [o * exp(i * k * x) + t * exp(-i * k * x) for o, t in zip(ones, twos)]
        wavefunctionsd = [wf.derivative(x) for wf in wavefunctions]

        equations = []

        wf0 , wf1 , wf2  = wavefunctions
        wf0d, wf1d, wf2d = wavefunctionsd

        equations.extend([
            f_inc(x=0) == wf0(x=0),
            wf0(x=1)   == f_out(x=0),

            f_inc(x=0) == wf1(x=0),
            wf1(x=1)   == wf2(x=0),
            wf2(x=L)   == f_out(x=0),

            -f_inc_d(x=0) + wf0d(x=0) + wf1d(x=0) == 0,
             f_out_d(x=0) - wf0d(x=1) - wf2d(x=L) == 0,
            -wf1d(x=1) + wf2d(x=0) == 0,
        ])

        print(len(equations))

        solutions = solve(
            equations,
            B, C, *(ones + twos),
            solution_dict=True
        )

        pprint(solutions[0])

        Bs = solutions[0][B].full_simplify()
        Cs = solutions[0][C].full_simplify()

        SM = asymmetric_Smatrix(Bs, Cs)
        return SM.det()

def icayley(x):
    return i * (1 + x) / (1 - x)

DPI = 200

# def plot_all(Sdet, suffix="", rrange=(10300, 10400), irange=(10.0, 10.2), points=1500):
#     print(n(abs(Sdet(rrange[0] + irange[0] * i)) ** 0.02))
#     # complex_plot(Sdet, rrange, irange, plot_points=points).save('plot{}.png'.format(suffix), figsize=[12, 2])
#     plot_abs = complex_plot(abs(Sdet) ** 0.04, rrange, irange, plot_points=points)
#     # l = line([(rrange[0], 5), (rrange[1], 5)], rgbcolor=(1, 0, 0))# complex_plot(lambda q: q.real() + 5 * i, rrange, irange, color='green')
#     l = plot(lambda x: 1.1 * ln(x), rrange, color='red')
#     t = var('t')
#     # l = parametric_plot((t, 5), (t, rrange[0], rrange[1]), color='red')
#     sum([plot_abs, l]).save('plot_abs{}.png'.format(suffix), figsize=[12, 2])
#     # complex_plot(ln(abs(Sdet) * 10), rrange, irange, plot_points=points).save('plot_ln{}.png'.format(suffix), figsize=[12, 2])
#     # unit_circle = circle((0, 0), 1)
#     # (complex_plot(ln(abs(Sdet(k=cayley(k)))), (-1, 1), (-1, 1), plot_points=points) + unit_circle).save('plot_circle.png', dpi=DPI)

def plot_all(Sdet, suffix="", rrange=(-2, 20), irange=(0, 2), points=1500):
    # print(n(abs(Sdet(rrange[0] + irange[0] * i)) ** 0.02))
    # complex_plot(Sdet, rrange, irange, plot_points=points).save('plot{}.png'.format(suffix), figsize=[12, 2])
    plot_abs = complex_plot(abs(Sdet), rrange, irange, plot_points=points)
    # l = line([(rrange[0], 5), (rrange[1], 5)], rgbcolor=(1, 0, 0))# complex_plot(lambda q: q.real() + 5 * i, rrange, irange, color='green')
    ll = [

        # plot(lambda x: ln(x) / 4, rrange, color='red'),
        # plot(lambda x: 1.05 / 4 * ln(x), rrange, color='green'),
        # plot(lambda x: 1.065 / 4 * ln(x), rrange, color='yellow'),
        # plot(lambda x: 1.10 / 4 * ln(x), rrange, color='blue'),
    ]
    t = var('t')
    # l = parametric_plot((t, 5), (t, rrange[0], rrange[1]), color='red')
    sum([plot_abs] + ll).save('abs{}.png'.format(suffix), figsize=[12, 2])
    # complex_plot(ln(abs(Sdet)), rrange, irange, plot_points=points).save('plot_ln{}.png'.format(suffix), figsize=[12, 2])
    # unit_circle = circle((0, 0), 1)
    # (complex_plot(ln(abs(Sdet(k=cayley(k)))), (-1, 1), (-1, 1), plot_points=points) + unit_circle).save('plot_circle.png', dpi=DPI)



# S = solve_popov_analytic(a_val=10)
# S = solve_double_loop(a_val=5, b_val=5)

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
            except RuntimeError as err:
                print(err)



a = var('a', domain='real')
b = var('b', domain='real')


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

# a = 1
a = rvar('a')
b = rvar('b')
L = rvar('L')


def contour_integral_analysis(expr):
    # denom(k) = cos(k) - 2 * i * sin(k)
    # denom(k) = (cos(k) + 2 * i * sin(k)) / (cos(k) - 2 * i * sin(k))
    ig(k) = ln(abs(expr(k=k))) / (k - 1) ** 2
    for R in [1, 5, 10, 20, 50, 100, 200, 400, 500, 1000, 10000]:
        pig(t) = ig(k = R * exp(i * t) + R * i) * R * i * exp(i * t)
        rpart, _  = numerical_integral(lambda q: pig(t=q).real(), 0, pi)
        impart, _ = numerical_integral(lambda q: pig(t=q).imag(), 0, pi)
        print("{} {}".format(rpart, impart))

def print_latex(expr):
    print("-------------------")
    print(latex(expr))
    print("-------------------")


def simple_case():
    W = var('W', domain='integer')
    u = rvar('u')
    v = rvar('v')

    a = 0
    b = 0
    L = 1
    num = (W**2 + 1) * sin(k) - 2 * i * W * cos(k)
    num *= num.conjugate()
    den = (W**2 + 1) * sin(k) + 2 * i * W * cos(k)
    den *= den.conjugate()
    num = (num(k = u + i * v)).full_simplify()
    den = den(k = u + i * v).full_simplify()
    Sdet = num / den
    print_latex(Sdet)
    view_later(Sdet)

# simple_case()
# view_all()
# sys.exit(0)

# a = 0
# b = 0
def test_aaa():
    L = 1
    W = var('W', domain='integer')
    solver = IntervalSolver(W, a=a, b=b, L=1)
    Sa = solver.solve_analytic()
    # view_later(Sa)
    print(print_latex(Sa))
    # Sa =
    u = rvar('u')
    v = rvar('v')
    Sa = Sa(k=u + v * i).full_simplify()
    print("Sa = " + str(Sa))
    # qqq = Sa.derivative(v)
    # print("Derivative = " + str())
    # for vv in range(0, 10):
    #     print(n(Sa(v=vv)))
    # plot(Sa, (v, 0, 10), plot_points=10).save('plot.png')
    view_later(Sa)
    print_latex(Sa)
    # view_later(solver.solve_analytic())

# test_aaa()
# view_all()
# sys.exit(0)

def try_same_length_analytic():
    L = 1
    solver = IntervalSolver([], a=a, b=b, L=L)
    Sa = solver.solve_analytic(W=var('W'))
    view_later(Sa)
    for W in [2, 3, 4, 5, 6]:
        solver = interval_solver_uniform(W, a, b)
        # S = solver.solve_symbolic()
        Sa = solver.solve_analytic()
        view_later(Sa)

    view_all()

# # a = 1
# # b = -a

def try_different_length():
    a = 0
    b = 0
    L = 1
    for W in [2, 3, 4]:
        view_later("Wires = " + str(W))
        # solver =  interval_solver_uniform(W, a, b)
        solver = interval_solver_nonuniform([i for i in range(1, W + 1)], a, b)

        S = solver.solve_symbolic_S()
        Sm = S.det() # .full_simplify()
        # view_later(Sm)

        Sa = solver.solve_analytic()
        plot_all(Sm, suffix="_nonuniform_" + str(W))

def try_different_length_2():
    a = 1
    b = 0
    L = 1
    for L2 in [1, 1.1, 1.3, 1.5, 2.0]:
        solver = interval_solver_nonuniform([L, L2], a, b)
        S = solver.solve_symbolic()
        plot_all(S, suffix="_two_wires_different_a1_" + "{:.1f}".format(float(L2)))


def try_different_length_a():
    # a = 1
    b = 0
    L = 1
    for a in arange(100, 500, 100):
        solver = interval_solver_nonuniform([1, 1], a, b)
        S = solver.solve_symbolic()
        plot_all(S, suffix="_two_wires__a_" + "{:.1f}".format(float(a)))
 
def try_different_length_symbolic():
    solver = interval_solver_nonuniform([1, 2], a, b)
    S = solver.solve_symbolic()
    view_later(S)
    view_all()

def try_same_length_fracional():
    a = 0
    b = 0
    L = 1
    for W in arange(-1, 1, 0.2):
        coeff = W**2 + 1
        # TODO 2 * k might be because the actual length is 2 * L instead of L :)
        ### CORRECT VERSION!!!!
        num = W * (a + b) * k * cos(L * k) + (a * b - coeff * k**2) * sin(L * k) + 2 * W * i * k**2 * cos(L * k) + i * (a + b) * k * sin(L * k)
        den = W * (a + b) * k * cos(L * k) + (a * b - coeff * k**2) * sin(L * k) - 2 * W * i * k**2 * cos(L * k) - i * (a + b) * k * sin(L * k)
        plot_all(num/den, suffix="_fractional_" + "{:.1f}".format(float(W)), rrange=(-2, 30), irange=(-4, 4), points=1500)

def two_wires_convergence():
    a = 0
    b = 0
    for L in [1, 0.99, 0.95, 0.9, 0.75, 0.5, 0.25, 0.10, 0.05, 0.01, 0.005]:
        solver = interval_solver_nonuniform([1, L], a, b)
        S = solver.solve_symbolic()
        plot_all(S, suffix="_convergence_" + "{:.3f}".format(float(L)), rrange=(-2, 60), irange=(-7, 7), points=1500)

def one_wire_popov():
    a = 0
    solver = PopovSolver(1, a=a, L=1)
    S = solver.solve_symbolic()
    Sa = solver.solve_analytic_2()
    Saa = solver.solve_analytic()
    # sanity_checks(S)
    test_matrices(S, Sa)
    test_matrices(Sa, Saa)
    view_later(Saa)
    view_all()
    # view_later(S)
    # view_all()
    # plot_all(Saa, suffix="_what", rrange=(-60, 60), irange=(-2, 2), points=1500)
    # Sa = solver.solve_analytic()
    # test_matrices(S, Sa)
    # view_later(Sa)
    # view_all()

def one_wire_loop():
    a = 2
    solver = LoopSolver(1, a=a, L=1)
    S = solver.solve_symbolic_S()
    sanity_checks(S)

    Sd = S.det()
    Sa = solver.solve_analytic()
    test_matrices(Sd, Sa)
    # plot_all(Sa, suffix="_loop_1", rrange=(-2, 100), irange=(-20, 20), points=1500)
    view_later(Sa)
    view_all()

def try_grid():
    a = 1
    for W in [1]:# , 2, 3, 4]:
        solver = GridSolver(W, a=a, L=1)
        S = solver.solve_symbolic()
        plot_all(S, suffix="_grid_1" + str(W), rrange=(-2, 40), irange=(-1.5, 1.5), points=1500)            

def try_triangle():
    for L in [0.11, 0.13, 0.15, 0.17, 0.20]:
        solver = TriangleSolver(a=0, L=L)
        S = solver.solve_symbolic()
        plot_all(S, suffix="_triangle_" + "{:.2f}".format(float(L)), rrange=(-2, 60), irange=(-7, 7), points=1500)

def try_interval_loop():
    for L in [1]: #.0, 0.1, 0.05, 0.01, 0.005, 0.001]:
        solver = IntervalLoopSolver(L=L)
        S = solver.solve_symbolic_S()
    # Sa = solver.solve_analytic()
    # sanity_checks(S)
        plot_all(S.det(), suffix="_interval_loop_{:.3f}".format(float(L)), rrange=(-2, 60), irange=(-20, 20), points=500)
    # test_matrices(S, Sa)
    # view_later(Sa)
    # view_all()

# try_same_length_fracional()
# two_wires_convergence()
# one_wire_loop()
# try_grid()
# try_triangle()
# try_same_length_analytic()
one_wire_loop()
# one_wire_popov()

# try_different_length_a()
    # denn = (a * b - (W + 1) * i * (a + b) * k - (W + 1)**2 * k**2) - (a * b + (W - 1) * i * (a + b) * k - (W - 1)**2 * k**2) * exp(2 * i * k)
    # nomm = denn(k=conjugate(k))
    # plot_all(nomm / denn, suffix="_intervam_" + str(W))

    # sanity_checks(Sa)
    # test_matrices(Sm, Sa)
    # contour_integral_analysis(Sm) #  * exp(2 * i * k))
# view_all()sa
