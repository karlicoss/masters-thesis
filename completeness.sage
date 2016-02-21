reset()

from sage.misc.viewer import viewer
viewer.pdf_viewer("atril")


_view_later = []
def view_later(formula):
    _view_later.append(formula)

def view_all():
    view(_view_later, tightpage=True)


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

    return SM(a=a_val, L=L_val)

def solve_interval():
    Q1, Q2 = var('C1 C2')
    a, b, L = var('a b L')
    assume(a, 'real')
    assume(b, 'real')
    assume(L, 'real')

    f2(x) = Q1 * exp(i * k * x) + Q2 * exp(-i * k * x)
    f2d = f2.derivative(x)

    solutions = solve(
        [
            f_inc(x=-L) == f2(x=-L), 
            f2(x=L) == f_out(x=L), 
            -f_inc_d(x=-L) + f2d(x=-L) == a * f_inc(x=-L), 
            -f2d(x=L) + f_out_d(x=L) == b * f_out(x=L)
        ],
        B, C, Q1, Q2, 
        solution_dict=True
    )


    Bs = solutions[0][B].full_simplify()
    Cs = solutions[0][C].full_simplify()
    SM = asymmetric_Smatrix(Bs, Cs)
    aa = -10
    bb = 10
    LL = 1
    return SM(a = aa, b = bb, L = LL)

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

def draw_double_loop():
    fontsize = 26
    length = 2
    radius = 0.5


    edge_inc = line([(-length, 0), (-radius,0)])
    edge_inc_inf = line([(-length - 0.2, 0), (-length, 0)], linestyle='dotted')
    edge_out = line([(radius, 0), (length, 0)])
    edge_out_inf = line([(length, 0), (length + 0.2, 0)], linestyle='dotted')
    edge_circle = circle((0, 0), radius)
    vertex_inc = point((-radius, 0), size=20)
    vertex_out = point((radius, 0), size=20)
    vertex_inc_label = text("$V_{inc}$", (-radius - 0.2, 0.1), fontsize=fontsize)
    vertex_out_label = text("$V_{out}$", (radius + 0.2, 0.1), fontsize=fontsize)
    edge_inc_label = text("$\Omega_{inc}$", (-(length - 0.3) - 0.2, 0.1), fontsize=fontsize)
    edge_out_label = text("$\Omega_{out}$", ((length - 0.3) + 0.2, 0.1), fontsize=fontsize)
    edge_up_label = text("$\Omega_1$", (0, radius + 0.1), fontsize=fontsize)
    edge_down_label = text("$\Omega_2$", (0, -radius - 0.15), fontsize=fontsize)
    objects = [
        edge_inc, edge_inc_inf, edge_inc_label,
        edge_circle, edge_up_label, edge_down_label,
        edge_out, edge_out_inf, edge_out_label,
        vertex_inc, vertex_inc_label,
        vertex_out, vertex_out_label,
    ]
    sum(objects).save('drawing.png', axes=False, figsize=12)

def solve_double_loop(a_val=1, b_val=-1, L_val=1):
    Q1, Q2, W1, W2 = var('Q1 Q2 W1 W2')
    a, b, L = var('a b L')
    assume(a, 'real')
    assume(b, 'real')
    assume(L, 'real')

    # TODO: for some reason, doesn't work, solve function stucks :(
    # f2(x) = Q1 * exp(-i * k * x) + Q2 * exp(-i * k * x)
    # f3(x) = W1 * exp(-i * k * x) + W2 * exp(-i * k * x)
    fQ(x) = Q1 * sin(k * x) + Q2 * cos(k * x)
    fW(x) = W1 * sin(k * x) + W2 * cos(k * x)
    
    fQd = fQ.derivative(x)
    fWd = fW.derivative(x)
    
    equations = [
        f_inc(x=0) == fQ(x=-L), 
        fQ(x=-L) == fW(x=-L), 
        fQ(x=L) == fW(x=L),
        fW(x=L) == f_out(x=0),
        -f_inc_d(x=0) + fQd(x=-L) + fWd(x=-L) == a * f_inc(x=0), 
        f_out_d(x=0) - fQd(x=L) - fWd(x=L) == a * f_out(x=0) # TODO FIXME TO b
    ]
    show(equations)

    solutions = solve(equations, B, C, Q1, Q2, W1, W2, solution_dict=True)
    Bs = solutions[0][B].full_simplify()
    Cs = solutions[0][C].full_simplify()

    view_later(Bs)
    view_later(Cs)
    SM = asymmetric_Smatrix(Bs, Cs)
    return SM(a=a_val, b=b_val, L=L_val)

def solve_double_loop_analytic(a_val=1):
    a = a_val
    b = a_val
    L = 1 # TODO L is ignored for now
    rp = (a**2 - 5 * k**2) * cos(k) * sin(k) - 4 * a * k * sin(k)**2 + 2 * a * k
    view_later(rp)
    ip = 2 * a * k * cos(k) * sin(k) - 4 * k**2 * sin(k)**2 + 2 * k**2
    view_later(ip)
    Sd = (rp + i * ip) / (rp - i * ip)
    return Sd


# interesting at a = 0
# S = solve_interval()


# sanity_checks(S)

rrange = (-10, 10)
irange = (-4, 4)
points = 500


def icayley(x):
    return i * (1 + x) / (1 - x)

DPI = 200

def plot_all(Sdet, suffix=""):
    complex_plot(Sdet, rrange, irange, plot_points=points).save('plot{}.png'.format(suffix), dpi=DPI)
    complex_plot(ln(abs(Sdet)), rrange, irange, plot_points=points).save('plot_ln{}.png'.format(suffix), dpi=DPI)
    # unit_circle = circle((0, 0), 1)
    # (complex_plot(ln(abs(Sdet(k=cayley(k)))), (-1, 1), (-1, 1), plot_points=points) + unit_circle).save('plot_circle.png', dpi=DPI)

# S = solve_double_loop(a_val=-1, b_val=-1)
# plot_all(S)

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
y = var('y')
# solve_double_loop(a_val=1)
# Sd = solve_double_loop_analytic(a_val=y)
# check_zero(Sd)
# plot_bound_energies()

# plot_all(Sd)

# view_later(Sd)
# view_all()
draw_double_loop()