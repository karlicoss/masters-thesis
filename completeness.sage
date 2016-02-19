reset()

R, T = var('R T')
k = var('k')
assume(k, 'imaginary')

f_inc(x) = exp(i * k * x) + R * exp(-i * k * x)
f_out(x) = T * exp(i * k * x)

f_inc_d = f_inc.derivative(x)
f_out_d = f_out.derivative(x)


def symmetric_Smatrix(R, T):
    return matrix([
        [R, T], 
        [T, R]
    ])


# Checks the matrix that it is a valid S-matrix
def sanity_checks(S): # variable: k
    for v in [-10, -1, 1, 10]:
        show("Unitary on {}:".format(v))
        show(n(S(k=v) * S(k=v).H))
        
    for v in [1 + i, 1 - i, -1 - i, -1 + i]:
        show("S-matrix property on {}:".format(v))
        show(n(S(k=v) * S(k=v.conjugate()).H))


# Simple form:
# Sdet = -(cos(k) + 2 * i * sin(k)) / (cos(k) - 2 * i * sin(k))
def solve_popov():
    P = var('P')
    a, L = var('a L')
    assume(a, 'real')
    assume(L, 'real')

    f2(x) = P * sin(k * x)
    f2d = f2.derivative(x)

    solutions = solve([f_inc(0) == f2(L), f2(L) == f_out(0), -f_inc_d(0) - f2d(L) + f_out_d(0) == a * f_inc(0)], R, T, P, solution_dict=True)

    Rs = solutions[0][R].full_simplify()
    Ts = solutions[0][T].full_simplify()

    aa = 0
    LL = 1
    return symmetric_Smatrix(Rs, Ts)(a = aa, L = LL)

def solve_interval():
    A, B = var('A B')
    a, b, L = var('a b L')
    assume(a, 'real')
    assume(b, 'real')
    assume(L, 'real')

    f2(x) = A * exp(i * k * x) + B * exp(-i * k * x)
    f2d = f2.derivative(x)

    solutions = solve(
        [
            f_inc(x=-L) == f2(x=-L), 
            f2(x=L) == f_out(x=L), 
            -f_inc_d(x=-L) + f2d(x=-L) == a * f_inc(x=-L), 
            -f2d(x=L) + f_out_d(x=L) == a * f_out(x=L)
        ],
        R, T, A, B, 
        solution_dict=True
    )

    Rs = solutions[0][R].full_simplify()
    Ts = solutions[0][T].full_simplify()

    aa = 0
    LL = 1
    return symmetric_Smatrix(Rs, Ts)(a = aa, L = LL)

def solve_loop():
    A, B = var('A B')
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

def solve_double_loop():
    A, B, C, D = var('A B C D')
    a, L = var('a L')
    assume(a, 'real')
    assume(L, 'real')

    # TODO: for some reason, doesn't work, solve function stucks :(
    # f2(x) = A * exp(-i * k * x) + B * exp(-i * k * x)
    # f3(x) = C * exp(-i * k * x) + D * exp(-i * k * x)
    f2(x) = A * sin(k * x) + B * cos(k * x)
    f3(x) = C * sin(k * x) + D * cos(k * x)
    
    f2d = f2.derivative(x)
    f3d = f3.derivative(x)
    
    equations = [
        f_inc(x=-L) == f2(x=-L), 
        f2(x=-L) == f3(x=-L), 
        f2(x=L) == f3(x=L),
        f3(x=L) == f_out(x=L),
        -f_inc_d(x=-L) + f2d(x=-L) + f3d(x=-L) == a * f_inc(x=-L), 
        f_out_d(x=L) - f2d(x=L) - f3d(x=L) == a * f_out(x=L)
    ]
    # show(equations)

    solutions = solve(equations, R, T, A, B, C, D, solution_dict=True)
    Rs = solutions[0][R].full_simplify()
    Ts = solutions[0][T].full_simplify()

    aa = 0
    LL = 1
    return Rs(a = aa, L = LL), Ts(a = aa, L = LL)


# interesting at a = 0
S = solve_popov()
sanity_checks(S)


show(S)
Sdet = S.det()

show("Det = ", Sdet)

rrange = (-30, 30)
irange = (-3, 3)
points = 500


def cayley(x):
    return i * (x + 1) / (x - 1)

DPI = 200

def plot_all(Sdet):
    complex_plot(abs(Sdet), rrange, irange, plot_points=points).save('plot.png', dpi=DPI)
    complex_plot(ln(abs(Sdet)), rrange, irange, plot_points=points).save('plot_ln.png', dpi=DPI)
    # unit_circle = circle((0, 0), 1)
    # (complex_plot(ln(abs(Sdet(k=cayley(k)))), (-1, 1), (-1, 1), plot_points=points) + unit_circle).save('plot_circle.png', dpi=DPI)

plot_all(Sdet)

# show(plot(abs(Sdet(z = x * i)), (x, -10, 10)), ymin=0, ymax=2)
