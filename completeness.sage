reset()

R, T = var('R T')
k = var('k')
assume(k, 'imaginary')

f_inc(x) = exp(i * k * x) + R * exp(-i * k * x)
f_out(x) = T * exp(i * k * x)

f_inc_d = f_inc.derivative(x)
f_out_d = f_out.derivative(x)

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

    return Rs(a = 2, L = 1), Ts(a = 2, L = 1)

def solve_interval():
    A, B = var('A B')
    a, L = var('a L')
    assume(a, 'real')
    assume(L, 'real')

    f2(x) = A * exp(i * k * x) + B * exp(-i * k * x)
    f2d = f2.derivative(x)

    solutions = solve([f_inc(-L) == f2(-L), f2(L) == f_out(L), -f_inc_d(-L) + f2d(-L) == a * f_inc(-L), -f2d(L) + f_out_d(L) == a * f_out(L)], R, T, A, B, solution_dict=True)

    Rs = solutions[0][R].full_simplify()
    Ts = solutions[0][T].full_simplify()

    aa = 1
    LL = 1
    return Rs(a = aa, L = LL), Ts(a = aa, L = LL)

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

    return Rs(a = 1, L = 1), Ts(a = 1, L = 1)

def solve_double_loop():
    A, B, C, D = var('A B C D')
    a, L = var('a L')
    assume(a, 'real')
    assume(L, 'real')
    # a = 1
    # L = 1

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
    show(equations)

    solutions = solve(equations, R, T, A, B, C, D, solution_dict=True)
    Rs = solutions[0][R].full_simplify()
    Ts = solutions[0][T].full_simplify()

    return Rs(a = 1, L = 1), Ts(a = 1, L = 1)


Rs, Ts = solve_double_loop()
show("R = ", Rs)
show("T = ", Ts)

S = matrix([
    [Ts, Rs], 
    [Rs, Ts]
])
Sdet = S.det()

show("Det = ", Sdet)

rrange = (-100, 100)
irange = (-5, 5)
points = 400


def cayley(x):
    return i * (x + 1) / (x - 1)

# Sdet = -(cos(k) + 2 * i * sin(k)) / (cos(k) - 2 * i * sin(k))

print(n(Sdet(k=0.1)))

DPI = 200

def plot_all(Sdet):
    complex_plot(abs(Sdet), rrange, irange, plot_points=points).save('plot.png', dpi=DPI)
    complex_plot(ln(abs(Sdet)), rrange, irange, plot_points=points).save('plot_ln.png', dpi=DPI)
    unit_circle = circle((0, 0), 1)
    (complex_plot(ln(Sdet(k=cayley(k))), (-1, 1), (-1, 1), plot_points=points) + unit_circle).save('plot_circle.png', dpi=DPI)
    # t = var('t')
    # R = 100
    # ft = ln(Sdet(k=R * exp(i * t) + R)) * R * i * exp(i * t)
    # print(n(ln(Sdet(k=1))))
    # print(n(ln(Sdet(k=i))))
    # print(n(ft(t=0)))
    # print(n(ft(t=pi/2)))
    # pieces = 100
    # print(sum([R * 2 * pi / pieces * n(ft(t = q / pieces * 2 * pi + 0.01)) for q in range(pieces)]))
    # plot(ft.real(), (t, 0, pi), plot_points=10).save('pppp.png')
    # show(ft.real().integrate(t, 0, pi))

plot_all(Sdet)


delta = 0.1

# ll = 0.4
# pieces = 50
# for zz in range(-pieces / 2, pieces / 2):
#     kk = pi + zz * ll / pieces
#     print(n(kk))
#     print(n(ln(Sdet(k=kk))))

# for qq in range(-100, 100):
#    print(n(ln(Sdet(k=qq * pi/2 + delta)) + ln(Sdet(k=qq * pi/2 - delta))))

# ff = ln(Sdet(k))
# ff2 = ln(Sdet(k)) * 1 / (k - 1) ** 2

# print(n(ff(k = 1000 * exp(i * 19 * pi / 20))))
# show(integrate())


# for qq in range(-1000, 1000):
#     print(n(ln(Sdet(k=100 * i + qq))))

# print(n(ff(k=(1 * i) * 100000 - 0.1))))
# print(n(ff(k=100 + 2 * i)))
# print(n(ff(k=-100 + 2 * i)))


# show((ln(S.det())).maxima_methods().residue(z, pi + 0.5 * i * ln(3)))
# show(n(S.det()(z = pi+ 0.5 * i * ln(3))))
# show(complex_plot((i * (exp(i * x) + 1) / (exp(i * x) - 1)).abs(), 0, 2 * pi), ymin=0, ymax=3)
# show(S(i * (k + 1) / (k - 1)).det().full_simplify())

# show(complex_plot(ln(S.det()) + 2 * pi * i, (-10, 10), (-3, 3), plot_points=points))
# show(complex_plot(ln(S(i * (k + 1) / (k - 1)).det()), rrange, irange, plot_points=points) + circle((0, 0), 1), aspect_ratio=1.0) #  + plot(10 * Ts.abs(), rrange))# , aspect_ratio=20.0)