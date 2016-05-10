def cvar(name):
    return var(name, domain='complex')


def rvar(name):
    return var(name, domain='real')


def pvar(name):
    return var(name, domain='positive')


def icayley(z):
    return i * (1 + z) / (1 - z)


def cayley(z):
    return (z - i) / (z + i)


def changevar(f, eqn, newvar):
    dx = diff(eqn.rhs(), newvar)
    # print(f)
    # print(f.subs(eqn))
    # print(dx)
    return f.subs(eqn) * dx


def complex_integral(expr, ff, tt):
    rpart, _ = numerical_integral(lambda q: expr(q).real().n(), ff, tt)
    impart, _ = numerical_integral(lambda q: expr(q).imag().n(), ff, tt)
    return rpart + i * impart


def estimate():
    k = cvar('k')
    t = rvar('t')
    R = pvar('R')
    # f = abs(ln(abs((cos(k) + (1 / (2 * k) + i) * sin(k) - 1) / (cos(k) + (1 / (2 * k) - i) * sin(k) - 1))))
    f = ln(k.imag())  # abs(ln(abs((cos(k) + (1 / (2 * k) + i) * sin(k) - 1) / (cos(k) + (1 / (2 * k) - i) * sin(k) - 1))))
    ig = f * diff(cayley(k), k)
    pig = changevar(ig, k == R * exp(i * t) + R * i, t).full_simplify()
    # for rr in [1, 2, 5, 10, 15, 20, 25, 30]:
    #     eee = pig(R=rr).simplify()
    #     print(complex_integral(eee, 0, 2 * pi))
    # print(eee(5 + i).real().n())
    print(pig.integrate(t, 0, 2 * pi))


estimate()

# contour_integral_analysis((k * cos(k) + a * sin(k) + 2 * i * k * sin(k)) / (k * cos(k) + a * sin(k) - 2 * i * k * sin(k)))

# 0 <= abs(f(k)) <= 1 ieverywhere in the upper complex plane, so ln(abs()) <= 0 everywhere
# therefore we can replace the function f(k) to f(k).imag(). which results in (3 - e^(2 * y)) (3 e^(2y) - 1)
# abs((1 + 2 i tg(k)) / (1 - 2 i tg(k))) < 1
# (1 + 4 tg(k) tg(k*) + 2 i (tg(k) - tg(k*))) < (1 + 4 tg(k) tg(k*) - 2 i (tg(k) - tg(k*))))
# true if (tg(k) - tg(k*)) is positive imaginary number (which indeed holds), expand as cosh and sinh


# looks same when a is not zero


# contour_integral_analysis(cos(k) + 2 * i * sin(k))
# contour_integral_analysis(cos(k) - 2 * i * sin(k))


# show(region_plot(lambda x, y: abs(x ** 2 + y ** 2) < 3, (x, -10, 10), (y, -10, 10)))

# region_plot(lambda x, y: abs(ln(abs(f(x + i * y))) / (x + i * y) ** 2) > 1, (x, -5, 5), (y, -5, 5))
