from numpy import arange

ff(k) = (cos(k) + 2 * i * sin(k)) / (cos(k) - 2 * i * sin(k))

def f(k):
    return ff(k)

x, y = var('x y')


def contour_integral_analysis(expr):
    # denom(k) = cos(k) - 2 * i * sin(k)
    # denom(k) = (cos(k) + 2 * i * sin(k)) / (cos(k) - 2 * i * sin(k))
    ig(k) = ln(abs(expr(k=k))) / (k - 1) ** 2
    for R in [1, 5, 10, 20, 50, 100, 200, 400, 500, 1000, 10000]:
        pig(t) = ig(k = R * exp(i * t) + R * i) * R * i * exp(i * t)
        show(numerical_integral(lambda q: pig(t=q).real(), 0, pi))
        show(numerical_integral(lambda q: pig(t=q).imag(), 0, pi))

a = 5
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