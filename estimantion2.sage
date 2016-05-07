W = 2
a = 0
b = -a
L = 1

def rvar(name):
    return var(name, domain='real')


W = var('W', domain='integer')
a = rvar('a')
b = rvar('b')
L = 1

k = var('k', domain='complex')

coeff = W**2 + 1
num = W * (a + b) * k * cos(L * k) + (a * b - coeff * k**2) * sin(L * k) + 2 * W * i * k**2 * cos(L * k) + i * (a + b) * k * sin(L * k)
den = W * (a + b) * k * cos(L * k) + (a * b - coeff * k**2) * sin(L * k) - 2 * W * i * k**2 * cos(L * k) - i * (a + b) * k * sin(L * k)


print(num)
print(den)
# res = abs(num / den) ** 2
# res2(k) = (num(k=k) * conjugate(num(k=k))) / (num(k=conjugate(k)) * conjugate(num(k=conjugate(k))))
# res = abs(num / den) ** 2
# res2(k) = abs(num(k=k) / (conjugate(num(k=k)))) ** 2
res = num / den
# res2(k) = num(k=k) / num(k=conjugate(k))

def plot_all(ff, suffix="", rrange=(0, 20), irange=(-3, 3), points=300):
    # print(n(abs(Sdet(rrange[0] + irange[0] * i)) ** 0.02))
    # complex_plot(Sdet, rrange, irange, plot_points=points).save('plot{}.png'.format(suffix), figsize=[12, 2])
    plot_abs = complex_plot(ff, rrange, irange, plot_points=points)
    # l = line([(rrange[0], 5), (rrange[1], 5)], rgbcolor=(1, 0, 0))# complex_plot(lambda q: q.real() + 5 * i, rrange, irange, color='green')
    ll = [

        # plot(lambda x: ln(x) / 4, rrange, color='red'),
        # plot(lambda x: 1.05 / 4 * ln(x), rrange, color='green'),
        # plot(lambda x: 1.065 / 4 * ln(x), rrange, color='yellow'),
        # plot(lambda x: 1.10 / 4 * ln(x), rrange, color='blue'),
    ]
    # t = var('t')
    # l = parametric_plot((t, 5), (t, rrange[0], rrange[1]), color='red')
    sum([plot_abs] + ll).save('estimation' + suffix + '.png'.format(suffix), figsize=[12, 5])
    # complex_plot(ln(abs(Sdet)), rrange, irange, plot_points=points).save('plot_ln{}.png'.format(suffix), figsize=[12, 2])
    # unit_circle = circle((0, 0), 1)
    # (complex_plot(ln(abs(Sdet(k=cayley(k)))), (-1, 1), (-1, 1), plot_points=points) + unit_circle).save('plot_circle.png', dpi=DPI)


def contour_integral_analysis(expr):
    # denom(k) = cos(k) - 2 * i * sin(k)
    # denom(k) = (cos(k) + 2 * i * sin(k)) / (cos(k) - 2 * i * sin(k))
    ig(k) = ln(abs(expr(k=k))) / (k - 1) ** 2
    for R in [1, 5, 10, 20, 50, 100, 200, 400, 500, 1000, 10000]:
        pig(t) = ig(k = R * exp(i * t) + R * i) * R * i * exp(i * t)
        rpart, _  = numerical_integral(lambda q: pig(t=q).real(), 0, pi)
        impart, _ = numerical_integral(lambda q: pig(t=q).imag(), 0, pi)
        print("{} {}".format(rpart, impart))

# numc(k) = num(k=conjugate(k))

from sage.misc.viewer import viewer
viewer.pdf_viewer("atril")

_view_later = []
def view_later(formula):
    _view_later.append(formula)

def view_all():
    view(_view_later, tightpage=True, debug=False)


def test_denom():
    W = 2
    v = rvar('v')
    u = rvar('u')
    f(u, v) = (sin(u)**2 + sinh(v)**2) * W**4 - 4 * W**3 * cosh(v) * sinh(v) - 2 * (sin(u)**2 - 3 * sinh(v)**2 - 2) * W**2 - 4 * W * cosh(v) * sinh(v) + sin(u)**2 + sinh(v)**2
    g(u, v) = (sin(u)**2 + sinh(v)**2) * W**4 + 4 * W**3 * cosh(v) * sinh(v) - 2 * (sin(u)**2 - 3 * sinh(v)**2 - 2) * W**2 + 4 * W * cosh(v) * sinh(v) + sin(u)**2 + sinh(v)**2
    res(u, v) = ln(f(u, v) / g(u, v))
    resc(k) = res(k.real(), k.imag())
    # plot_all(resc, suffix='_test_denom')
    plot(res(u=pi + 0.0), (v, 0, 2)).save('test_denom_1.png', figsize=[12, 5])
    plot(res(u=pi + 0.1), (v, 0, 2)).save('test_denom_2.png', figsize=[12, 5])
    plot(res(u=pi + 0.2), (v, 0, 2)).save('test_denom_3.png', figsize=[12, 5])
    plot(res(u=pi + 0.3), (v, 0, 2)).save('test_denom_4.png', figsize=[12, 5])

k = var('k')
S(k) = (k * exp(i * k) + 0 * sin(k) - k) / (k * exp(-i * k) + 0 * sin(k) - k)
contour_integral_analysis(S)


# test_denom()

# u = rvar('u')
# v = rvar('v')
# aa = (num(W=2, a=1, b=2)(k=u + v * i) / cosh(v)).full_simplify()
# view_later(aa)
# view_all()

# print(n(ln(abs(res(k=1000 * (1 + i))))))
# print(latex(res))
# print(res)
# # print(n(abs(res(k=1))))
# for rr in [1, 2, 4, 8, 16, 32, 64]:
#     print("Radius = " + str(rr))
#     print("Num")
#     print(n(num(k=rr*exp(i * pi/4))))
#     print("Denom")
#     print(n(den(k=rr*exp(i * pi/4))))
#     print("Res")
#     print(n(res(k=rr*exp(i * pi/4))))
#     print("=============")
# # contour_integral_analysis(res)
# plot_all(num, suffix='num')
# plot_all(den, suffix='den')
# plot_all(res, suffix='all')
#

# for qq in [1, 0.5, 0.2, 0.1, 0.05, 0.01, 0.001]:
#     print(n(abs(res(k=qq + 0.01 * i))))
#
#