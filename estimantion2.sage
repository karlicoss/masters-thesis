W = 2
a = 1
b = 1
L = 1

k = var('k', domain='complex')

coeff = W**2 + 1
num = W * k * ((a + b) + 2 * i * k) * cos(k * L) * exp(-i * L * k) + (a * b + i * (a + b) * k - coeff * k**2) * exp(-i * L * k) * sin(k * L)
den = W * k * ((a + b) - 2 * i * k) * cos(k * L) * exp( i * L * k) + (a * b - i * (a + b) * k - coeff * k**2) * exp( i * L * k) * sin(k * L)

# num = exp(-i * k) * (cos(k) + 2 * i * sin(k))
# den = exp(i * k) *  (cos(k) - 2 * i * sin(k))

# num = num # * exp(i * L * k)
# den = den # * exp(-i * L * k)

# num = 9 * (exp(2 * i * k) - 1) ** 2 - 64 * exp(2 * i * k)
# den = (exp(4 * i * k) - 9 * exp(2 * i * k)) * (exp(2 * i * k) - 9)
den = (a * b - (W + 1) * i * (a + b) * k - (W + 1)**2 * k**2) - (a * b + (W - 1) * i * (a + b) * k - (W - 1)**2 * k**2) * exp(2 * i * k)
num = den(k=conjugate(k))

print(num)
print(den)
# res = abs(num / den) ** 2
# res2(k) = (num(k=k) * conjugate(num(k=k))) / (num(k=conjugate(k)) * conjugate(num(k=conjugate(k))))
# res = abs(num / den) ** 2
# res2(k) = abs(num(k=k) / (conjugate(num(k=k)))) ** 2
res = num / den
# res2(k) = num(k=k) / num(k=conjugate(k))

def plot_all(ff, suffix="", rrange=(-10, 30), irange=(-2, 2), points=1000):
    # print(n(abs(Sdet(rrange[0] + irange[0] * i)) ** 0.02))
    # complex_plot(Sdet, rrange, irange, plot_points=points).save('plot{}.png'.format(suffix), figsize=[12, 2])
    plot_abs = complex_plot(abs(ff), rrange, irange, plot_points=points)
    # l = line([(rrange[0], 5), (rrange[1], 5)], rgbcolor=(1, 0, 0))# complex_plot(lambda q: q.real() + 5 * i, rrange, irange, color='green')
    ll = [

        # plot(lambda x: ln(x) / 4, rrange, color='red'),
        # plot(lambda x: 1.05 / 4 * ln(x), rrange, color='green'),
        # plot(lambda x: 1.065 / 4 * ln(x), rrange, color='yellow'),
        # plot(lambda x: 1.10 / 4 * ln(x), rrange, color='blue'),
    ]
    # t = var('t')
    # l = parametric_plot((t, 5), (t, rrange[0], rrange[1]), color='red')
    sum([plot_abs] + ll).save('estimation' + suffix + '.png'.format(suffix), figsize=[12, 2])
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

# print(n(ln(abs(res(k=1000 * (1 + i))))))
# print(latex(res))
print(res)
# f(k) = exp(-2 * I * L * k)
# contour_integral_analysis(f)
# contour_integral_analysis(res2)
# plot_all(num, suffix='num')
plot_all(den, suffix='den')
# plot_all(res, suffix='all')
# plot_all(res2, suffix='allc')

# plot_all(numc, suffix='num_conj')


# for qq in [1, 0.5, 0.2, 0.1, 0.05, 0.01, 0.001]:
#     print(n(abs(res(k=qq + 0.01 * i))))
#
#