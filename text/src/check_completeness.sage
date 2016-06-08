from numpy import arange

def cvar(name):
    return var(name, domain=CC)

def rvar(name):
    return var(name, domain=RR)

def icayley(z):
    return i * (1 + z) / (1 - z)

def cayley(z):
    return (z - i) / (z + i)

def my_numerical_integral(expr, x, ff, tt):
    f = fast_callable(expr, vars=[x], domain=CC)
    res, err = numerical_integral(f, ff, tt)
    if err > 0.1:
        raise RuntimeError("Error is too large")
    return res 


k = cvar('k')
R, t = rvar('R'), rvar('t')
x, y = rvar('x'), rvar('y')


CCC(r) = ((icayley(r) + icayley(-r)) / 2).imag()
RRR(r) = ((icayley(r) - icayley(-r)) / 2).imag()


def plot_all(lines):
    p = sum(lines)
    p.set_legend_options(handlelength=5)
    p.axes_labels(['$q$', ''])
    return p

def test_convergence(s, color='blue', style='-', label=''):
    rs     = []
    ints   = []
    for q in range(3, 30, 1):
        r = 1 - 2 ^ (-q)

        rr = RRR(r=r).n()
        cc = CCC(r=r).n()

        f = ln(s(k=rr * exp(i * t) + cc * i)) * sqrt(rr) / (rr * exp(i * t) + cc * i + i)
        g = sqrt(rr) / (rr * exp(i * t) + cc * i + i)

        int_f = my_numerical_integral(norm(f), t, 0, 2 * pi)
        int_g = my_numerical_integral(norm(g), t, 0, 2 * pi)

        rs.append(q)
        ints.append(sqrt(int_f * int_g))
        print("F = " + str(int_f))

    return line(zip(rs, ints), rgbcolor=color, linestyle=style, legend_label=label)

a_values = [100, 10, 1, 0.1]
a_zero   = 0.0

a = rvar('a')
s = abs((cos(k) + (a / (2 * k) + i) * sin(k) - 1) / (cos(k) + (a / (2 * k) - i) * sin(k) - 1))

plots = []
for aa, color, style in zip(a_values, rainbow(len(a_values)), ['-', '--', '-.', ':']):
    label = "a = {:.2f}".format(float(aa))
    plots.append(test_convergence(s(a=aa), color=color, style=style, label=label))


plot_all(plots).save('plots_converging.png', ymin=0, ymax=6, dpi=150)
plot_all([test_convergence(s(a=a_zero))]).save('plots_diverging.png', ymin=0, dpi=150)