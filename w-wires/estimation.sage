
def cvar(name):
    return var(name, domain=CC)


def rvar(name):
    return var(name, domain=RR)


def pvar(name):
    return var(name, domain='positive')


def icayley(z):
    return i * (1 + z) / (1 - z)


def cayley(z):
    return (z - i) / (z + i)


def complex_integral(expr, x, ff, tt):
    rp = fast_callable(expr.real(), vars=[x], domain=CC)
    ip = fast_callable(expr.imag(), vars=[x], domain=CC)
    eps = 1e-02
    rpart, _ = numerical_integral(rp, ff, tt, eps_rel=eps, rule=1)
    impart, _ = numerical_integral(ip, ff, tt, eps_rel=eps, rule=1)
    return rpart + i * impart


def changevar(f, eqn, newvar):
    dx = diff(eqn.rhs(), newvar)
    return f.subs(eqn) * dx


k = cvar('k')
R, t = rvar('R'), rvar('t')
x, y = rvar('x'), rvar('y')



W = 2
coeff = W ** 2 + 1
a = 0
b = 0
L = 1
num = W * (a + b) * k * cos(L * k) + (a * b - coeff * k ** 2) * sin(L * k) + 2 * W * i * k ** 2 * cos(L * k) + i * (a + b) * k * sin(L * k)
den = W * (a + b) * k * cos(L * k) + (a * b - coeff * k ** 2) * sin(L * k) - 2 * W * i * k ** 2 * cos(L * k) - i * (a + b) * k * sin(L * k)
S = abs(num / den)

def compute_integral(f):
    ig = f * diff(cayley(k), k)
    return ig
    # print(ig)
    # pig = changevar(ig, k == R * exp(i * t) + R * i, t).full_simplify()
    # print(pig(t=0.4999 * pi).simplify()(R=1000).n())  # .limit(R=oo))
    # for rr in [1, 2, 5, 10, 15, 20, 25]:
    #     print(complex_integral(pig(R=rr), t, 0, 2 * pi))
    # print(pig.integrate(t, 0, 2 * pi))
    # pig = f(k=)

    # for R in [2, 4, 8, 16, 32, 64]:

# compute_integral(S)
g(k) = sin(k) / (k - i)

left = -10
right = 10

i1 = complex_integral(g(k=x), x, left, right)
print(i1)

gg(z) = changevar(g, k == icayley(z), z) # TODO differential?
ggg(t) = changevar(gg, z == exp(i * t), t)
i2 = complex_integral(ggg, t, arg(cayley(left)), 2 * pi + arg(cayley(right)))
print(i2)