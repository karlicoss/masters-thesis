from numpy import arange


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


# f(z): Cayley domain,   z in D
# F(k): standard domain, k in C_+

def from_standard_to_cayley():
    F(k) = sin(k) / (k - i)

    Left = -10
    Right = 10

    i1 = complex_integral(F(k=x), x, Left, Right)
    print(i1)

    f(z) = changevar(F, k == icayley(z), z)
    ff(t) = changevar(f, z == exp(i * t), t)


    left = arg(cayley(Left))
    right = arg(cayley(Right))
    if right < left:
        right += 2 * pi
    i2 = complex_integral(ff, t, left, right)
    print(i2)


def from_cayley_to_standard():
    f(z) = sin(z) * z
    left = pi/2
    right = 3 * pi/2
    ff(z) = changevar(f, z == exp(i * t), t)
    i1 = complex_integral(ff, t, left, right)
    print(i1)

    F(k) = changevar(f, z == cayley(k), k)
    Left = icayley(exp(i * left))
    Right = icayley(exp(i * right))
    i2 = complex_integral(F, k, Left, Right)
    print(i2)

# given F(k), f(z) = F(icayley(z))
# intagral in Cayley space is   f(z) dz
# integral in standard space is f(z=cayley(k)) * cayley'(k) = F(k) * cayley'(k) dk = F(k) * (2 * i / (k + i)^2)


C(r) = (icayley(r) + icayley(-r)) / 2
R(r) = ((icayley(r) - icayley(-r)) / (2 * i)).real() # TODO divide by i

# line([(icayley(R * exp(i * q)) - C(r=R))(R=0.1).n() for q in arange(0, 2 * pi, 0.01 * pi)]).show(aspect_ratio='equal')

# F(k) = ln(S(k=k)) # sin(k) / (k - i)

# Left = -10
# Right = 10

# i1 = complex_integral(F(k=x), x, Left, Right)
# print(i1)

# f(z) = changevar(F, k == icayley(z), z)
# ff(t) = changevar(f, z == exp(i * t), t)


# left = arg(cayley(Left))
# right = arg(cayley(Right))
# if right < left:
#     right += 2 * pi
# i2 = complex_integral(ff, t, left, right)
# print(i2)

def test_integral():
    ig(k) = ln(S(k=k)) * cayley(k).diff()

    for rr in [1, 5, 10, 20, 40, 100, 200, 400, 1000]:
        tig(t) = changevar(ig, k == rr * exp(i * t) + rr * i, t)
        print(complex_integral(tig, t, 0, 2 * pi))
    # integrate it now!

def test_first_segment(f):
    Z = atanh(4/5)
    for rr in [1, 5, 10, 20, 50, 100, 1000, 2000, 5000]:
        phi = acos((rr - Z) / rr)
        cv(t) = changevar(f * 1 / (k + i)^2, k == rr * exp(i * t) + rr * i, t)
        print(cv)
        print(complex_integral(cv, t, -pi/2 - phi, -pi/2 + phi))

def test_second_segment(f):
    Z = atanh(4/5)
    for rr in [1, 5, 10, 20, 50, 100, 1000, 2000, 5000]:
        phi = acos((rr - Z) / rr)
        cv(t) = changevar(f * 1 / (k + i)^2, k == rr * exp(i * t) + rr * i, t)
        print(complex_integral(cv, t, -pi/2 + phi, -pi/2 - phi + 2 * pi))

# integrate 1/(k + i)^2 for k from 0 to oo, apparently converges. 
# estimate by integrate 1/(k^2 + 1) for k from 0 to oo ?

def test_all(S):
    # upper bound
    # test_first_segment(S - 1)
    test_first_segment(ln(S))
    # lower bound 
    test_first_segment(1 - 1 / S)

    # test_second_segment(S - 1)
    test_second_segment(ln(S))
    test_second_segment(1 - 1 / S) # 1 - 1/S ? TODO 1 converges, so sufficient condition is convergence of 1/S ?



# cosh(atanh(x)) - 2 * sinh(atanh(x)) TODO????

# tanh doesn't help, still divergent....

# C = 0.5
# sage: f(x) = 1 - 1/(C - tanh(x))
# sage: f.integrate(x, 0, atanh(C))
# -1.3333333333333333*log(2) - 46.95090622996493
# sage: f.integrate(x)
# x |--> 3*x + 4/3*log(3*e^(-2*x) - 1)
# clearly divergent!