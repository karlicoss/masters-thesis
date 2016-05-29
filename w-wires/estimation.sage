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


C(r) = ((icayley(r) + icayley(-r)) / 2).imag()
R(r) = ((icayley(r) - icayley(-r)) / 2).imag()

CCC(r) = ((icayley(r) + icayley(-r)) / 2).imag()
RRR(r) = ((icayley(r) - icayley(-r)) / 2).imag()


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


# cv(t) = changevar(f * 1 / (k + i)^2, k == rr * exp(i * t) + rr * i, t)
def test_schwarz(f):
    Z = atanh(4/5)
    for rr in [2 ** q for q in range(1, 30)]:
        phi = n(acos((rr - Z) / rr))
        f1 = f(k=rr * exp(i * t) + rr * i) * sqrt(rr) / (rr * exp(i * t) + rr * i + i) # ???
        f2 = 1 / (rr * exp(i * t) + rr * i + i) * i * sqrt(rr) * exp(i * t) # looks like it converges to pi
        p  = f1 * f2
        p1 = fast_callable(norm(f1)     , vars=[t], domain=CC)
        p2 = fast_callable(norm(f2)     , vars=[t], domain=CC)
        print("Phi: " + str(phi))
        intervals = [
            (-pi/2 - phi, -pi/2 + phi),
            (-pi/2 + phi, -pi/2 + 2 * pi - phi),
            (0, 2 * pi),
        ]
        for tl, tr in intervals:
            res = norm(complex_integral(p, t, tl, tr))
            print("Complex: " + str(res))

            res1, _ = numerical_integral(p1, tl, tr)
            print("Logarithmic: " + str(res1))

            res2, _ = numerical_integral(p2, tl, tr)
            print("Jacobian: " + str(res2))

            print("Upper bound: " + str(res1 * res2))
            print("---------------")            
        print("================================")


Z = atanh(4/5)
# test_schwarz(ln(abs(1/Z * (Z - k.imag()))))
# test_schwarz(ln(S))

# test_schwarz(k / k)


def test_schwarz_paper():
    W = 2
    coeff = W ** 2 + 1
    a = 0
    b = 0
    L = 1
    num = W * (a + b) * k * cos(L * k) + (a * b - coeff * k ** 2) * sin(L * k) + 2 * W * i * k ** 2 * cos(L * k) + i * (a + b) * k * sin(L * k)
    den = W * (a + b) * k * cos(L * k) + (a * b - coeff * k ** 2) * sin(L * k) - 2 * W * i * k ** 2 * cos(L * k) - i * (a + b) * k * sin(L * k)
    S = abs(num / den)



    V = 2 * W / (W^2 + 1)
    Z = atanh(V)
    for q in [3, 4, 5, 6, 7, 8, 9, 10, 11]:
        r_cayley = 1 - 2 ** (-q)
        rr = R(r=r_cayley).n()
        cc = C(r=r_cayley).n()
        # rr = 2 ** q
        # cc = rr

        Z0 = cc + rr

        # rr = 5 ** q
        print("R = " + str(rr))
        phi(y) = acos((cc - y) / rr)

        # TODO one more estimate?
        s1(y) = 1 - 2 * V / (tanh(y) + V)
        functions = [
            2/V * (0.5 * V - k.imag()),
            (W^2 - 1)^2 / (4 * (W^2 + 1) * W) * (Z - k.imag()),
            s1(y=Z0) / (Z0 - Z) * (k.imag() - Z),
            s1(y=Z0),
        ]

        ZZ = solve(functions[0] == functions[1], k)[0].rhs()

        # lol intervals depend on the estimations
        intervalsY = [
            (cc - rr, ZZ),
            (ZZ, Z),
            (Z, Z0),
            # (Z0, cc + rr),
        ]

        intervalsT = [(-pi/2 + phi(y=yf), -pi/2 + phi(y=yt)) for yf, yt in intervalsY]

        w = var('w', domain=RR)
        # plots = sum([plot(ff(k=rr * exp(i * t) + rr * i), w, tl, tr) for (tl, tr), ff, color in zip(intervals2, functions, rainbow(4))])
        plots = [plot(ff(k=w * i), w, tl, tr, color=color) for (tl, tr), ff, color in zip(intervalsY, functions, rainbow(4))]
        plots.append(plot(S(k=1 + w * i), w, intervalsY[0][0], intervalsY[-1][1], color='black'))
        sum(plots).save('plots_{}.png'.format(q))



        sum_logarithmic = 0.0
        sum_jacobian = 0.0
        sum_sup = 0.0
        for (tl, tr), ff in zip(intervalsT, functions):
            print("Interval: {:.2f}-{:.2f}".format(float(tl.real()), float(tr.real())))
            f1 = ln(ff(k=rr * exp(i * t) + cc * i)) * sqrt(rr) / (rr * exp(i * t) + cc * i)
            f2 = 1 / (rr * exp(i * t) + cc * i + i) * i * sqrt(rr) * exp(i * t) # looks like it converges to pi
            p  = f1 * f2
            p1 = fast_callable(norm(f1)     , vars=[t], domain=CC)
            p2 = fast_callable(norm(f2)     , vars=[t], domain=CC)
        
            res = norm(complex_integral(p, t, tl, tr))
            print("Complex: " + str(res))

            res1, _ = numerical_integral(p1, tl, tr)
            sum_logarithmic += res1
            print("Logarithmic: " + str(res1))

            res2, _ = numerical_integral(p2, tl, tr)
            sum_jacobian += res2
            print("Jacobian: " + str(res2))

            up = res1 * res2
            sum_sup += up
            print("Upper bound: " + str(up))
            print("---------------")            
        print("Total logarithmic: " + str(sum_logarithmic))
        print("Total Jacobian   : " + str(sum_jacobian))
        print("Total upper bound: " + str(sum_sup))
        print("================================")

# test_schwarz_paper()

def test_schwarz_paper_x():
    W = 2
    coeff = W ** 2 + 1
    a = 0
    b = 0
    L = 1
    num = W * (a + b) * k * cos(L * k) + (a * b - coeff * k ** 2) * sin(L * k) + 2 * W * i * k ** 2 * cos(L * k) + i * (a + b) * k * sin(L * k)
    den = W * (a + b) * k * cos(L * k) + (a * b - coeff * k ** 2) * sin(L * k) - 2 * W * i * k ** 2 * cos(L * k) - i * (a + b) * k * sin(L * k)
    S = abs(num / den)



    V = 2 * W / (W^2 + 1)
    Z = atanh(V)
    for q in range(3, 30, 2):
        r_cayley = 1 - 2 ** (-q)
        rr = R(r=r_cayley).n()
        cc = C(r=r_cayley).n()
        # rr = 2 ** q
        # cc = rr

        print("R = " + str(rr))


        f1 = ln(1/(2 * rr - Z) * (x - Z)) ** 2
        f2 = 1 / (x * sqrt(x * (2 * rr - x)))
        f = f1 * f2

        g1 = f1.diff(x)
        g2 = f2.integrate(x)

        
        print(g1 * g2)

        tl = Z
        tr = 2 * rr

        # mv, mx = find_local_minimum(, tl, tr)

        pf = fast_callable(f            , vars=[x], domain=CC)
        # p1 = fast_callable(norm(f1)     , vars=[x], domain=CC)
        # p2 = fast_callable(norm(f2)     , vars=[x], domain=CC)
        pg = fast_callable(g1 * g2      , vars=[x], domain=CC)

        resf , _ = numerical_integral(pf, tl, tr)
        resg , _ = numerical_integral(pg, tl, tr)
        # res1, _ = numerical_integral(p1, tl, tr)
        # res2, _ = numerical_integral(p2, tl, tr)


        print("Resf: " + str(resf))
        print("Resg: " + str(resg))
        # print("Est: " + str(res1 * res2))


def test_schwarz_paper_w():
    Z = 1
    R = var('R', domain='positive')
    # assume(R > 1)
    # assume(-4*R^2/(4*R^2-4*R+1)+4*R/(4*R^2-4*R+1)+3/(4*R^2-4*R+1)-4*R^2/(2*R-1)+6*R/(2*R-1)+4/(2*R-1)+2*R+1 > 0)
    # assume(-4*R^2/(4*R^2-4*R+1)<0)

    f1 = ln(0.1 / (2 * R - Z) * (x - Z)) ** 2 
    f2 = 1 / (x * sqrt(x * (2 * R - x)))
    f = f1 * f2

    Rv = 100
    plot(f(R=Rv) , Z, 2 * Rv).save('plot.png')
    plot(f1(R=Rv), Z, 2 * Rv).save('plot1.png')
    plot(f2(R=Rv), Z, 2 * Rv).save('plot2.png')
    it = lambda rr: numerical_integral(f1(R=rr), Z, 2 * rr)
    for q in range(1, 5):
        print(it(10**q))

    
    # print(f1.diff(x))
    # print(f2.integrate(x))

    # print(f)


    # for q in range(3, 30, 2):
    #     r_cayley = 1 - 2 ** (-q)
    #     rr = R(r=r_cayley).n()
    #     cc = C(r=r_cayley).n()
    #     # rr = 2 ** q
    #     # cc = rr

    #     print("R = " + str(rr))


    #     f1 = ln(1/(2 * rr - Z) * (x - Z)) ** 2
    #     f2 = 1 / (x * sqrt(x * (2 * rr - x)))
    #     f = f1 * f2

    #     g1 = f1.diff(x)
    #     g2 = f2.integrate(x)

        
    #     print(g1 * g2)

    #     tl = Z
    #     tr = 2 * rr

    #     # mv, mx = find_local_minimum(, tl, tr)

    #     pf = fast_callable(f            , vars=[x], domain=CC)
    #     # p1 = fast_callable(norm(f1)     , vars=[x], domain=CC)
    #     # p2 = fast_callable(norm(f2)     , vars=[x], domain=CC)
    #     pg = fast_callable(g1 * g2      , vars=[x], domain=CC)

    #     resf , _ = numerical_integral(pf, tl, tr)
    #     resg , _ = numerical_integral(pg, tl, tr)
    #     # res1, _ = numerical_integral(p1, tl, tr)
    #     # res2, _ = numerical_integral(p2, tl, tr)


    #     print("Resf: " + str(resf))
    #     print("Resg: " + str(resg))
    #     # print("Est: " + str(res1 * res2))

# Wolfram: integrate ln(0.5/(2 * 50000 - 1)(x - 1))/(x * sqrt(x * (2 * 50000 - x))) for x from 1 to 2 * 50000


def test_schwarz_paper_3():
    Z = 1
    R = var('R', domain='positive')
    C = var('C', domain='positive')
    e(y) = 0.5 / (C + R - Z) * (y - Z) # TODO!!!!
    phi = asin((C - Z) / R)
    f  = ln(e(y=R * sin(t) + C))^2 * 1 / (R + C^2/R + 2 * C * sin(t))
    g  = ln(e(y=x))^2 * R / (-C^2 + R^2 + 2 * C * x) * 1 / sqrt((R + C - x) * (R - C + x))
    # min at x = C!
    qq = 1 / sqrt(R + C - x)
    g1 = ln(e(y=x))^2 * R / (-C^2 + R^2 + 2 * C * x) * qq(x = C) # * 1 / sqrt(R - C + x) # denominator is > 1, we are safe to throw it away, overestimation
    # print((ln(e(y=x))^2 * R / (-C^2 + R^2 + 2 * C * x)).integrate(x))
    # g1: x < C
    # Ok, integrate ln^2(a x + b) / (c x + d), and then estimate?
    g2 = ln(e(y=C))^2 * R / (-C^2 + R^2 + 2 * C * x) * 1 / sqrt(-C^2 + R^2 + 2 * C * x - x^2)
    

    for q in range(3, 20, 2):
        print("=======")
        r_cayley = 1 - 2 ** (-q)
        rr = RRR(r=r_cayley).n()
        cc = CCC(r=r_cayley).n()
        print("R = " + str(rr))
        print("C = " + str(cc))

        print(numerical_integral(f(R=rr, C=cc), -phi(R=rr, C=cc), pi/2)[0])
        print(numerical_integral(g(R=rr, C=cc), Z, rr + cc)[0])
        print(numerical_integral(g1(R=rr, C=cc), Z, cc)[0])
        print(numerical_integral(g2(R=rr, C=cc), cc, rr + cc)[0])



def test_schwarz_paper_4():
    Z = 1
    R = var('R', domain='positive')
    C = var('C', domain='positive')
    e(y) = 0.5 / (2 * R - Z) * (y - Z)

    f  = ln(e(y=x))^2 * 1 / (x * sqrt(x * (2 * R - x)))
    # f1 = ln(e(y=x))^2 * 1 / (x * sqrt(2 * R - x)) # TODO R - x?
    qq = 1 / sqrt(2 * R - x)
    # (y - y1) / (x - x1) = (y2 - y1) / (x2 - x1)
    # x1 = Z
    # x2 = R
    # ww = (qq(x=R) - qq(x=Z)) / (R - Z) * (x - Z) + qq(x = Z)
    ww = qq(x=R)
    f1 = ln(e(y=x))^2 * 1 / x * ww
    f2 = ln(e(y=R))^2 * 1 / (x * sqrt(x * (2 * R - x)))

    # f21 = norm(ln(e(y=x)) * R / (-(C - x)^2 + R^2 + x^2) * sqrt(sqrt(-(C - x)^2 + R^2)))
    # f22 = norm(ln(e(y=x)) * 1 / sqrt(sqrt(-(C - x)^2 + R^2)))

    # print(numerical_integral(f(R=rr), Z, 2 * rr)[0])

    for q in range(3, 30, 3):
        print("=======")
        r_cayley = 1 - 2 ** (-q)
        rr = RRR(r=r_cayley).n()
        cc = CCC(r=r_cayley).n()
        print("R = " + str(rr))

        print(numerical_integral(f(R=rr), Z, 2 * rr)[0])
        print(numerical_integral(f1(R=rr), Z, rr)[0])
        print(numerical_integral(f2(R=rr), rr, 2 * rr)[0])

test_schwarz_paper_3()