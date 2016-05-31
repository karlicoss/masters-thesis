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
    for q in [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]:
        r_cayley = 1 - 2 ** (-q)
        rr = R(r=r_cayley).n()
        cc = C(r=r_cayley).n()

        Z0 = cc + rr # TODO Ymax

        print("R = " + str(rr))
        print("C = " + str(cc))

        phi(y) = acos((cc - y) / rr)

    
        s1(y) = 1 - 2 * V / (tanh(y) + V)
        functions = [
            -2 / V * y + 1,
            1 / (2 * V) * (V^2 - 1) * (y - Z),
            s1(y=Z0) / (Z0 - Z) * (y - Z),
        ]

        # ZZ = solve(functions[0] == functions[1], y)[0].rhs()
        # ZZ = (Z + (cc - rr)) / 2
        # ZZ = solve(functions[0] == 0, y)[0].rhs()
        ZZ = (cc - rr) + V/10 # ???? wtf???

        # lol intervals depend on the estimations
        intervalsY = [
            (cc - rr, ZZ),
            (ZZ, Z),
            (Z, Z0),
        ]

        intervalsT = [(-pi/2 + phi(y=yf), -pi/2 + phi(y=yt)) for yf, yt in intervalsY]

        w = var('w', domain=RR)
        # plots = sum([plot(ff(k=rr * exp(i * t) + rr * i), w, tl, tr) for (tl, tr), ff, color in zip(intervals2, functions, rainbow(4))])
        plots = [plot(ff(y=w), w, tl, tr, color=color) for (tl, tr), ff, color in zip(intervalsY, functions, rainbow(len(functions)))]
        plots.append(plot(S(k=i * w), w, intervalsY[0][0], intervalsY[-1][1], color='black'))
        sum(plots).save('plots_{}.png'.format(q))

        sum_logarithmic = 0.0
        sum_jacobian = 0.0
        sum_sup = 0.0
        for (tl, tr), ff in zip(intervalsT, functions):
            print("Interval: {:.2f}-{:.2f}".format(float(tl.real()), float(tr.real())))
            f1 = ln(ff(y=imag(rr * exp(i * t) + cc * i))) * sqrt(rr) / (rr * exp(i * t) + cc * i) # TODO apparently adding i here doesn't help a lot?
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
    g1 = ln(e(y=x))^2 * R / (-C^2 + R^2 + 2 * C * x) * 1 / sqrt(R + C - C) * 1 / sqrt(R - C + Z) # denominator is > 1, we are safe to throw it away, overestimation
    # print((ln(e(y=x))^2 * R / (-C^2 + R^2 + 2 * C * x)).integrate(x))
    # g1: x < C
    # Ok, integrate ln^2(a x + b) / (c x + d), and then estimate?
    # g2 = ln(e(y=C))^2 * R / (-C^2 + R^2 + 2 * C * x) * 1 / sqrt(R + C - x) # sqrt((R + C - x) * (R - C + x))
    g2 = ln(e(y=C))^2 * R / (-C^2 + R^2 + 2 * C * C) * 1 / sqrt(R + C - x) * 1  / sqrt((R - C + C))
    

    for q in range(3, 30, 3):
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
    e(y) = 0.5 / (C + R - Z) * (y - Z) # TODO!!!!
    phi = asin((C - Z) / R)
    f  = ln(e(y=R * sin(t) + C))^2 * 1 / (R + C^2/R + 2 * C * sin(t))
    g  = ln(e(y=x))^2 * R / (-C^2 + R^2 + 2 * C * x) * 1 / sqrt((R + C - x) * (R - C + x))
    # min at x = C!
    g1 = ln(e(y=x))^2 * R / (-C^2 + R^2 + 2 * C * x) * 1 / sqrt(R + C - C) * 1 / sqrt(R - C + Z) # denominator is > 1, we are safe to throw it away, overestimation
    # print((ln(e(y=x))^2 * R / (-C^2 + R^2 + 2 * C * x)).integrate(x))
    # g1: x < C
    # Ok, integrate ln^2(a x + b) / (c x + d), and then estimate?
    # g2 = ln(e(y=C))^2 * R / (-C^2 + R^2 + 2 * C * x) * 1 / sqrt(R + C - x) # sqrt((R + C - x) * (R - C + x))
    g2 = ln(e(y=C))^2 * R / (-C^2 + R^2 + 2 * C * C) * 1 / sqrt(R + C - x) * 1  / sqrt((R - C + C))
    

    for q in range(3, 30, 3):
        print("=======")
        r_cayley = 1 - 2 ** (-q)
        rr = RRR(r=r_cayley).n()
        cc = CCC(r=r_cayley).n()
        print("R = " + str(rr))
        print("C = " + str(cc))

        # print(numerical_integral(f(R=rr, C=cc), -phi(R=rr, C=cc), pi/2)[0])
        print(numerical_integral(g(R=rr, C=cc), 0, Z)[0])
        # print(numerical_integral(g1(R=rr, C=cc), Z, cc)[0])
        # print(numerical_integral(g2(R=rr, C=cc), cc, rr + cc)[0])

def test_schwarz_paper_5():
    Z = 1
    R = var('R', domain='positive')
    C = var('C', domain='positive')



    # tl = C - R
    # tr = Z/2

    # f  = ln(e(y=R * sin(t) + C))^2 * 1 / (R + C^2/R + 2 * C * sin(t))
    e(y) = 1 - y
    g  = ln(e(y=x))^2 * R / (-C^2 + R^2 + 2 * C * x) * 1 / sqrt((R + C - x) * (R - C + x))
    # gg = ln(e(y=x))^2 * R / (-C^2 + R^2 + 2 * C * x) * 1 / sqrt((R + C - x) * (R - C + x))
    gg = ln(e(y=x))^2 * R / (-C^2 + R^2 + 2 * C * x) * 1 / sqrt((R + C - Z/2) * (R - C + x))
    gg = ((-x)/(1 + (-x)))^2 * R / (-C^2 + R^2 + 2 * C * x) * 1 / sqrt((R + C - Z/2) * (R - C + x)) # TODO !!!! estimate!
    # should be easily integrable after?

    ggg = ln(e(y=x))^2 * R / (-C^2 + R^2 + 2 * C * Z/2) * 1 / sqrt((R + C - Z) * (R - C + Z/2)) # easy lol

    # TODO second part: similar to third?....

    # TODO prove about arbitrary linear function?...


    # f  = ln(e(y=x))^2 * R / (2 * R * x) * 1 / sqrt((2 * R - x) * (x))
    # ff = ln(e(y=tr))^2 * R / (2 * R * x) * 1 / sqrt((2 * R - tr) * x)
    # min at x = C!
    # g1 = ln(e(y=x))^2 * R / (-C^2 + R^2 + 2 * C * x) * 1 / sqrt(R + C - C) * 1 / sqrt(R - C + Z) # denominator is > 1, we are safe to throw it away, overestimation
    # print((ln(e(y=x))^2 * R / (-C^2 + R^2 + 2 * C * x)).integrate(x))
    # g1: x < C
    # Ok, integrate ln^2(a x + b) / (c x + d), and then estimate?
    # g2 = ln(e(y=C))^2 * R / (-C^2 + R^2 + 2 * C * x) * 1 / sqrt(R + C - x) # sqrt((R + C - x) * (R - C + x))
    # g2 = ln(e(y=C))^2 * R / (-C^2 + R^2 + 2 * C * C) * 1 / sqrt(R + C - x) * 1  / sqrt((R - C + C))
    

    for q in range(3, 50, 3):
        print("=======")
        r_cayley = 1 - 2 ** (-q)
        rr = RRR(r=r_cayley).n()
        cc = CCC(r=r_cayley).n()
        print("R = " + str(rr))
        print("C = " + str(cc))

        # print(numerical_integral(f(R=rr, C=cc), -phi(R=rr, C=cc), pi/2)[0])
        # print(numerical_integral(g(R=rr, C=cc), cc - rr, Z)[0])
        # print(numerical_integral(gg(R=rr, C=cc), cc - rr, Z)[0])
        print(numerical_integral(g(R=rr, C=cc), cc - rr, Z/2)[0])
        print(numerical_integral(gg(R=rr, C=cc), cc - rr, Z/2)[0])

        print(numerical_integral(g(R=rr, C=cc), Z/2, Z)[0])
        print(numerical_integral(ggg(R=rr, C=cc), Z/2, Z)[0])

# from 0 to Z
def test_schwarz_paper_first():
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


    R = var('R', domain='positive')
    C = var('C', domain='positive')

    ff = C - R
    mm = V / 4 - R + R# C - R + 0.01 # TODO ???? second estimate depends on this variable
    tt = Z - R + R

    o(y)   = ln((V - tanh(y)) / (tanh(y) + V))^2 * R / (-C^2 + R^2 + 2 * C * (y + 1) + 1) * 1 / sqrt((R + C - y ) * (R - C + y))

    g(y)   = ln(1 - 2 / V * y)^2 * R / (-C^2 + R^2 + 2 * C * (y + 1) + 1) * 1 / sqrt((R + C - y ) * (R - C + y))
    g1(y)  = ln(1 - 2 / V * y)^2 * R / (-C^2 + R^2 + 2 * C * (y + 1) + 1) * 1 / sqrt((R + C - tt) * (R - C + y))
    g2(y) =  ((- 2 / V * y) / (1 + - 2 / V * y))^2 * R / (-C^2 + R^2 + 2 * C * (y + 1) + 1) * 1 / sqrt((R + C - mm) * (R - C + y))
    g3(y) =  ((- 2 / V * y) / (1 + - 2 / V * mm))^2 * R / (-C^2 + R^2 + 2 * C * (C - R + 1) + 1) * 1 / sqrt((R + C - mm) * (R - C + y))

    f(y)   = ln(1 / (2 * V) * (V^2 - 1) * (y - Z))^2 * R / (-C^2 + R^2 + 2 * C * (y + 1) + 1) * 1 / sqrt((R + C - y) * (R - C + y))
    f1(y)  = ln(1 / (2 * V) * (V^2 - 1) * (y - Z))^2 * R / (-C^2 + R^2 + 2 * C * ((C - R) + 1) + 1) * 1 / sqrt((R + C - y) * (R - C + y))
    f2(y)  = ln(1 / (2 * V) * (V^2 - 1) * (y - Z))^2 * R / (-C^2 + R^2 + 2 * C * ((C - R) + 1) + 1) * 1 / sqrt((R + C - y) * (R - C + mm))
    f3(y)  = ln(1 / (2 * V) * (V^2 - 1) * (y - Z))^2 * R / (-C^2 + R^2 + 2 * C * ((C - R) + 1) + 1) * 1 / sqrt((R + C - tt) * (R - C + mm))

    for q in range(3, 30, 3):
        print("=======")
        r_cayley = 1 - 2 ** (-q)
        rr = RRR(r=r_cayley).n()
        cc = CCC(r=r_cayley).n()
        # print((-C^2 + R^2 + 2 * (C + 1) * y)(R=rr, C=cc)(y=cc-rr))
        # print(g(R=rr, C=cc)(y=cc - rr))

        print("R = " + str(rr))
        print("C = " + str(cc))

        fff = ff(R=rr, C=cc)
        mmm = mm(R=rr, C=cc)
        ttt = tt(R=rr, C=cc)


        ### plots
        # TODO plot original function as well
        (
            plot( o(R=rr, C=cc), fff, ttt, color='black') + 
            plot( g(R=rr, C=cc), fff, mmm) + 
            plot( f(R=rr, C=cc), mmm, ttt) +
            plot(g1(R=rr, C=cc), fff, mmm, color='magenta') +
            plot(f1(R=rr, C=cc), mmm, ttt, color='magenta') +
            plot(g2(R=rr, C=cc), fff, mmm, color='green') +
            plot(f2(R=rr, C=cc), mmm, ttt, color='green') +
            plot(g3(R=rr, C=cc), fff, mmm, color='red') +
            plot(f3(R=rr, C=cc), mmm, ttt, color='red')
        ).save("plot_%d.png" % q, ymin=0, ymax=1)

        ###

        # plot(g(R=rr, C=cc), ff(R=rr, C=cc), tt(R=rr, C=cc)).save("plotttt.png")
        # raise RuntimeError
        print("!!! First part:")
        for gi in [g, g1, g2, g3]:
            print(numerical_integral(gi(R=rr, C=cc), fff, mmm)[0])

        print("!!! Second part:")
        for fi in [f, f1, f2, f3]:
            # plot(ff(R=rr, C=cc), mm(R=rr, C=cc), tt(R=rr, C=cc)).save("plotttt.png")
            print(numerical_integral(fi(R=rr, C=cc), mmm, ttt)[0])

# test_schwarz_paper()
# test_schwarz_paper_5()
test_schwarz_paper_first()