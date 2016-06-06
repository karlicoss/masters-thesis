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

def my_numerical_integral(expr, x, ff, tt):
    f = fast_callable(expr, vars=[x], domain=CC)
    res, err = numerical_integral(f, ff, tt)
    if err > 0.1: # meh
        raise RuntimeError("Error is too large")
    return res 

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

        Z0 = cc + rr

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
            f1 = ln(ff(y=imag(rr * exp(i * t) + cc * i))) * sqrt(rr) / (rr * exp(i * t) + cc * i) # NOTE apparently adding i here doesn't help a lot?
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


def test_logarithmic_paper():
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

    w(y) = abs(S(k=i * y))
    s(y)   = ln(w(y=y))^2 * R / (-C^2 + R^2 + 2 * (C + 1) * y + 1) * 1 / sqrt((R + C - y ) * (R - C + y))
    # o(y)   = ln(abs((V - tanh(y)) / (tanh(y) + V)))^2 * R / (-C^2 + R^2 + 2 * (C + 1) * y + 1) * 1 / sqrt((R + C - y ) * (R - C + y))
    o = s

    # from C - R to Z
    def first(r_val, c_val, tag=''):
        ff = C - R
        mm = V / 4 - R + R
        tt = Z - R + R

        est_base = [
            o,
            ln(w(y=y))^2 * R / (-C^2 + R^2 + 2 * (C + 1) * y + 1) * 1 / (sqrt(R - C + y) * sqrt(R + C - Z)),
        ]

        # we can't use ln(w(y=mm)) here, integral grows too fast, and we get constant
        # we can't replace y in denominator to ff, still too fast
        e(y) = 1 - 2 / V * y
        est_left = est_base + [
            # ln(w(y=mm))^2 * R / (-C^2 + R^2 + 2 * (C + 1) * mm + 1) * 1 / (sqrt(R - C + y) * sqrt(R + C - Z)), ugh, hard to prove that function increases
            ln(e(y=y))^2                          * R / (-C^2 + R^2 + 2 * (C + 1) * y + 1)  * 1 / (sqrt(R - C + y) * sqrt(R + C - Z)),
            # -x / (1 - x) <= log(1 - x)
            ((- 2 / V * y) / (1 - 2 / V * y ))^2  * R / (-C^2 + R^2 + 2 * (C + 1) * y + 1)  * 1 / (sqrt(R - C + y) * sqrt(R + C - Z)),
            # since x / sqrt(C + x) is always increasing,
            ((- 2 / V * y) / (1 - 2 / V * mm))^2  * R / (-C^2 + R^2 + 2 * (C + 1) * y + 1)  * 1 / (sqrt(R - C + y) * sqrt(R + C - Z)),
            # notice that function y^2/(a + b * y) increases
            ((- 2 / V * mm) / (1 - 2 / V * mm))^2 * R / (-C^2 + R^2 + 2 * (C + 1) * mm + 1) * 1 / (sqrt(R - C + y) * sqrt(R + C - Z)),
            # yay!!!
        ]

        ee(y) = 1 / (2 * V) * (V^2 - 1) * (y - Z)
        est_right = est_base + [
            # ln(w(y=y))^2  * R / (-C^2 + R^2 + 2 * (C + 1) * mm + 1) * 1 / (sqrt(R - C + mm) * sqrt(R + C - Z)),
            # ln(ee(y=y))^2 * R / (-C^2 + R^2 + 2 * (C + 1) * mm + 1) * 1 / (sqrt(R - C + mm) * sqrt(R + C - Z)),
        ]
        # f3: singularity of type ln^2(x) around x = 0 is integrable and bounded by constant since V/4 is independent of R and C

        fff = ff(R=rr, C=cc)
        mmm = mm(R=rr, C=cc)
        ttt = tt(R=rr, C=cc)

        ### plots
        colors = [
            'black',
            'blue',
            'magenta',
            'green',
            'red'
        ]

        # plot1 = sum(plot(fi(R=rr, C=cc), fff, mmm, color=col) for fi, col in zip(est_left, colors))
        # plot2 = sum(plot(fi(R=rr, C=cc), mmm, ttt, color=col) for fi, col in zip(est_right, colors))

        # (plot1 + plot2).save("plot_%s.png" % tag, ymin=0, ymax=1)

        ### integrals
        print("!!! First part:")
        for fi in est_left:
            print(my_numerical_integral(fi(R=rr, C=cc), y, fff, mmm))

        print("!!! Second part:")
        for fi in est_right:
            # plot(ff(R=rr, C=cc), mm(R=rr, C=cc), tt(R=rr, C=cc)).save("plotttt.png")
            print(my_numerical_integral(fi(R=rr, C=cc), y, mmm, ttt))
        ##

    # from Z to C + R
    def second(r_val, c_val, tag=''):
        ff = Z - R + R
        ll = Z + 1 / R # huh?
        mm = C
        tt = C + R

        e(y) = w(y=tt) / (C + R - Z) * (y - Z)
        # e(y) = w(y=tt) / (C + R) * (y - Z) # NOTE ok for big enough C + R !!!!

        fff = ff(R=rr, C=cc)
        lll = ll(R=rr, C=cc)
        mmm = mm(R=rr, C=cc)
        ttt = tt(R=rr, C=cc)


        est_base = [
            o,
            ln(e(y=y))^2 * R / (-C^2 + R^2 + 2 * (C + 1) * y + 1) * 1 / sqrt((R + C - y) * (R - C + y)),
            ln(e(y=y))^2 * R / (-C^2 + R^2 + 2 * (C + 1) * y + 1) * 1 / sqrt((R + C - y) * (R - C + ff)),
        ]

        # ff to ll
        est_left = est_base + [
            ln(e(y=y))^2 * R / (-C^2 + R^2 + 2 * (C + 1) * ff + 1) * 1 / sqrt((R + C - ll) * (R - C + ff)),
        ]

        # ll to mm
        est_mid = est_base + [
            ln(e(y=ll))^2 * R / (-C^2 + R^2 + 2 * (C + 1) * y + 1) * 1 / sqrt((R + C - mm) * (R - C + ff)),
        ]

        # mm to tt
        est_right = est_base + [
            ln(e(y=mm))^2 * R / (-C^2 + R^2 + 2 * (C + 1) * mm + 1) * 1 / sqrt((R + C - y) * (R - C + ff)),
        ]

        colors = [
            'black',
            # 'blue',
            'magenta',
            # 'green',
            'red'
        ]

        # print(n(fff))
        # print(n(lll))
        # print(n(mmm))
        # print(n(ttt))

        # plot1 = sum(plot(fi(R=rr, C=cc), fff, lll, color=col) for fi, col in zip(est_left, colors))
        # plot2 = sum(plot(fi(R=rr, C=cc), lll, mmm, color=col) for fi, col in zip(est_mid, colors))
        # plot3 = sum(plot(fi(R=rr, C=cc), mmm, ttt, color=col) for fi, col in zip(est_right, colors))
        # (plot1 + plot2 + plot3).save("plot_%s.png" % tag, ymin=0, ymax=4)




        print("Left:")
        for fi in est_left:
            print(my_numerical_integral(fi(R=rr, C=cc), y, fff, lll))

        print("Mid:")
        for fi in est_mid:
            print(my_numerical_integral(fi(R=rr, C=cc), y, lll, mmm))

        print("Right:")
        for fi in est_right:
            print(my_numerical_integral(fi(R=rr, C=cc), y, mmm, ttt))


    for q in range(3, 30, 3):
        print("=======")
        r_cayley = 1 - 3 ** (-q)
        # r_cayley = 0.69
        rr = RRR(r=r_cayley).n()
        cc = CCC(r=r_cayley).n()

        print("R = " + str(rr))
        print("C = " + str(cc))

        first(rr, cc, tag='f_' + str(q))
        # second(rr, cc, tag='second_' + str(q))


test_logarithmic_paper()