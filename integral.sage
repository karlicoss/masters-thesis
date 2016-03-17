sage: region_plot(lambda x, y: abs(f(x, y)) < 1, (x, -10, 10), (y, -10, 10))
---------------------------------------------------------------------------
NameError                                 Traceback (most recent call last)
<ipython-input-1-e64d593ef8fe> in <module>()
----> 1 region_plot(lambda x, y: abs(f(x, y)) < Integer(1), (x, -Integer(10), Integer(10)), (y, -Integer(10), Integer(10)))

NameError: name 'y' is not defined
sage: def f(k):
....:     return (cos(k) + 2 * i * sin(k)) / (cos(k) - 2 * i * sin(k))
....: 
sage: region_plot(lambda x, y: abs(f(x, y)) < 1, (x, -10, 10), (y, -10, 10))
---------------------------------------------------------------------------
NameError                                 Traceback (most recent call last)
<ipython-input-3-e64d593ef8fe> in <module>()
----> 1 region_plot(lambda x, y: abs(f(x, y)) < Integer(1), (x, -Integer(10), Integer(10)), (y, -Integer(10), Integer(10)))

NameError: name 'y' is not defined
sage: x, y = var('x y')
sage: region_plot(lambda x, y: f(x, y) < 0, (x, -10, 10), (y, -10, 10))
---------------------------------------------------------------------------
TypeError                                 Traceback (most recent call last)
<ipython-input-5-4ca383134c9a> in <module>()
----> 1 region_plot(lambda x, y: f(x, y) < Integer(0), (x, -Integer(10), Integer(10)), (y, -Integer(10), Integer(10)))

/L/soft/SageMath/local/lib/python2.7/site-packages/sage/misc/decorators.pyc in wrapper(*args, **kwds)
    548                 options['__original_opts'] = kwds
    549             options.update(kwds)
--> 550             return func(*args, **options)
    551 
    552         #Add the options specified by @options to the signature of the wrapped

/L/soft/SageMath/local/lib/python2.7/site-packages/sage/plot/contour_plot.pyc in region_plot(f, xrange, yrange, plot_points, incol, outcol, bordercol, borderstyle, borderwidth, alpha, **options)
    944     xy_data_arrays = numpy.asarray([[[func(x, y) for x in xsrange(*ranges[0], include_endpoint=True)]
    945                                      for y in xsrange(*ranges[1], include_endpoint=True)]
--> 946                                     for func in f_all[neqs::]],dtype=float)
    947     xy_data_array=numpy.abs(xy_data_arrays.prod(axis=0))
    948     # Now we need to set entries to negative iff all

/L/soft/SageMath/local/lib/python2.7/site-packages/sage/plot/contour_plot.pyc in <lambda>(x, y)
   1021     from sage.symbolic.expression import is_Expression
   1022     if not is_Expression(f):
-> 1023         return lambda x,y: -1 if f(x,y) else 1
   1024 
   1025     op = f.operator()

<ipython-input-5-4ca383134c9a> in <lambda>(x, y)
----> 1 region_plot(lambda x, y: f(x, y) < Integer(0), (x, -Integer(10), Integer(10)), (y, -Integer(10), Integer(10)))

TypeError: f() takes exactly 1 argument (2 given)
sage: region_plot(lambda x, y: f(x + i * y) < 0, (x, -10, 10), (y, -10, 10))
Launched png viewer for Graphics object consisting of 1 graphics primitive
sage: region_plot(lambda x, y: abs(f(x + i * y)) < 0.1, (x, -10, 10), (y, -10, 10))
Launched png viewer for Graphics object consisting of 1 graphics primitive
sage: region_plot(lambda x, y: abs(f(x + i * y)) < 1.0, (x, -10, 10), (y, -10, 10))
Launched png viewer for Graphics object consisting of 1 graphics primitive
sage: region_plot(lambda x, y: abs(f(x + i * y)) < 0.5, (x, -10, 10), (y, -10, 10))
Launched png viewer for Graphics object consisting of 1 graphics primitive
sage: region_plot(lambda x, y: abs(x ** 2 + y ** 2) < 1, (x, -10, 10), (y, -10, 10))
Launched png viewer for Graphics object consisting of 1 graphics primitive
sage: region_plot(lambda x, y: abs(f(x + i * y)) < 0.5, (x, -10, 10), (y, -10, 10))
Launched png viewer for Graphics object consisting of 1 graphics primitive
sage: region_plot(lambda x, y: abs(f(x + i * y)) < 0.3, (x, -10, 10), (y, -10, 10))
Launched png viewer for Graphics object consisting of 1 graphics primitive
sage: region_plot(lambda x, y: abs(f(x + i * y)) < 0.3, (x, -10, 10), (y, -10, 10))
Launched png viewer for Graphics object consisting of 1 graphics primitive
sage: region_plot(lambda x, y: abs(f(x + i * y)) < 0.1, (x, -10, 10), (y, -10, 10))
Launched png viewer for Graphics object consisting of 1 graphics primitive
sage: region_plot(lambda x, y: abs(f(x + i * y)) < 0.1, (x, -10, 10), (y, -10, 10), plot_points=10)
Launched png viewer for Graphics object consisting of 1 graphics primitive
sage: region_plot(lambda x, y: abs(f(x + i * y)) < 0.1, (x, -10, 10), (y, -10, 10), plot_points=100)
Launched png viewer for Graphics object consisting of 1 graphics primitive
sage: region_plot(lambda x, y: abs(f(x + i * y)) < 0.1, (x, -10, 10), (y, -10, 10), plot_points=100)
KeyboardInterrupt
sage: def ff(k):
....:     return cos(k) - 2 * i * sin(k)
....: 
sage: region_plot(lambda x, y: abs(ff(x + i * y)) < 0.1, (x, -10, 10), (y, -10, 10), plot_points=100)
Launched png viewer for Graphics object consisting of 1 graphics primitive
sage: region_plot(lambda x, y: abs(ff(x + i * y)) < 1, (x, -10, 10), (y, -10, 10), plot_points=100)
Launched png viewer for Graphics object consisting of 1 graphics primitive
sage: region_plot(lambda x, y: abs(ff(x + i * y)) < 10, (x, -10, 10), (y, -10, 10), plot_points=100)
Launched png viewer for Graphics object consisting of 1 graphics primitive
sage: region_plot(lambda x, y: 1 < abs(ff(x + i * y)) < 10, (x, -10, 10), (y, -10, 10), plot_points=100)
Launched png viewer for Graphics object consisting of 1 graphics primitive
sage: complex_plot(cos, (-10, 10), (-10, 10))
Launched png viewer for Graphics object consisting of 1 graphics primitive
sage: region_plot(lambda x, y: abs(f(x + i * y)) < 0.5, (x, -10, 10), (y, -10, 10))
KeyboardInterrupt
sage: def f(k):
    return (cos(k) + 2 * i * sin(k)) / (cos(k) - 2 * i * sin(k))
KeyboardInterrupt
sage: region_plot(lambda x, y: abs(ff(x + i * y)) < 10, (x, -10, 10), (y, -10, 10), plot_points=100)
KeyboardInterrupt
sage: region_plot(lambda x, y: abs(ff(x + i * y)) < 10, (x, -10, 10), (y, -10, 10), plot_points=100)
KeyboardInterrupt
sage: def fff(k):
....:     return cos(k) + 2 * i * sin(k)
....: 
sage: region_plot(lambda x, y: abs(fff(x + i * y)) < 0.1, (x, -10, 10), (y, -10, 10), plot_points=10)
Launched png viewer for Graphics object consisting of 1 graphics primitive
sage: region_plot(lambda x, y: abs(fff(x + i * y)) < 0.1, (x, -10, 10), (y, -10, 10), plot_points=50)
Launched png viewer for Graphics object consisting of 1 graphics primitive
sage: region_plot(lambda x, y: abs(fff(x + i * y)) < 1, (x, -10, 10), (y, -10, 10), plot_points=50)
Launched png viewer for Graphics object consisting of 1 graphics primitive
sage: region_plot(lambda x, y: abs(fff(x + i * y)) < 0.1, (x, -10, 10), (y, -10, 10), plot_points=100)
Launched png viewer for Graphics object consisting of 1 graphics primitive
