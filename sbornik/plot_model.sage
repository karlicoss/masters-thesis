fontsize = 26
length = 2
radius = 0.3

color='black'

edge_L = line([(-length, 0), (0, 0)], rgbcolor=color)
edge_L_inf = line([(-length - 0.2, 0), (-length, 0)], linestyle='dotted', rgbcolor=color)
edge_R = line([(0, 0), (length, 0)], rgbcolor=color)
edge_R_inf = line([(length, 0), (length + 0.2, 0)], linestyle='dotted', rgbcolor=color)
edge_circle = circle((0, radius), radius, rgbcolor=color)
vertex = point((0, 0), size=20, rgbcolor=color)
vertex_label = text("$V$", (0 - 0.11, -0.1), fontsize=fontsize, rgbcolor=color)
edge_L_label = text("$\Omega_L$", (-(length - 0.3) - 0.2, 0.1), fontsize=fontsize, rgbcolor=color)
edge_R_label = text("$\Omega_R$", ((length - 0.3) + 0.2, 0.1), fontsize=fontsize, rgbcolor=color)
edge_circle_label = text("$\Omega$", (0, 2 * radius - 0.1), fontsize=fontsize, rgbcolor=color)
objects = [
    edge_L, edge_L_inf, edge_L_label,
    edge_circle, edge_circle_label,
    edge_R, edge_R_inf, edge_R_label,
    vertex, vertex_label,
]
sum(objects).save('graph.eps', axes=False, figsize=12)