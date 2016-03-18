
fontsize = 26
length = 2
radius = 0.3

edge_inc = line([(-length, 0), (-radius,0)])
edge_inc_inf = line([(-length - 0.2, 0), (-length, 0)], linestyle='dotted')
edge_out = line([(radius, 0), (length, 0)])
edge_out_inf = line([(length, 0), (length + 0.2, 0)], linestyle='dotted')
edge_circle = circle((0, 0), radius)
vertex_inc = point((-radius, 0), size=20)
vertex_out = point((radius, 0), size=20)
vertex_inc_label = text("$V_{inc}$", (-radius - 0.11, 0.1), fontsize=fontsize)
vertex_out_label = text("$V_{out}$", (radius + 0.1, 0.1), fontsize=fontsize)
edge_inc_label = text("$\Omega_{inc}$", (-(length - 0.3) - 0.2, 0.1), fontsize=fontsize)
edge_out_label = text("$\Omega_{out}$", ((length - 0.3) + 0.2, 0.1), fontsize=fontsize)
edge_up_label = text("$\Omega_1$", (0, radius - 0.1), fontsize=fontsize)
edge_down_label = text("$\Omega_2$", (0, -radius + 0.1), fontsize=fontsize)
objects = [
    edge_inc, edge_inc_inf, edge_inc_label,
    edge_circle, edge_up_label, edge_down_label,
    edge_out, edge_out_inf, edge_out_label,
    vertex_inc, vertex_inc_label,
    vertex_out, vertex_out_label,
]
sum(objects).save('graph.eps', axes=False, figsize=12)

