import sympy as sp
import numpy as np


a = sp.Symbol("a")
b = sp.Symbol("b")
c = sp.Symbol("c")
d = sp.Symbol("d")

x = sp.Symbol("x")
y = sp.Symbol("y")
z = sp.Symbol("z")

b2 = sp.Matrix([1, 0])
b0 = sp.Matrix([-1/2, np.sqrt(3)/2])
b1 = sp.Matrix([-1/2, -np.sqrt(3)/2])

phi0 = sp.Matrix([a + c*x, b + c*y])
phi1 = sp.Matrix([a + c*x, b + c*y])
phi2 = sp.Matrix([a + c*x, b + c*y])

def integrate_against_basis(phi):
    print(phi[0])
    print(b2[0])
    phib2 = sp.integrate((phi[0]*b2[1] - phi[1]*b2[0]).subs({y: -np.sqrt(3) / 3}), (x, -1, 1))
    phib1 = sp.integrate((phi[0]*b1[1] - phi[1]*b1[0]).subs({x: (- z - 1) / 2, y: np.sqrt(3) * (3 * -z + 1) / 6}), (z, -1, 1))
    phib0 = sp.integrate((phi[0]*b0[1] - phi[1]*b0[0]).subs({x: (1 - z) / 2, y: np.sqrt(3) * (3 * z + 1) / 6}), (z, -1, 1))

    return phib0, phib1, phib2


# (phi1b0, phi1b1, phi1b2) = integrate_against_basis(phi1)
# print(phi1b2)
# print(phi1b0)
# print(phi1b1)

# soln1, = sp.linsolve([phi1b2 - 1, phi1b0, phi1b1], (a, b, c))
# soln2, = sp.linsolve([phi1b2, phi1b0 - 1, phi1b1], (a, b, c))
# soln3, = sp.linsolve([phi1b2, phi1b0, phi1b1 - 1], (a, b, c))

# print(soln1)
# print(soln2)
# print(soln3)



# g2 = (z, -np.sqrt(3) / 3)
# g1 = ((- z - 1) / 2, np.sqrt(3) * (3 * -z + 1) / 6)
# g0 = ((1 - z) / 2, np.sqrt(3) * (3 * z + 1) / 6)

# aa, bb, cc = soln1
# v2 = phi0.subs({a: aa, b: bb, c: cc})
# attached = v2.subs({x: g2[0], y: g2[1]})
# print(sp.integrate(attached.dot(b2), (z, -1, 1)))
# print(sp.integrate(v2.dot(b2).subs({y: -np.sqrt(3) / 3}), (x, -1, 1)))
# print(sp.integrate(attached.dot(b1), (z, -1, 1)))
# print(sp.integrate(attached.dot(b0), (z, -1, 1)))
# print(integrate_against_basis(v2))

# print("2")
# d, e, f = soln2
# v1 = phi0.subs({a: d, b: e, c: f})
# attached = v1.subs({x: g1[0], y: g1[1]})
# # print(g1[0].subs({z: -1}))
# # print(g1[0].subs({z: 1}))
# # print(g1[1].subs({z: -1}))
# # print(g1[1].subs({z: 1}))

# print(attached.dot(b1))
# print(sp.integrate(attached.dot(b1), (z, -1, 1)))
# print(sp.integrate(v1.dot(b1), (x, 0, -1), (y, 2*np.sqrt(3)/3, -np.sqrt(3)/3)))

# print("3")
# g, h, i = soln3
# v0 = phi0.subs({a: g, b: h, c: i})
# attached = v0.subs({x: g0[0], y: g0[1]})
# # print(g0[0].subs({z: -1}))
# # print(g0[0].subs({z: 1}))
# # print(g0[1].subs({z: -1}))
# # print(g0[1].subs({z: 1}))
# print(sp.integrate(attached.dot(b0), (z, -1, 1)))
# print(sp.integrate(v0.dot(b0), (x, 1, 0), (y, -np.sqrt(3)/3, 2*np.sqrt(3)/3)))

# P = sp.MatrixSymbol("P", 10, 6)
# x = sp.Symbol("x")
# y = sp.Symbol("y")

# coords = [(-1, -np.sqrt(3)/3), (1, -np.sqrt(3)/3), (0, 2*np.sqrt(3)/3),
#           (0, -np.sqrt(3)/3), (-0.5, np.sqrt(3)/6), (-0.5, np.sqrt(3)/6),
#           (0, 0)]

# vander = sp.Matrix([[((x**i)*(y**j)).subs({x: x1, y: y1})
#                     for i in range(4) for j in range(4) if i + j <= 3] for (x1, y1) in coords])
# print(vander.shape)
# b = sp.Matrix([2, 0])
# print(b.shape)

# V = vander
# # V = np.array([vander.subs({x: x1, y: y1})[0].T for (x1, y1) in coords])
# I_1 = sp.Matrix(np.eye(10)[0])
# print(I_1)
# print(V.shape)
# # print(V.inv())
# # print(sp.MatMul(V, P))
# print(V.solve(I_1))
# # print(sp.solve(I = sp.MatMul(P, V)))

twod = np.array([[1, -1, -np.sqrt(3)/3],
                 [1, 1, -np.sqrt(3)/3],
                 [1, 0, 2*np.sqrt(3)/3],])

# A = [-1, -np.sqrt(3)/3, -np.sqrt(3)/3]
# B = [1, -np.sqrt(3)/3, -np.sqrt(3)/3]
# C = [0, 0, 3*np.sqrt(2)/4]
# D = [0, 2*np.sqrt(3)/3, -np.sqrt(3)/3]

A = [-1, 0, -1/np.sqrt(2)]
B = [1, 0, -1/np.sqrt(2)]
C = [0, 1, 1/np.sqrt(2)]
D = [0, -1, 1/np.sqrt(2)]

x = sp.Symbol("x")
y = sp.Symbol("y")
xy = sp.Matrix([1, x, y])

# threed1a = np.array([[0, 0, 3*np.sqrt(2)/4],
#                    [0, 2*np.sqrt(3)/3, -np.sqrt(3)/3],
#                    [-1, -np.sqrt(3)/3, -np.sqrt(3)/3]])
threed1 = np.array([D, A, C])
threed2 = np.array([A, B, D])
threed3 = np.array([A, B, C])
threed4 = np.array([B, D, C])
# print(threed1a)
# print(threed1)
# print(threed2)
# print(threed3)
# print(threed4)

res1 = np.linalg.solve(twod, threed1)
res2 = np.linalg.solve(twod, threed2)
res3 = np.linalg.solve(twod, threed3)
res4 = np.linalg.solve(twod, threed4)

print(xy)
print(xy.T * res1)
print(xy.T * res2)
print(xy.T * res3)
print(xy.T * res4)


func1 =lambda x, y: [-0.577350269189626*y - 1/3,
                                     -x,
                                     0.235702260395516 - 0.816496580927726*y]
i = 0 
for a, x, y in twod:
    print(threed1[i])
    print(func1(x, y))
    i += 1 

# func2 = lambda x, y: [-0.5*x - 0.866025403784439*y,
#                       0.866025403784439*x - 0.5*y ,
#                     -0.577350269189626]

# i = 0 
# for a, x, y in twod:
#     print(threed2[i])
#     print(func2(x, y))
#     i += 1 