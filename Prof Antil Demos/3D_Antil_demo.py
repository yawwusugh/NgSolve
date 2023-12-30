from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from ngsolve import *
# mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))
mesh = Mesh(unit_cube.GenerateMesh(maxh=0.2))
Draw (mesh)

Sigma = HDiv(mesh, order=2)
V = L2(mesh, order=1)
X = Sigma*V

sigma,u = X.TrialFunction()
tau,v = X.TestFunction()

a = BilinearForm(X)
a += (sigma*tau+div(sigma)*v+div(tau)*u)*dx

f = LinearForm(X)
f += -1*v*dx

a.Assemble()
f.Assemble()

gfu = GridFunction(X)
gfu.vec.data = a.mat.Inverse() * f.vec

Draw (gfu.components[0], mesh, "flux")
#Draw (gfu.components[1], mesh, "u")

n = specialcf.normal(mesh.dim)
Integrate (gfu.components[0]*n, mesh, BND)

