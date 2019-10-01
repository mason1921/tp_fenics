import dolfin
import fenics

DIM = 2

def strain_displacement(u):
    return dolfin.sym(dolfin.grad(u))

def stress_strain(eps, lambda_, mu):
    I = dolfin.Identity(DIM)
    return lambda_*dolfin.tr(eps)*I+2*mu*eps

class Cantilever:
    def __init__(self, d, h, ell=1., degree=1, nu=0.3, E=1.0):
        self.d = d
        self.ell = ell
        self.degree = degree
        mesh = fenics.RectangleMesh(fenics.Point(0., -0.5*d),
                                    fenics.Point(ell, 0.5*d),
                                    int(ell/h), int(d/h))
        self.left = dolfin.CompiledSubDomain('on_boundary && x[0] <= atol',
                                             atol=1e-5*h)
        self.right = dolfin.CompiledSubDomain('on_boundary && x[0] >= ell-atol',
                                              ell=ell,
                                              atol=1e-5*h)
        element = dolfin.VectorElement('P',
                                       cell=mesh.ufl_cell(),
                                       degree=degree,
                                       dim=DIM)
        self.V = dolfin.FunctionSpace(mesh, element)
        self.lambda_ = dolfin.Constant(E*nu/(1.+nu)/(1.-2.*nu))
        self.mu = dolfin.Constant(E/2./(1.+nu))

    def bilinear_form(self, u, v):
        return dolfin.inner(stress_strain(strain_displacement(u), self.lambda_, self.mu),
                            strain_displacement(v))*dolfin.dx

    def linear_form(self, v):
        # expr = 'x[0] < x_right ? 0.0 : -12.*M/d/d/d*x[1]'
        # traction = dolfin.Expression((expr, '0.0'),
        #                              x_right=self.ell-self.left.get_property('atol'),
        #                              M=M, d=self.d, degree=1)
        # return dolfin.dot(traction, v)*dolfin.ds
        g = dolfin.Constant((0., 1.))
        return dolfin.dot(g, v)*dolfin.dx


    def variational_problem(self):
        u = dolfin.TrialFunction(self.V)
        v = dolfin.TestFunction(self.V)
        lhs = self.bilinear_form(u, v)
        rhs = self.linear_form(v)
        u_prescribed = dolfin.Constant((0., 0.))
        bcs = [dolfin.DirichletBC(self.V, u_prescribed, self.left)]
        u = dolfin.Function(self.V)
        return dolfin.LinearVariationalProblem(lhs, rhs, u, bcs)

    def solution(self):
        problem = self.variational_problem()
        solver = dolfin.LinearVariationalSolver(problem)
        solver.solve()
        return problem.u_ufl
