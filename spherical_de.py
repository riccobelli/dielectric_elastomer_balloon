# Copyright (C) 2023 Davide Riccobelli
#
# This file is part of dielectric_elastomer_baloon library for FEniCS.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 2.1 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

from bifenics import (
    BifenicsProblem,
    ArclengthContinuation,
    parallel_max,
    parallel_min,
)

from dolfin import (
    FiniteElement,
    VectorElement,
    MixedElement,
    FunctionSpace,
    Function,
    conditional,
    split,
    det,
    inv,
    tr,
    sqrt,
    inner,
    Constant,
    Measure,
    interpolate,
    assemble,
    as_tensor,
    as_vector,
    dx,
    ln,
    TestFunction,
    derivative,
    DirichletBC,
    Expression,
    MeshFunction,
    SubDomain,
    near,
    Mesh,
    SpatialCoordinate,
    DOLFIN_PI,
    parameters,
    File,
    project,
    FacetNormal,
)
import numpy as np
import scipy.special as spl
import scipy.io
import os
import ufl

os.environ["OMP_NUM_THREADS"] = " 1"
parameters["form_compiler"]["quadrature_degree"] = 10


class ElectroactiveSphere(BifenicsProblem):
    class Outside(SubDomain):
        def inside(self, x, on_boundary):
            TOL = 1e-2
            return on_boundary and near(x[0] ** 2 + x[1] ** 2, 1, TOL)

    class Inside(SubDomain):
        def __init__(self, Ri):
            self.Ri = Ri
            SubDomain.__init__(self)

        def inside(self, x, on_boundary):
            TOL = 1e-2
            return on_boundary and near(x[0] ** 2 + x[1] ** 2, self.Ri**2, TOL)

    class Inner_Layer(SubDomain):
        def __init__(self, thickness_out):
            self.thickness_out = thickness_out
            SubDomain.__init__(self)

        def inside(self, x, on_boundary):
            TOL = 1e-4
            return x[0] ** 2 + x[1] ** 2 < (1 - self.thickness_out) ** 2 + TOL

    class Outer_Layer(SubDomain):
        def __init__(self, thickness_out):
            self.thickness_out = thickness_out
            SubDomain.__init__(self)

        def inside(self, x, on_boundary):
            TOL = 1e-4
            return x[0] ** 2 + x[1] ** 2 > (1 - self.thickness_out) ** 2 - TOL

    class Left(SubDomain):
        def inside(self, x, on_boundary):
            TOL = 1e-4
            return on_boundary and near(x[0], 0, TOL)

    def __init__(
        self,
        thickness_in=0.1,
        thickness_out=0.1,
        eps_out=1,
        Jm=10,
        mui=1,
        muo=1,
        K=1000,
        m=4,
        dh=1e-4,
        dielectric_out=True,
    ):
        self.thickness_in = thickness_in
        self.thickness_out = thickness_out
        self.Jm = Constant(Jm)
        self.G = conditional(
            Expression("x[0]*x[0]+x[1]*x[1]", degree=2) > (1 - thickness_out) ** 2,
            muo,
            mui,
        )
        self.K = Constant(K)
        self.eps_out = Constant(eps_out)
        self.m = m
        self.counter = 0
        self.dh = dh
        self.dielectric_out = dielectric_out
        self.export_data = {
            "deltaphi": np.array([]),
            "ro_min": np.array([]),
            "ro_max": np.array([]),
            "ri_max": np.array([]),
            "ri_min": np.array([]),
            "p0": np.array([]),
            "Vint": np.array([]),
            "energy": np.array([]),
            "lambda_i_min": np.array([]),
            "lambda_i_max": np.array([]),
            "u_max_o": np.array([]),
            "u_min_o": np.array([]),
        }

    def apply_perturbation(self, x, y):
        Pm = spl.legendre(self.m)
        dPm = np.polyder(Pm)
        THETA = np.arctan2(x, y)
        if self.m != 1:
            condition = True

        else:
            condition = (
                x**2 + y**2 <= (1.03 - self.thickness_in - self.thickness_out) ** 2
            )
        uR = self.dh * Pm(np.cos(THETA))
        uT = (
            -self.dh
            * dPm(np.cos(THETA))
            * np.sin(THETA)
            / (self.m * (self.m + 1)) ** 0.5
        )

        x_bar = np.where(condition, x + uR * np.sin(THETA) + uT * np.cos(THETA), x)
        y_bar = np.where(condition, y + uR * np.cos(THETA) - uT * np.sin(THETA), y)
        return [x_bar, y_bar]

    def mesh(self):
        mesh = Mesh("mesh.xml")

        self.boundaries = MeshFunction("size_t", mesh, 1)
        self.boundaries.set_all(0)
        inside = self.Inside(1 - self.thickness_in - self.thickness_out)
        outside = self.Outside()
        inner_layer = self.Inner_Layer(self.thickness_out)
        outer_layer = self.Outer_Layer(self.thickness_out)
        left = self.Left()
        inside.mark(self.boundaries, 1)
        outside.mark(self.boundaries, 2)
        left.mark(self.boundaries, 3)

        self.no_diel_layer_dof = MeshFunction("size_t", mesh, 1)
        self.no_diel_layer_dof.set_all(0)
        if self.dielectric_out is True:
            inner_layer.mark(self.no_diel_layer_dof, 1)
        else:
            outer_layer.mark(self.no_diel_layer_dof, 1)
        file_bound = File("output/boundaries.pvd")
        file_bound << self.no_diel_layer_dof

        self.domains = MeshFunction("size_t", mesh, mesh.topology().dim())
        self.domains.set_all(0)
        if self.dielectric_out is True:
            inner_layer.mark(self.domains, 1)
        else:
            outer_layer.mark(self.domains, 1)

        x = mesh.coordinates()[:, 0]
        y = mesh.coordinates()[:, 1]

        x_bar, y_bar = self.apply_perturbation(x, y)
        xy_bar_coor = np.array([x_bar, y_bar]).transpose()
        mesh.coordinates()[:] = xy_bar_coor

        self.my_mesh = mesh
        return self.my_mesh

    def function_space(self, mesh):
        VP2elem = VectorElement("CG", mesh.ufl_cell(), 2)
        P1elem = FiniteElement("CG", mesh.ufl_cell(), 1)
        Relem = FiniteElement("R", mesh.ufl_cell(), 0)
        elem = MixedElement([VP2elem, P1elem, Relem])
        V = FunctionSpace(mesh, elem)
        return V

    def parameters(self):
        deltaphi = Constant(0)
        p0 = Constant(0)
        return {"deltaphi": deltaphi, "p0": p0}

    def boundary_conditions(self, mesh, V):
        bciphi = DirichletBC(V.sub(1), Constant(0), self.boundaries, 1)
        bcophi = DirichletBC(V.sub(1), Constant(0), self.boundaries, 2)
        bcl = DirichletBC(V.sub(0).sub(0), Constant(0), self.boundaries, 3)
        bcnodielectric = DirichletBC(V.sub(1), Constant(0), self.no_diel_layer_dof, 1)
        bcs = [bcl, bciphi, bcophi, bcnodielectric]
        return bcs

    def grad_cyl(self, u, X):
        grad_u = as_tensor(
            [
                [u[0].dx(0), u[0].dx(1), 0],
                [u[1].dx(0), u[1].dx(1), 0],
                [0, 0, u[0] / X[0]],
            ]
        )
        return grad_u

    def residual(self, uphihalpha, vpsibeta, parameters):
        deltaphi = parameters["deltaphi"]
        p0 = parameters["p0"]
        u, phih, alpha = split(uphihalpha)
        v, psi, beta = split(vpsibeta)
        X = SpatialCoordinate(self.my_mesh)
        if self.dielectric_out is True:
            linear_ramp_phi = (
                deltaphi
                * ((X[0] * X[0] + X[1] * X[1]) ** 0.5 - (1 - self.thickness_out))
                / self.thickness_out
            )
            self.phibc = conditional(
                X[0] * X[0] + X[1] * X[1] > (1 - self.thickness_out) ** 2,
                linear_ramp_phi,
                0,
            )
        else:
            linear_ramp_phi = (
                (X[0] * X[0] + X[1] * X[1]) ** 0.5
                - (1 - self.thickness_out - self.thickness_in)
            ) / self.thickness_in
            self.phibc = conditional(
                X[0] * X[0] + X[1] * X[1] < (1 - self.thickness_out) ** 2,
                deltaphi * linear_ramp_phi,
                deltaphi,
            )

        phi = self.phibc + phih

        R, Z = X[0], X[1]
        r, z = R + u[0], Z + u[1]

        Jgeo = 2 * DOLFIN_PI * R

        F = as_tensor(
            [[r.dx(0), r.dx(1), 0.0], [z.dx(0), z.dx(1), 0.0], [0.0, 0.0, r / R]]
        )

        self.deformation_gradient = F

        invF = inv(F)
        J = det(F)
        C = F.T * F
        invC = inv(C)
        Fbar = F * J ** (-1 / 3)
        Cbar = Fbar.T * Fbar
        I1bar = tr(Cbar)
        I1 = tr(C)

        # Electric
        Er = -as_vector([phi.dx(0), phi.dx(1), 0])
        E = inv(F.T) * Er

        # Stress
        elastic_energy_dens = (
            -1 / 2 * self.G * self.Jm * ln(1 - (I1bar - 3) / self.Jm)
            + 1 / 2 * self.G * self.K * ln(J) ** 2
        )

        electrical_energy_dens = -J * self.eps_out / 2 * inner(Er, invC * Er)

        zero_mean_disp = alpha * u[1]

        my_dx = Measure("dx", domain=self.my_mesh, subdomain_data=self.domains)
        self.energy = (
            Jgeo * elastic_energy_dens * dx
            + Jgeo * zero_mean_disp * dx
            + Jgeo * electrical_energy_dens * my_dx(0)
        )

        der = derivative(self.energy, uphihalpha, vpsibeta)

        N = FacetNormal(self.my_mesh)
        ext = -(
            Jgeo
            * p0
            * inner(
                as_vector([v[0], v[1], 0]), ufl.cofac(F) * as_vector([N[0], N[1], 0])
            )
        )
        my_ds = Measure("ds", domain=self.my_mesh, subdomain_data=self.boundaries)

        FF = der - ext * my_ds(1)

        return FF

    def monitor(self, uphihalpha, parameters, xdmf_file):
        deltaphi = parameters["deltaphi"]
        p0 = parameters["p0"]
        uexp = Function(uphihalpha, 0, name="displacement")
        phiexp = Function(uphihalpha, 1, name="phi")
        Vtmp = FunctionSpace(self.my_mesh, "CG", 1)
        realphi = project(phiexp + self.phibc, Vtmp)
        phibc = project(self.phibc, Vtmp)
        phiexp.assign(realphi)
        phibc.rename("phibc", "phibc")

        t = float(deltaphi)
        t = round(abs(t), 10)
        self.counter += 1

        xdmf_file.write(uexp, self.counter)
        xdmf_file.write(phiexp, self.counter)
        xdmf_file.write(phibc, self.counter)

        X = SpatialCoordinate(self.my_mesh)
        x = X + uexp
        r = (x[0] * x[0] + x[1] * x[1]) ** 0.5
        rfunc = project(r, Vtmp)
        my_ds = Measure("ds", domain=self.my_mesh, subdomain_data=self.boundaries)
        a_function = interpolate(Constant(1), Vtmp)
        dof_surface_out = assemble(a_function * TestFunction(Vtmp) * my_ds(2))
        dof_surface_out_vect = dof_surface_out[:]
        r_surface_out = rfunc.vector()[dof_surface_out_vect != 0]
        romax = parallel_max(r_surface_out, self.my_mesh.mpi_comm())
        romin = parallel_min(r_surface_out, self.my_mesh.mpi_comm())
        dof_surface_in = assemble(a_function * TestFunction(Vtmp) * my_ds(1))
        dof_surface_in_vect = dof_surface_in[:]
        r_surface_in = rfunc.vector()[dof_surface_in_vect != 0]
        rimax = parallel_max(r_surface_in, self.my_mesh.mpi_comm())
        rimin = parallel_min(r_surface_in, self.my_mesh.mpi_comm())

        normu = sqrt(inner(uexp, uexp))
        normu_func = project(normu, Vtmp)
        normu_surface_out = normu_func.vector()[dof_surface_out_vect != 0]
        umaxo = parallel_max(normu_surface_out, self.my_mesh.mpi_comm())
        umino = parallel_min(normu_surface_out, self.my_mesh.mpi_comm())

        F = self.deformation_gradient
        lambda3 = F[2, 2]
        reduced_F = as_tensor([[F[0, 0], F[0, 1]], [F[1, 0], F[1, 1]]])
        C = reduced_F.T * reduced_F
        edelta = tr(C) * tr(C) - 4 * det(C)
        e1 = 1.0 / 2.0 * (tr(C) - sqrt(edelta))
        e2 = 1.0 / 2.0 * (tr(C) + sqrt(edelta))
        lambda1 = sqrt(e1)
        lambda2 = sqrt(e2)
        eigenvalues = as_vector([lambda1, lambda2, lambda3])

        DGtmp = FunctionSpace(self.my_mesh, "DG", 0)
        for j in range(3):
            expeigen = project(eigenvalues[j], DGtmp)
            expeigen.rename(f"lambda{j+1}", "")
            xdmf_file.write(expeigen, self.counter)

        lambda1fun = project(lambda1, DGtmp)
        lambda2fun = project(lambda2, DGtmp)
        lambda3fun = project(lambda3, DGtmp)
        a_function_DG = interpolate(Constant(1), DGtmp)
        dof_surface_in_DG = assemble(a_function_DG * TestFunction(DGtmp) * my_ds(1))
        dof_surface_in_vect_DG = dof_surface_in_DG[:]
        v1 = lambda1fun.vector()[dof_surface_in_vect_DG != 0]
        v2 = lambda2fun.vector()[dof_surface_in_vect_DG != 0]
        v3 = lambda3fun.vector()[dof_surface_in_vect_DG != 0]
        v = np.concatenate((v1, v2, v3))
        lam_min = parallel_min(v, self.my_mesh.mpi_comm())
        lam_max = parallel_max(v, self.my_mesh.mpi_comm())

        # Update arrays with results
        self.export_data["deltaphi"] = np.append(
            self.export_data["deltaphi"], float(deltaphi)
        )
        self.export_data["p0"] = np.append(self.export_data["p0"], float(p0))
        self.export_data["ro_min"] = np.append(self.export_data["ro_min"], romin)
        self.export_data["ro_max"] = np.append(self.export_data["ro_max"], romax)
        self.export_data["ri_min"] = np.append(self.export_data["ri_min"], rimin)
        self.export_data["ri_max"] = np.append(self.export_data["ri_max"], rimax)
        self.export_data["lambda_i_min"] = np.append(
            self.export_data["lambda_i_min"], lam_min
        )
        self.export_data["lambda_i_max"] = np.append(
            self.export_data["lambda_i_max"], lam_max
        )
        self.export_data["u_max_o"] = np.append(self.export_data["u_max_o"], umaxo)
        self.export_data["u_min_o"] = np.append(self.export_data["u_min_o"], umino)

        N = FacetNormal(self.my_mesh)
        N = as_vector([N[0], N[1], 0])

        Jgeo = 2 * DOLFIN_PI * X[0]
        internalVol = assemble(
            -1
            / 3
            * Jgeo
            * det(F)
            * inner(as_vector([x[0], x[1], 0]), inv(F.T) * N)
            * my_ds(1)
        )
        Ri = 1 - self.thickness_in - self.thickness_out
        initial_internal_volume = 4 / 3 * DOLFIN_PI * Ri**3
        deltaV = internalVol - initial_internal_volume
        total_energy = assemble(self.energy) - float(p0) * deltaV
        self.export_data["Vint"] = np.append(self.export_data["Vint"], deltaV)
        self.export_data["energy"] = np.append(self.export_data["energy"], total_energy)
        # Save to mat file
        scipy.io.savemat("output/data.mat", self.export_data)

    def ac_monitor(self, uphihalpha, arclength_param, parameters, xdmf_file):
        tmp_parameters = dict(parameters)
        value_type = [str(type(k)) for k in parameters.values()]
        idx = value_type.index("<class 'ufl.indexed.Indexed'>")
        key = list(parameters.keys())[idx]
        tmp_parameters[key] = arclength_param
        assert tmp_parameters[key] != parameters[key]
        self.monitor(uphihalpha, tmp_parameters, xdmf_file)

    def solver_parameters(self):
        parameters = {
            "nonlinear_solver": "snes",
            "snes_solver": {
                "linear_solver": "mumps",
                "absolute_tolerance": 1e-8,
                "relative_tolerance": 1e-8,
                "maximum_iterations": 6,
            },
        }
        return parameters
