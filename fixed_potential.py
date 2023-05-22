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

from spherical_de import ElectroactiveSphere
from bifenics import ArclengthContinuation, ParameterContinuation
from dolfin import Function, project, parameters, Constant

parameters["allow_extrapolation"] = False
XDMF_options = {
    "flush_output": True,
    "functions_share_mesh": True,
    "rewrite_function_mesh": False,
}

td = 0.9
te = 0.9
Rm = te
Ri = Rm * td
t1 = 1 - Rm
t2 = Rm - Ri
H = t1 + t2
dh = 1e-5
K = 1000


sphere = ElectroactiveSphere(
    thickness_in=t2,
    thickness_out=t1,
    mui=1,
    muo=10,
    eps_out=1,
    dh=dh,
    m=4,
    K=K,
    Jm=97.2,
    dielectric_out=False,
)

analysis = ParameterContinuation(
    sphere,
    "deltaphi",
    start=0,
    end=0.3 * t2,
    dt=0.25,
    saving_file_parameters=XDMF_options,
    save_output=False,
)
actual_state = analysis.run()


class ElectroactiveSphereRestart(ElectroactiveSphere):
    def mesh(self):
        self.my_mesh = sphere.my_mesh
        self.boundaries = sphere.boundaries
        self.no_diel_layer_dof = sphere.no_diel_layer_dof
        self.domains = sphere.domains
        return self.my_mesh

    def initial_guess(self, V):
        initial_guess = Function(V)
        projected_actual_state = project(actual_state[0], V)
        initial_guess.assign(projected_actual_state)
        return initial_guess

    def parameters(self):
        deltaphi = actual_state[1]["deltaphi"]
        p0 = actual_state[1]["p0"]
        return {"deltaphi": Constant(float(deltaphi)), "p0": Constant(float(p0))}


restartsphere = ElectroactiveSphereRestart(
    thickness_in=t2,
    thickness_out=t1,
    mui=1,
    muo=10,
    dh=dh,
    m=4,
    K=K,
    Jm=97.2,
    dielectric_out=False,
)

analysis = ArclengthContinuation(
    restartsphere,
    "p0",
    start=0,
    end=-1,
    ds=0.01,
    initial_direction=-1,
    saving_file_parameters=XDMF_options,
    save_output=False,
    first_step_with_parameter_continuation=False,
    max_steps=1000,
    predictor_type="secant",
    n_step_for_doubling=2,
)
analysis.run()
