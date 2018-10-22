"""Derivation of variable `clhtkisccp`."""


import iris
from iris import Constraint

from ._derived_variable_base import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `clhtkisccp`."""

    def get_required(self, frequency):
        """Get variable `short_name` and `field` pairs required for derivation.

        Parameters
        ----------
        frequency : str
            Frequency of the desired derived variable.

        Returns
        -------
        list of tuples
            List of tuples (`short_name`, `field`) of all variables required
            for derivation.

        """
        return [('clisccp', 'T4' + frequency)]

    def calculate(self, cubes, fx_files=None):
        """Compute ISCCP high level thick cloud area fraction.

        Parameters
        ----------
        cubes : iris.cube.CubeList
            `CubeList` containing `clisccp` (`isccp_cloud_area_fraction`).
        fx_files : dict, optional
            If required, dictionary containing fx files  with `short_name`
            (key) and path (value) of the fx variable.

        Returns
        -------
        iris.cube.Cube
            `Cube` containing ISCCP high level thick cloud area fraction.

        """
        clisccp_cube = cubes.extract_strict(
            Constraint(name='isccp_cloud_area_fraction'))

        tau = iris.Constraint(
            atmosphere_optical_thickness_due_to_cloud=lambda t: t > 23.)
        plev = iris.Constraint(air_pressure=lambda p: p <= 44000.)
        clhtkisccp_cube = clisccp_cube
        clhtkisccp_cube = clhtkisccp_cube.extract(tau & plev)
        coord_names = [
            coord.standard_name for coord in clhtkisccp_cube.coords()
            if len(coord.points) > 1
        ]
        if 'atmosphere_optical_thickness_due_to_cloud' in coord_names:
            clhtkisccp_cube = clhtkisccp_cube.collapsed(
                'atmosphere_optical_thickness_due_to_cloud', iris.analysis.SUM)
        if 'air_pressure' in coord_names:
            clhtkisccp_cube = clhtkisccp_cube.collapsed('air_pressure',
                                                        iris.analysis.SUM)

        return clhtkisccp_cube