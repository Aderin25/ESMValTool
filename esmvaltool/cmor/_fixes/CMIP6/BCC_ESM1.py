# pylint: disable=invalid-name, no-self-use, too-few-public-methods
"""Fixes for BCC-ESM1."""
import numpy as np
import iris

from ..fix import Fix


class allvars(Fix):
    """Common fixes to all vars"""

    def fix_metadata(self, cubes):
        """
        Fix metadata.

        Fixes error in time coordinate, sometimes contains trailing zeros

        Parameters
        ----------
        cube: iris.cube.CubeList

        Returns
        -------
        iris.cube.CubeList

        """
        # time
        for cube in cubes:
            try:
                old_time = cube.coord('time')
                if old_time.is_monotonic():
                    pass

                time_units = old_time.units
                time_data = old_time.points

                idx_zeros = np.where(time_data == 0.0)[0]
                time_diff = time_units.num2date(time_data[1]) \
                    - time_units.num2date(time_data[0])
                days = time_diff.days

                for idx in idx_zeros:
                    if idx == 0:
                        continue
                    correct_time = time_units.num2date(time_data[idx - 1])
                    if days <= 31 and days >= 28:  # assume monthly time steps
                        new_time = \
                            correct_time.replace(month=correct_time.month + 1)
                    else:  # use "time[1] - time[0]" as step
                        new_time = correct_time + time_diff
                    old_time.points[idx] = time_units.date2num(new_time)

                # create new time bounds
                old_time.bounds = None
                old_time.guess_bounds()

                # replace time coordinate with "repaired" values
                new_time = iris.coords.DimCoord.from_coord(old_time)
                time_idx = cube.coord_dims(old_time)
                cube.remove_coord('time')
                cube.add_dim_coord(new_time, time_idx)

            except iris.exceptions.CoordinateNotFoundError:
                pass

        # Remove 1D lat and lon as they are incorrect, this is an irregular
        # 2D grid and lat/lon need to be 2d.
        for cube in cubes:
            coords_to_remove = []
            for coord in cube.coords():
                if coord.var_name in ['time', ]: continue
                if coord.var_name in ['lat', 'lon']:
                    coords_to_remove.append(coord)

            for coord in coords_to_remove:
                cube.remove_coord(coord)

            latitudes = cube.coord('latitude')
            longitudes = cube.coord('longitude')
            latitudes.var_name = 'lat'
            longitudes.var_name = 'lon'

        return cubes
