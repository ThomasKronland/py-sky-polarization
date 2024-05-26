import numpy as np

"""
Methods
-------
- sph2cart(azimuth, elevation, r)
    Transform spherical coordinates to cartesian coordinates
"""


def sph2cart(azimuth, elevation, r):
    """
    Transform spherical coordinates to cartesian coordinates

    Parameters
    ----------
    azimuth   : float number
        The azimuth coordinate in radians
    elevation : float number
        The elevation coordinate in radians
    r         : float number
        The distance coordinate

    Returns
    -------
    x, y, z : float triplet
        The cartesian coordinates

    Raises
    ------
    TypeError
        if azimuth, elevation or r is not a float number

    Example
    --------
    az = np.pi/4
    el = np.pi/6
    [OMx,OMy,OMz]=sph2cart(az,el,1)
    """
    if not isinstance(azimuth, (int, float)) or not isinstance(elevation, (int, float)) \
            or not isinstance(r, (int, float)):
        raise TypeError("the inputs should all be a float objects")
    x = r * np.cos(elevation) * np.cos(azimuth)
    y = r * np.cos(elevation) * np.sin(azimuth)
    z = r * np.sin(elevation)
    return x, y, z
