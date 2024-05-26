import datetime
from astropy.coordinates import get_sun, AltAz, EarthLocation, HADec, PrecessedGeocentric
import astropy.coordinates as coord
from astropy.time import Time
import astropy.units as u


class SunEphemeris:
    """
    A class used to compute the position of the Sun

    Attributes
    ----------
    date : string
        date in format "YYYY-MM-DD"
    hour_utc : int
        UTC hour
    minute_utc :int
        UTC minute
    second_utc :int
        UTC second
    latitude : float
        latitude, in degrees
    longitude : float
        longitude, in degrees

    Methods
    -------
    get_sun_horizontal_coordinates()
        Computes the position of the Sun in the horizontal frame at a given moment in azimuth and elevation/altitude
    get_sun_equatorial_coordinates_ha_dec()
        Computes the position of the Sun in the equatorial frame at a given moment in hour angle and declination
    get_sun_equatorial_coordinates_ra_dec()
        Computes the position of the Sun in the equatorial frame at a given moment in right ascension and declination

    """

    def __init__(self, date, hour_utc, minute_utc, second_utc, latitude=None, longitude=None):
        """
        Parameters
        ----------
        date : string
            date in format "YYYY-MM-DD"
        hour_utc : int
            UTC hour
        minute_utc :int
            UTC minute
        second_utc :int
            UTC second
        latitude : float
            latitude, in degrees
        longitude : float
            longitude, in degrees

        Raises
        ------
        ValueError
            date format should be "YYYY-MM-DD"
        TypeError
            hour_utc, minute_utc and second_utc should be integer
            latitude and longitude should either be None or float
        """
        if date is not None:
            try:
                datetime.date.fromisoformat(date)
            except ValueError:
                raise ValueError("Incorrect date format, should be YYYY-MM-DD")
        if not isinstance(hour_utc, int):
            raise TypeError("hour_utc attribute should be an integer")
        if not isinstance(minute_utc, int):
            raise TypeError("minute_utc attribute should be an integer")
        if not isinstance(second_utc, int):
            raise TypeError("second_utc attribute should be an integer")
        if not isinstance(latitude, (int, float)) and latitude is not None:
            raise TypeError("camera_latitude attribute should be a float object")
        if not isinstance(longitude, (int, float)) and longitude is not None:
            raise TypeError("camera_longitude attribute should be a float object")
        self.date = date
        self.hour_utc = hour_utc
        self.minute_utc = minute_utc
        self.second_utc = second_utc
        self.latitude = latitude
        self.longitude = longitude
        self.year = datetime.date.fromisoformat(date).year
        self.month = datetime.date.fromisoformat(date).month
        self.day = datetime.date.fromisoformat(date).day
        self.datetime = datetime.datetime(self.year, self.month, self.day, self.hour_utc, self.minute_utc,
                                          self.second_utc)

    def get_sun_horizontal_coordinates(self):
        """
        Computes the position of the Sun in the horizontal frame at a given moment in azimuth and elevation/altitude

        Parameters
        ----------

        Returns
        -------
        sun_azimuth : float
           sun azimuth in the horizontal frame in degrees
        sun_altitude : float
            sun altitude/elevation in the horizontal frame in degrees

        Raises
        ------
        AttributeError
            longitude and latitude should be provided
        """
        if self.longitude is None or self.latitude is None:
            raise AttributeError("Longitude and latitude should be provided")
        # verifier si besoin de prendre en compte l'utc offset
        loc = coord.EarthLocation(lon=self.longitude * u.deg, lat=self.latitude * u.deg)
        time_date = Time(self.datetime.isoformat(), format='isot', scale='utc')
        altaz = coord.AltAz(location=loc, obstime=time_date)
        sun_altitude = get_sun(time_date).transform_to(altaz).alt.degree
        sun_azimuth = get_sun(time_date).transform_to(altaz).az.degree
        return sun_azimuth, sun_altitude

    def get_sun_equatorial_coordinates_ha_dec(self):
        """
        Computes the position of the Sun in the equatorial frame at a given moment in hour angle and declination

        Parameters
        ----------

        Returns
        -------
        sun_hour_angle : float
           sun hour angle in the equatorial frame in degrees
        sun_declination : float
            sun declination in the equatorial frame in degrees

        Raises
        ------
        AttributeError
            longitude and latitude should be provided
        """
        if self.longitude is None or self.latitude is None:
            raise TypeError("Longitude and latitude should be provided")
        # verifier si besoin de prendre en compte l'utc offset
        loc = coord.EarthLocation(lon=self.longitude * u.deg, lat=self.latitude * u.deg)
        time_date = Time(self.datetime.isoformat(), format='isot', scale='utc')
        hadec = coord.HADec(location=loc, obstime=time_date)
        sun_hour_angle = get_sun(time_date).transform_to(hadec).ha.degree
        sun_declination = get_sun(time_date).transform_to(hadec).dec.degree
        return sun_hour_angle, sun_declination

    def get_sun_equatorial_coordinates_ra_dec(self):
        """
        Computes the position of the Sun in the equatorial frame at a given moment in right ascension and declination

        Parameters
        ----------

        Returns
        -------
        sun_right_ascension : float
           sun right ascension in the equatorial frame in degrees
        sun_declination : float
            sun declination in the equatorial frame in degrees
        """
        time_date = Time(self.datetime.isoformat(), format='isot', scale='utc')
        # we use PrecessedGeocentric frame, accounting for Earth precession (not nutation) for better results
        # (cf : https://github.com/astropy/astropy/issues/10713)
        pre_geo = PrecessedGeocentric(equinox=time_date)
        sun_right_ascension = get_sun(time_date).transform_to(pre_geo).ra.degree
        sun_declination = get_sun(time_date).transform_to(pre_geo).dec.degree
        return sun_right_ascension, sun_declination


if __name__ == "__main__":
    sun = SunEphemeris("2024-03-18", 19, 39, 36, 46.67747, 13.90172)
    print(sun.get_sun_horizontal_coordinates())
    print(sun.get_sun_equatorial_coordinates_ra_dec())
    print(sun.get_sun_equatorial_coordinates_ha_dec())
