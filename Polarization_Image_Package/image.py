import numpy as np
import datetime
from Sun_Ephemeris_Package.sunEphemeris import SunEphemeris
from Camera_Package.camera import Camera


class Image:
    """
    A class used to represent an image

    Attributes
    ----------
    image : float matrix
        the pixel image
    camera : Camera
        the camera used to take the image
    date : str
        date of the image in format "YYYY-MM-DD"
    hour_utc : int
        UTC hour of the image
    minute_utc :int
        UTC minute of the image
    second_utc :int
        UTC second of the image
    camera_exposure_time : float
        camera exposure time used to take the image, in microseconds
    camera_latitude : float
        latitude of the camera when the image was taken, in degrees
    camera_longitude : float
        longitude of the camera when the image was taken, in degrees

    Methods
    -------
    get_image_size()
        Computes the image size
    get_sun_horizontal_coordinates()
        Computes the position of the Sun in the horizontal frame at the moment when the image was taken
    get_sun_equatorial_coordinates_ha_dec()
        Computes the position of the Sun in the equatorial frame at the moment when the image was taken (declination
        and hour angle)
    get_sun_equatorial_coordinates_ra_dec()
        Computes the position of the Sun in the equatorial frame at the moment when the image was taken (declination
        and right ascension)
    has_same_camera(image)
        Test if two images are taken with the same camera
    has_same_position(image)
        Test if two images are taken at the same position
    """

    def __init__(self, image, camera=None, date=None, hour_utc=None, minute_utc=None, second_utc=None,
                 camera_exposure_time=None, camera_latitude=None, camera_longitude=None):
        """
        Parameters
        ----------
        image : ndarray
            the pixel image
        camera : Camera
            the camera used to take the image
        date : str
            date of the image in format "YYYY-MM-DD"
        hour_utc : int
            UTC hour of the image
        minute_utc :int
            UTC minute of the image
        second_utc :int
            UTC second of the image
        camera_exposure_time : float
            camera exposure time used to take the image in microseconds
        camera_latitude : float
            latitude of the camera when the image was taken, in degrees
        camera_longitude : float
            longitude of the camera when the image was taken, in degrees

        Raises
        ------
        TypeError
            if image is not a numpy ndarray or a list, or if camera is not a Camera, or if hour, minutes and seconds are
            not integers, or if camera exposure time, camera latitude and longitude are not float numbers
        ValueError
            if date has incorrect format, different from "YYYY-MM-DD"
        """
        if not isinstance(image, (list, np.ndarray)):
            raise TypeError("image attribute should be a numpy.ndarray or list object")
        if not isinstance(camera, Camera) and camera is not None:
            raise TypeError("camera attribute should be a Camera object")
        if date is not None:
            try:
                datetime.date.fromisoformat(date)
            except ValueError:
                raise ValueError("Incorrect date format, should be YYYY-MM-DD")
        if not isinstance(hour_utc, int) and hour_utc is not None:
            raise TypeError("hour_utc attribute should be an integer")
        if not isinstance(minute_utc, int) and minute_utc is not None:
            raise TypeError("minute_utc attribute should be an integer")
        if not isinstance(second_utc, int) and second_utc is not None:
            raise TypeError("second_utc attribute should be an integer")
        if not isinstance(camera_exposure_time, (int, float)) and camera_exposure_time is not None:
            raise TypeError("camera_exposure_time attribute should be a float object")
        if not isinstance(camera_latitude, (int, float)) and camera_latitude is not None:
            raise TypeError("camera_latitude attribute should be a float object")
        if not isinstance(camera_longitude, (int, float)) and camera_longitude is not None:
            raise TypeError("camera_longitude attribute should be a float object")
        self.image = image
        self.camera = camera
        self.date = date
        self.hour_utc = hour_utc
        self.minute_utc = minute_utc
        self.second_utc = second_utc
        self.camera_exposure_time = camera_exposure_time
        self.camera_latitude = camera_latitude
        self.camera_longitude = camera_longitude

    def get_image_size(self):
        """
        Computes the image size

        Parameters
        ----------

        Returns
        -------
        tuple of int
            the height and width of the image
        """
        return len(self.image), len(self.image[0])

    def get_sun_horizontal_coordinates(self):
        """
        Computes the position of the Sun in the horizontal frame when the image was taken

        Parameters
        ----------

        Returns
        -------
        tuple of float
            the sun coordinates in the horizontal frame in degrees (sun azimuth, sun altitude)

        Raises
        ------
        AttributeError
            the date, hour, minutes, seconds, of the image and latitude, and longitude of the camera should be provided
        """
        if self.date is None or self.hour_utc is None or self.minute_utc is None or self.second_utc is None or \
                self.camera_latitude is None or self.camera_longitude is None:
            raise AttributeError("the date, hour, minutes, seconds, of the image and latitude, and longitude of the "
                                 "camera should be provided")

        sun = SunEphemeris(self.date, self.hour_utc, self.minute_utc, self.second_utc, self.camera_latitude,
                           self.camera_longitude)

        return sun.get_sun_horizontal_coordinates()

    def get_sun_equatorial_coordinates_ha_dec(self):
        """
        Computes the position of the Sun in the equatorial coordinate frame when the image was taken, in declination and
        hour angle

        Parameters
        ----------

        Returns
        -------
        tuple of float
           the sun coordinates in the equatorial frame in degrees (sun hour angle, sun declination)

        Raises
        ------
        AttributeError
            the date, hour, minutes, seconds, of the image should be provided
        """
        if self.date is None or self.hour_utc is None or self.minute_utc is None or self.second_utc is None or \
                self.camera_latitude is None or self.camera_longitude is None:
            raise AttributeError("the date, hour, minutes, seconds, of the image and latitude, and longitude of the "
                                 "camera should be provided")

        sun = SunEphemeris(self.date, self.hour_utc, self.minute_utc, self.second_utc, self.camera_latitude,
                           self.camera_longitude)
        return sun.get_sun_equatorial_coordinates_ha_dec()

    def get_sun_equatorial_coordinates_ra_dec(self):
        """
        Computes the position of the Sun in the equatorial coordinate frame when the image was taken, in declination and
        right ascension

        Parameters
        ----------

        Returns
        -------
        tuple of float
           the sun coordinates in the equatorial frame in degrees (sun right ascension, sun declination)

        Raises
        ------
        AttributeError
            the date, hour, minutes, seconds, of the image should be provided
        """
        if self.date is None or self.hour_utc is None or self.minute_utc is None or self.second_utc is None:
            raise AttributeError("the date, hour, minutes, seconds, of the image should be provided")

        sun = SunEphemeris(self.date, self.hour_utc, self.minute_utc, self.second_utc, self.camera_latitude,
                           self.camera_longitude)
        return sun.get_sun_equatorial_coordinates_ra_dec()

    def has_same_camera(self, image):
        """
        Test if two images are taken with the same camera

        Parameters
        ----------
        image : Image
            an image to compare the corresponding camera with

        Returns
        -------
        bool
            True if both images were taken with the same camera, else False

        Raises
        ------
        TypeError
            the input image should be an Image object
        AttributeError
            the camera of the images should be provided
        """
        if not isinstance(image, Image):
            raise TypeError("the input image is not an Image object")
        if self.camera is None or image.camera is None:
            raise AttributeError("Camera should be provided")
        if image.camera != self.camera:
            return False
        return True

    def has_same_position(self, image):
        """
        Test if two images are taken at the same position

        Parameters
        ----------
        image : Image
            an image to compare the corresponding position with

        Returns
        -------
        bool
            True if both images were taken at the same position, else False

        Raises
        ------
        TypeError
            the input image should be an Image object
        AttributeError
            the camera longitude and latitude of the images should be provided
        """
        if not isinstance(image, Image):
            raise TypeError("the input image is not an Image object")
        if self.camera_latitude is None or self.camera_longitude is None or image.camera_latitude is None \
                or image.camera_longitude is None:
            raise AttributeError("Camera longitude and latitude should be provided")
        if image.camera_latitude != self.camera_latitude or image.camera_longitude != self.camera_longitude:
            return False
        return True


if __name__ == "__main__":
    from Camera_Package.sensor import Sensor
    from Camera_Package.dofpSensor import DoFPSensor
    from Camera_Package.lens import Lens
    from Camera_Package.camera import Camera
    import numpy as np

    sensor_pixel_size1 = 3.45
    sensor_pixel_size2 = 3.2
    sensor_dimensions_row_column1 = (2048, 2448)
    sensor_type1 = "monochrome"

    sensor1 = Sensor(sensor_pixel_size1, sensor_dimensions_row_column1, sensor_type1)
    sensor2 = Sensor(sensor_pixel_size2, sensor_dimensions_row_column1, sensor_type1)

    lens_focal_length1 = 1.8
    lens_conjugation_model1 = 'r2'

    lens1 = Lens(lens_focal_length1, lens_conjugation_model1)

    camera_center_coordinates1 = (1024, 1224)

    # Create a camera
    cam1 = Camera(sensor1, lens1, camera_center_coordinates1)
    cam2 = Camera(sensor2, lens1, camera_center_coordinates1)

    cam_lon1 = 10
    cam_lon2 = 20
    cam_lat1 = 30
    cam_lat2 = 40

    image1 = Image(image=np.zeros((10, 10)), camera=cam1, camera_latitude=cam_lat1, camera_longitude=cam_lon1)
    image2 = Image(image=np.ones((10, 10)), camera=cam2, camera_latitude=cam_lat1, camera_longitude=cam_lon1)
    image3 = Image(image=np.ones((10, 10)), camera=cam2, camera_latitude=cam_lat2, camera_longitude=cam_lon1)
    image4 = Image(image=np.ones((10, 10)), camera=cam1, camera_latitude=cam_lat1, camera_longitude=cam_lon2)
    image5 = Image(image=np.ones((10, 10)), camera=cam1, camera_latitude=cam_lat2, camera_longitude=cam_lon2)

    print(image1.has_same_position(image2))
    print(image1.has_same_position(image3))
    print(image1.has_same_position(image4))
    print(image1.has_same_camera(image2))
    print(image1.has_same_camera(image3))
    print(image1.has_same_camera(image4))
    print(image1.has_same_position(image5))
