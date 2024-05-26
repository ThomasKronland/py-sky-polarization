from Sun_Ephemeris_Package.sunEphemeris import SunEphemeris
from Camera_Package.dofpCamera import DoFPCamera
from Polarization_Image_Package.image import Image


class DoFPProcessedPolarizationImage(Image):
    """
    A class used to represent a processed polarization image (angle of polarization or degree of polarization) from a
    division of focal plane polarimetric camera

    Attributes
    ----------
    image : float matrix
        the pixel image
    polarization_attribute : str
        the polarization attribute represented by the image ("aop_image_frame", "aop_meridian_frame", or "dolp")
    method : str
        method used to process the raw polarization image. If None, image should be the same size as the camera sensor
    camera : DoFPCamera
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
    color_chanel : str
        color channel used to process the polarization image, only for bayer (RGGB) sensor.
        We consider the camera with the following architecture (may vary depending on the camera)
        Camera color-polarization pattern
        - color order of the pixel
         R  | Gr
         ---|---
         Gb | B
        - polarizer orientation
         0°   |  45°
         -----|-----
         135° | 90°
    dofp_processed : bool
        true if the image is processd from a division of focal plane polarimetric camera

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

    def __init__(self, image, polarization_attribute, method="adjacent_super_pixel", camera=None, date=None,
                 hour_utc=None, minute_utc=None, second_utc=None, camera_exposure_time=None, camera_latitude=None,
                 camera_longitude=None, color_chanel=None):
        """
        Parameters
        ----------
        image : float matrix
            the pixel image
        polarization_attribute : str
            the polarization attribute represented by the image ("aop_image_frame", "aop_meridian_frame", or "dolp")
        method : str
            method used to process the raw polarization image
        camera : DoFPCamera
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
        color_chanel : str
            color channel used to process the polarization image, only for bayer (RGGB) sensor.
            We consider the camera with the following architecture (may vary depending on the camera)
            Camera color-polarization pattern
            - color order of the pixel
             R  | Gr
             ---|---
             Gb | B
            - polarizer orientation
             0°   |  45°
             -----|-----
             135° | 90°

        Raises
        ------
        TypeError
            if polarization_attribute or method is not a string, or camera is provided and is not a dofp camera
        ValueError
            if polarization_attribute is not one of the following: aop_image_frame,
            aop_meridian_frame, or dolp

        """
        if not isinstance(polarization_attribute, str):
            raise TypeError("polarization_attribute attribute should be an string object")
        if (polarization_attribute != "aop_image_frame") and (polarization_attribute != "aop_meridian_frame") \
                and (polarization_attribute != "dolp"):
            raise ValueError("polarization_attribute should be aop_image_frame, aop_meridian_frame or dolp")
        if not isinstance(method, str):
            raise TypeError("method attribute should be a string object")
        if isinstance(camera, DoFPCamera) or camera is None:
            super().__init__(image, camera, date, hour_utc, minute_utc, second_utc, camera_exposure_time,
                             camera_latitude, camera_longitude)
            self.dofp_processed = True
            self.polarization_attribute = polarization_attribute
            self.method = method
            self.color_chanel = color_chanel
        else:
            raise TypeError("The camera is not a DoFP camera")

