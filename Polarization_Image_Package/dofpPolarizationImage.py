from Polarization_Image_Package.image import Image
from Camera_Package.dofpCamera import DoFPCamera
from Polarization_Image_Package.dofpProcessedPolarizationImage import DoFPProcessedPolarizationImage
import numpy as np


class DoFPPolarizationImage(Image):
    """
    A class used to represent DoFP polarization image

    Attributes
    ----------
    image : float matrix
        the pixel image
    camera : DoFPCamera
        the camera used to take the image. If not all the properties of the camera are known, a camera can be created
        with random parameters but by providing the camera type to compute the AoLP and DoLP
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
    is_polarized : bool
        true if the image is taken from a DoFP camera

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
    get_polarization_processed_image(method="adjacent_super_pixel", color_chanel=None)
        Computes, from the image attribute, the angle of polarization image in the image and meridian plane,
        and the degree of polarization image. If the camera is bayer type, the images are computed from pixel
        corresponding to the input color chanel.
   """

    def __init__(self, image, camera=None, date=None, hour_utc=None, minute_utc=None, second_utc=None,
                 camera_exposure_time=None, camera_latitude=None, camera_longitude=None):
        """
        Parameters
        ----------
        image : float matrix
            the pixel image
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

        Raises
        ------
        AttributeError
            an instance of the class DoFPPolarizationImage can only be created if the camera is from the class
            DoFPCamera

        """

        if isinstance(camera, DoFPCamera) or camera is None:
            super().__init__(image, camera, date, hour_utc, minute_utc, second_utc, camera_exposure_time,
                             camera_latitude, camera_longitude)
            self.is_polarized = True
        else:
            raise AttributeError("The camera is not a DoFP camera")

    def get_polarization_processed_image(self, method="adjacent_super_pixel", color_chanel=None):
        """
        Computes, from the image attribute, the angle of polarization image in the image and meridian plane,
        and the degree of polarization image

        Parameters
        ----------
        method : str
            method used to process the polarization image ("adjacent_super_pixel" or "superimposed_super_pixel")
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

        Returns
        -------
        aop_data_processing_image_frame : DoFPProcessedPolarizationImage
            Image with all attributes being the same as the input DoFPPolarizationImage, except the image attribute
            being the angle of polarization measured from the input image in the image frame.
            The size of the output image is :
            - (self.get_image_size[0]/2, self.get_image_size[1]/2) if camera is monochrome and the method used for
            processing is the adjacent super pixel method
            - (self.get_image_size[0]/4, self.get_image_size[1]/4) if camera is Bayer RGB type and the method used for
            processing is the adjacent super pixel method
        aop_data_processing_meridian_frame : DoFPProcessedPolarizationImage
            Image with all attributes being the same as the input DoFPPolarizationImage, except the image attribute
            being the angle of polarization measured from the input image in the meridian frame.
            The size of the output image is :
            - (self.get_image_size[0]/2, self.get_image_size[1]/2) if camera is monochrome and the method used for
            processing is the adjacent super pixel method
            - (self.get_image_size[0]/4, self.get_image_size[1]/4) if camera is Bayer RGB type and the method used for
            processing is the adjacent super pixel method
        dolp_data_processing : DoFPProcessedPolarizationImage
            Image with all attributes being the same as the input DoFPPolarizationImage, except the image attribute
            being the degree of polarization measured from the input image in the image frame.
            The size of the output image is :
            - (self.get_image_size[0]/2, self.get_image_size[1]/2) if camera is monochrome and the method used for
            processing is the adjacent super pixel method
            - (self.get_image_size[0]/4, self.get_image_size[1]/4) if camera is Bayer RGB type and the method used for
            processing is the adjacent super pixel method

        Raises
        ------
        AttributeError
            the method can only be called if the attribute camera is provided. The type of the camera should be known
            (monochrome or bayer). If the precise properties of the camera are unknown, create a new camera with random
            parameters and its real type
        TypeError
            if method is not a string, or if camera is bayer and color chanel is not a string
        ValueError
            if the camera is bayer type, the color chanel should be provided
        NotImplementedError
            this function is only implemented for default adjacent_super_pixel method on bayer or monochrome camera
        """
        if self.camera is None:
            raise AttributeError("The camera attribute should be provided for its type (monochrome or bayer)")
        if not isinstance(method, str):
            raise TypeError("method should be a string")
        if method == "adjacent_super_pixel":
            if self.camera.sensor.sensor_type == "monochrome":
                rows, cols = self.get_image_size()
                rows_final_image = int(rows / 2)
                cols_final_image = int(cols / 2)
                # Make output matrix, half the high and half the width of input matrix :
                aop_data_processing_image_frame = np.zeros((rows_final_image, cols_final_image))
                aop_data_processing_meridian_frame = np.zeros((rows_final_image, cols_final_image))
                dolp_data_processing = np.zeros((rows_final_image, cols_final_image))
                for k in (range(rows_final_image)):
                    for j in (range(cols_final_image)):
                        # Compute Stokes parameters :
                        s0 = 0.5 * (self.image[2 * k, 2 * j] + self.image[2 * k, 2 * j + 1] + self.image[
                            2 * k + 1, 2 * j] + self.image[2 * k + 1, 2 * j + 1])
                        s1 = self.image[2 * k + 1, 2 * j + 1] - self.image[2 * k, 2 * j]
                        s2 = self.image[2 * k, 2 * j + 1] - self.image[2 * k + 1, 2 * j]
                        dolp_data_processing[k, j] = np.sqrt((s1 ** 2) + (s2 ** 2)) / s0
                        aop_data_processing_image_frame[k, j] = 0.5 * np.arctan2(s2, s1)
                        # Compute pixels' coordinates:
                        x_super_pixel = + j - (cols_final_image / 2) - 0.5
                        y_super_pixel = - k + (rows_final_image / 2) + 0.5
                        aop_data_processing_meridian_frame[k, j] = 0.5 * np.angle(
                            (s1 + s2 * 1j) / (x_super_pixel + y_super_pixel * 1j) / (
                                    x_super_pixel + y_super_pixel * 1j))
            elif self.camera.sensor.sensor_type == "bayer":
                if color_chanel is None:
                    raise ValueError("the color chanel on which to compute the algorithm should be provided")
                if not isinstance(color_chanel, str):
                    raise TypeError("color_chanel should be a string")
                rows, cols = self.get_image_size()
                rows_final_image = int(rows / 4)
                cols_final_image = int(cols / 4)
                # Make output matrix, half the high and half the width of input matrix :
                aop_data_processing_image_frame = np.zeros((rows_final_image, cols_final_image))
                aop_data_processing_meridian_frame = np.zeros((rows_final_image, cols_final_image))
                dolp_data_processing = np.zeros((rows_final_image, cols_final_image))
                if color_chanel == "R":
                    color = 0
                elif color_chanel == "Gr":
                    color = 1
                elif color_chanel == "Gb":
                    color = 2
                elif color_chanel == "B":
                    color = 3
                else:
                    raise ValueError("color_chanel should be one of the following R, Gr, Gb or B")
                # Processing AOP and DOP for the selected color
                # AOP values calculation
                for k in (range(rows_final_image)):
                    for j in (range(cols_final_image)):
                        # Compute Stokes parameters :
                        s0 = 0.5 * (self.image[4 * k + 2 * (color // 2), 4 * j + 2 * (color % 2)]
                                    + self.image[4 * k + 2 * (color // 2), 4 * j + 1 + 2 * (color % 2)]
                                    + self.image[4 * k + 1 + 2 * (color // 2), 4 * j + 2 * (color % 2)]
                                    + self.image[4 * k + 1 + 2 * (color // 2), 4 * j + 1 + 2 * (color % 2)])
                        s1 = self.image[4 * k + 1 + 2 * (color // 2), 4 * j + 1 + 2 * (color % 2)] - self.image[
                            4 * k + 2 * (color // 2), 4 * j + 2 * (color % 2)]
                        s2 = self.image[4 * k + 2 * (color // 2), 4 * j + 1 + 2 * (color % 2)] - self.image[
                            4 * k + 1 + 2 * (color // 2), 4 * j + 2 * (color % 2)]
                        dolp_data_processing[k, j] = np.sqrt((s1 ** 2) + (s2 ** 2)) / s0
                        aop_data_processing_image_frame[k, j] = 0.5 * np.arctan2(s2, s1)
                        # Compute pixels' coordinates:
                        x_super_pixel = + j - (cols_final_image / 2) - 0.5
                        y_super_pixel = - k + (rows_final_image / 2) + 0.5
                        aop_data_processing_meridian_frame[k, j] = 0.5 * np.angle(
                            (s1 + s2 * 1j) / (x_super_pixel + y_super_pixel * 1j) / (
                                    x_super_pixel + y_super_pixel * 1j))
            else:
                raise NotImplementedError("this method is not implemented for this type of sensor")
        if method == "superimposed_super_pixel":
            if self.camera.sensor.sensor_type == "monochrome":
                raise NotImplementedError("superimposed_super_pixel method is not implemented for monochrome sensor")
            elif self.camera.sensor.sensor_type == "bayer":
                raise NotImplementedError("superimposed_super_pixel method is not implemented for bayer sensor")
            else:
                raise NotImplementedError("this method is not implemented for this type of sensor")

        aop_data_processing_image_frame_class_image = DoFPProcessedPolarizationImage(
            image=aop_data_processing_image_frame,
            polarization_attribute="aop_image_frame",
            method=method,
            camera=self.camera,
            date=self.date,
            hour_utc=self.hour_utc,
            minute_utc=self.minute_utc,
            second_utc=self.second_utc,
            camera_exposure_time=self.camera_exposure_time,
            camera_latitude=self.camera_latitude,
            camera_longitude=self.camera_longitude,
            color_chanel=color_chanel
        )
        aop_data_processing_meridian_frame_class_image = DoFPProcessedPolarizationImage(
            image=aop_data_processing_meridian_frame,
            polarization_attribute="aop_meridian_frame",
            method=method,
            camera=self.camera,
            date=self.date,
            hour_utc=self.hour_utc,
            minute_utc=self.minute_utc,
            second_utc=self.second_utc,
            camera_exposure_time=self.camera_exposure_time,
            camera_latitude=self.camera_latitude,
            camera_longitude=self.camera_longitude,
            color_chanel=color_chanel
        )
        dolp_data_processing_class_image = DoFPProcessedPolarizationImage(
            image=dolp_data_processing,
            polarization_attribute="dolp",
            method=method,
            camera=self.camera, date=self.date,
            hour_utc=self.hour_utc,
            minute_utc=self.minute_utc,
            second_utc=self.second_utc,
            camera_exposure_time=self.camera_exposure_time,
            camera_latitude=self.camera_latitude,
            camera_longitude=self.camera_longitude,
            color_chanel=color_chanel
        )

        return aop_data_processing_image_frame_class_image, aop_data_processing_meridian_frame_class_image, \
               dolp_data_processing_class_image


if __name__ == "__main__":
    from Camera_Package.sensor import Sensor
    from Camera_Package.dofpSensor import DoFPSensor
    from Camera_Package.lens import Lens
    from Camera_Package.dofpCamera import DoFPCamera
    import numpy as np
    import matplotlib.pyplot as plt
    import Polarization_Simulation_Package.polarizationSimulation as pola_simu
    from Camera_Image_Display_Package.imageDisplay import ImageDisplay

    sensor_pixel_size1 = 3.45
    sensor_dimensions_row_column1 = (2048, 2448)
    sensor_type1 = "bayer"

    sensor_dimensions_row_column2 = (2048 // 2, 2448 // 2)
    sensor_type2 = "monochrome"

    sensor1 = DoFPSensor(sensor_pixel_size1, sensor_dimensions_row_column1, sensor_type1)
    sensor2 = DoFPSensor(sensor_pixel_size1, sensor_dimensions_row_column2, sensor_type2)

    lens_focal_length1 = 1.8
    lens_focal_length2 = 1.8 / 2
    lens_conjugation_model1 = 'r2'

    lens1 = Lens(lens_focal_length1, lens_conjugation_model1)
    lens2 = Lens(lens_focal_length2, lens_conjugation_model1)

    camera_center_coordinates1 = (1024, 1224)
    camera_center_coordinates2 = (1024 // 2, 1224 // 2)

    # Create a camera
    cam1 = DoFPCamera(sensor1, lens1, camera_center_coordinates1)
    cam2 = DoFPCamera(sensor2, lens2, camera_center_coordinates2)
    bits_image = pola_simu.simu_micro_polarizers(cam1, "Rayleigh", 50, 20)
    bits_image2 = pola_simu.simu_micro_polarizers(cam2, "Rayleigh", 50, 20)
    aop_data_processing_image_frame_class_image, aop_data_processing_meridian_frame_class_image, \
        dolp_data_processing_class_image = bits_image.get_polarization_processed_image(color_chanel="B")

    aop_data_processing_image_frame_class_image2, aop_data_processing_meridian_frame_class_image2, \
        dolp_data_processing_class_image2 = bits_image2.get_polarization_processed_image()

    aolp_image_display = ImageDisplay(aop_data_processing_image_frame_class_image, "aolp image bayer")
    aolp_meridian_display = ImageDisplay(aop_data_processing_meridian_frame_class_image, "aolp meridian bayer")
    dolp_im_display = ImageDisplay(dolp_data_processing_class_image, "dolp bayer")

    aolp_image_display2 = ImageDisplay(aop_data_processing_image_frame_class_image2, "aolp image monochrome")
    aolp_meridian_display2 = ImageDisplay(aop_data_processing_meridian_frame_class_image2, "aolp meridian monochrome")
    dolp_im_display2 = ImageDisplay(dolp_data_processing_class_image2, "dolp bayer monochrome")

    aolp_image_display.display_image()
    aolp_meridian_display.display_image()
    dolp_im_display.display_image()

    aolp_image_display2.display_image()
    aolp_meridian_display2.display_image()
    dolp_im_display2.display_image()

    plt.show()
