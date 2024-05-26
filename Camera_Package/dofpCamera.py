import numpy as np
from Camera_Package.lens import Lens
from Camera_Package.dofpSensor import DoFPSensor
from Camera_Package.camera import Camera
import skimage.measure
from scipy.stats import circmean


class DoFPCamera(Camera):
    """
    A class used to represent a division of focal plane polarimetric camera

    Attributes
    ----------
    dofp_sensor : DoFPSensor
        the camera division of focal plane polarimetric sensor
    lens : Lens
        the camera lens
    center_coordinates : tuple of int
        the center coordinate of the lens in the camera image plane
    dofp : bool
        true if the camera is a dofp camera

    Methods
    -------
    get_camera_azimuth_elevation_matrix()
        Computes the transformation matrices from pixel coordinates to azimuth and elevation
    camera_azimuth_elevation_matrix_pooling(method='adjacent_super_pixel', color_chanel=None)
        Reduces the size of the transformation matrices from pixel coordinates to azimuth and elevation to make
        it correspond to AoLP and DoLP images computed from separated super_pixel
    """

    def __init__(self, dofp_sensor, lens, center_coordinates):
        """
        Parameters
        ----------
        dofp_sensor : DoFPSensor
            the division of focal plane polarimetric sensor
        lens : Lens
            the camera lens
        center_coordinates : tuple of int
            the center coordinate of the lens in the camera image plane

        Raises
        ------
        TypeError
            an instance of the class DoFPCamera can only be created if the sensor is from the class DoFPSensor

        """
        if isinstance(dofp_sensor, DoFPSensor):
            super().__init__(dofp_sensor, lens, center_coordinates)
            self.dofp = True
        else:
            raise TypeError("The sensor is not a DoFP sensor")

    def __eq__(self, other):
        """
        Test if two dofp cameras are the same. Meaning same sensor, lens and center_coordinates attributes

        Parameters
        ----------
        other : DoFPCamera
            the camera to compare with

        Returns
        -------
        bool
            True if the cameras are the same, else False

        """
        if isinstance(other, DoFPCamera):
            if self.sensor == other.sensor and self.lens == other.lens and \
                    self.center_coordinates == other.center_coordinates:
                return True
        return False

    def get_camera_azimuth_elevation_matrix_pooled(self, method='adjacent_super_pixel', color_chanel=None):
        """
        Reduces the size of the transformation matrices from pixel coordinates to azimuth and elevation to make
        it correspond to AoLP and DoLP images computed from super_pixel

        Parameters
        ----------
        method : string
            the method used to compute the AoLP and/or DoLP from a raw image
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
        reduced_incident_rays_azimuth_matrix
           a reduced matrix of transformation matrices from pixel coordinates to azimuth
        reduced_incident_rays_elevation_matrix
            a reduced matrix of transformation matrices from pixel coordinates to azimuth
        Raises
        ------
        NotImplementedError
           the method is only implemented for the adjacent super pixel method, and for monochrome or bayer sensors
        """
        if method == 'adjacent_super_pixel':
            if self.sensor.sensor_type == "monochrome":
                incident_rays_azimuth_matrix, incident_rays_elevation_matrix, _, _ = \
                    self.get_camera_azimuth_elevation_matrix()
                reduced_incident_rays_azimuth_matrix = skimage.measure.block_reduce(incident_rays_azimuth_matrix,
                                                                                    (2, 2), circmean,
                                                                                    func_kwargs={'high': np.pi,
                                                                                                 'low': -np.pi})
                reduced_incident_rays_elevation_matrix = skimage.measure.block_reduce(incident_rays_elevation_matrix,
                                                                                      (2, 2), circmean,
                                                                                      func_kwargs={'high': np.pi,
                                                                                                   'low': -np.pi})
            elif self.sensor.sensor_type == "bayer" and color_chanel is not None:
                incident_rays_azimuth_matrix, incident_rays_elevation_matrix, _, _ = \
                    self.get_camera_azimuth_elevation_matrix()
                reduced_incident_rays_azimuth_matrix_first_step = skimage.measure.block_reduce(
                                                                                        incident_rays_azimuth_matrix,
                                                                                        (2, 2), circmean,
                                                                                        func_kwargs={'high': np.pi,
                                                                                                     'low': -np.pi})
                reduced_incident_rays_elevation_matrix_first_step = skimage.measure.block_reduce(
                                                                                        incident_rays_elevation_matrix,
                                                                                        (2, 2), circmean,
                                                                                        func_kwargs={'high': np.pi,
                                                                                                     'low': -np.pi})
                if color_chanel == "R":
                    color = 0
                if color_chanel == "Gr":
                    color = 1
                if color_chanel == "Gb":
                    color = 2
                if color_chanel == "B":
                    color = 3
                reduced_incident_rays_azimuth_matrix = reduced_incident_rays_azimuth_matrix_first_step[
                                                       (color // 2)::2, (color % 2)::2]
                reduced_incident_rays_elevation_matrix = reduced_incident_rays_elevation_matrix_first_step[
                                                         (color // 2)::2, (color % 2)::2]
            elif self.sensor.sensor_type == "bayer":
                # ok if we consider all the color chanel, else we may have to create different matrices for each color
                # chanel, for that we can do a first pooling with (2,2) kernel, and then a subsampling deepending on
                # the color chanel we choose
                incident_rays_azimuth_matrix, incident_rays_elevation_matrix, _, _ = \
                    self.get_camera_azimuth_elevation_matrix()
                reduced_incident_rays_azimuth_matrix = skimage.measure.block_reduce(incident_rays_azimuth_matrix,
                                                                                    (4, 4),
                                                                                    circmean,
                                                                                    func_kwargs={'high': np.pi,
                                                                                                 'low': -np.pi})
                reduced_incident_rays_elevation_matrix = skimage.measure.block_reduce(incident_rays_elevation_matrix,
                                                                                      (4, 4), circmean,
                                                                                      func_kwargs={'high': np.pi,
                                                                                                   'low': -np.pi})
            else:
                raise NotImplementedError("not implemented for other camera sensor types")
        else:
            raise NotImplementedError("not implemented for other methods")
        return reduced_incident_rays_azimuth_matrix, reduced_incident_rays_elevation_matrix


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from Polarization_Image_Package.image import Image
    from Camera_Package.dofpSensor import DoFPSensor
    from Camera_Package.lens import Lens
    from Camera_Image_Display_Package.imageDisplay import ImageDisplay
    from Camera_Package.camera import Camera

    sensor_pixel_size = 3.45
    sensor_dimensions_row_column = (2048, 2448)
    sensor_type = "bayer"
    lens_focal_length = 1.8
    lens_conjugation_model = 'r2'
    camera_center_coordinates = (1024, 1224)
    # Create a camera
    sensor = DoFPSensor(sensor_pixel_size, sensor_dimensions_row_column, sensor_type)
    lens = Lens(lens_focal_length, lens_conjugation_model)
    cam = DoFPCamera(sensor, lens, camera_center_coordinates)

    incident_rays_azimuth_matrix, incident_rays_elevation_matrix, _, _ = \
        cam.get_camera_azimuth_elevation_matrix()
    reduced_incident_rays_azimuth_matrix, reduced_incident_rays_elevation_matrix = \
        cam.get_camera_azimuth_elevation_matrix_pooled(color_chanel="Gr")
    incident_rays_azimuth_matrix_image = Image(incident_rays_azimuth_matrix)
    incident_rays_elevation_matrix_image = Image(incident_rays_elevation_matrix)
    reduced_incident_rays_azimuth_matrix_image = Image(reduced_incident_rays_azimuth_matrix)
    reduced_incident_rays_elevation_matrix_image = Image(reduced_incident_rays_elevation_matrix)
    azimuth_matrix_image_display = ImageDisplay(incident_rays_azimuth_matrix_image, "azimuth matrix")
    elevation_matrix_image_display = ImageDisplay(incident_rays_elevation_matrix_image, "elevation matrix")
    reduced_azimuth_matrix_image_display = ImageDisplay(reduced_incident_rays_azimuth_matrix_image,
                                                        "azimuth matrix reduced")
    reduced_elevation_matrix_image_display = ImageDisplay(reduced_incident_rays_elevation_matrix_image,
                                                          "elevation matrix reduced")
    azimuth_matrix_image_display.display_image()
    elevation_matrix_image_display.display_image()
    reduced_azimuth_matrix_image_display.display_image()
    reduced_elevation_matrix_image_display.display_image()
    plt.show()

    # test equality
    sensor_pixel_size1 = 3.45
    sensor_dimensions_row_column1 = (2048, 2448)
    sensor_type1 = "monochrome"

    sensor_pixel_size2 = 3.2
    sensor_dimensions_row_column2 = (2000, 2448)
    sensor_type2 = "bayer"

    sensor1 = DoFPSensor(sensor_pixel_size1, sensor_dimensions_row_column1, sensor_type1)
    sensor2 = DoFPSensor(sensor_pixel_size2, sensor_dimensions_row_column1, sensor_type1)
    sensor3 = DoFPSensor(sensor_pixel_size1, sensor_dimensions_row_column1, sensor_type1)

    lens_focal_length1 = 1.8
    lens_conjugation_model1 = 'r2'

    lens_focal_length2 = 1.2
    lens_conjugation_model2 = 'r0'

    lens1 = Lens(lens_focal_length1, lens_conjugation_model1)
    lens2 = Lens(lens_focal_length2, lens_conjugation_model1)
    lens3 = Lens(lens_focal_length1, lens_conjugation_model2)

    camera_center_coordinates1 = (1024, 1224)
    camera_center_coordinates2 = (1000, 1224)

    # Create a camera
    cam1 = DoFPCamera(sensor1, lens1, camera_center_coordinates1)
    cam2 = DoFPCamera(sensor2, lens1, camera_center_coordinates1)
    cam3 = DoFPCamera(sensor3, lens1, camera_center_coordinates1)
    cam4 = DoFPCamera(sensor1, lens2, camera_center_coordinates1)
    cam5 = DoFPCamera(sensor1, lens1, camera_center_coordinates2)
    cam6 = Camera(sensor3, lens1, camera_center_coordinates1)

    print(cam1 == cam1)
    print(cam1 == cam2)
    print(cam1 == cam3)
    print(cam1 == cam4)
    print(cam1 == cam5)
    print(cam3 == cam6)
