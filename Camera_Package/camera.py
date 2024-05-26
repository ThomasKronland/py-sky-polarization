from Camera_Package.lens import Lens
from Camera_Package.sensor import Sensor
import numpy as np


class Camera:
    """
    A class used to represent a camera

    Attributes
    ----------
    sensor : Sensor
        the camera sensor
    lens : Lens
        the camera lens
    center_coordinates : tuple of int
        the center coordinate of the lens in the camera image plane (height center, width center)

    Methods
    -------
    get_camera_azimuth_elevation_matrix()
        Computes the transformation matrices from pixel coordinates to azimuth and elevation in radians
    """

    def __init__(self, sensor, lens, center_coordinates):
        """
        Parameters
        ----------
        sensor : Sensor
            the camera sensor
        lens : Lens
            the camera lens
        center_coordinates : tuple of int
            the center coordinate of the lens in the camera image plane (height center, width center)

        Raises
        ------
        TypeError
            if sensor is not a Sensor, or lens is not a Lens, or center_coordinates is not a tuple of int
        ValueError
            if center_coordinates has not a length of 2
        """
        if not isinstance(sensor, Sensor):
            raise TypeError("sensor attribute should be a Sensor object")
        if not isinstance(lens, Lens):
            raise TypeError("lens attribute should be a Lens object")
        if not isinstance(center_coordinates, tuple):
            raise TypeError("center_coordinates attribute should be a tuple object")
        if len(center_coordinates) != 2:
            raise ValueError("center_coordinates attribute must have a length of 2")
        if not isinstance(center_coordinates[0], int) or not isinstance(center_coordinates[1], int):
            raise TypeError("center_coordinates attribute must be a tuple of int")
        self.sensor = sensor
        self.lens = lens
        self.center_coordinates = center_coordinates
        self.incident_rays_azimuth_matrix = None
        self.incident_rays_elevation_matrix = None
        self.x_coordinate_pixels = None
        self.y_coordinate_pixels = None

    def __eq__(self, other):
        """
        Test if two cameras are the same. Meaning same sensor, lens and center_coordinates attributes

        Parameters
        ----------
        other : Camera
            the camera to compare with

        Returns
        -------
        bool
            True if the cameras are the same, else False

        """
        if isinstance(other, Camera):
            if self.sensor == other.sensor and self.lens == other.lens and \
                    self.center_coordinates == other.center_coordinates:
                return True
        return False

    def get_camera_azimuth_elevation_matrix(self):
        """
        Computes the transformation matrices from pixel coordinates to azimuth and elevation

        Parameters
        ----------

        Returns
        -------
        incident_rays_azimuth_matrix : ndarray
            a matrix of dimension (self.sensor.sensor_dimensions[0], self.sensor.sensor_dimensions[1]) with each value
            corresponding to the azimuth each pixel sees in skydome in radians
        incident_rays_elevation_matrix : ndarray
            a matrix of dimension (self.sensor.sensor_dimensions[0], self.sensor.sensor_dimensions[1]) with each value
            corresponding to the elevation each pixel sees in skydome in radians
        x_coordinate_pixels : ndarray
        y_coordinate_pixels : ndarray

         Raises
        ------
        NotImplementedError
            if lens conjugation model is not "r0", "r1", "r2", "r3", or "r4"
        """
        if self.incident_rays_azimuth_matrix is None:
            # put focal length in micrometers
            focal_length = self.lens.focal_length * 1000
            # pixels coordinates in pixels
            sensor_height = self.sensor.sensor_dimensions[0]
            sensor_width = self.sensor.sensor_dimensions[1]

            # Attention : verifier que les coordonnees du centre doivent bien être corrigees de cette maniere
            x_coordinate_pixels = np.linspace((sensor_width - 1) / 2, - (sensor_width - 1) / 2,
                                              sensor_width) - (self.center_coordinates[1] - sensor_width // 2)
            y_coordinate_pixels = np.linspace((sensor_height - 1) / 2, - (sensor_height - 1) / 2,
                                              sensor_height) - (self.center_coordinates[0] - sensor_height // 2)

            # pixels coordinates in micrometers:
            sensor_pixel_size = self.sensor.pixel_size
            x_coordinate_micrometers = sensor_pixel_size * x_coordinate_pixels
            y_coordinate_micrometers = sensor_pixel_size * y_coordinate_pixels
            complex_sensor_plane = np.ones((sensor_height, 1)) * x_coordinate_micrometers + 1j * (
                np.transpose(y_coordinate_micrometers)[:, np.newaxis]) * np.ones((1, sensor_width))
            incident_rays_azimuth_matrix = np.angle(complex_sensor_plane)
            conjugation_model = self.lens.conjugation_model
            if conjugation_model == 'r0':
                incident_rays_elevation_matrix = (np.pi / 2) - np.arctan(np.abs(complex_sensor_plane) / focal_length)

            elif conjugation_model == 'r1':
                incident_rays_elevation_matrix = (np.pi / 2) - 2 * np.arctan(
                    np.abs(complex_sensor_plane) / focal_length / 2)

            elif conjugation_model == 'r2':
                incident_rays_elevation_matrix = (np.pi / 2) - (np.abs(complex_sensor_plane) / focal_length)

            elif conjugation_model == 'r3':
                incident_rays_elevation_matrix = (np.pi / 2) - 2 * np.arcsin(
                    np.abs(complex_sensor_plane) / focal_length / 2)

            elif conjugation_model == 'r4':
                incident_rays_elevation_matrix = (np.pi / 2) - np.arcsin(np.abs(complex_sensor_plane) / focal_length)

            else:
                raise NotImplementedError('Elevation matrix can not be computed yet for this lens conjugation')
            self.incident_rays_azimuth_matrix = incident_rays_azimuth_matrix
            self.incident_rays_elevation_matrix = incident_rays_elevation_matrix
            self.x_coordinate_pixels = x_coordinate_pixels
            self.y_coordinate_pixels = y_coordinate_pixels
        return self.incident_rays_azimuth_matrix, self.incident_rays_elevation_matrix, self.x_coordinate_pixels, \
            self.y_coordinate_pixels


if __name__ == "__main__":
    from Camera_Package.sensor import Sensor
    from Camera_Package.dofpSensor import DoFPSensor
    from Camera_Package.lens import Lens
    from Camera_Package.dofpCamera import DoFPCamera

    sensor_pixel_size1 = 3.45
    sensor_dimensions_row_column1 = (2048, 2448)
    sensor_type1 = "monochrome"

    sensor_pixel_size2 = 3.2
    sensor_dimensions_row_column2 = (2000, 2448)
    sensor_type2 = "bayer"

    sensor1 = Sensor(sensor_pixel_size1, sensor_dimensions_row_column1, sensor_type1)
    sensor2 = Sensor(sensor_pixel_size2, sensor_dimensions_row_column1, sensor_type1)
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
    cam1 = Camera(sensor1, lens1, camera_center_coordinates1)
    cam2 = Camera(sensor2, lens1, camera_center_coordinates1)
    cam3 = Camera(sensor3, lens1, camera_center_coordinates1)
    cam4 = Camera(sensor1, lens2, camera_center_coordinates1)
    cam5 = Camera(sensor1, lens1, camera_center_coordinates2)
    cam6 = DoFPCamera(sensor3, lens1, camera_center_coordinates1)

    print(cam1 == cam1)
    print(cam1 == cam2)
    print(cam1 == cam3)
    print(cam1 == cam4)
    print(cam1 == cam5)
    print(cam3 == cam6)
