from Camera_Package.sensor import Sensor


class DoFPSensor(Sensor):
    """
    A class used to represent a division of focal plane polarimetric sensor

    Attributes
    ----------
    pixel_size : float
        the size of a pixel of the camera in micrometer
    sensor_dimensions : tuple of float
        the number of pixels in height and width such as (height, width)
    sensor_type : str
        type of the sensor ("monochrome", "bayer")
    dofp : bool
        true if the sensor is a division of focal plane sensor

    Methods
    -------

    """

    def __init__(self, pixel_size, sensor_dimensions, sensor_type):
        """
        Parameters
        ----------
        pixel_size : float
            the size of a pixel of the camera in micrometer
        sensor_dimensions : tuple of float
            the number of pixels in height and width such as (height, width)
        sensor_type : str
            type of the sensor ("monochrome", "bayer")
        """
        super().__init__(pixel_size, sensor_dimensions, sensor_type)
        self.dofp = True

    def __eq__(self, other):
        """
        Test if two dofp sensors are the same. Meaning same pixel size, dimensions and type.

        Parameters
        ----------
        other : DoFPSensor
            the sensor to compare with

        Returns
        -------
        bool
            True if the sensors are the same, else False

        """
        if isinstance(other, DoFPSensor):
            if self.pixel_size == other.pixel_size and self.sensor_dimensions == other.sensor_dimensions and \
                    self.sensor_type == other.sensor_type:
                return True
        return False


if __name__ == "__main__":
    sensor_pixel_size1 = 3.45
    sensor_dimensions_row_column1 = (2048, 2448)
    sensor_type1 = "monochrome"

    sensor_pixel_size2 = 3.2
    sensor_dimensions_row_column2 = (2000, 2448)
    sensor_type2 = "bayer"

    sensor1 = DoFPSensor(sensor_pixel_size1, sensor_dimensions_row_column1, sensor_type1)
    sensor2 = DoFPSensor(sensor_pixel_size2, sensor_dimensions_row_column1, sensor_type1)
    sensor3 = DoFPSensor(sensor_pixel_size1, sensor_dimensions_row_column2, sensor_type1)
    sensor4 = DoFPSensor(sensor_pixel_size1, sensor_dimensions_row_column1, sensor_type2)
    sensor5 = Sensor(sensor_pixel_size1, sensor_dimensions_row_column1, sensor_type1)

    print(sensor1 == sensor1)
    print(sensor1 == sensor2)
    print(sensor1 == sensor3)
    print(sensor1 == sensor4)
    print(sensor1 == sensor5)
