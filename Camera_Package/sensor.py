class Sensor:
    """
    A class used to represent a sensor

    Attributes
    ----------
    pixel_size : float
        the size of a pixel of the camera in micrometer
    sensor_dimensions : tuple of float
        the number of pixels in height and width such as (height, width)
    sensor_type : str
        type of the sensor ("monochrome", "bayer")

    Methods
    -------

    """

    def __init__(self, pixel_size, sensor_dimensions, sensor_type):
        """
        Parameters
        ----------
        pixel_size : float
            the size of a pixel of the camera in micrometer
        sensor_dimensions : tuple of int
            the number of pixels in height and width such as (height, width)
        sensor_type : str
            type of the sensor ("monochrome", "bayer")
        """
        if not isinstance(pixel_size, (int, float)):
            raise TypeError("pixel_size attribute should be a float object")
        if not isinstance(sensor_dimensions, tuple):
            raise TypeError("sensor_dimensions attribute should be a tuple object")
        if len(sensor_dimensions) != 2:
            raise ValueError("sensor_dimensions attribute must have a length of 2")
        if not isinstance(sensor_dimensions[0], int) or not isinstance(sensor_dimensions[1], int):
            raise ValueError("sensor_dimensions attribute must be a tuple of int")
        if not isinstance(sensor_type, str):
            raise TypeError("sensor_type attribute should be a string object")
        self.pixel_size = pixel_size
        self.sensor_dimensions = sensor_dimensions
        self.sensor_type = sensor_type

    def __eq__(self, other):
        """
        Test if two sensors are the same. Meaning same pixel size, dimensions and type.

        Parameters
        ----------
        other : Sensor
            the sensor to compare with

        Returns
        -------
        bool
            True if the sensors are the same, else False

        """
        if isinstance(other, Sensor):
            if self.pixel_size == other.pixel_size and self.sensor_dimensions == other.sensor_dimensions and \
                    self.sensor_type == other.sensor_type:
                return True
        return False


if __name__ == "__main__":
    from Camera_Package.dofpSensor import DoFPSensor

    sensor_pixel_size1 = 3.45
    sensor_dimensions_row_column1 = (2048, 2448)
    sensor_type1 = "monochrome"

    sensor_pixel_size2 = 3.2
    sensor_dimensions_row_column2 = (2000, 2448)
    sensor_type2 = "bayer"

    sensor1 = Sensor(sensor_pixel_size1, sensor_dimensions_row_column1, sensor_type1)
    sensor2 = Sensor(sensor_pixel_size2, sensor_dimensions_row_column1, sensor_type1)
    sensor3 = Sensor(sensor_pixel_size1, sensor_dimensions_row_column2, sensor_type1)
    sensor4 = Sensor(sensor_pixel_size1, sensor_dimensions_row_column1, sensor_type2)
    sensor5 = DoFPSensor(sensor_pixel_size1, sensor_dimensions_row_column1, sensor_type1)

    print(sensor1 == sensor1)
    print(sensor1 == sensor2)
    print(sensor1 == sensor3)
    print(sensor1 == sensor4)
    print(sensor1 == sensor5)


