import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
from Polarization_Image_Package.image import Image


class ImageDisplay:
    """
    A class used to display an image

    Attributes
    ----------
    image : Image
        the image to display
    figure_title : str
        title of the figure
    figure_number : int
        number of the figure
    cmap : str
        colormap used to display the image

    Methods
    -------
    display_image(center_image = False)
        method used to display the image
    """

    def __init__(self, image, figure_title, figure_number=None, cmap='viridis'):
        """
        Parameters
        ----------
        image : Image
            the image to display
        figure_title : str
            title of the figure
        figure_number : int
            number of the figure
        cmap : str
            colormap used to display the image

        Raises
        ------
        TypeError
            if image is not an Image, or figure_number is provided and not an integer, or figure_title is not a string,
            or cmap is not a string
        ValueError
            if cmap is not a valid color map
        """
        if not isinstance(image, Image):
            raise TypeError("image attribute should be an Image object")
        if not isinstance(figure_number, int) and figure_number is not None:
            raise TypeError("figure_number attribute should be an int object")
        if not isinstance(figure_title, str):
            raise TypeError("figure_title attribute should be a string object")
        if not isinstance(cmap, str):
            raise TypeError("cmap attribute should be a string object")
        if cmap not in list(colormaps):
            raise ValueError("unknown cmap value, to check the list of possible values, check "
                             "list(matplotlib.colormaps)")
        self.image = image
        self.figure_number = figure_number
        self.cmap = cmap
        self.figure_title = figure_title

    def display_image(self, center_image=False):
        """
        displays an image from class Image

        Parameters
        ----------
        center_image : bool
            a bool to consider, if True, the center of the image, corresponding to the optical center of the camera, as
            the center of the displayed image

        Returns
        -------

        Raises
        ------
        AttributeError
            if center_image is True, the camera used to take the image should provided

        """
        image_height = self.image.get_image_size()[0]
        image_width = self.image.get_image_size()[1]

        center_coordinates = (image_height // 2, image_width // 2)
        if center_image:
            if self.image.camera is not None:
                center_coordinates = self.image.camera.center_coordinates
            else:
                raise AttributeError("The camera used to take the image should be provided")

        x_coordinate_pixels = np.linspace((image_width - 1) / 2, - (image_width - 1) / 2,
                                          image_width) - (center_coordinates[1] - image_width // 2)
        y_coordinate_pixels = np.linspace((image_height - 1) / 2, - (image_height - 1) / 2,
                                          image_height) - (center_coordinates[0] - image_height // 2)
        x_pixel_mesh = (np.ones(y_coordinate_pixels.size)[:, np.newaxis]) * x_coordinate_pixels[np.newaxis, :]
        y_pixel_mesh = np.ones(x_coordinate_pixels.size)[np.newaxis, :] * y_coordinate_pixels[:, np.newaxis]

        plt.figure(self.figure_number)
        h = plt.pcolormesh(x_pixel_mesh, y_pixel_mesh, self.image.image, cmap=self.cmap)
        plt.colorbar()
        plt.title(self.figure_title)


if __name__ == "__main__":
    import numpy as np
    from Polarization_Image_Package.image import Image
    from Camera_Package.camera import Camera
    from Camera_Package.sensor import Sensor
    from Camera_Package.lens import Lens

    sensor_pixel_size = 3.45
    sensor_dimensions_row_column = (2048, 2448)
    sensor_type = "monochrome"
    lens_focal_length = 1.8
    lens_conjugation_model = 'r2'
    camera_center_coordinates = (1024, 1224)
    # Create a camera
    sensor = Sensor(sensor_pixel_size, sensor_dimensions_row_column, sensor_type)
    lens = Lens(lens_focal_length, lens_conjugation_model)
    cam = Camera(sensor, lens, camera_center_coordinates)
    # Create an image
    im = np.random.rand(30, 20)
    image1 = Image(im)
    # display image
    dispIm = ImageDisplay(image1, "test")
    dispIm.display_image()
    plt.show()
