from Polarization_Image_Package.image import Image


class ImageDataset:
    """
    A class used to create a dataset of image

    Attributes
    ----------
    images : list of Image
        a list containing images

    Methods
    -------
    add_image(image)
       adds an image to the dataset
    dataset_has_same_camera():
        Test if all images of the dataset are taken with the same camera
    dataset_has_same_position():
        Test if all images of the dataset are taken at the same position
    """

    def __init__(self):
        """
        Parameters
        ----------
        """
        self.images = []

    def add_image(self, image):
        """
        Adds an image to the dataset

        Parameters
        ----------
        image : DoFPPolarizationImage
            a division of focal plane polarization image

        Returns
        -------

        Raises
        ------
        AttributeError
            only Image objects can be added

        """
        if not isinstance(image, Image):
            raise TypeError("image should be an Image object")
        self.images.append(image)

    def __len__(self):
        """
        Computes the length of the dataset

        Parameters
        ----------

        Returns
        -------
        int
            the number of images in the dataset
        """
        return len(self.images)

    def __getitem__(self, index):
        """
        Returns an image in the dataset given an index

        Parameters
        ----------
        index
            the index of the image in the dataset images list

        Returns
        -------
        Image
           the image corresponding the input index
        """
        return self.images[index]

    def dataset_has_same_camera(self):
        """
        Test if all images of the dataset are taken with the same camera

        Parameters
        ----------

        Returns
        -------
        bool
          True if all images of the dataset were taken with the same camera, else False

        Raises
        ------
        ValueError
            the dataset should not be empty
        """
        if len(self) == 0:
            raise ValueError("Dataset is empty")
        first_image = self[0]
        for image in self.images:
            if not image.has_same_camera(first_image):
                return False
        return True

    def dataset_has_same_position(self):
        """
        Test if all images of the dataset are taken at the same position

        Parameters
        ----------

        Returns
        -------
        bool
            True if all images of the dataset were taken at the same position, else False

        Raises
        ------
        ValueError
            the dataset should not be empty
        """
        if len(self) == 0:
            raise ValueError("Dataset is empty")
        first_image = self[0]
        for image in self.images:
            if not image.has_same_position(first_image):
                return False
        return True


if __name__ == "__main__":
    from Camera_Package.sensor import Sensor
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

    image1 = Image(image=np.zeros((10,10)), camera=cam1, camera_latitude=cam_lat1, camera_longitude=cam_lon1)
    image2 = Image(image=np.ones((10,10)), camera=cam2, camera_latitude=cam_lat1, camera_longitude=cam_lon1)
    image3 = Image(image=np.ones((10,10)), camera=cam2, camera_latitude=cam_lat2, camera_longitude=cam_lon1)
    image4 = Image(image=np.ones((10,10)), camera=cam1, camera_latitude=cam_lat1, camera_longitude=cam_lon2)
    image5 = Image(image=np.ones((10,10)), camera=cam1, camera_latitude=cam_lat2, camera_longitude=cam_lon2)

    dataset1 = ImageDataset()
    for i in range(5):
        dataset1.add_image(image1)

    dataset2 = ImageDataset()
    dataset2.add_image(image2)
    dataset2.add_image(image3)

    dataset3 = ImageDataset()
    dataset3.add_image(image4)
    dataset3.add_image(image5)
    dataset3.add_image(image1)

    dataset4 = ImageDataset()
    dataset4.add_image(image1)
    dataset4.add_image(image2)

    print(dataset1.dataset_has_same_camera())
    print(dataset2.dataset_has_same_camera())
    print(dataset3.dataset_has_same_camera())
    print(dataset4.dataset_has_same_camera())
    print(dataset1.dataset_has_same_position())
    print(dataset2.dataset_has_same_position())
    print(dataset3.dataset_has_same_position())
    print(dataset4.dataset_has_same_position())



