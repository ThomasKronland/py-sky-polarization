from Polarization_Image_Package.dofpPolarizationImage import DoFPPolarizationImage
from Polarization_Image_Package.imageDataset import ImageDataset


class DoFPPolarizationImageDataset(ImageDataset):
    """
    A class used to create a dataset of DoFP polarization image

    Attributes
    ----------
    images : list of DoFPPolarizationImage
        a list containing division of focal plane polarization images

    Methods
    -------
    add_image(image)
       adds an image to the dataset
    dataset_has_same_camera():
        Test if all images of the dataset are taken with the same camera
    dataset_has_same_position():
        Test if all images of the dataset are taken at the same position
    """

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
        TypeError
            only division of focal plane polarization images can be added

        """
        if isinstance(image, DoFPPolarizationImage):
            self.images.append(image)
        else:
            raise TypeError("image input is not a DoFPPolarizationImage object")