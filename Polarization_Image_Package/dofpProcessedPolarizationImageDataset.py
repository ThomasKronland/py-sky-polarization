from Polarization_Image_Package.dofpProcessedPolarizationImage import DoFPProcessedPolarizationImage
from Polarization_Image_Package.imageDataset import ImageDataset


class DoFPProcessedPolarizationImageDataset(ImageDataset):
    """
    A class used to create a dataset of DoFP processed polarization image

    Attributes
    ----------
    images : list of DoFPProcessedPolarizationImage
        a list containing division of focal plane processed polarization images

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
        image : DoFPProcessedPolarizationImage
            a division of focal plane polarization processsed image

        Returns
        -------

        Raises
        ------
        TypeError
            only division of focal plane polarization processed images can be added

        """
        if isinstance(image, DoFPProcessedPolarizationImage):
            self.images.append(image)
        else:
            raise TypeError("image input is not a DoFPPolarizationImage object")