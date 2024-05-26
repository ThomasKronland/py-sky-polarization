class Lens:
    """
    A class used to represent a camera lens

    Attributes
    ----------
    focal_length : float
        the focal length of the lens in millimeters
    conjugation model : str
        the conjugation model of the lens ('r0','r1', 'r2','r3', or 'r4')

    Methods
    -------

    """

    def __init__(self, focal_length, conjugation_model):
        """
        Parameters
        ----------
        focal_length : float
            the focal length of the lens in millimeters
        conjugation model : str
            the conjugation model of the lens ('r0','r1', 'r2','r3', or 'r4')
            r0 : perspective imaging -> thin lens model
            r1 : stereographic imaging
            r2 : equidistant imaging -> twice the incident angle, twice the image length
            r3 : equi-solid angle imaging -> twice the incident solid angle, twice the image area
            r4 : vertical imaging

        Raises
        ----------
        TypeError
            if focal_length is not a float or int and conjugation_model is not a string
        """
        if not isinstance(focal_length, (int, float)):
            raise TypeError("focal_length attribute should be a float object")
        if not isinstance(conjugation_model, str):
            raise TypeError("conjugation_model attribute should be a string object")
        self.focal_length = focal_length
        self.conjugation_model = conjugation_model

    def __eq__(self, other):
        """
        Test if two lens are the same. Meaning same focal_length and conjugation model.

        Parameters
        ----------
        other : Lens
            the lens to compare with

        Returns
        -------
        bool
            True if the lens are the same, else False

        """
        if isinstance(other, Lens):
            if self.focal_length == other.focal_length and self.conjugation_model == other.conjugation_model:
                return True
        return False


if __name__ == "__main__":
    from Camera_Package.dofpSensor import DoFPSensor

    sensor_pixel_size1 = 3.45
    sensor_dimensions_row_column1 = (2048, 2448)
    sensor_type1 = "monochrome"

    lens_focal_length1 = 1.8
    lens_conjugation_model1 = 'r2'
    lens_focal_length2 = 1.2
    lens_conjugation_model2 = 'r0'
    lens1 = Lens(lens_focal_length1, lens_conjugation_model1)
    lens2 = Lens(lens_focal_length2, lens_conjugation_model1)
    lens3 = Lens(lens_focal_length1, lens_conjugation_model2)
    sensor = DoFPSensor(sensor_pixel_size1, sensor_dimensions_row_column1, sensor_type1)

    print(lens1 == lens1)
    print(lens1 == lens2)
    print(lens1 == lens3)
    print(lens1 == sensor)

