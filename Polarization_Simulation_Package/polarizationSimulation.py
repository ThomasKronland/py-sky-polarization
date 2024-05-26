import numpy as np
from Camera_Package.dofpCamera import DoFPCamera
from Camera_Package.camera import Camera
from Polarization_Image_Package import DoFPPolarizationImage

"""
Methods
-------
- simulate_rayleigh_sky(camera, sun_elevation, sun_azimuth, dolp_max)
    Simulates angle of linear polarization and degree of polarization images based on Rayleigh single scattering model, 
    for each pixel of a given camera, for a given sun position and maximum degree of linear polarization
- simulate_berry_sky(camera, sun_elevation, sun_azimuth, dolp_max, neutral_point_angular_distance_delta)
    Simulates angle of linear polarization and degree of polarization images based on Berry model, for each pixel of a 
    given camera, for a given sun position, maximum degree of linear polarization and neutral points angular distance
- simu_sky_intensity_cie(camera, sun_elevation, sun_azimuth, cie_sky_number)
    Simulates an image of relative sky intensity, for each pixel of a given camera, for a given sun position and a given 
    CIE model.
- simu_micro_polarizers(camera, simu_model, sun_elevation, sun_azimuth, dolp_max=1, tolerance_deg=0,
                          extinction_ratio=1, cie_sky_number=5, pixel_saturation_ratio=1.1, number_of_bits=16,
                          snr_gauss=10 ** 9, neutral_point_angular_distance_delta=0)
    Simulates a pixel polarization image using an input sky polarization model, as seen by a polarimetric dofp camera,
    for a given sun position relative to the camera, maximum degree of linear polarization (and neutral points angular
    distance, if the simulation is based on Berry model)
"""


def simulate_rayleigh_sky(camera, sun_elevation, sun_azimuth, dolp_max):
    """
    Simulates angle of linear polarization and degree of polarization images based on Rayleigh single scattering model,
    for each pixel of a given camera, for a given sun position and maximum degree of linear polarization

    Parameters
    ----------
    camera : Camera
        camera used to simulate the polarization pattern
    sun_elevation : float
        sun elevation in degrees
    sun_azimuth : float
        sun azimuth in degrees
    dolp_max : float
        maximum degree of polarization in the sky

    Returns
    -------
    aop_matrix_global : ndarray
       image of angle of polarization in the image frame as seen by the camera
    dolp_matrix : ndarray
        image of degree of linear polarization as seen by the camera

    Raises
    ------
    TypeError
        camera should be an instance of the Camera class, and sun_elevation, sun_azimuth and dolp_max should be float
        objects

    """
    if not isinstance(camera, Camera):
        raise TypeError("camera should be an instance of the Camera class")
    if not isinstance(sun_elevation, (int, float)) or not isinstance(sun_azimuth, (int, float)) or not isinstance(
            dolp_max, (int, float)):
        raise TypeError("sun_elevation, sun_azimuth and dolp_max should be float objects")

    sun_elevation_rad = sun_elevation * np.pi / 180
    sun_azimuth_rad = sun_azimuth * np.pi / 180

    # Sun vector in XYZ coordinates
    x_sun = np.cos(sun_elevation_rad) * np.cos(sun_azimuth_rad)
    y_sun = np.cos(sun_elevation_rad) * np.sin(sun_azimuth_rad)
    z_sun = np.sin(sun_elevation_rad)

    sky_particle_azimuth_matrix, sky_particle_elevation_matrix, x_coordinate_pixels, y_coordinate_pixels \
        = camera.get_camera_azimuth_elevation_matrix()

    # Particle vectors in XYZ coordinates
    x_particle = np.multiply(np.cos(sky_particle_elevation_matrix), np.cos(sky_particle_azimuth_matrix))
    y_particle = np.multiply(np.cos(sky_particle_elevation_matrix), np.sin(sky_particle_azimuth_matrix))
    z_particle = np.sin(sky_particle_elevation_matrix)

    # AoP_g Matrix :
    # the local AoP in 3D frame is the Angle between the E vector (OP x OS) and
    # the plane which contains OP and Z.

    tan_aop_l = (np.multiply(np.sin(sky_particle_elevation_matrix) * np.cos(sun_elevation_rad),
                             np.cos(sun_azimuth_rad - sky_particle_azimuth_matrix)) - np.sin(
        sun_elevation_rad) * np.cos(
        sky_particle_elevation_matrix)) / (
                        np.cos(sun_elevation_rad) * np.sin(- sun_azimuth_rad + sky_particle_azimuth_matrix))
    aop_matrix_global = np.arctan(np.tan(np.arctan(tan_aop_l) + sky_particle_azimuth_matrix))
    # here we use "AoP_g=atan(tan (AoP_L - alpha_p))" to have AoP_g between -pi/2 and +pi/2

    # DoLP Array :

    cos_diffusion_angle = x_sun * x_particle + y_sun * y_particle + z_sun * z_particle
    dolp_matrix = (1 - np.multiply(cos_diffusion_angle, cos_diffusion_angle)) / (
            1 + np.multiply(cos_diffusion_angle, cos_diffusion_angle))
    dolp_matrix = np.maximum(0, np.minimum(dolp_max, 1)) * dolp_matrix

    # pixels pointing towards directions under the horizon are not interesting for us
    aop_matrix_global[sky_particle_elevation_matrix < 0] = 10 ** (-10)
    dolp_matrix[sky_particle_elevation_matrix < 0] = 10 ** (-10)

    return aop_matrix_global, dolp_matrix


def simulate_berry_sky(camera, sun_elevation, sun_azimuth, dolp_max, neutral_point_angular_distance_delta):
    """
    Simulates angle of linear polarization and degree of polarization images based on Berry model, for each pixel of a
    given camera, for a given sun position, maximum degree of linear polarization and neutral points angular distance

    Parameters
    ----------
    camera : Camera
        camera used to simulate the polarization pattern
    sun_elevation : float
        sun elevation in degrees
    sun_azimuth : float
        sun azimuth in degrees
    dolp_max : float
        maximum degree of polarization in the sky
    neutral_point_angular_distance_delta : float
        the angular distance between the neutral points in degrees

    Returns
    -------
    aop_matrix_global : ndarray
       image of angle of polarization in the image frame as seen by the camera
    dolp_matrix : ndarray
        image of degree of linear polarization as seen by the camera

    Raises
    ------
    TypeError
        camera should be an instance of the Camera class, and sun_elevation, sun_azimuth, dolp_max and
        neutral_point_angular_distance_delta should be float objects
    """
    if not isinstance(camera, Camera):
        raise TypeError("camera should be an instance of the Camera class")
    if not isinstance(sun_elevation, (int, float)) or not isinstance(sun_azimuth, (int, float)) or not isinstance(
            dolp_max, (int, float)) or not isinstance(neutral_point_angular_distance_delta, (int, float)):
        raise TypeError("sun_elevation, sun_azimuth, dolp_max and neutral_point_angular_distance_delta should be float "
                        "objects")
    sun_elevation_rad = sun_elevation * np.pi / 180
    sun_azimuth_rad = sun_azimuth * np.pi / 180
    neutral_point_angular_distance_delta_rad = neutral_point_angular_distance_delta * np.pi / 180

    sky_particle_azimuth_matrix, sky_particle_elevation_matrix, x_coordinate_pixels, y_coordinate_pixels \
        = camera.get_camera_azimuth_elevation_matrix()

    # Neutral points' complex projection:
    z_babinet = np.tan((np.pi / 2 - sun_elevation_rad - neutral_point_angular_distance_delta_rad / 2) / 2) * np.exp(
        1j * sun_azimuth_rad)
    z_brewster = np.tan(
        (np.pi / 2 - sun_elevation_rad + neutral_point_angular_distance_delta_rad / 2) / 2) * np.exp(
        1j * sun_azimuth_rad)
    z_arago = - 1 / np.conj(z_brewster)
    z_fourth = - 1 / np.conj(z_babinet)
    # Particles' complex projection :
    z_particles = np.multiply(np.tan((np.pi / 2 - sky_particle_elevation_matrix) / 2),
                              np.exp(1j * sky_particle_azimuth_matrix))
    # Berry's function "w" :
    num_w_particules = np.multiply(
        np.multiply(np.multiply(- 4 * (z_particles - z_babinet), (z_particles - z_brewster)),
                    (z_particles - z_arago)), (z_particles - z_fourth))
    den_w_particules = ((1 + np.abs(z_particles) ** 2) ** 2) * np.abs(z_babinet - z_fourth) * np.abs(
        z_brewster - z_arago)
    w = num_w_particules / den_w_particules
    # AoP_g Array
    aop_matrix_global = 0.5 * np.angle(np.multiply(w, np.exp(- 2 * 1j * sun_azimuth_rad)))
    # DoP
    dolp_matrix = dolp_max * (np.abs(w) / (2 - np.abs(w)))
    # See "Polarization singularities in the clear sky" by M V Berry, M R Dennis
    # and R L Lee Jr, in 2004.

    # pixels pointing towards directions under the horizon are not interesting for us
    aop_matrix_global[sky_particle_elevation_matrix < 0] = 10 ** (-10)
    dolp_matrix[sky_particle_elevation_matrix < 0] = 10 ** (-10)

    return aop_matrix_global, dolp_matrix


def simu_sky_intensity_cie(camera, sun_elevation, sun_azimuth, cie_sky_number):
    """
    Simulates an image of relative sky intensity, for each pixel of a given camera, for a given sun position and a given
    CIE model.

    Parameters
    ----------
    camera : Camera
        camera used to simulate the polarization pattern
    sun_elevation : float
        sun elevation in degrees
    sun_azimuth : float
        sun azimuth in degrees
    cie_sky_number : int
        cie sky number between 1 and 15 such as
        1  : CIE standard overcast sky, steep luminance gradation towards zenith azimuthal uniformity
        2  : Overcast, with steep luminance gradation and slight brightening towards the sun
        3  : Overcast, moderately graded with azimuthal uniformity
        4  : Overcast, moderately graded and slight brightening towards the sun
        5  : Sky of uniform luminance
        6  : Partly cloudy sky, no gradation towards zenith, slight brightening towards the sun
        7  : Partly cloudy sky, no gradation towards zenith, brighter circumsolar region
        8  : Partly cloudy sky, no gradation towards zenith, distinct solar corona
        9  : Partly cloudy, with the obscured sun
        10 : Partly cloudy, with brighter circumsolar region
        11 : White-blue sky with distinct solar corona
        12 : CIE standard clear sky, low luminance turbidity
        13 : CIE standard clear sky, polluted atmosphere
        14 : Cloudless turbid sky with broad solar corona
        15 : White blue turbid sky with broad solar corona
    Returns
    -------
    skylight_relative_intensity : ndarray
       image of relative sky intensity for a given CIE model as seen by the camera
    Raises
    ------
    TypeError
        camera should be an instance of the Camera class, sun_elevation and sun_azimuth should be float objects and
        cie_sky_number should be an int objects
    ValueError
        cie_sky_number should be comprised between 1 and 15
    """
    if not isinstance(camera, Camera):
        raise TypeError("camera should be an instance of the Camera class")
    if not isinstance(sun_elevation, (int, float)) or not isinstance(sun_azimuth, (int, float)):
        raise TypeError("sun_elevation and sun_azimuth should be float objects")
    if not isinstance(cie_sky_number, int):
        raise TypeError("cie_sky_number should be an int objects")
    if cie_sky_number not in range(1, 16):
        raise ValueError("cie_sky_number should be comprised between 1 and 15")
    sun_elevation_rad = sun_elevation * np.pi / 180
    sun_azimuth_rad = sun_azimuth * np.pi / 180

    sky_particle_azimuth_matrix, sky_particle_elevation_matrix, x_coordinate_pixels, y_coordinate_pixels \
        = camera.get_camera_azimuth_elevation_matrix()

    cie_coef = np.array([[4.0, - 0.7, 0, - 1, 0], [4.0, - 0.7, 2, - 1.5, 0.15], [1.1, - 0.8, 0, - 1, 0],
                         [1.1, - 0.8, 2, - 1.5, 0.15], [0, - 1, 0, - 1, 0], [0, - 1, 2, - 1.5, 0.15],
                         [0, - 1, 5, - 2.5, 0.3], [0, - 1, 10, - 3, 0.45], [- 1, - 0.55, 2, - 1.5, 0.15],
                         [- 1, - 0.55, 5, - 2.5, 0.3], [- 1, - 0.55, 10, - 3, 0.45], [- 1, - 0.32, 10, - 3, 0.45],
                         [- 1, - 0.32, 16, - 3, 0.3], [- 1, - 0.15, 16, - 3, 0.3], [- 1, - 0.15, 24, - 2.8, 0.15]])
    # here some intermediate outcome
    x_sun = np.cos(sun_elevation_rad) * np.cos(sun_azimuth_rad)
    y_sun = np.cos(sun_elevation_rad) * np.sin(sun_azimuth_rad)
    z_sun = np.sin(sun_elevation_rad)
    sun_zenith_angle_rad = np.pi / 2 - sun_elevation_rad
    x_sky_particles = np.multiply(np.cos(sky_particle_elevation_matrix), np.cos(sky_particle_azimuth_matrix))
    y_sky_particles = np.multiply(np.cos(sky_particle_elevation_matrix), np.sin(sky_particle_azimuth_matrix))
    z_sky_particles = np.sin(sky_particle_elevation_matrix)
    sky_particles_zenith_angle_rad = np.pi / 2 - sky_particle_elevation_matrix
    # Calculus does not work for negative values of sky_particle_elevation_matrix
    sky_particles_zenith_angle_rad[sky_particles_zenith_angle_rad >= np.pi / 2] = np.pi / 2 - 10 ** (-5)
    a_coef = cie_coef[cie_sky_number - 1, 0]
    b_coef = cie_coef[cie_sky_number - 1, 1]
    c_coef = cie_coef[cie_sky_number - 1, 2]
    d_coef = cie_coef[cie_sky_number - 1, 3]
    e_coef = cie_coef[cie_sky_number - 1, 4]
    cos_scatt_angle = (x_sun * x_sky_particles + y_sun * y_sky_particles + z_sun * z_sky_particles)
    cos_square_scatt_angle = cos_scatt_angle ** 2
    scatt_angle_rad = np.arccos(cos_scatt_angle)
    # here phy ratio and f ratio in CIE model (see "Analysis of vertical sky
    # components under various CIE standard general skies" by D.H.W. Li1,
    # C. Li1, S.W. Lou1, E.K.W. Tsang2 and J.C. Lam, in 2015 ):
    phy_cie_ratio = (1 + a_coef * np.exp(b_coef / np.cos(sky_particles_zenith_angle_rad))) / (
            1 + a_coef * np.exp(b_coef))
    f_cie_ratio = (1 + c_coef * (np.exp(d_coef * scatt_angle_rad) - np.exp(
        d_coef * np.pi / 2)) + e_coef * cos_square_scatt_angle) / (1 + c_coef * (
            np.exp(d_coef * sun_zenith_angle_rad) - np.exp(d_coef * np.pi / 2)) + e_coef * (
                                                                           np.cos(sun_zenith_angle_rad) ** 2))
    # end of calculation
    skylight_relative_intensity = np.multiply(phy_cie_ratio, f_cie_ratio)

    return skylight_relative_intensity


def simu_micro_polarizers(camera, simu_model, sun_elevation, sun_azimuth, dolp_max=1, tolerance_deg=0,
                          extinction_ratio=1, cie_sky_number=5, pixel_saturation_ratio=1.1, number_of_bits=16,
                          snr_gauss=10 ** 9, neutral_point_angular_distance_delta=0):
    """
    Simulates a pixel polarization image using an input sky polarization model, as seen by a polarimetric camera,
    for a given sun position relative to the camera, maximum degree of linear polarization (and neutral points angular
    distance, if the simulation is based on Berry model)

    Parameters
    ----------
    camera : DoFPCamera
        camera used to simulate the polarization pattern
    simu_model : str
        model used to simulate the sky polarization pattern ("Berry" or "Rayleigh")
    sun_elevation : float
        sun elevation in degrees
    sun_azimuth : float
        sun azimuth in degrees
    dolp_max : float
        maximum degree of polarization in the sky
    tolerance_deg : float
        mechanical tolerance in polarizer orientation in degrees
    extinction_ratio : float
        If T1 is the intensity transmittance for an incident ray totally linearly polarized along transmission axis and
        T2 is the intensity transmittance for an incident ray totally linearly polarized at 90 degrees from transmission
        axis, then the extinction ratio is (T1-T2)/(T1+T2). Because we work with relative intensity there is no need to
        take into account absorbance, so we suppose null absorbance, which means T1+T2=1.
    cie_sky_number : int
        cie sky number between 1 and 15 such as
        1  : CIE standard overcast sky, steep luminance gradation towards zenith azimuthal uniformity
        2  : Overcast, with steep luminance gradation and slight brightening towards the sun
        3  : Overcast, moderately graded with azimuthal uniformity
        4  : Overcast, moderately graded and slight brightening towards the sun
        5  : Sky of uniform luminance
        6  : Partly cloudy sky, no gradation towards zenith, slight brightening towards the sun
        7  : Partly cloudy sky, no gradation towards zenith, brighter circumsolar region
        8  : Partly cloudy sky, no gradation towards zenith, distinct solar corona
        9  : Partly cloudy, with the obscured sun
        10 : Partly cloudy, with brighter circumsolar region
        11 : White-blue sky with distinct solar corona
        12 : CIE standard clear sky, low luminance turbidity
        13 : CIE standard clear sky, polluted atmosphere
        14 : Cloudless turbid sky with broad solar corona
        15 : White blue turbid sky with broad solar corona
    pixel_saturation_ratio : float
        ratio between the sensor irradiance saturation value and the maximum relative irradiance coming on sensor's
        pixels
    number_of_bits : int
        number of output bits to encode pixel values
    snr_gauss : float
        per pixel gaussian noise Signal to Noise Ratio
    neutral_point_angular_distance_delta : float
        the angular distance between the neutral points in degrees

    Returns
    -------
    bits_matrix_image : DoFPPolarizationImage
       pixel image of the sky as seen by the camera

    Raises
    ------
    TypeError
        Camera should be an instance of DoFPCamera
        simu_model should be a string object
        sun_elevation, sun_azimuth, dolp_max, neutral_point_angular_distance_delta, tolerance_deg, extinction_ratio,
        pixel_saturation_ratio, and snr_gauss should be float objects
        cie_sky_number and number_of_bits should be int objects
    ValueError
        cie_sky_number should be comprised between 1 and 15
    NotImplementedError
        The simulation model should be 'Rayleigh' or 'Berry'
    """
    if not isinstance(camera, DoFPCamera):
        raise TypeError("Camera should be an instance of DoFPCamera")
    if not isinstance(simu_model, str):
        raise TypeError("simu_model should be a string object")
    if not isinstance(sun_elevation, (int, float)) or not isinstance(sun_azimuth, (int, float)) or not isinstance(
            dolp_max, (int, float)) or not isinstance(neutral_point_angular_distance_delta, (int, float)) or not \
            isinstance(tolerance_deg, (int, float)) or not isinstance(extinction_ratio, (int, float)) or not isinstance(
            pixel_saturation_ratio, (int, float)) or not isinstance(snr_gauss, (int, float)):
        raise TypeError("sun_elevation, sun_azimuth, dolp_max, neutral_point_angular_distance_delta, tolerance_deg, "
                        "extinction_ratio, pixel_saturation_ratio, and snr_gauss should be float objects")
    if not isinstance(cie_sky_number, int) or not isinstance(number_of_bits, int):
        raise TypeError("cie_sky_number and number_of_bits should be int objects")
    if cie_sky_number not in range(1, 16):
        raise ValueError("cie_sky_number should be comprised between 1 and 15")

    # create micro-polarizer directions matrix
    tolerance_rad = tolerance_deg * np.pi / 180

    rows, cols = camera.sensor.sensor_dimensions
    sky_radiance = simu_sky_intensity_cie(camera, sun_elevation, sun_azimuth, cie_sky_number)
    if simu_model == "Rayleigh":
        aop_matrix_global_rad, dolp_matrix = simulate_rayleigh_sky(camera, sun_elevation, sun_azimuth,
                                                                   dolp_max)
    elif simu_model == "Berry":
        aop_matrix_global_rad, dolp_matrix = simulate_berry_sky(camera, sun_elevation, sun_azimuth, dolp_max,
                                                                neutral_point_angular_distance_delta)
    else:
        raise NotImplementedError("The simulation model should be 'Rayleigh' or 'Berry'")
    # Matrix of sensor size with only written the 135deg polarizers :
    mat135 = np.abs(((3 * np.pi / 4) * np.sin(np.pi / 2 * (np.arange(1, rows + 1, 1))))[:, np.newaxis] * np.cos(
        np.pi / 2 * (np.arange(1, cols + 1, 1)))[np.newaxis, :])
    # Matrix of sensor size with only written the 90deg polarizers :
    mat90 = np.abs(((np.pi / 2) * np.sin(np.pi / 2 * (np.arange(1, rows + 1, 1))))[:, np.newaxis] * np.sin(
        np.pi / 2 * (np.arange(1, cols + 1, 1)))[np.newaxis, :])
    # Matrix of sensor size with only written the 45deg polarizers :
    mat45 = np.abs(((np.pi / 4) * np.cos(np.pi / 2 * (np.arange(1, rows + 1, 1))))[:, np.newaxis] * np.sin(
        np.pi / 2 * (np.arange(1, cols + 1, 1)))[np.newaxis, :])
    # Matrix of sensor size with angular defects in polarizers' direction
    defects_matrix = tolerance_rad * (1 - 2 * np.random.rand(dolp_matrix.shape[0], dolp_matrix.shape[1]))
    # Matrix of sensor size with micro-polarizer directions
    mat_polarizer_angle = mat90 + mat45 + mat135 + defects_matrix
    intensity_on_pixels_matrix = np.multiply(0.5 * sky_radiance, (1 + np.multiply(
        extinction_ratio * dolp_matrix, np.cos(2 * (aop_matrix_global_rad - mat_polarizer_angle)))))

    # add properties of camera
    with_noise_matrix = np.multiply((1 / snr_gauss) * intensity_on_pixels_matrix,
                                    np.random.randn(intensity_on_pixels_matrix.shape[0],
                                                    intensity_on_pixels_matrix.shape[
                                                        1])) + intensity_on_pixels_matrix
    max_pixel_value = (2 ** number_of_bits) - 1
    max_intensity_comming_on_pixels = np.amax(np.amax(intensity_on_pixels_matrix))
    bits_matrix = np.minimum((np.floor(
        (max_pixel_value / (max_intensity_comming_on_pixels * pixel_saturation_ratio)) * with_noise_matrix)).astype(
        int), max_pixel_value)
    bits_matrix_image = DoFPPolarizationImage(bits_matrix, camera)
    return bits_matrix_image


if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt
    from Polarization_Image_Package.image import Image
    from Camera_Package.dofpCamera import DoFPCamera
    from Camera_Package.dofpSensor import DoFPSensor
    from Camera_Package.lens import Lens
    from Camera_Image_Display_Package.imageDisplay import ImageDisplay

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
    # simulate an image
    im = simu_micro_polarizers(cam, "Berry", 50, 0, dolp_max=1, tolerance_deg=0,
                               extinction_ratio=1, cie_sky_number=5, pixel_saturation_ratio=1.1, number_of_bits=16,
                               snr_gauss=10 ** 9, neutral_point_angular_distance_delta=0)
    aop_image_frame, aop_meridian_frame, dolp = im.get_polarization_processed_image(color_chanel="Gr")
    # display image
    dispIm = ImageDisplay(dolp, "test")
    dispIm.display_image()
    plt.show()
