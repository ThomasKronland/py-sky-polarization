import numpy as np
import matplotlib.pyplot as plt
import Polarization_Simulation_Package.polarizationSimulation as pola_simu
from Camera_Package import DoFPCamera, DoFPSensor, Lens
from Polarization_Image_Package import DoFPPolarizationImageDataset, Image, DoFPPolarizationImage
from Camera_Image_Display_Package import ImageDisplay
from Sun_Ephemeris_Package import SunEphemeris


def main():
    # Camera attributes
    sensor_pixel_size = 3.45
    sensor_dimensions_row_column = (2048, 2448)
    sensor_type = "bayer"
    lens_focal_length = 1.8
    lens_conjugation_model = 'r2'
    camera_center_coordinates = (1024, 1224)

    # extrinsic parameters
    camera_latitude = 43
    camera_longitude = 20
    date = "2024-06-12"
    hours = [15, 15, 16, 16, 16]
    minutes = [20, 40, 0, 20, 40]
    seconds = [0, 0, 0, 0, 0]
    num_images = len(hours)

    # Create a camera
    dofp_sensor = DoFPSensor(sensor_pixel_size, sensor_dimensions_row_column, sensor_type)
    lens = Lens(lens_focal_length, lens_conjugation_model)
    cam = DoFPCamera(dofp_sensor, lens, camera_center_coordinates)

    # Create a simulated polarization image dataset
    polarization_dataset = DoFPPolarizationImageDataset()
    for i in range(num_images):
        sun = SunEphemeris(date=date, hour_utc=hours[i], minute_utc=minutes[i], second_utc=seconds[i],
                           latitude=camera_latitude, longitude=camera_longitude)
        sun_azimuth, sun_altitude = sun.get_sun_horizontal_coordinates()
        polarization_image = pola_simu.simu_micro_polarizers(cam, "Rayleigh", sun_altitude, sun_azimuth)
        polarization_image.camera_longitude = camera_longitude
        polarization_image.camera_latitude = camera_latitude
        polarization_dataset.add_image(polarization_image)

    # Load an image and process it
    sky_pola_array = np.load('image.npy')
    sky_pola_image = DoFPPolarizationImage(image=sky_pola_array, camera=cam)
    _, _, dolp_image = sky_pola_image.get_polarization_processed_image(color_chanel="R")
    sky_pola_image_display = ImageDisplay(image=dolp_image, figure_title="Sky dolp image")
    sky_pola_image_display.display_image()
    plt.show()


if __name__ == "__main__":
    main()
