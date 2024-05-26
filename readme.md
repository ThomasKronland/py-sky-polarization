# py-sky-polarization
A python library to process skylight polarization images

## Description

This project is composed of the following packages:
* Camera_Package \
This package is made to create a camera object based on specifications on the lens 
and on the sensor. \
This package is, for now, made specifically for division of focal plane (DoFP) polarimetric cameras, but can easily be adapted 
to other types of sensors. \
This package contains the following classes :
  * Lens : to create a lens object from its focal length and conjugation model
  * Sensor : to create a sensor object from its pixel size, number of pixels in row and column, type (monochrome or bayer)
    * DoFPSensor : subclass of Sensor, to create a division of focal plane polarimetric sensor with the same attributes as input as a Sensor object
  * Camera : to create a camera object from a sensor and a lens objects
    * DoFPCamera : subclass of Camera, to create a division of focal plane polarimetric camera from a lens and a dofpsensor object \
  
    For a given camera, the azimuth and elevation corresponding to the direction each pixel points at can be obtained from 
the function "get_camera_azimuth_elevation_matrix()".


* Sun_Ephemeris_Package \
This package is made to obtain the position of the Sun for a given date and position on Earth.\
It contains the following class:
  * SunEphemeris : to create a sun object from the date, the time in utc and the location in latitude and longitude \
  The Sun position can then be obtained either in azimuth and elevation, hour angle and declination or right ascension and 
  declination (last not requiring location) (based on AstroPy library)


* Polarization_Image_Package \
This package is made to associate an image with some metadata such as the camera used to take the image, the date, time 
and location at which the image was taken.\
This package contains the following classes:
  * Image : to create an image object from an image array, a camera object, the date, the time in utc and, the location 
  in latitude and longitude of the camera and the camera's exposure time
    * DoFPPolarizationImage : subclass of Image, same but requires a dofp camera 
    * DoFPProcessedPolarizationImage : subclass of Image, for processed dofp image such as angle or degree of polarization images.
  * ImageDataset : to create a dataset of image objects
    * DoFPPolarizationImageDataset : subclass of ImageDataset, to create a dataset of dofp images
    * DoFPProcessedPolarizationImageDataset : subclass of ImageDataset, to create a dataset of dofp processed images 


* Polarization_Simulation_Package \
This package is made to simulate skylight polarization images based essentialy on the OpenSky simulator (cite). \
It contains only functions.


* Camera_Image_Display_Package \
This package is made to display an image. \
This package contains the following class:
  * ImageDisplay : to create a figure from an image object

* Mathematical_Tools_Package \
Math functions for calculus

## Getting Started

### Requirements

Cf. requirements.txt

### Installing

This code can be used alone, for simple processing of polarization images, but is intended to simplify tests of skylight 
polarization image processing algorithms.\
Functionalities may of course be added depending on the users needs.

### Executing program

* Start by importing some packages
  ```
  import numpy as np
  import matplotlib.pyplot as plt
  import Polarization_Simulation_Package.polarizationSimulation as pola_simu
  from Camera_Package import DoFPCamera, DoFPSensor, Lens
  from Polarization_Image_Package import DoFPPolarizationImageDataset, Image, DoFPPolarizationImage
  from Camera_Image_Display_Package import ImageDisplay
  from Sun_Ephemeris_Package import SunEphemeris
  ```
* Then create a camera object, for instance
  ```
  sensor_pixel_size = 3.45
  sensor_dimensions_row_column = (2048, 2448)
  sensor_type = "bayer"
  lens_focal_length = 1.8
  lens_conjugation_model = 'r2'
  camera_center_coordinates = (1024, 1224)
  
  # Create a camera
  dofp_sensor = DoFPSensor(sensor_pixel_size, sensor_dimensions_row_column, sensor_type)
  lens = Lens(lens_focal_length, lens_conjugation_model)
  cam = DoFPCamera(dofp_sensor, lens, camera_center_coordinates)
  ```
* Then you can, for instance, simulate a polarization image
  * by first finding the position of the Sun at a given position and a given time
    ```
    camera_latitude = 43
    camera_longitude = 20
    date = "2024-06-12"
    hour_utc = 15
    minute_utc = 20
    second_utc = 0
    sun = SunEphemeris(date=date, hour_utc=hour_utc, minute_utc=minute_utc, second_utc=second_utc,
                           latitude=camera_latitude, longitude=camera_longitude)
    sun_azimuth, sun_altitude = sun.get_sun_horizontal_coordinates()
    ```
  * and then simulate the polarization image for your camera, based on Rayleigh single scattering or Berry model
    ```
    polarization_image = pola_simu.simu_micro_polarizers(cam, "Rayleigh", sun_altitude, sun_azimuth)
    polarization_image.camera_longitude = camera_longitude
    polarization_image.camera_latitude = camera_latitude
    ```
* Or you can upload an image taken with your own camera
  ```
  sky_pola_array = np.load('image.npy')
  sky_pola_image = DoFPPolarizationImage(image=sky_pola_array, camera=cam)
  ```
* Then you can process the image (simulated or real) to obtain the angle or degree of polarization (aop and dolp respectively) in the sky (by choosing here the color chanel you want for a bayer camera, else not required)
  ```
  aop_g, aop_l, dolp_image = polarization_image.get_polarization_processed_image(color_chanel="R")
  ``` 
* You can also gather you images in a dataset for future work
  ```
  polarization_dataset = DoFPPolarizationImageDataset()
  polarization_dataset.add_image(polarization_image)
  polarization_dataset.add_image(sky_pola_image)
  ```
* Finally, you can display your images 
  ```
  dolp_image_display = ImageDisplay(image=dolp_image, figure_title="Sky dolp image", cmap='jet')
  sky_pola_image_display.display_image()
  plt.show()
  ```



## Help

Any advise for common problems or issues, see the documentation with the code or contact me

## Authors

Thomas Kronland-Martinet  
mail : thomas.kronland-martinet@univ-amu.fr

## Version History

* 0.1
    * Initial Release

## License

This project is licensed under the MIT License - see the LICENSE.md file for details

## Acknowledgments
This works contains code taken from Antoine Moutenet and Léo Poughon's published work:
* [skylight polarization simulation code](https://github.com/MoutenetA/OpenSky)
* [skylight polarization database you may consider using](https://github.com/mol-1/Long-Term-Skylight-Polarization-Measurement-Device)
