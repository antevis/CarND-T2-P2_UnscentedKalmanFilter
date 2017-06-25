[//]: # (Image References)
[lidarNIS]: ./data/LidarNIS.png
[radarNIS]: ./data/LidarNIS.png

# Unscented Kalman Filter Project
**Self-Driving Car Engineer Nanodegree Program**

`main.cpp` has been modified to save **NIS** data to file `data/nis.txt` (lines 111 - 118).

With the chosen values for `std_a_` and `std_yawdd_` the **NIS** charts look like this:

![alt text][LidarNIS]

![alt text][radarNIS]

This project involves the [Term 2 Simulator](https://github.com/udacity/self-driving-car-sim/releases).

The project requires [uWebSocketIO](https://github.com/uWebSockets/uWebSockets) to be installed for either Linux or Mac systems.

Once the install for **uWebSocketIO** is complete, the main program can be built and ran by doing the following from the project top directory.

1. `mkdir build`
2. `cd build`
3. `cmake ..`
4. `make`
5. `./UnscentedKF`
