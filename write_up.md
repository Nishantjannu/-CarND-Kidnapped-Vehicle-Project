# **Kidnapped Vehicle** 

---

**Particle Filter (Markov Localization) Project**

The goals / steps of this project are the following:
* Implement a 2 dimensional particle filter in C++.
* Particle filter will be given a map and some initial localization information (analogous to what a GPS would provide).
* At each time step filter will also get observation and control data.

[//]: # (Image References)

[image1]: /success.png 
[image2]: /pf.png 

## Rubric Points
### The [rubric points](https://review.udacity.com/#!/rubrics/747/view) for this project are: 
* Particle filter localizes the vehicle to within the desired accuracy
* Particle filter runs within the specified time of 100 seconds
---

### Reading in the Data

The code that will read in and parse the data files is in the main.cpp file. 

First, it reads in the map data into a map object using the read_map_data helper function

The main.cpp file then creates an instance of a ParticleFilter. 

The code reads in the data file line by line. If the filter is uninitialised, it parses the data into GPS co-ordinates and calls the init function. Otherwise, it parses in the control data and calls the prediction function.

Next, the landmark obervation data is pushed onto a noisy_observations list and the updateWeights, resample functions are called.


### File Structure

#### Files in the Github src Folder
* main.cpp - communicates with the Term 2 Simulator receiving data measurements, calls functions to run the Particle filter and calculates Average, Highest weights at a given time step.
* particle_filter.cpp - defines the init, prediction, updateWeights, resample functions
* map.h - defines the map class
* helper_functions.h - defines a few helper functions which can be inserted inline and a few useful structs.

### Final Results

## Simulator
 
Below is an image of what it looks like when the simulator successfully is able to track the car to a particle. 

Notice that the green laser sensors from the car nearly overlap the blue laser sensors from the particle,

this means that the particle transition calculations were done correctly.

![Simulator][image1]

The error between the ground truth and particle filter data into the simulator is within predefined limits.

![Simulator][image2]
