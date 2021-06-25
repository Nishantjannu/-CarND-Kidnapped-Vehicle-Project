/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;
using std::uniform_int_distribution;
using std::uniform_real_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle. (SAME as Above?)
   // Doesn't make sense otherwise (superimposed noises?????)
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles
  
  std::default_random_engine gen;

  // create a normal (Gaussian) distribution for x, y, theta
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  
  for (int i = 0; i < num_particles; ++i) {
    
    Particle particle_temp; // temporary particle variable
    
    // sample from the normal distributions
    //  "gen" is the random engine initialized earlier
    particle_temp.id = i;
    particle_temp.x = dist_x(gen);
    particle_temp.y = dist_y(gen);
    particle_temp.theta = dist_theta(gen);
    particle_temp.weight = 1.0;
    
    // pushback to particles vector of the ParticleFilter class
    particles.push_back(particle_temp);
    
  }
  
  is_initialized = true; 
  
    
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  
  std::default_random_engine gen;

  // create a normal (Gaussian) distribution for x, y, theta
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);
  
  // in a given time-step, the same noise will be added to all particles
  // we do not have uncertainties in velocity, yaw_rate 
  // so instead, we add uncertainty to the prediction using a gaussian noise
  double noise_x = dist_x(gen);
  double noise_y = dist_y(gen);
  double noise_theta = dist_theta(gen);
  
  
  if(fabs(yaw_rate) < 0.00001){
    for (int i = 0; i < num_particles; ++i) {
      particles[i].x += velocity*delta_t*cos(particles[i].theta) + noise_x;
      particles[i].y += velocity*delta_t*sin(particles[i].theta) + noise_y;
    }    
  }
  
  else{
    double coeff = velocity/yaw_rate ;

    for (int i = 0; i < num_particles; ++i) {
      // using non-zero yaw-rate equations
      particles[i].x += coeff*( sin(particles[i].theta + yaw_rate*delta_t) -  sin(particles[i].theta) ) + noise_x;
      particles[i].y -= coeff*( cos(particles[i].theta + yaw_rate*delta_t) -  cos(particles[i].theta) ) + noise_y;
      particles[i].theta += yaw_rate*delta_t + noise_theta;
    }
  }
  
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  for (int i = 0; i < observations.size(); ++i) {
    
    double min_dist = dist(observations[i].x, observations[i].y, predicted[0].x, predicted[0].y); 
    observations[i].id = predicted[0].id;
    observations[i].landmark_x = predicted[0].x;
    observations[i].landmark_y = predicted[0].y;
    
    for (int j = 1; j < predicted.size(); ++j) {
      double current_dist = dist(observations[i].x, observations[i].y, 
                                 predicted[j].x, predicted[j].y);
      
      if(current_dist < min_dist){
        min_dist = current_dist;
        observations[i].id = predicted[j].id;
        observations[i].landmark_x = predicted[j].x;
        observations[i].landmark_y = predicted[j].y;
      }  
    }
    
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  
  /**
   * For each particle,
   *
   * STEP1: Transform sensor measurements to map co-ordinates
   * STEP2: Predict measurements to all landmarks within sensor range       
   * STEP3: Associate sensor measurements to map landmarks
   * STEP4: Calculate new weight
   * STEP5: Normalise the weights for resampling step ahead
   */
  
  
  // predicted measurements between one particular particle 
  // and all the map landmarks within sensor range
  
  vector<LandmarkObs> predicted;
  vector<LandmarkObs> observations_map;
  double sum_weights = 0;
  
  // For each particle
  for (int i = 0; i < num_particles; ++i) {
    
    // STEP1
    for (int j = 0; j < observations.size(); ++j) {

      // transform to map x coordinate
      LandmarkObs obs_temp;
      obs_temp.x = particles[i].x + (cos(particles[i].theta) * observations[j].x) 
                   - (sin(particles[i].theta) * observations[j].y);
      // transform to map y coordinate
      obs_temp.y = particles[i].x + (sin(particles[i].theta) * observations[j].x) 
                + (cos(particles[i].theta) * observations[j].y);
    
      observations_map.push_back(obs_temp);  
    }
    
    // STEP2  
    for (int j = 0; j < map_landmarks.landmark_list.size(); ++j) {
      
      LandmarkObs pred_temp;
      
      pred_temp.x = map_landmarks.landmark_list[j].x_f;
      pred_temp.y = map_landmarks.landmark_list[j].y_f;
      pred_temp.id = map_landmarks.landmark_list[j].id_i;
      
      if( dist(particles[i].x, particles[i].y, pred_temp.x, pred_temp.y) <= sensor_range) { 
        predicted.push_back(pred_temp);
      }
    }
    
    // STEP3
    dataAssociation(predicted, observations_map);
    
    // STEP4
    double weight_particle = 1.0;
    
    for (int j = 0; j < observations_map.size(); ++j) {
      
      double weight_j = multiv_prob(std_landmark[0], std_landmark[1], observations_map[j].x, observations_map[j].y,
                           observations_map[j].landmark_x, observations_map[j].landmark_y);
      
      weight_particle *= weight_j;
    }
    
    particles[i].weight = weight_particle;
    sum_weights += particles[i].weight;
           
    predicted.clear();
    observations_map.clear();
  }
  
  // STEP5
  // For each particle
  for (int i = 0; i < num_particles; ++i) {
    particles[i].weight = (particles[i].weight)/sum_weights;
  }  
  
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  
  // Resampling Wheel
  
  vector<Particle> particles_next;
  
  std::default_random_engine gen;
  uniform_int_distribution<int> dist_i(0,num_particles-1);
  int index = dist_i(gen);  
  
  double beta = 0.0;
  double w_max = 0;
  
  for (int i = 0; i < num_particles; ++i) {
    if(particles[i].weight > w_max) {
      w_max = particles[i].weight;
    }  
  }
  
  uniform_real_distribution<double> dist_r(0.0,2*w_max);
  
  for (int i = 0; i < num_particles; ++i) {
    beta += dist_r(gen);
    while(particles[index].weight < beta) {
      beta -= particles[index].weight;
      index = (index + 1)%num_particles;
    }
    particles_next.push_back(particles[index]);
  }
  
  particles.clear();
  particles = particles_next;
  particles_next.clear();

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}