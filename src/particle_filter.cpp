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
   * Sets the number of particles. Initializes all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   */
  
  std::cout << "Initialising PF with GPS( x, y, theta ) "<< x << " " << y << " " << theta << std::endl;
  
  num_particles = 100;  // Set the number of particles
  
  std::default_random_engine gen;

  // create a normal (Gaussian) distribution for x, y, theta
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  
  for (int i = 0; i < num_particles; ++i) {
    
    Particle particle_temp; // temporary particle variable
    
    // sample from the normal distributions
    // "gen" is the random engine initialized earlier
    particle_temp.id = i;
    particle_temp.x = dist_x(gen);
    particle_temp.y = dist_y(gen);
    particle_temp.theta = dist_theta(gen);
    particle_temp.weight = 1.0;
    
    // pushback to particles vector of the ParticleFilter class
    particles.push_back(particle_temp);
    
  }
  
  is_initialized = true; 
  
  std::cout << "PF Initialisation Complete" << std::endl;
}


void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  
  // Adds measurements to each particle and adds random Gaussian noise.
  
  std::default_random_engine gen;

  // create a normal (Gaussian) distribution for x, y, theta
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);
  
  /**
  * we do not have uncertainties in velocity, yaw_rate 
  * so instead, we add uncertainty to the prediction using a gaussian noise
  */
  
  for (int i = 0; i < num_particles; ++i) {
    
    double theta = particles[i].theta;
    
    if(fabs(yaw_rate) < 0.00001) {
      // using zero yaw-rate equations
      particles[i].x += velocity * delta_t * cos(theta) + dist_x(gen);
      particles[i].y += velocity * delta_t * sin(theta) + dist_y(gen);
      particles[i].theta += dist_theta(gen);
    }    
    else {
      double coeff = velocity/yaw_rate ;
      // using non-zero yaw-rate equations
      particles[i].x += coeff * (sin(theta + yaw_rate*delta_t) - sin(theta)) + dist_x(gen);
      particles[i].y += coeff * (cos(theta) - cos(theta + yaw_rate*delta_t)) + dist_y(gen);
      particles[i].theta += yaw_rate*delta_t + dist_theta(gen);
    }
  }  
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations_map) {
  /**
   * Finds the predicted measurement that is closest to each 
   * observed measurement and assign the observed measurement to this 
   * particular landmark.
   * NOTE: this method is used as a helper function during the update weights phase 
   */
  
  // For each observation
  for (unsigned int i = 0; i < observations_map.size(); ++i) {
    
    double min_dist = std::numeric_limits<double>::max(); 
    
    // Loop though each landmark within sensor range
    for (unsigned int j = 0; j < predicted.size(); ++j) {
      double current_dist = dist(observations_map[i].x, observations_map[i].y, 
                                 predicted[j].x, predicted[j].y);
      
      if(current_dist < min_dist){
        min_dist = current_dist;
        observations_map[i].id = predicted[j].id;
        // added landmark_x, landmark_y variables to LandmarkObs struct
        // makes it easier to find probability in next step
        // instead of searching through map for particular id again
        observations_map[i].landmark_x = predicted[j].x;
        observations_map[i].landmark_y = predicted[j].y;
      }  
    }
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * Updates the weights of each particle using a mult-variate Gaussian distribution
   * For each particle,
   * STEP1: Predict measurements to all landmarks within sensor range  
   * STEP2: Transform sensor measurements to map co-ordinates
   * STEP3: Associate sensor measurements to map landmarks
   * STEP4: Calculate new weight
   * STEP5: Normalise the weights for resampling step ahead
   */
  

  double sum_weights = 0;
  
  // For each particle
  for (int i = 0; i < num_particles; ++i) {
    
    vector<LandmarkObs> predicted;
    // vector of predicted measurements between the particle 
    // and all the map landmarks within sensor range
    
    // STEP1: Find Landmarks in sensor range
    for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); ++j) {
      
      LandmarkObs pred_temp;
      double dist_temp;
      
      pred_temp.x = map_landmarks.landmark_list[j].x_f;
      pred_temp.y = map_landmarks.landmark_list[j].y_f;
      pred_temp.id = map_landmarks.landmark_list[j].id_i;
      
      dist_temp = dist(particles[i].x, particles[i].y, pred_temp.x, pred_temp.y);
      
      if( dist_temp <= sensor_range) { 
        predicted.push_back(pred_temp);
      }
    }
    
        
    /**
    * NOTE: if there is no landmark near the particle 
    * ignore the measurement and move onto next particle
    * but we do not discard this particle (we keep it with its previous weight)
    */
    
    if(predicted.empty()) {
      continue;
    }
    
    // Reset weight of particle
    particles[i].weight = 1.0;
    
    // STEP2: Transform Co-ordinates
    vector<LandmarkObs> observations_map;
    
    for (unsigned int j = 0; j < observations.size(); ++j) {
      
      double theta = particles[i].theta;
      LandmarkObs obs_temp;
      // transform to map x coordinate
      obs_temp.x = particles[i].x + (cos(theta) * observations[j].x) - (sin(theta) * observations[j].y);
      // transform to map y coordinate
      obs_temp.y = particles[i].y + (sin(theta) * observations[j].x) + (cos(theta) * observations[j].y);
      // push to transformed observations vector
      observations_map.push_back(obs_temp);  
      
      
      // STEP3
      
      /** 
      * NOTE: not using the below function, because it creates another for loop
      * that inadvertently slow the program down
      * its utlity is directly incorporated below
      */
      
      // dataAssociation(predicted, observations_map);
      
      double min_dist = std::numeric_limits<double>::max(); 
    
      // Loop though each landmark within sensor range
      for (unsigned int k = 0; k < predicted.size(); ++k) {
        double current_dist = dist(observations_map[j].x, observations_map[j].y, 
                                   predicted[k].x, predicted[k].y);

        if(current_dist < min_dist){
          min_dist = current_dist;
          observations_map[j].id = predicted[k].id;
          // added landmark_x, landmark_y variables to LandmarkObs struct
          // makes it easier to find probability in next step
          // instead of searching through map for particular id again
          observations_map[j].landmark_x = predicted[k].x;
          observations_map[j].landmark_y = predicted[k].y;
        }  
      }

      // STEP4
      double weight_j = multiv_prob(std_landmark[0], std_landmark[1], observations_map[j].x, observations_map[j].y,
                                    observations_map[j].landmark_x, observations_map[j].landmark_y);

      // take product for all pairs
      particles[i].weight *= weight_j;
    }
    sum_weights += particles[i].weight;
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
  uniform_int_distribution<int> dist_i(0,num_particles -1);
  // random index to start at
  int index = dist_i(gen);  
  
  double beta = 0.0;
  double w_max = 0;
  
  // find w_max
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