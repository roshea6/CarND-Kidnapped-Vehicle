/**
 * particle_filter.cpp
 *
 * Author: Ryan O'Shea
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
#include <limits>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;

// Create random generator
// Might be better to make this a class variable
static std::default_random_engine gen;

// Std is an array containing the noise values for each of the measurements in the initial reading
void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 500;  // TODO: Set the number of particles

  // Create Gaussian distribution for each of the measurement variables
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  // Loop through the number of particles and create new particles using the Gaussian distributions
  for(int i = 0; i < num_particles; i++)
  {
    // Create a new particle
    Particle new_part;

    // Set ID to particle number
    new_part.id = i;

    // Draw from the corresponding Gaussian distributions for it's state variables
    new_part.x = dist_x(gen);
    new_part.y = dist_y(gen);
    new_part.theta = dist_theta(gen);
    new_part.weight = 1.0;

    // Add the new particle to the vector of Particles
    particles.push_back(new_part);

    // Add a new weight of 1 to the weights vector
    // weights.push_back(double(1.0));
  }

  // Set the status to initialized
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
  // Create random generator
  // std::default_random_engine gen;

  // Loop through all the particles and add the effects of motion to them
  for(int i = 0; i < num_particles; i++)
  {
    Particle part = particles[i];

    // Extract the initial values from the particle
    double x_0 = part.x;
    double y_0 = part.y;
    double theta_0 = part.theta;

    // Calculate the final positions of the particles after motion is applied
    double x_f, y_f, theta_f;
    // Check to make sure yaw rate isn't near 0 to avoid a divide by 0 error which was occasionally happening
    if(abs(yaw_rate) < .00005)
    {
      // If it is then use the equations for straight line driving
      x_f = x_0 + velocity * delta_t * cos(theta_0);
      y_f = y_0 + velocity * delta_t * sin(theta_0);
      theta_f = theta_0;
    }
    else
    {
      // Otherwise use the normal equations
      x_f = x_0 + (velocity/yaw_rate) * (sin(theta_0 + yaw_rate*delta_t) - sin(theta_0));
      y_f = y_0 + (velocity/yaw_rate) * (cos(theta_0) - cos(theta_0 + yaw_rate*delta_t));
      theta_f = theta_0 + yaw_rate*delta_t;
    }
    

    // Create Gaussian distribution for each of updated values
    normal_distribution<double> dist_x(x_f, std_pos[0]);
    normal_distribution<double> dist_y(y_f, std_pos[1]);
    normal_distribution<double> dist_theta(theta_f, std_pos[2]);

    // Draw the final values from each of the distributions
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);
  }

  

}

// Use the nearest neighbor method to figure out which landmark matches with each observation
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

  // Loop through all the observations and find their closest landmarks
  for(int i = 0; i < observations.size(); i ++)
  {
    // Set the best distance to the highest possible starting value
    double best_dist = std::numeric_limits<double>::max();

    // Extract the state variables from the observation for easier referencing
    double obs_x = observations[i].x;
    double obs_y = observations[i].y;

    // Loop through all the landmarks to find the closest one to the point
    for(int j = 0; j < predicted.size(); j++)
    {
      // Extract the state variables from the predicted landmarks for easier referencing
      double pred_x = predicted[j].x;
      double pred_y = predicted[j].y;

      // Get the distance between the predicted point and the landmark
      double distance = sqrt((pow(obs_x - pred_x, 2)) + (pow(obs_y - pred_y, 2)));

      // If a closer landmark is found then update distance and id
      if(distance < best_dist)
      {
        best_dist = distance;
        observations[i].id = predicted[j].id; 
      }
    }
    
  }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a multi-variate Gaussian 
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

  // Loop through all the particles and update their weights based on the observations
  for(int i = 0; i < num_particles; i++)
  {
    double part_x = particles[i].x;
    double part_y = particles[i].y;
    double part_theta = particles[i].theta;

    // Get all the landmarks within sensor range of the particle
    vector<LandmarkObs> landmarks_in_range;
    for(int j = 0; j < map_landmarks.landmark_list.size(); j++)
    {
      // Calculate distance from particle to landmark
      double dist = sqrt((pow(part_x - map_landmarks.landmark_list[j].x_f, 2)) + (pow(part_y - map_landmarks.landmark_list[j].y_f, 2)));

      // If the landmark is in range then add it to the vector
      if(dist < sensor_range)
      {
        // Create a new landmark struct to append to the vector
        LandmarkObs landmark;
        landmark.id = map_landmarks.landmark_list[j].id_i;
        landmark.x = map_landmarks.landmark_list[j].x_f;
        landmark.y = map_landmarks.landmark_list[j].y_f;

        // Add the landmark
        landmarks_in_range.push_back(landmark);
      }
    }

    // Transform the observations into the maps frame
    vector<LandmarkObs> map_frame_obs;
    for(int k = 0; k < observations.size(); k++)
    {
      double obs_x = observations[k].x;
      double obs_y = observations[k].y;
      
      // Transform x and y
      double trans_x = part_x + (cos(part_theta) * obs_x) - (sin(part_theta) * obs_y);
      double trans_y = part_y + (sin(part_theta) * obs_x) + (cos(part_theta) * obs_y);

      // Create a new landmark and populate it with the relevant data
      LandmarkObs trans_obs;
      trans_obs.x = trans_x;
      trans_obs.y = trans_y;
      // trans_obs.id = observations[k].id;

      // Add the transformed observation to the vector
      map_frame_obs.push_back(trans_obs);
    }

    // Use data association to match each observation to its nearest landmark
    dataAssociation(landmarks_in_range, map_frame_obs);

    // Calculate the particle probability based on the observations
    double particle_prob = 1.0;

    for(int j = 0; j < map_frame_obs.size(); j++)
    {
      double obs_x = map_frame_obs[j].x;
      double obs_y = map_frame_obs[j].y;
      double obs_id = map_frame_obs[j].id;

      // Use the observation id to find the coordinates of the best fit landmark
      double mean_x;
      double mean_y;
      for(int k = 0; k < landmarks_in_range.size(); k++)
      {
        // Once the landmark is found extract it's location
        if(obs_id == landmarks_in_range[k].id)
        {
          mean_x = landmarks_in_range[k].x;
          mean_y = landmarks_in_range[k].y;
          break;
        }
      }

      double std_x = std_landmark[0];
      double std_y = std_landmark[1];

      // Use the landmark and observation locations to calculate the probability with a multivariate Gaussian
      double divisor = 2 * M_PI * std_x * std_y;
      double prob = exp(-((pow(obs_x - mean_x, 2)/(2*std_x*std_x)) + (pow(obs_y - mean_y, 2)/(2*std_y*std_y))));

      particle_prob *= prob/divisor;

    }

    // Update the weight of the particle
    particles[i].weight = particle_prob;
    // Add the weight to the weights vector as well so we can easily reference them later
    // weights.push_back(particle_prob);

  }

  // Calculate total weight so we can normalize the probabilites
  double total_weight = 0.0;
  for(int i = 0; i < num_particles; i++)
  {
    total_weight += particles[i].weight;
  }

  // Normalize each weight
  for(int i = 0; i < num_particles; i++)
  {
    // Add the smallest possible values in case the total weight happens to be 0
    particles[i].weight /= total_weight + std::numeric_limits<double>::epsilon();
  }

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  // Create random generator
  // std::default_random_engine gen;

  // Add the all the particle weights to a single vector
  for(int i = 0; i < num_particles; i++)
  {
    weights.push_back(particles[i].weight);
  }

  // Create a weighted distribution with our particle weights as the weights
  std::discrete_distribution<int> weighted_dist(weights.begin(), weights.end());

  // Create a new vector of particles to hold our new samples
  vector<Particle> resampled_parts;

  // Take new samples
  for(int i = 0; i < num_particles; i++)
  {
    // Get the index of the chosen particle
    int idx = weighted_dist(gen);

    // Add the particle to the new set of particles
    resampled_parts.push_back(particles[idx]);
  }

  // Update our set of particles to match the newly generated on
  particles = resampled_parts;

  // Reset particle weights
  for(int i = 0; i < num_particles; i++)
  {
    particles[i].weight = 1.0;
  }
  weights.clear();

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