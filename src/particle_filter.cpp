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
#include "particle_filter.h"

using std::string;
using std::vector;
using namespace std;
using std::string;
using std::vector;
using std::normal_distribution;
using std::discrete_distribution;
using std::uniform_real_distribution;
using std::uniform_int_distribution;
static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles
  
  default_random_engine gen;
  
  // Create a normal (Gaussian) distribution for x, y, and theta
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
  
  for(int i = 0; i < num_particles; i++) {
		Particle particle;

		particle.id = i;
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);
		particle.weight = 1.0;

		particles.push_back(particle);

		weights.push_back(1.0);

		cout << "Particle ID#" << particle.id << " initialized: "
				 << particle.x << ", " << particle.y << ", " << particle.theta
				 << endl;
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
  // Create Gaussian distribution for x, y, and theta
  
  double x;
  double y;
  double theta;
	default_random_engine gen;
	normal_distribution<double> noise_x(0.0, std_pos[0]);
	normal_distribution<double> noise_y(0.0, std_pos[1]);
	normal_distribution<double> noise_theta(0.0, std_pos[2]);
  
    for (int i = 0; i < num_particles; ++i) {
      
    theta = particles[i].theta;
      x = particles[i].x;
      y = particles[i].y;
      
      
    if (fabs(yaw_rate)<0.0001)
      
    {
    x += velocity * delta_t * cos(theta);
    y += velocity * delta_t * sin(theta);
    //theta = particles[i].theta;
    }
      
    else
      
    {    
    x += (velocity/yaw_rate)*(sin(theta +yaw_rate*delta_t)-sin(theta));
    y += (velocity/yaw_rate)*(cos(theta)-cos(theta + yaw_rate*delta_t));
    theta += yaw_rate*delta_t;
    }
      
       //ADD Noise
particles[i].x = x + noise_x(gen);
particles[i].y = y + noise_y(gen);
particles[i].theta = theta + noise_theta(gen);

}
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations, Particle &particle) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  
  	std::vector<int> associations;
	std::vector<double> sense_x;
	std::vector<double> sense_y;
    
  
for (unsigned int i = 0; i < observations.size(); ++i)
  {
  double x_obs = observations[i].x;
  double y_obs = observations[i].y;
  
    double min_dist = std::numeric_limits<double>::max(); 
    //int id_in_map = -1; 
    for (unsigned int j = 0; j < predicted.size(); ++j)
    {
      double x_pred = predicted[j].x;
      double y_pred = predicted[j].y;
      double dx = x_obs-x_pred;
      double dy = y_obs-y_pred;
      double distance = sqrt(dx*dx+dy*dy);
      double diff = distance;
      if (diff < min_dist)
      {
        min_dist = diff;
        observations[i].id = predicted[j].id;
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
  // for each particle...
  double sig_x = std_landmark[0];
  double weight_sum = 0;
  
  double sig_y = std_landmark[1];
  double gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);
  
  for (int i = 0; i < num_particles; ++i) {

    // get the particle x, y coordinates
    double p_x = particles[i].x;
    double p_y = particles[i].y;
    double p_theta = particles[i].theta;
        // reinit weight
    particles[i].weight = 1.0;
    
    //Loop over observations and transform
    vector<LandmarkObs> Transformed_OBS;
    for (unsigned int j = 0; j < observations.size(); ++j)
    {
      double x_map = (observations[j].x * cos(p_theta)) - (observations[j].y * sin(p_theta)) + p_x;
      double y_map = (observations[j].x * sin(p_theta)) + (observations[j].y * cos(p_theta)) + p_y;
      Transformed_OBS.push_back (LandmarkObs{observations[j].id, x_map, y_map}); 

      //Loop over landmarks
      vector<LandmarkObs> predicted; 
    for (unsigned int k = 0; k < map_landmarks.landmark_list.size(); ++k)
    {
      
      int land_ID = map_landmarks.landmark_list[k].id_i;
      float land_x = map_landmarks.landmark_list[k].x_f;
      float land_y = map_landmarks.landmark_list[k].y_f;
      double pred_distance = dist(p_x, p_y, land_x, land_y);
      if (pred_distance <= sensor_range)
      {
        predicted.push_back(LandmarkObs{land_ID, land_x, land_y});
      }
    }
      
      //Calculate weight contribution of single observation and multiple to weight
      
      // Step3. Nearest Neighbor Data Association
    dataAssociation(predicted, Transformed_OBS, particles[i]);
    // Step4. Compute WEIGHT of particle
    vector<int> association;
    vector<double> sense_x;
	vector<double> sense_y;
    
    particles[i].weight = 1.0;
    double map_x, map_y, mu_x, mu_y;
    for (unsigned int t = 0; t < Transformed_OBS.size(); ++t)
    {
      map_x =  Transformed_OBS[t].x;
      map_y =  Transformed_OBS[t].y;
      for (unsigned int p = 0; p < predicted.size(); ++p)
      {
        // Associate prediction with transformed observation
        if (predicted[p].id == Transformed_OBS[t].id)
        {
          mu_x = predicted[p].x;
          mu_y = predicted[p].y;
        }
      }
      
      // Compute exponent
      double exponent = (0.5*pow( (map_x - mu_x)/sig_x, 2.0 )+0.5*pow( (map_y - mu_y)/sig_y, 2.0 ));
      // Compute weight using normalization terms and exponent
      double p_weight = gauss_norm * exp(- exponent);
      particles[i].weight *= p_weight;
    
            // Append particle associations
      association.push_back(Transformed_OBS[t].id);
      sense_x.push_back(map_x);
      sense_y.push_back(map_y);
    }
weights[i] = particles[i].weight;
      weight_sum += weights[i];
  }
    
    
}
  
  //Normalize weights
    if (fabs(weight_sum)> 0.0)
    {
      for (unsigned int m = 0; m < weights.size() ; ++m)
      {
        weights[m] = weights[m]/weight_sum;
      }
    }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
vector<Particle> new_particles;
  vector<double> weights;
  default_random_engine gen;

  for (int i = 0; i < num_particles; i++) {
    weights.push_back(particles[i].weight);
  }

  uniform_int_distribution<int> particle_idx(0, num_particles-1);
  int index = particle_idx(gen);

  double max_weight = *max_element(weights.begin(), weights.end());

  double beta = 0.0;

  uniform_real_distribution<double> rand_weight(0.0, max_weight);

  for (int i = 0; i < num_particles; i++) {
    
    beta += rand_weight(gen) * 2.0;
    while (beta > weights[index]) {
      beta -= weights[index];
      index = (index + 1) % num_particles;
    }
    new_particles.push_back(particles[index]);

  }

  particles = new_particles;

}

Particle ParticleFilter::SetAssociations(Particle& particle, 
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
  
  return particle;
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