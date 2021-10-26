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
using namespace std;

//random engine to use across different method calls
static default_random_engine gen;

using std::string;
using std::vector;
using std::normal_distribution;
using std::discrete_distribution;
using std::uniform_real_distribution;
using std::uniform_int_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  if (is_initialized) {
    return;
  }
  num_particles = 150;  // TODO: Set the number of particles

  double std_x = std[0];
  double std_y = std[1];
  double std_theta = std[2];
  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_theta(theta, std_theta);
  
  for (int i = 0; i < num_particles; ++i) {
    
    
    Particle particle;
    particle.id = i;
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1.0;
    
    particles.push_back(particle);
    // Print your samples to the terminal.
    //std::cout << "Sample " << i + 1 << " " << sample_x << " " << sample_y << " " 
            //  << sample_theta << std::endl;
  }
  // Filter is initialized.
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
  double std_x = std_pos[0];
  double std_y = std_pos[1];
  double std_theta = std_pos[2];
  normal_distribution<double> dist_x(0, std_x);
  normal_distribution<double> dist_y(0, std_y);
  normal_distribution<double> dist_theta(0, std_theta);
  for (int i = 0; i < num_particles; ++i) {
    double theta = particles[i].theta;
    if (fabs(yaw_rate)<0.0001)
    {
    particles[i].x = particles[i].x + (velocity)*delta_t*cos(theta);
    particles[i].y = particles[i].y + (velocity)*delta_t*sin(theta);
    particles[i].theta = particles[i].theta;
    }
    else
    {    
    particles[i].x = particles[i].x + (velocity/yaw_rate)*(sin(theta + (yaw_rate*delta_t))-sin(theta));
    particles[i].y = particles[i].y + (velocity/yaw_rate)*(cos(theta)-cos(theta + (yaw_rate*delta_t)));
    particles[i].theta = particles[i].theta + yaw_rate*delta_t;
    }
    
    //ADD Noise
    particles[i].x = particles[i].x + dist_x(gen);
    particles[i].y = particles[i].y + dist_y(gen);
	particles[i].theta = particles[i].theta + dist_theta(gen);

  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations,Particle &particle) {
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
  int num_predicted = predicted.size();
  int num_observed = observations.size();
  
  
  for (unsigned int i = 0; i<num_observed; i++){
    int map_id = -1;
    double min = numeric_limits<double>::max();
    for (unsigned int j=0; j<num_predicted; j++){
      double dist = sqrt((pow(predicted[j].x - observations[i].x, 2)) + (pow(predicted[j].y - observations[i].y, 2)));
      if (dist < min) {
            min=dist;
            map_id = predicted[j].id;;
        }
    }
      observations[i].id = map_id;
    
    associations.push_back(map_id);
	sense_x.push_back(observations[i].x);
	sense_y.push_back(observations[i].y);
    }
  
  // Set assocations for visualization purpose only
  particle = SetAssociations(particle, associations, sense_x, sense_y);

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
  //obs from veh coordinates to landmark coordinates
  
  for (int i = 0; i < num_particles; i++) {

    // get the particle x, y coordinates
    double p_x = particles[i].x;
    double p_y = particles[i].y;
    double theta = particles[i].theta;
    

    // create a vector to hold the map landmark locations predicted to be within sensor range of the particle
    vector<LandmarkObs> predictions;
    
    //Each map landmark for loop
	for(unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {

			//Get id,x,y coordinates
			float lm_x = map_landmarks.landmark_list[j].x_f;
			float lm_y = map_landmarks.landmark_list[j].y_f;
			int lm_id = map_landmarks.landmark_list[j].id_i;
			
			//Consider landmarks within sensor range of particles
			if(fabs(lm_x - p_x) <= sensor_range && fabs(lm_y - p_y) <= sensor_range) {
			predictions.push_back(LandmarkObs{ lm_id, lm_x, lm_y });
			}
		}
    
    //Transform Observed coordinates
    vector<LandmarkObs> mappedObs;
    for(unsigned int j =0; j<=observations.size(); j++){
      
     double mapped_x = p_x + (cos(theta) * observations[j].x) - (sin(theta) * observations[j].y);
     double mapped_y = p_y + (sin(theta) * observations[j].x) + (cos(theta) * observations[j].y);
     mappedObs.push_back(LandmarkObs{ observations[j].id, mapped_x, mapped_y });
    }
    
    //Associate transformed obs to landmark. Find closest landmark to each obs
    dataAssociation(predictions, mappedObs, particles[i]);
    // Set assocations for visualization purpose only 
    //sense_x = getSenseX;
    //sense_y = getSenseY;
    //particle = SetAssociations(particle, associations, sense_x, sense_y);
    
    //Initial weight to 1
    particles[i].weight = 1;
    
    for (unsigned int j=0; j<= mappedObs.size(); j++){
      double obs_x = mappedObs[j].x;
      double obs_y = mappedObs[j].y;
      
      unsigned int LandmarkID = mappedObs[j].id;
      double mu_x, mu_y;
      for (unsigned int k = 0; k<= predictions.size(); j++){
        if(k == LandmarkID){
           mu_x = predictions[k].x;
           mu_y = predictions[k].y;
          break;
        }
      }
      double s_x = std_landmark[0];
      double s_y = std_landmark[1];
      double x_diff = obs_x - mu_x;
      double y_diff = obs_y- mu_y;
      double weight = (1/(2*M_PI*s_x*s_y))*exp(-(((x_diff*x_diff)/(2*s_x*s_x))+ ((y_diff*y_diff)/(2*s_y*s_y))));
      particles[i].weight *= weight;  //Product of this obs weight with total obs weight           
      cout << "> Particle ID#" << particles[i].id << " weight updated"  << ": "<< particles[i].weight << "\n"
				 << endl;               
      
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
  //Weights and Max weight
  vector<double> weights;
  double maxWeight = numeric_limits<double>::min();
  double beta = 0;
  for(int i = 0; i < num_particles; i++) {
		weights.push_back(particles[i].weight);
		if(particles[i].weight > maxWeight) {
			maxWeight = particles[i].weight;
		}
	}
  uniform_real_distribution<double> distDouble(0.0, maxWeight);
  uniform_int_distribution<int> distInt(0, num_particles - 1);

  // Generating index.
  int index = distInt(gen);


  // Resample
  vector<Particle> resampledParticles;
  for(int i = 0; i < num_particles; i++) {
    beta += distDouble(gen) * 2.0;
    while( beta > weights[index]) {
      beta -= weights[index];
      index = (index + 1) % num_particles;
    }
    resampledParticles.push_back(particles[index]);
  }

  particles = resampledParticles;

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