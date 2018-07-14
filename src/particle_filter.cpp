#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <limits>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	// TODO Complete

	/*Number is particles are chosen in order to run the algorithm in almost real time and 
	introduce lowest possible error in localization. This is a tunable parameter.*/
	num_particles = 20;
	default_random_engine gen;
	weights.resize(num_particles);
  	particles.resize(num_particles);

	// Standard deviations for GPS measurements
	double std_x = std[0];
	double std_y = std[1];
	double std_theta = std[2];
	

	// Normal distributions for coordinates and orientation
	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);

	
	for (int i = 0; i < num_particles; ++i) {

		// Sample from normal distribution
		double sample_x = dist_x(gen);
		double sample_y = dist_y(gen);
		double sample_theta = dist_theta(gen);
		particles[i].x=sample_x;
		particles[i].y=sample_y;
		particles[i].theta=sample_theta;
		particles[i].weight=1;
		weights[i]=1;

		//cout << "Init: Particle " << i << "\tX:"<<sample_x<< "\tY:"<<sample_y<< "\tThe:"<<sample_theta<<endl;
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	//TODO complete

	default_random_engine gen;
	double std_x = std_pos[0];
	double std_y = std_pos[1];
	double std_theta = std_pos[2];
	
	for (unsigned int j = 0; j < particles.size(); ++j){
		double curX, curY, curTheta;
		curX = particles[j].x;
		curY = particles[j].y;
		curTheta = particles[j].theta;

		if (fabs(yaw_rate) < 0.0001) {
			curX = curX + velocity * cos(curTheta) * delta_t;
			curY = curY + velocity * sin(curTheta) * delta_t;
		} else {
			curX = curX + (velocity/yaw_rate) * (sin(curTheta + (yaw_rate * delta_t)) - sin(curTheta));
			curY = curY + (velocity/yaw_rate) * (cos(curTheta) - cos(curTheta + (yaw_rate * delta_t)));
			curTheta = curTheta + (yaw_rate * delta_t);
		}
		
		normal_distribution<double> dist_x(curX, std_x);
		normal_distribution<double> dist_y(curY, std_y);
		normal_distribution<double> dist_theta(curTheta, std_theta);
		
		particles[j].x = dist_x(gen);
		particles[j].y = dist_y(gen);
		particles[j].theta = dist_theta(gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// Implements nearest neighbour alg between the observations and predicted measurements
	for (unsigned int i = 0; i < observations.size(); i++) {
		double lowest_dist = numeric_limits<double>::max();
		int closest_landmark_id = -1;
		double obs_x = observations[i].x;
		double obs_y = observations[i].y;

		for (unsigned int j = 0; j < predicted.size(); j++) {
		  double pred_x = predicted[j].x;
		  double pred_y = predicted[j].y;
		  int pred_id = predicted[j].id;
		  double current_dist = dist(obs_x, obs_y, pred_x, pred_y);
		  if (current_dist < lowest_dist) {
		    lowest_dist = current_dist;
		    closest_landmark_id = pred_id;
		  }
		}
		observations[i].id = closest_landmark_id;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		 const std::vector<LandmarkObs> &observations,
		const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	//   TODO complete

  

  double weights_total = 0;
  const double std_x = std_landmark[0];
  const double std_y = std_landmark[1];
  const double std_x_2 = pow(std_x, 2);
  const double std_y_2 = pow(std_y, 2);
  const double gauss_norm = (1.0/(2.0 * M_PI * std_x * std_y));

  for (unsigned int i = 0; i < num_particles; i++) {
    double particle_x = particles[i].x;
    double particle_y = particles[i].y;
    double particle_theta = particles[i].theta;

    /*Step 1: Transform observations from vehicle co-ordinates to map co-ordinates.*/
    //Vector containing observations transformed to map co-ordinates w.r.t. current particle.
    vector<LandmarkObs> particle_observations;
	vector<LandmarkObs> predicted_landmarks;

    //Transform observations from vehicle's co-ordinates to map co-ordinates.
    for (unsigned int j = 0; j < observations.size(); j++) {
      LandmarkObs particle_obs;
      particle_obs.id = j;
      particle_obs.x = particle_x + (cos(particle_theta) * observations[j].x) - (sin(particle_theta) * observations[j].y);
      particle_obs.y = particle_y + (sin(particle_theta) * observations[j].x) + (cos(particle_theta) * observations[j].y);
      particle_observations.push_back(particle_obs);
    }

    /* Retrieve the predicted landmark measurements */
    
    for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {
      Map::single_landmark_s current_landmark = map_landmarks.landmark_list[j];
      if ((fabs((particle_x - current_landmark.x_f)) <= sensor_range) && (fabs((particle_y - current_landmark.y_f)) <= sensor_range)) {
        predicted_landmarks.push_back(LandmarkObs {current_landmark.id_i, current_landmark.x_f, current_landmark.y_f});
      }
    }

    /*Step 3: Associate observations to lpredicted andmarks using nearest neighbor algorithm.*/
    //Associate observations with predicted landmarks
    dataAssociation(predicted_landmarks, particle_observations);

    /* Weight calculation through multivariate gaussian */
    particles[i].weight = 1.0;

    particles[i].associations.clear();
	particles[i].sense_x.clear();
	particles[i].sense_y.clear();

    for (unsigned int k = 0; k < particle_observations.size(); k++) {
      double obs_x = particle_observations[k].x;
      double obs_y = particle_observations[k].y;
      int obs_id = particle_observations[k].id;
      double multi_prob = 1.0;

	  
	  particles[i].associations.push_back(obs_id);
	  particles[i].sense_x.push_back(obs_x);
	  particles[i].sense_y.push_back(obs_y);
      
	  for (unsigned int l = 0; l < predicted_landmarks.size(); l++) {
        double pred_x = predicted_landmarks[l].x;
        double pred_y = predicted_landmarks[l].y;
        int pred_id = predicted_landmarks[l].id;

        if (obs_id == pred_id) {
          multi_prob = gauss_norm * exp(-1.0 * ((pow((obs_x - pred_x), 2)/(2.0 * std_x_2)) + (pow((obs_y - pred_y), 2)/(2.0 * std_y_2))));
          particles[i].weight *= multi_prob;
        }
      }
    }
    weights_total += particles[i].weight;
  }

	// Weights normalization
  for (int i = 0; i < particles.size(); i++) {
    particles[i].weight /= weights_total;
    weights[i] = particles[i].weight;
  }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	//   TODO complete

	vector<Particle> resampled_particles;

	// Create a generator to be used for generating random particle index and beta value
	default_random_engine gen;
	
	//Generate random particle index
	uniform_int_distribution<int> particle_index(0, num_particles - 1);
	
	int current_index = particle_index(gen);
	
	double beta = 0.0;
	
	double max_weight_2 = 2.0 * *max_element(weights.begin(), weights.end());
	
	for (int i = 0; i < particles.size(); i++) {
		uniform_real_distribution<double> random_weight(0.0, max_weight_2);
		beta += random_weight(gen);

	  while (beta > weights[current_index]) {
	    beta -= weights[current_index];
	    current_index = (current_index + 1) % num_particles;
	  }
	  resampled_particles.push_back(particles[current_index]);
	}
	particles = resampled_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
		                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

  return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}