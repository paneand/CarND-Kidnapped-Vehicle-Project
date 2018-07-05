/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

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
	default_random_engine gen;

	// Set particles number and resize accordingly weights and particles vectors

	num_particles = 100;			
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

		// Create particle with sampled coordinates and orientation and put it in the particles vector with weight equals to 1
		 Particle particle = {};
		 particle.x=sample_x;
		 particle.y=sample_y;
		 particle.theta=sample_theta;
		 particle.weight=1;
		 particles.push_back(particle);
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	vector<Particle>::iterator it;
	default_random_engine gen;
	double std_x = std_pos[0];
	double std_y = std_pos[1];
	double std_theta = std_pos[2];
	for (it = particles.begin(); it != particles.end(); ++it){
        Particle curPar = *it;
		double curX, curY, curTheta;
		curX = curPar.x;
		curY = curPar.y;
		curTheta = curPar.theta;
		if(fabs(yaw_rate)>0.0001){			// yaw_rate !=0 
		curX = curX + velocity/yaw_rate * (sin(curTheta+(yaw_rate*delta_t))-sin(curTheta));
		curY = curY + velocity/yaw_rate * (cos(curTheta)-cos(curTheta+(yaw_rate*delta_t)));
		curTheta = curTheta + yaw_rate*delta_t;
		}	// yaw_rate = 0
		else{
			curX = curX + cos(curTheta)*velocity*delta_t;					// TO BE CHECKED
			curY = curY + sin(curTheta)*velocity*delta_t;
		}
		// Add gaussian noise
		normal_distribution<double> dist_x(curX, std_x);
		normal_distribution<double> dist_y(curY, std_y);
		normal_distribution<double> dist_theta(curTheta, std_theta);
		curPar.x = dist_x(gen);
		curPar.y = dist_y(gen);
		curPar.theta = dist_theta(gen);
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	vector<LandmarkObs>::iterator itLand;
	for (itLand = observations.begin(); itLand != observations.end(); ++itLand){
        LandmarkObs curObs = *itLand;
		double minDistance = numeric_limits<double>::max();
		for (unsigned int l=0; l< predicted.size(); ++l) {
			LandmarkObs curPred = predicted[l];
			double diff_x, diff_y, distance;
			diff_x = curObs.x - curPred.x;
			diff_y = curObs.y - curPred.y;
			distance = sqrt((diff_x*diff_x)+(diff_y*diff_y));
			if(distance < minDistance){
				curObs.id=curPred.id;
				minDistance = distance;
			}
		}
	}



}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {

	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html


	// Gaussian multivariate constant terms 
	const double gauss_norm = 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]);
  	const double x_denom = 2 * std_landmark[0] * std_landmark[0];
  	const double y_denom = 2 * std_landmark[1] * std_landmark[1];


	// Iterators for particles and LandmarksObs
	vector<Particle>::iterator itPart;
	vector<LandmarkObs>::const_iterator itLand;

	double totalWeights = 0; // Used for weight normalization




	for (itPart = particles.begin(); itPart != particles.end(); ++itPart){
        Particle curPar = *itPart;
		vector<LandmarkObs> particleObservations;
		vector<LandmarkObs> particlePredictedObservations;			// Pseudo-ranges

		// Transform observation in particle coordinate system
		for (itLand = observations.begin(); itLand != observations.end(); ++itLand){
			// Transform current observation from car reference system to current particle one
			// Using homogeneous transform
			LandmarkObs carObs = *itLand;
			LandmarkObs actualObsFromParticle;

			double x_part = curPar.x;
			double y_part = curPar.y;
			double theta = curPar.theta;
			
			double x_obs = carObs.x;
			double y_obs = carObs.y;

			// transform to map x coordinate
			actualObsFromParticle.x = x_part + (cos(theta) * x_obs) - (sin(theta) * y_obs);

			// transform to map y coordinate
			actualObsFromParticle.y = y_part + (sin(theta) * x_obs) + (cos(theta) * y_obs);

			particleObservations.push_back(actualObsFromParticle);
		}


		// Estimates pseudo ranges from particle perspective
        for (unsigned int l=0; l< map_landmarks.landmark_list.size(); ++l) {
			LandmarkObs predictedObs;
			predictedObs.x = map_landmarks.landmark_list[l].x_f - curPar.x;
			predictedObs.y = map_landmarks.landmark_list[l].y_f - curPar.y;
			predictedObs.id = map_landmarks.landmark_list[l].id_i;
			double obsDistance = sqrt((predictedObs.x*predictedObs.x)+(predictedObs.y*predictedObs.y));
			if(obsDistance<=sensor_range){
				particlePredictedObservations.push_back(predictedObs);
			}
        }

		// Data association:
		// Foreach measurement, associate the corresponding nearest landmark
		dataAssociation(particlePredictedObservations,particleObservations);


		
		// Update weight and stores into the particle the perceived landmarks

		double currentWeight = 1;

		vector<int> associations;
		vector<double> sense_x,sense_y;

		for (unsigned int l=0; l< particleObservations.size(); ++l) {
			int id = particleObservations[l].id;
			double x_obs = particleObservations[l].x;
			double y_obs = particleObservations[l].y;

			associations.push_back(id);
			sense_x.push_back(x_obs);
			sense_y.push_back(y_obs);

			double mu_x, mu_y;
			for (unsigned int l=0; l< map_landmarks.landmark_list.size(); ++l) {
				if(map_landmarks.landmark_list[l].id_i==id){
					mu_x = map_landmarks.landmark_list[l].x_f; 
					mu_y = map_landmarks.landmark_list[l].y_f; 
					break;
				}
			}
			double xd2 = (x_obs - mu_x)*(x_obs - mu_x);
			double yd2 = (y_obs - mu_y)*(y_obs - mu_y);

			double exponent= xd2/x_denom + yd2/y_denom;
			double multivariate = pow(gauss_norm,-exponent);
			
			currentWeight*=multivariate;
		}

		SetAssociations(curPar,associations,sense_x,sense_y);
		
		curPar.weight = currentWeight;
		totalWeights+=currentWeight;
	}

	// Weight normalization
	for (itPart = particles.begin(); itPart != particles.end(); ++itPart){
		Particle curPar = *itPart;
		curPar.weight=curPar.weight/totalWeights;
		currentWeightsArray.push_back(curPar.weight);
	}

}

void ParticleFilter::resample() {
	// Vector for new particles
  	vector<Particle> new_particles (num_particles);
  	default_random_engine gen;
  	for (int i = 0; i < num_particles; ++i) {
		
		// Discrete distribution between 0 and num_particles-1 where each index of the particles vector can be extracted
		// with the corresponding particle's weight probability 
    	discrete_distribution<int> index(weights.begin(), weights.end());
    	
		// The same particle could be pick up multiple times
		new_particles[i] = particles[index(gen)];
    
  }
  
  // Replace old particles with the resampled particles
  particles = new_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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
