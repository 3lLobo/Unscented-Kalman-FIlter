#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.2;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  // initalize state dimension and augmented state estimation
  n_x_ = 5;
  n_aug_ = 7
  lambda_ = 3 - n_aug_
  is_initialized_ = false;
  
  weights_ = VectorXd(2*n_aug_+1);
  weights_(0) = lambda_/(lambda_+n_aug_);
  for (int i=1; i<2*n_aug_+1; i++) { 
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
 
  if (!is_initialized_){
	time_us_ = meas_package.timestamp_;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
	  float px = rho * cos(phi)
	  float py = rho * sin(phi)
      x_ << px, py, 0, 0, 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
      P_(0,0) = 0.6;
      P_(1,1) = 0.6;
    }
    is_initialized_ = true;
    return; 
  }
  double dt = (double)(meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
  time_us_ = meas_package.timestamp_;
  Prediction(dt);
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
  	UpdateRadar(meas_package);
  }
  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
  	UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  // aug mean state and aug covariance matrix
  VectorXd x_aug = VectorXd(7);
  MatrixXd P_aug = MatrixXd(7, 7);
  x_aug(5) = 0;
  x_aug(6) = 0;
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;  
  MatrixXd L = P_aug.llt().matrixL();
  MatrixXd Xsig_aug = MatrixXd(7, 15);
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++){
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+7) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  for (int i = 0; i< 2*n_aug_+1; i++){
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);
    double px_p, py_p;
    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v/yawd * ( sin (yaw_p) - sin(yaw));
      py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw_p) );
    }
    else {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }
	//add noise
    double dt2 = delta_t * delta_t;
    px_p += 0.5*nu_a*dt2 * cos(yaw);
    py_p += 0.5*nu_a*dt2 * sin(yaw);
    v_p  += nu_a*delta_t;

    yaw_p  += 0.5*nu_yawdd*dt2;
    yawd_p += nu_yawdd*delta_t;

	// write into prediction vector
	Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p; 
    Xsig_pred_(3,i) = yaw_p; 
    Xsig_pred_(4,i) = yawd_p;
	
	VectorXd x_pred_ = Xsig_pred_ * weights_ ;
	MatrixXd P_pred_ = MatrixXd(5, 5);
    P_pred_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++){
		
	    VectorXd x_diff = Xsig_pred_.col(i) - x_pred_;
		// normalization
		while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
		P_pred_ = P_pred_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
  P_ = P_pred_;
  x_ = x_pred_;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  	int n_z = 2;  
	VectorXd z = meas_package.raw_measurements_;
	// sigma point matrix
	MatrixXd Z_sig = Xsig_pred_.block(0,0, n_z, 2 * n_aug_ + 1);
	//Covariance matrix
  	MatrixXd S = MatrixXd(n_z,n_z);	
	S.fill(0.0);
  	VectorXd z_pred_ = VectorXd(n_z);
	z_pred_.fill(0.0);
	for (int i=0; i < 2*n_aug_+1; i++) {
		z_pred_ = z_pred_ + weights_(i) * Z_sig.col(i);
	}
  	for (int i=0; i < n_sig_; i++){
  		//residual 
      	VectorXd z_diff = Zsig.col(i) - z_pred_;
      
      	// normalization between -pi and pi
      	while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
      	while (z_diff(1) <-M_PI) z_diff(1) += 2.*M_PI;
      
      	S = S + weights_(i) * z_diff * z_diff.transpose();
  	}
	//cross correlation Matrix
	MatrixXd T_c = MatrixXd(n_x_, n_z);
	T_c.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) { 
		VectorXd z_diff = Z_sig.col(i) - z_pred_;
			while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
			while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

		VectorXd x_diff = Xsig_pred_.col(i) - x_;
	    	while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
			while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
		T_c = T_c + weights_(i) * x_diff * z_diff.transpose();
	}
	
	S = S + R_lidar_;
  //Kalman gain
  MatrixXd K = T_c * S.inverse();
  VectorXd z_diff = z - z_pred_;
  // Update step
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  // calculate lidar NIS
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;  
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  // dimension for radar measurement
  int n_z = 3;
  MatrixXd Z_sig = MatrixXd(n_z, 2 * n_aug_ + 1);
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_pred_ = VectorXd(n_z);
  z_pred_.fill(0.0);
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  MatrixXd T_c = MatrixXd(n_x_, n_z);
  T_c.fill(0.0);
  //transformation into measrement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  

    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v   = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double help_1 = cos(yaw)*v;
    double help_2 = sin(yaw)*v;
    double radius_ = sqrt(p_x*p_x + p_y*p_y);
    Z_sig(0,i) = radius_;                   
    Z_sig(1,i) = atan2(p_y,p_x);          
    Z_sig(2,i) = (p_x*help_1 + p_y*help_2)/radius_;      
  }

  for (int i=0; i < 2*n_aug_+1; i++) {
    z_pred_ = z_pred_ + weights_(i) * Z_sig.col(i);
  }  
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Z_sig.col(i) - z_pred_;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();

  }
  S = S + R_radar_;
  //Cross correlation matrix

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    VectorXd z_diff = Z_sig.col(i) - z_pred_;
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    VectorXd x_diff = X_sig_pred_.col(i) - x_;
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    T_c = T_c + weights_(i) * x_diff * z_diff.transpose();
  }
  // get measurements
  VectorXd z= meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred_;
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  
  MatrixXd K = T_c * S.inverse();
  
  //Radar update step
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  NIS_radar_ = z.transpose() * S.inverse() * z;


}
