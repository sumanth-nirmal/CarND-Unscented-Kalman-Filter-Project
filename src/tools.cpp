#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

    VectorXd rmse = VectorXd(4);
    rmse << 0, 0, 0, 0;

    // Check the validity of the following inputs:
    // The estimation vector size should not be zero
    if(estimations.size() == 0){
        cout << "Input is empty" << endl;
    }
        // The estimation vector size should equal ground truth vector size
    else if(estimations.size() != ground_truth.size()){
        cout << "Invalid estimation or ground_truth. Data should have the same size" << endl;
    }
    else {
        //accumulate squared residuals
        for (unsigned int i = 0; i < estimations.size(); ++i) {

            VectorXd residual = estimations[i] - ground_truth[i];

            //coefficient-wise multiplication
            residual = residual.array() * residual.array();
            rmse += residual;
        }

        //calculate the mean
        rmse = rmse / estimations.size();

        //calculate the squared root
        rmse = rmse.array().sqrt();
    }
    //return the result
    return rmse;
}