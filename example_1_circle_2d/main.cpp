#include <iostream>
#include <vector>
#include <cmath>
using namespace std;
/*
This program implements the Gauss Newton Algorithm to solve a specific problem

problem statement:
given a circle in the 2D plane,
find the point of the circle that is closer to the origin of the axes.

the algorithm takes a given point on the circle and iterates, each iteration gives a point on
that is closer to the optimal minimum point.
 */

/*
 * structure to represents a circle in the 2D plane.
 */
struct Circle {
    double xc;
    double yc;
    double radius;
};
/*
 * Return a specific point on a circle.
 *
 * \param[in] alpha [radiants].
 * \param[in] circle [Circle struct].
 * \param[out]: point in the 2D belonging to the circle.
 */
vector<double> r_fcn (double alpha, Circle circle) {
    vector<double> r = {0,0};
    double xc = circle.xc;
    double yc = circle.yc;
    double rho = circle.radius;
    r[0] = xc+rho*cos(alpha);
    r[1] = yc+rho*sin(alpha);
    return r;
}

/*
 * Return the Jacobian of the circle function at a given point.
 *
 * \param[in] alpha [radiants].
 * \param[in] circle [Circle struct].
 * output: jacobian of the circle function at the specific alpha
 */
vector<double> J_fcn (double alpha, Circle circle) {
    vector<double> J = {0,0};
    double rho = circle.radius;
    J[0] = -rho*sin(alpha);
    J[1] = rho*cos(alpha);
    return J;
}


double JTrk_fcn(vector<double> JT, vector<double> rk) {
    return JT[0]*rk[0]+JT[1]*rk[1];
}

int main() {
    cout << "This program performs minimization of the circle coordinate distance from the origin" << endl;
    double xc = 1;
    double yc = 1;
    double rho = 1;
    // a Circle structure is defined with the desired circle parameters
    // in this way, the structure is given as input to the circle and jacobian functions.
    Circle circle;
    circle.xc = xc;
    circle.yc = yc;
    circle.radius = rho;


    cout << "The circle has center in (" << xc << ", " << yc << ") and radius " << rho << "." << endl;
    cout << endl;
    double JTJ = pow(rho ,2 );
    double JTJ1 = pow(JTJ ,-1 );

    double alpha_0 = 1.5*M_PI;
    double alpha_k = alpha_0;

    cout << "The input for the algorithm is the point at " << 180*alpha_0/M_PI << " degrees." << endl;
    cout << endl;

    int N = 9;

    vector<double> alpha_vec;

    for (int i = 0; i <N; i++)
    {
        vector<double> rk = r_fcn(alpha_k, circle); // compute the point on the circle
        vector<double> J = J_fcn(alpha_k, circle); // compute the Jacobian of the circle function
        double JTrk = JTrk_fcn(J,rk); // perform vector by vector multiplication
        double h = -JTJ1*JTrk; // perform scalar by vector multiplication
        alpha_k += h; // update alpha
        cout << "At iteration " << i+1 << " alpha is " << 180*alpha_k/M_PI << " degrees." << endl; // print updated alpha

        alpha_vec.push_back(alpha_k); // add alpha to the iteration history vector
    }

    cout << "Final iteration is " << 180*alpha_k/M_PI << " degrees." << endl;
    /*
    cose da implementare
    - scrivere su file i risultati
    - scrivere su foglio di testo i risultati
    - plottare i risultati richiamandoli da un file
     */


    return 0;
}



