#pragma once
#include <cmath>
#ifndef M_PI
namespace
{
	const double M_PI = acos(-1.0);
}
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>

using namespace std;

//Data structures
vector<double> x;
vector<double> uold;
vector<double> unew;

//Variables that we'll need to use
int N;
double t, tfinal, dt, dx;
double a, nu; //advection and diffusion parameters;

//For writing data to disk
void writeFields() {

	ofstream write;

	write.open("numVsExact.dat");

	write << "TITLE = \"Exact solution vs numerical solution\"" << endl;
	write << "VARIABLES = \"x\", \"exact solution\", \"numerical solution\"" << endl;
	write << "ZONE T=\"Only Zone\", I=" << N << ", F=POINT" << endl;

	for (int i = 0; i < N; i++){
		write << x[i] << " " << exp(-nu*t)*sin(x[i] - a*t) << " " << unew[i] << endl;
	}

	write.close();

}

//To calculate conserved quantity
double calcConserve() {

	//Calculate conservation
	double In = 0.0;

	for (int i = 0; i < N; i++) In += dx*unew[i];

	cout << "I( t = " << t << " ) = " << In << endl;

	return In;

}

//Write a function that will update unew and uold at each timestep
void update() {

	double dxi = 1.0 / dx; //Inverse of dx

	//First update the stencil

	//Interior points
	for (int i = 1; i < N - 1; i++) {

		double Fi_p = a*uold[i] - nu*(uold[i + 1] - uold[i]) / dx;
		double Fi_m = a*uold[i - 1] - nu*(uold[i] - uold[i - 1]) / dx;

		//V*ddtu = -(F_+ - F_-)
		unew[i] = uold[i] - (dt / dx)*(Fi_p - Fi_m);

	}
	//Left-most cell
  {
	  double Fi_p = a*uold[0] - nu*(uold[1] - uold[0]) / dx;
	  double Fi_m = a*uold[N - 1] - nu*(uold[0] - uold[N - 1]) / dx;

	  //V*ddtu = -(F_+ - F_-)
	  unew[0] = uold[0] - (dt / dx)*(Fi_p - Fi_m);
  }
  //Right-most cell
  {
	  double Fi_p = a*uold[N - 1] - nu*(uold[0] - uold[N - 1]) / dx;
	  double Fi_m = a*uold[N - 2] - nu*(uold[N - 1] - uold[N - 2]) / dx;

	  //V*ddtu = -(F_+ - F_-)
	  unew[N - 1] = uold[N - 1] - (dt / dx)*(Fi_p - Fi_m);
  }

  //Now kick everything up one timestep
  for (int i = 0; i < N; i++){
	  uold[i] = unew[i];
  }

}

int main() {

	cout << "Beginning finite volume code for 1D advection-diffusion equation with periodic BCs." << endl;

	cout << "Pre-processing: Initialize parameters and data structures." << endl;
	cout << "Input N:";
	cin >> N;

	a = 1.0;
	nu = 1.0;

	dx = 2.0*M_PI / double(N);
	dt = __min(0.5*dx / a, 0.25*dx*dx / nu);
	tfinal = 0.5;

	//Initialize data structures
	x.resize(N);
	uold.resize(N);
	unew.resize(N);
	for (int i = 0; i < N; i++){
		x[i] = (0.5 + double(i))*dx;
		uold[i] = sin(x[i]);
		unew[i] = sin(x[i]);
	}

	cout << "Initial integral:" << endl;
	calcConserve();

	cout << "Beginning timestepping." << endl;
	t = 0;
	while (t < tfinal){
		update();
		t += dt;
	}
	cout << "Final integral:" << endl;
	calcConserve();

	cout << "Timestepping complete - outputing final state for post-processing." << endl;
	writeFields();

	cout << "Compute error compared to exact solution." << endl;
	double error = 0.0;
	for (int i = 0; i < N; i++){
		error += dx*abs(unew[i] - exp(-nu*t)*sin(x[i] - a*t));
	}
	cout << "(dt,dx,error) = (" << dt << "," << dx << "," << error << ")" << endl;

	cout << "Finished." << endl;

	return 0;
}

