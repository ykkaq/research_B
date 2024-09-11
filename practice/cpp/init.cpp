#include "init.h"

#include <iostream>
#include <string>

using namespace std;

void Student::setter(Info input_info) {  // Change this line
  info = input_info;                     // Add this line
}

Info Student::getter() { return info; }

void Student::showInformation() {
  cout << "Name: " << info.name << endl;    // Change this line
  cout << "Age: " << info.age << endl;      // Change this line
  cout << "Grade: " << info.grade << endl;  // Change this line
  cout << "Class: " << info.cls << endl;    // Change this line
}

Car::Car() : speed(0.0), migration(0.0) {
  cout << "Generate instance." << endl;
}
Car::~Car() { cout << "Delete instance." << endl; }

void Car::setSpeed(double input_speed) { speed = input_speed; }

double Car::getMigration() { return speed; }

void Car::drive(double hour) {
  cout << speed << " km/h, " << hour << " hour(s)" << endl;
  cout << speed * hour << "km moving" << endl;
  migration += speed * hour;
}
