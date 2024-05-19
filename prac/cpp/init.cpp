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