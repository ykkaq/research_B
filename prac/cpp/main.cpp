#include <iostream>
#include <string>

#include "init.cpp"

using namespace std;

void input(Info* student_info);

int main(int argc, char** argv) {
  Car car;
  car.setSpeed(60);
  car.drive(1.5);
  car.setSpeed(100);
  car.drive(2.5);
  cout << "sum: " << car.getMigration() << " km" << endl;

  return 0;
}

void input(Info* student_info) {
  cout << "Enter name: ";
  cin >> student_info->name;
  cout << "Enter age: ";
  cin >> student_info->age;
  cout << "Enter grade: ";
  cin >> student_info->grade;
  cout << "Enter class: ";
  cin >> student_info->cls;
}