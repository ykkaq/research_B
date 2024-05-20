#ifndef _INIT_H_
#define _INIT_H_

#include <string>  // Add this line

using namespace std;

struct Info {
  string name;
  int grade;
  string cls;
  int age;
};

class Student {
 public:
  void showInformation();
  void setter(Info input_info);  // Change this line
  Info getter();

 private:
  Info info;
};

class Car {
 public:
  Car();
  ~Car();
  void setSpeed(double speed);
  double getSpeed();
  double getMigration();
  void drive(double hour);

 private:
  double speed;
  double migration;
};

#endif  // _INIT_H_