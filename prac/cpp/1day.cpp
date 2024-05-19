#include <iostream>
#include <string>

#include "init.cpp"

using namespace std;

void input(Info* student_info);

int main(int argc, char** argv) {
  int num = 2;
  Student student[num];

  for (int i = 0; i < num; i++) {
    Info student_info;  // Add this line

    input(&student_info);

    student[i].setter(student_info);  // Change this line
  }

  for (int i = 0; i < num; i++) {
    student[i].showInformation();
  }

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