// lets define a functor
#ifndef __FUNCTOR_HPP
#define __FUNCTOR_HPP

#include <iostream>

// A Functor
#if defined(DEFAULT)

class __attribute__((visibility("default"))) Increment
{
private:
    int num;
public:
    Increment() = default;
    Increment(int n) : num(n) {  }
  
    // This operator overloading enables calling
    // operator function () on objects of increment
    int operator () (int arr_num) const {
        return num + arr_num;
    }
    void visibility_state() {std::cout<<" default visibility "<<std::endl;}
};
extern __attribute__((visibility("default"))) Increment myincrementor;

#elif defined(NODEFAULTVISIBILITY)

// Was a Functor but now just class with default constructor
class Increment
{
private:
    int num;
public:
    Increment() = default;
    Increment(int n) : num(n) {  }
  
    // This operator overloading enables calling
    // operator function () on objects of increment
    int operator () (int arr_num) const {
        return num + arr_num;
    }
    void visibility_state() {std::cout<<" no visibility attribute set "<<std::endl;}
};
extern Increment myincrementor;

#else 
class __attribute__((visibility("hidden"))) Increment
{
private:
    int num;
public:
    Increment() = default;
    Increment(int n) : num(n) {  }
  
    // This operator overloading enables calling
    // operator function () on objects of increment
    int operator () (int arr_num) const {
        return num + arr_num;
    }
    void visibility_state() {std::cout<<" hidden visibility "<<std::endl;}
};
extern __attribute__((visibility("hidden"))) Increment myincrementor;
#endif

#endif
