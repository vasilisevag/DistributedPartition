#pragma once

#include <iostream>
#include <initializer_list>
#include <cassert>
using namespace std;

template <typename data_type>
class Point{
   public:
      Point() : dimension(0), data(nullptr) {}
      Point(unsigned int dimension) : dimension(dimension), data(new data_type[dimension]) {}
      Point(const initializer_list<data_type>& values) : dimension(values.size()), data(new data_type[dimension]){
         copy(begin(values), end(values), data);
      }
      Point(const Point<data_type>& point) : dimension(point.dimension), data(new data_type[dimension]) {
         for(int i = 0; i < dimension; i++)
            data[i] = point.data[i];
      }
      ~Point(){delete[] data;}
      Point<data_type>& operator=(const Point<data_type>& point){
         if(dimension != point.dimension){
            delete[] data;
            dimension = point.dimension;
            data = new data_type[dimension];   
         }
         for(int i = 0; i < dimension; i++)
            data[i] = point.data[i];
         return *this;
      }
      data_type& operator[](int i){
         assert(i >= 0 && i < dimension);
         return data[i];
      }
      const data_type& operator[](int i) const {
         assert(i >= 0 && i < dimension);
         return data[i];
      }
      int dim() const {
         return dimension;
      }
      data_type* dataptr() {
         return data;
      } 
   private:
      unsigned int dimension;
      data_type* data;
};

template <typename data_type>
data_type operator-(const Point<data_type>& lhs_point, const Point<data_type>& rhs_point){
   data_type square_diff = 0;
   for(int i = 0; i < lhs_point.dim(); i++)
      square_diff += (lhs_point[i] - rhs_point[i])*(lhs_point[i] - rhs_point[i]);
   return square_diff;
}

template <typename data_type>
istream& operator>>(istream& in, Point<data_type>& point){
   for(int i = 0; i < point.dim(); i++)
      in >> point[i];
   return in;
} 