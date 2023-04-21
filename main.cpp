#include <iostream>
#include <cassert>
#include <initializer_list> //for debugging purposes
#include "Point.h"
#include <ctime>
#include <vector>
#include <algorithm>
#include <utility>
#include <chrono>
#include <ratio>
#include <fstream>
#include "Mpi_Vector.h"
using namespace std;

template<typename data_type>                         //for debugging purposes
void display(Mpi_Vector<data_type>& mpi_vector){
   for(int i = 0 ;i < mpi_vector.getsize(); i++){
      for(int j = 0; j < mpi_vector[0].dim(); j++)
         cout << mpi_vector[i][j] << ' ';
      cout << endl;
   }
}

template<typename data_type>
data_type* calculate_distances(Mpi_Vector<data_type>& mpivec, int left, int right, Point<data_type>& pivot){
   data_type* distances;
   distances = (mpivec.rank() != left) ? new data_type[mpivec.getsize()] : new data_type[(right - left)*mpivec.getsize()];
   for(int i = 0; i < mpivec.getsize(); i++)
      distances[i] = mpivec[i] - pivot;
   return distances;
}

template<typename data_type>
double sendrecv_distances_calculate_and_distribute_median(Mpi_Vector<data_type>& mpivec, int left, int right, Point<data_type>& pivot, data_type* distances){
   double median;
   data_type *distances_part;
   if(mpivec.rank() != left){
      mpivec.send(left, distances, mpivec.getsize());
      mpivec.recv(left, &median, 1, 'd');
   } 
   else{
      distances_part = new data_type[mpivec.getsize()];
      for(int i = 0; i < mpivec.getsize(); i++) distances_part[i] = distances[i];
      for(int i = left + 1; i < right; i++)
         mpivec.recv(i, distances + (i-left)*mpivec.getsize(), mpivec.getsize());
      nth_element(&distances[0], &distances[(right - left)*mpivec.getsize()/2], &distances[(right - left)*mpivec.getsize()]);
      nth_element(&distances[0], &distances[((right - left)*mpivec.getsize()-1)/2], &distances[(right - left)*mpivec.getsize()]);
      median = (distances[((right - left)*mpivec.getsize()-1)/2] + distances[(right - left)*mpivec.getsize()/2])/2.0;
      for(int i = left + 1; i < right; i++)
         mpivec.send(i, &median, 1, 'd');
      for(int i = 0; i < mpivec.getsize(); i++) distances[i] = distances_part[i];
      delete[] distances_part;  
   }
   return median;
}

template<typename data_type>
void swap_elements(Mpi_Vector<data_type>& mpivec, data_type* distances, int leftptr, int rightptr){
   Point<data_type> point_temp = mpivec[leftptr];
   mpivec[leftptr] = mpivec[rightptr];
   mpivec[rightptr] = point_temp;
   data_type dist_temp = distances[leftptr];
   distances[leftptr] = distances[rightptr];
   distances[rightptr] = dist_temp;
}

template <typename data_type>
int partition_points_and_ret_num_of_elems_less_than_median(Mpi_Vector<data_type>& mpivec, Point<data_type>& pivot, data_type* distances, double median){
   int num_of_elems_less_than_median, leftptr, rightptr;
   leftptr = -1; rightptr = mpivec.getsize();
   while(true){
      while(leftptr < mpivec.getsize()-1 && distances[++leftptr] < median);
      while(rightptr > 0 && distances[--rightptr] > median);
      if(leftptr >= rightptr) break;
      else swap_elements(mpivec, distances, leftptr, rightptr);
   }
   num_of_elems_less_than_median = leftptr;
   if(leftptr == mpivec.getsize() - 1 && distances[mpivec.getsize() - 1] < median) num_of_elems_less_than_median++;
   delete[] distances;
   return num_of_elems_less_than_median;
}

template<typename data_type>
int* build_and_distribute_exchange_table(Mpi_Vector<data_type>& mpivec, int left, int right, int num_of_elems_less_than_median){
   int* exchange_table = new int[right-left];
   exchange_table[mpivec.rank()-left] = num_of_elems_less_than_median;
   for(int i = left; i < right; i++)
      if(i != mpivec.rank()) mpivec.send(i, exchange_table + mpivec.rank() - left, 1, 'i');
   for(int i = left; i < right; i++)
      if(i != mpivec.rank()) mpivec.recv(i, exchange_table + i - left, 1, 'i');
   return exchange_table;
   // this behaves like a barrier!
}

template<typename data_type>
int* calculate_prefix_sum(Mpi_Vector<data_type>& mpivec, int left, int right, int* exchange_table){
   int num_of_points = mpivec.getsize();
   int num_of_proc = right - left;
   int* prefix_sum = new int[2*num_of_proc];
   prefix_sum[0] = exchange_table[0];
   for(int i = 1; i < 2*num_of_proc; i++) 
      if(i < num_of_proc) prefix_sum[i] = prefix_sum[i - 1] + exchange_table[i];
      else prefix_sum[i] = prefix_sum[i - 1] + num_of_points - exchange_table[i - num_of_proc];
   return prefix_sum;
}

template<typename data_type>
pair<MPI_Request**, int> send_points_less_than_median(Mpi_Vector<data_type>& mpivec, int left, int right, int* exchange_table, int* prefix_sum, vector<Point<data_type>>& keep_for_me){
   int num_of_points = mpivec.getsize(), request_num = 0;
   MPI_Request** mpireq = new MPI_Request*[exchange_table[mpivec.rank() - left]];
   for(int i = 0; i < exchange_table[mpivec.rank() - left]; i++){
      if(mpivec.rank() != left){
         if((prefix_sum[mpivec.rank() - left - 1]+i)/num_of_points == mpivec.rank() - left){
            keep_for_me.push_back(mpivec[i]);
            continue;
         }
         mpireq[request_num++] = mpivec.sendp(left + (prefix_sum[mpivec.rank() - left - 1]+i)/num_of_points, i, false, 10);
      }
      else{
         if(i/num_of_points == 0){
            keep_for_me.push_back(mpivec[i]);
            continue;
         }
         mpireq[request_num++] = mpivec.sendp(left + i/num_of_points, i, false, 10);
      }
   }
   return {mpireq, request_num};
}

template<typename data_type>
pair<MPI_Request**, int> send_points_greater_than_median(Mpi_Vector<data_type>& mpivec, int left, int right, int* exchange_table, int* prefix_sum, vector<Point<data_type>>& keep_for_me){
   int num_of_points = mpivec.getsize(), request_num = 0;
   MPI_Request** mpireq = new MPI_Request*[num_of_points - exchange_table[mpivec.rank() - left]];
   for(int i = 0; i < num_of_points - exchange_table[mpivec.rank() - left]; i++){
      if(mpivec.rank() != left){
         if((prefix_sum[mpivec.rank() + right - 2*left - 1]+i)/num_of_points == mpivec.rank() - left){
            keep_for_me.push_back(mpivec[i + exchange_table[mpivec.rank() - left]]);
            continue;
         }
         mpireq[request_num++] = mpivec.sendp(left + (prefix_sum[mpivec.rank() + right - 2*left - 1]+i)/num_of_points, exchange_table[mpivec.rank() - left] + i, false, 10);
      }
      else{
         if((prefix_sum[right - left - 1] + i)/num_of_points == 0){
            keep_for_me.push_back(mpivec[i + exchange_table[0]]);
            continue;
         }
         mpireq[request_num++] = mpivec.sendp(left + (prefix_sum[right - left - 1] + i)/num_of_points, exchange_table[0] + i, false, 10);
      }
   }
   return {mpireq, request_num};
}

void wait_buffers_to_be_ready_to_receive(pair<MPI_Request**, int>& less_reqs, pair<MPI_Request**, int>& greater_reqs){
   MPI_Status mpistat;
   for(int i = 0; i < less_reqs.second; i++)
      MPI_Wait(less_reqs.first[i], &mpistat);
   for(int i = 0; i < greater_reqs.second; i++)
      MPI_Wait(greater_reqs.first[i], &mpistat);
   delete[] less_reqs.first;
   delete[] greater_reqs.first;
}

template<typename data_type>
void exchange_elements(Mpi_Vector<data_type>& mpivec, int left, int right, int* exchange_table){
   int num_of_points = mpivec.getsize();
   int* prefix_sum = calculate_prefix_sum(mpivec, left, right, exchange_table);
   
   vector<Point<data_type>> keep_for_me;
   pair<MPI_Request**, int> less_reqs = send_points_less_than_median(mpivec, left, right, exchange_table, prefix_sum, keep_for_me);
   pair<MPI_Request**, int> greater_reqs = send_points_greater_than_median(mpivec, left, right, exchange_table, prefix_sum, keep_for_me);
   
   wait_buffers_to_be_ready_to_receive(less_reqs, greater_reqs);

   delete[] prefix_sum;
   delete[] exchange_table;

   for(int i = 0; i < num_of_points - keep_for_me.size(); i++) mpivec.recvp(-1, i, 10); // from any source!
   for(int i = 0, r = mpivec.getsize() - 1; i < keep_for_me.size(); i++, r--) mpivec[r] = keep_for_me[i];
}

template <typename data_type>
void distribute_by_median(Mpi_Vector<data_type>& mpivec, int left, int right, Point<data_type>& pivot){
   if(right - left == 1) return;
   data_type* distances = calculate_distances(mpivec, left, right, pivot);
   double median = sendrecv_distances_calculate_and_distribute_median(mpivec, left, right, pivot, distances);
   int num_of_elems_less_than_median = partition_points_and_ret_num_of_elems_less_than_median(mpivec, pivot, distances, median);
   int* exchange_table = build_and_distribute_exchange_table(mpivec, left, right, num_of_elems_less_than_median); 
   exchange_elements(mpivec, left, right, exchange_table);
   int mid = (left + right)/2;
   if(mpivec.rank() >= left && mpivec.rank() < mid)
      distribute_by_median(mpivec, left, mid, pivot);
   if(mpivec.rank() >= mid && mpivec.rank() < right)
      distribute_by_median(mpivec, mid, right, pivot);
}

template<typename data_type>
double find_min(Mpi_Vector<data_type>& mpivec, Point<data_type>& pivot){
   double mindist, dist;
   mindist = mpivec[0] - pivot;
   for(int i = 1; i < mpivec.getsize(); i++){
      dist = mpivec[i] - pivot;
      if(mindist > dist) mindist = dist;
   }
   return mindist;
}

template<typename data_type>
double find_max(Mpi_Vector<data_type>& mpivec, Point<data_type>& pivot){
   double maxdist, dist;
   maxdist = mpivec[0] - pivot;
   for(int i = 1; i < mpivec.getsize(); i++){
      dist = mpivec[i] - pivot;
      if(maxdist < dist) maxdist = dist;
   }
   return maxdist;
}

template<typename data_type>
bool check_solution(Mpi_Vector<data_type>& mpivec, int comm_size, Point<data_type>& pivot){
   double rhsmindist;
   double maxdist, mindist = find_min(mpivec, pivot);
   if(mpivec.rank()) 
      mpivec.send(mpivec.rank() - 1, &mindist, 1, 'd');
   if(mpivec.rank() != comm_size - 1){
      maxdist = find_max(mpivec, pivot);
      mpivec.recv(mpivec.rank() + 1, &rhsmindist, 1, 'd');
      if(rhsmindist < maxdist)
         return false;
   }
   return true;
}

template<typename data_type>
void sendrecv_pivot(Mpi_Vector<data_type>& mpivec, int left, int right, Point<data_type>& pivot){
   if(mpivec.rank() == left)
      for(int i = left + 1; i < right; i++)
         mpivec.send(i, pivot.dataptr(), pivot.dim());
   else mpivec.recv(left, pivot.dataptr(), pivot.dim());
}

template <typename data_type>
Point<data_type> randomly_select_pivot(Mpi_Vector<data_type>& mpivec, int comm_size){
   Point<data_type> pivot(mpivec[0].dim());
   if(!mpivec.rank())
      pivot = mpivec[rand() % mpivec.getsize()];
   sendrecv_pivot(mpivec, 0, comm_size, pivot);
   return pivot;
}

template <typename data_type>
bool partition_points(Mpi_Vector<data_type>& mpivec, int comm_size){
   Point<data_type> pivot = randomly_select_pivot(mpivec, comm_size); 
   distribute_by_median(mpivec, 0, comm_size, pivot);
   //display(mpivec); // debugging!
   return check_solution(mpivec, comm_size, pivot);
}

void execution_time(int comm_rank, int comm_size, double exec_time){
   double max_time = exec_time;
   if(comm_rank) MPI_Send(&exec_time, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
   else{
      double* exec_times = new double[comm_size];
      exec_times[0] = exec_time;
      for(int i = 1; i < comm_size; i++) MPI_Recv(&exec_times[i], 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      for(int i = 1; i < comm_size; i++) if(max_time < exec_times[i]) max_time = exec_times[i];
      cout << max_time << endl;
      delete exec_times;
   } 
}

int main(int argc, char** argv) {
   srand(time_t(time(NULL)));
   chrono::high_resolution_clock::time_point start, end;
   int comm_rank, comm_size;
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
   MPI_Comm_size(MPI_COMM_WORLD, &comm_size);  
   Mpi_Vector<double> mpivec(comm_rank, to_string(comm_rank) + ".txt");
   MPI_Barrier(MPI_COMM_WORLD); // start at the same time!
   start = chrono::high_resolution_clock::now();
   if(partition_points(mpivec, comm_size) == false) cout << "partition failed!";
   end = chrono::high_resolution_clock::now();
   double exex_time = chrono::duration_cast<chrono::duration<double>>(end - start).count();
   execution_time(comm_rank, comm_size, exex_time);
   MPI_Finalize();
}