#pragma once
#include <string>
#include <fstream>
#include <stdlib.h>
#include "Point.h"
#include <mpi.h>
using namespace std;

template <typename data_type>
class Mpi_Vector{
   public:
      Mpi_Vector(int comm_rank, const string& file_name) : comm_rank(comm_rank) {
        ifstream fin(file_name);
        if(fin.fail()) {cout << "couldn't open the file"; exit(1);}
        int dimension;
        fin >> dimension >> size;
        data = new Point<data_type>[size];
        for(int i = 0; i < size; i++){
            data[i] = Point<data_type>(dimension);
            fin >> data[i];
        }
      }
      ~Mpi_Vector() {delete[] data;}
      Point<data_type>& operator[](int i){
          assert(i >= 0 && i < size);
          return data[i];
      }
      int getsize() const {
          return size;
      }
      int rank() const {
          return comm_rank;
      }
      MPI_Request* sendp(int receiver_id, int ptr, bool destroy_req = true, int tag = 42){
            MPI_Request* mpireq = new MPI_Request;
            if(typeid(data_type).name()[0] == 'i')
                MPI_Isend(data[ptr].dataptr(), data[0].dim(), MPI_INT, receiver_id, tag, MPI_COMM_WORLD, mpireq);
            else
                MPI_Isend(data[ptr].dataptr(), data[0].dim(), MPI_DOUBLE, receiver_id, tag, MPI_COMM_WORLD, mpireq);
            if(destroy_req){
                delete mpireq;
                return nullptr;
            }
            return mpireq;
      }
      void recvp(int sender_id, int ptr, int tag = 42){
            if(typeid(data_type).name()[0] == 'i'){
                if(sender_id == -1)
                    MPI_Recv(data[ptr].dataptr(), data[0].dim(), MPI_INT, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                else
                    MPI_Recv(data[ptr].dataptr(), data[0].dim(), MPI_INT, sender_id, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
            }
            else{
                if(sender_id == -1)
                    MPI_Recv(data[ptr].dataptr(), data[0].dim(), MPI_DOUBLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                else
                    MPI_Recv(data[ptr].dataptr(), data[0].dim(), MPI_DOUBLE, sender_id, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
            }
      }
      MPI_Request* send(int receiver_id, void* addr, int size, char type = 0, bool destroy_req = true){
        MPI_Request* mpireq = new MPI_Request;
            if(type){
                if(type == 'i')
                    MPI_Isend(addr, size, MPI_INT, receiver_id, 42, MPI_COMM_WORLD, mpireq);
                else
                    MPI_Isend(addr, size, MPI_DOUBLE, receiver_id, 42, MPI_COMM_WORLD, mpireq);   // care with that!
            }
            else{
                if(typeid(data_type).name()[0] == 'i')
                    MPI_Isend(addr, size, MPI_INT, receiver_id, 42, MPI_COMM_WORLD, mpireq);
                else
                    MPI_Isend(addr, size, MPI_DOUBLE, receiver_id, 42, MPI_COMM_WORLD, mpireq);   
            }
        if(destroy_req){
            delete mpireq;
            return nullptr;
        }
        return mpireq;
    }
    void recv(int sender_id, void* addr, int size, char type = 0){
        if(type){
            if(type == 'i')
                MPI_Recv(addr, size, MPI_INT, sender_id, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
            else
                MPI_Recv(addr, size, MPI_DOUBLE, sender_id, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else{
            if(typeid(data_type).name()[0] == 'i')
                MPI_Recv(addr, size, MPI_INT, sender_id, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
            else
                MPI_Recv(addr, size, MPI_DOUBLE, sender_id, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
   private:
      int size;
      Point<data_type>* data;
      int comm_rank;
};
