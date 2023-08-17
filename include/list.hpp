#ifndef LIST_HPP
#define LIST_HPP

// Libraries
#include "basic_lib.hpp"
#include "utilities.hpp"

namespace basic_lib{
// Definitions

template <typename T>
struct Linked_List{
    public:
    T data;
    Linked_List* next=NULL;
};

template <typename T>
class List{
private:
    int N_elem=0;
    Linked_List<T> *head=NULL;
    T get_elem (const int index) const {
        check_error(index>=this->N_elem||index<0, "list index out of range!");
        Linked_List<T> *next=this->head->next;
        for (int i=0; i<index; i++){
            next = next->next;
        }
        return next->data;
    }
public:
    List(){
        this->head = (Linked_List<T>*) calloc(1, sizeof(Linked_List<T>));
		assert(this->head!=NULL);
    }
    ~List(){
        clear();
        free(this->head);
    }
    void clear(){
        const int N=this->N_elem;
        for (int i=0; i<N; i++){
            remove(0);
        }
    }
    void append(const T elem){
        Linked_List<T> *next=this->head;
        while(next->next!=NULL){
            next = next->next;
        }
        next->next = (Linked_List<T>*) calloc(1, sizeof(Linked_List<T>));
        assert(next->next!=NULL);
        next->next->data = elem;
        this->N_elem++;
    }
    // for debugging only!
    // void print(){
    //     printf("List of %d elements:\n", this->N_elem);
    //     if (this->head!=NULL){
    //         Linked_List<T> *next=this->head;
    //         while(next->next!=NULL){
    //             next = next->next;
    //             std::cout << next->data << " at " << next << " next pointer is " << next->next << std::endl;
    //         }
    //     }
    // }
    void remove (const int index){
        check_error(index>=this->N_elem||index<0, "list index out of range!");
        Linked_List<T> *next=this->head->next;
        Linked_List<T> *prev=this->head;
        for (int i=0; i<index; i++){
            next = next->next;
            prev = prev->next;
        }
        prev->next = next->next;
        free(next);
        this->N_elem--;
    }
    T operator() (const int index) const {
        check_error(index>=this->N_elem||index<0, "list index out of range!");
        Linked_List<T> *next=this->head->next;
        for (int i=0; i<index; i++){
            next = next->next;
        }
        return next->data;
    }
    int len(){
        return this->N_elem;
    }
    void sort(int (*compar)(const void*, const void*)){
        T *data=(T*)calloc(this->N_elem, sizeof(T));
		assert(data!=NULL);
        for (int i=0; i<this->N_elem; i++){
            data[i] = this->get_elem(i);
        }
        qsort(data, this->N_elem, sizeof(T), compar);
        const int N=this->N_elem;
        this->clear();
        for (int i=0; i<N; i++){
            this->append(data[i]);
        }
        free(data);
    }
};

// Functions

}

#endif