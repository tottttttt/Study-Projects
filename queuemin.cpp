#include <fstream>
#include <iostream>


using namespace std;

class Stack
{
private:
    long long size = 0;
    long long top_point = -1;
    long long* s_array = nullptr;
    long long* s_array_min = nullptr;

    void Add_Memory(int size) {

        long long* temp_s_array = new long long[size];
        long long* temp_s_array_min = new long long[size];

        for (long long i = 0; i < top_point + 1; i++) {
            temp_s_array[i] = s_array[i];
            temp_s_array_min[i] = s_array_min[i];
        }
        delete[] s_array;
        delete[] s_array_min;

        s_array = temp_s_array;
        s_array_min = temp_s_array_min;
    }

    long long Add_MIN(long long a) {
        if (a < s_array_min[top_point - 1]) {
            return a;
        }
        else {
            return s_array_min[top_point - 1];
        }
    }

public:


    void push(long long elem) {

        if (!size) {
            s_array = new long long[1];
            s_array_min = new long long[1];
            top_point++;
            size++;
            s_array[top_point] = elem;
            s_array_min[top_point] = elem;
        }
        else if (size - top_point - 1) {
            top_point++;
            s_array[top_point] = elem;
            s_array_min[top_point] = Add_MIN(elem);
        }
        else {

            size *= 2;
            Add_Memory(size);

            top_point++;
            s_array[top_point] = elem;
            s_array_min[top_point] = Add_MIN(elem);

        }
    }



    long long pop() {
        top_point--;
        return s_array[top_point + 1];
    }

    long long get_min() {
        return s_array_min[top_point];
    }

    bool empty() {
        if (top_point + 1) {
            return false;
        }
        return true;
    }

    ~Stack() {
        delete[] s_array;
        delete[] s_array_min;
    }
};


class Queue {
private:
    Stack Push_Stack;
    Stack Pop_Stack;
public:

    void push(long long elem) {
        Push_Stack.push(elem);
    }

    long long pop() {
        if (Pop_Stack.empty()) {
            while (!Push_Stack.empty()) {
                Pop_Stack.push(Push_Stack.pop());
            }
        }
        return Pop_Stack.pop();
    }

    long long get_min() {
        if (Pop_Stack.empty()) return Push_Stack.get_min();
        if (Push_Stack.empty()) return Pop_Stack.get_min();
        return min(Push_Stack.get_min(), Pop_Stack.get_min());
    }

};


int main(int argc, char* argv[]) {
    if (argc < 3) {
        return 1;
    }
    ifstream inFile(argv[1]);
    if (!inFile) {
        return 2;
    }
    ofstream outFile(argv[2]);
    if (!outFile) {
        return 3;
    }


    Queue my_queue;
    int n;
    inFile >> n;
    string task;
    long long elem;

    for (int i = 0; i < n; i++) {
        inFile >> task;
        switch (task[0])
        {
        case '+':
            inFile >> elem;
            my_queue.push(elem);
            break;
        case '-':
            my_queue.pop();
            break;
        case '?':
            outFile << my_queue.get_min() << endl;
            break;
        }

    }

    inFile.close();
    outFile.close();
    return 0;
}