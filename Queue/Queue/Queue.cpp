#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>

using namespace std;

const int DELETED = -1e9 - 7, NOT_PUSH = -1e9 - 9;


class binaryHeap
{
public:
    vector<pair<int, int>> array;
    //vector<int> position;
    int position[100000];
    binaryHeap() {
        array.push_back({ 0,0 });
        //position.push_back(0);
    }
    
    void siftUp(int x)
    {
        //int prev = x, i = (x -1)/ 2;
        //while (prev > 0 && array[prev] < array[i])
        //{
        //    swap(array[prev], array[i]);
        //    swap(position[array[prev].second], position[array[i].second]);
        //    prev = i;
        //    i = (i - 1) / 2;
        //}

        if (x == 1) {
            return;
        }

        if (array[x / 2].first > array[x].first) {
            position[array[x].second] = x / 2;
            position[array[x / 2].second] = x;
            swap(array[x], array[x / 2]);
            siftUp(x / 2);
        }
    }

    void siftDown(int i)
    {   
        if (i * 2 >= array.size()) return;
            
        int min_ind = i;
        if (i * 2 + 1 == array.size()) min_ind = i * 2;
        else if (array[i * 2] <= array[i * 2 + 1]) min_ind = i * 2;
        else min_ind = i * 2 + 1;
            
        if (array[i].first > array[min_ind].first) {
            position[array[i].second] = min_ind;
            position[array[min_ind].second] = i;
            swap(array[i], array[min_ind]);
            siftDown(min_ind);
        }
    }

   

    void push(pair<int, int> x)
    {
        array.push_back(x);
        //position.push_back(array.size() - 1);
        position[x.second] = array.size() - 1;
        siftUp(array.size() - 1);
    }

    pair<int, int> extractMin()
    {
        if (array.size()>1) {
            pair<int, int> ans = array[1];
            position[array.back().second] = 1;
            array[1] = array.back();
            array.pop_back();
            position[ans.second] = -1;
            siftDown(1);
            return ans;
        }
        else {
            return {1e9 + 1,0};
        }
    }

    void decreaseKey(int ind, int x)
    {
       /*if (position[ind] < 0 || array[position[ind]].first < x)
        {
            return;
        }*/
        array[position[ind]].first = x;
        siftUp(position[ind]);

        /*array[position[ind]].first = x;
        siftUp(position[ind]);*/
    }
};


int main(/*int argc, char* argv[]*/) {
    /*if (argc<3) {
      return -1;
    }
    std::ifstream inFile(argv[1]);
    if (!inFile){
      return -2;
    }
    std::ofstream outFile(argv[2]);
    if (!outFile){
      return -3;
    }
    */
    ifstream inFile("C:/Users/User/Documents/GitHub/task9-priorityqueue-tottttttt/tests/in/input6.txt");
    ofstream outFile("C:/Users/User/Desktop/result.txt");
    
    int count = 1;
    binaryHeap Heap;
    while (inFile) {
        string command;
        int x;
        inFile >> command;
        
        if (command == "push") {
            inFile >> x;
            pair<int, int> X(x, count);
            Heap.push(X);
           
        }
        else if (command == "extract-min") {
            pair<int,int> result = Heap.extractMin();
            if (result.first != 1e9+1) cout << result.first << endl;
            else cout << '*' << '\n';
            
        }
        else if (command == "decrease-key") {
            int index;
            inFile >> index;
            inFile >> x;
           
            Heap.decreaseKey(index,x);
        }
        count++;
    }

    inFile.close();
    outFile.close();
    return 0;
}

