// Hash-Table.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#include <iostream>
#include<fstream>
#include<string>
using namespace std;

class Node {
public:
    int64_t data;
    Node* next;
    Node() {
        next = nullptr;
    }
};

bool Check_Simple(int64_t n) {
    for (int64_t i = 2; i <= sqrt(n); i++) {
        if (n % i == 0) {
            return false;
        }
    }
    return true;
}

int64_t Find_Simple(int64_t n) {
    while (!Check_Simple(n + 1)) {
        n++;
    }
    return ++n;
}

class SL_List {
public:
    Node* head;
    SL_List(){
        head = nullptr;
    }
    void Add(int64_t elem) {
        Node* new_node = new Node;
        new_node->data = elem;
        new_node->next = head;
        head = new_node;
    }
    Node* Search(int64_t elem) {
        Node* current_node = head;
        Node* prev_node = nullptr;
        if (current_node) {
            if (head->data == elem) return head;
            while (current_node->data != elem) {
                prev_node = current_node;
                current_node = current_node->next;
                if (!current_node) return nullptr;
            }
            return prev_node;
        }
        else return nullptr;
    }
    bool Delete(int64_t elem) {
        Node* current_node = Search(elem);
        
        if (current_node) {
            if (current_node == head) {
                Node* deleted_node = head;
                head = head->next;
                delete deleted_node;
                return true;
            }
            else {
                Node* deleted_node = current_node->next;
                current_node->next = current_node->next->next;
                delete deleted_node;
                return true;
            }
        }
        else {
            return false;
        }
    }
};

class Hash_Table {
public:
    SL_List* table_array;
    int64_t size;
    Hash_Table(int n) {
        n = Find_Simple(n/3);
        table_array = new SL_List[n];
        size = n;
    }

    int H_F(int64_t elem) {
        return elem % size;
    }

    void Add(int64_t elem) {
        table_array[H_F(elem)].Add(elem);
    }

    void Delete(int64_t elem) {
        table_array[H_F(elem)].Delete(elem);
    }
    bool Search(int64_t elem) {
        if (table_array[H_F(elem)].Search(elem) != nullptr) return true;
        return false;
    }
    void print() {
        system("cls");
        for (int i = 0; i < size; i++) {
            Node* current_node = table_array[i].head;
            while (current_node) {
                cout << current_node->data << " ";
                current_node = current_node->next;
            }
            cout << endl;
        }
    }
};




int main(/*int argc, char* argv[]*/) {
    /*if (argc < 3) {
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
    */

    ifstream inFile("input5.txt");
    ofstream outFile("result.txt");

    int n;
    
    inFile >> n;
    Hash_Table HT(n);
    for (int i = 0; i < n; i++) {
        string command;
        long long int x;
        inFile >> command;

        switch (command[0])
        {
        case '+':
            inFile >> x;
            HT.Add(x);
            //HT.print();
            break;
        case '-':
            inFile >> x;
            HT.Delete(x);
            //HT.print();
            break;
        case '?':
            inFile >> x;
            if (HT.Search(x)) outFile << "true" << endl;
            else outFile << "false" << endl;
            break;
        }

    }
    inFile.close();
    outFile.close();
    return 0;
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
