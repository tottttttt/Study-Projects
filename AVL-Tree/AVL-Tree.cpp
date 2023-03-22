// AVL-Tree.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#include <iostream>
#include<fstream>
#include<string>
using namespace std;

class Node {
public:
    long long int data;
    Node* left, * right, *parent;
    
    int d_balance;
    Node(long long int x);
};
Node::Node(long long int x) {
    data = x;
    left = nullptr;
    parent = nullptr;
    right = nullptr;
    d_balance = 0;
}

class AVL_Tree {
public:
    Node* root;
    AVL_Tree() {
        root = nullptr;
    }

    void Add(long long int elem) {
        Node* current_node = root;
        Node* parent_node = nullptr;
        bool Left_Right;

        while (current_node != nullptr) {
            parent_node = current_node;
            if (current_node->data < elem) {
                current_node = current_node->right;
                Left_Right = false;
            }
            else {
                current_node = current_node->left;
                Left_Right = true;
            }
        }

        current_node = new Node(elem);
        current_node->parent = parent_node;
        if (parent_node == nullptr) {
            root = current_node;
        }
        else {
            if (Left_Right) {
                parent_node->left = current_node;
                
                Balance(parent_node->left);
            }
            else {
                parent_node->right = current_node;
                
                Balance(parent_node->right);

            }
        }



    };

    

   

    void Left_Rotation(Node* node) {
        Node* a = node->parent;
        Node* b = node;
        Node* c = node->left;

        
        node->parent->right = node->left;
        if (node->left != nullptr) node->left->parent = node->parent; 
        node->left = node->parent;              
        if (node->left != root) {
            node->parent = node->left->parent;
            if (node->parent->left == node->left) node->parent->left = node;
            else node->parent->right = node;
        }
        else {
            root = node;
            node->parent = nullptr;
        }
        node->left->parent = node;

    }

    void Right_Rotation(Node* node) {
       
        node->parent->left = node->right;
        if (node->right != nullptr) node->right->parent = node->parent;
        node->right = node->parent;
        
        
        if (node->right != root) {
            node->parent = node->right->parent;
            
            if (node->parent->left == node->right) node->parent->left = node;
            else node->parent->right = node;
        }
        else {
            root = node;
            node->parent = nullptr;
        }
        node->right->parent = node;
        
    }

    void Double_Left_Rotation(Node* node) {
        Right_Rotation(node);
        Left_Rotation(node);
    }

    void Double_Right_Rotation(Node* node) {
        Left_Rotation(node);
        Right_Rotation(node);
    }

    void Balance(Node* node) {
        Node* parent_node = node->parent;
        Node* current_node = node;
        bool left_right; // left = true, right = false
        do {
             
            if (current_node != root) {
                if (parent_node->left == current_node) parent_node->d_balance--;
                else parent_node->d_balance++;
            }
            
            current_node = parent_node;
            parent_node = current_node->parent;

            switch (current_node->d_balance)
            {
            case 2:
                switch (current_node->right->d_balance)
                {
                case 1:
                    Left_Rotation(current_node->right);
                    //current_node = current_node->parent;
                    /*if(current_node != root) current_node->parent->d_balance--;*/ //высота поддерева уменьшилась на 1, значит надо об этом сообщить родителю
                    current_node->parent->left->d_balance = 0;
                    current_node->parent->d_balance = 0;
                    
                    break;
                case -1:
                    switch (current_node->right->left->d_balance)
                    {
                    case 1:
                        Double_Left_Rotation(current_node->right->left);
                        //current_node = current_node->parent;
                        current_node->parent->left->d_balance = -1;
                        current_node->parent->right->d_balance = 0;
                        current_node->parent->d_balance = 0;
                        /*if (current_node != root) current_node->parent->d_balance--;*/

                        
                        break;
                    case 0:
                        Double_Left_Rotation(current_node->right->left);
                        //current_node = current_node->parent;
                        current_node->parent->left->d_balance = 0;
                        current_node->parent->right->d_balance = 0;
                        current_node->parent->d_balance = 0;
                        /*if (current_node != root) current_node->parent->d_balance--;*/

                        
                        break;
                    case -1:
                        Double_Left_Rotation(current_node->right->left);
                        //current_node = current_node->parent;
                        current_node->parent->left->d_balance = 0;
                        current_node->parent->right->d_balance = 1;
                        current_node->parent->d_balance = 0;
                        /*if (current_node != root) current_node->parent->d_balance--;*/

                        
                        break;
                    }
                    break;
                case 0:
                    Left_Rotation(current_node->right);
                    //current_node = current_node->parent;
                    current_node->parent->left->d_balance = 1;
                    current_node->parent->d_balance = -1;
                    
                    break;
                }
                current_node = current_node->parent;
                break;
            case -2:
                switch (current_node->left->d_balance)
                {
                case 1:
                    switch (current_node->left->right->d_balance)
                    {
                    case -1:
                        Double_Right_Rotation(current_node->left->right);
                        
                        current_node->parent->left->d_balance = 0;
                        current_node->parent->right->d_balance = 1;
                        current_node->parent->d_balance = 0;
                       /* if (current_node != root) current_node->parent->d_balance--;*/ // высота поддерева уменьшилась на 1 

                        
                        break;
                    case 0:
                        Double_Right_Rotation(current_node->left->right);
                        //current_node = current_node->parent;
                        current_node->parent->left->d_balance = 0;
                        current_node->parent->right->d_balance = 0;
                        current_node->parent->d_balance = 0;
                        /*if (current_node != root) current_node->parent->d_balance--;*/

                        break;
                    case 1:
                        Double_Right_Rotation(current_node->left->right);
                        //current_node = current_node->parent;
                        current_node->parent->left->d_balance = -1;
                        current_node->parent->right->d_balance = 0;
                        current_node->parent->d_balance = 0;
                       /* if (current_node != root) current_node->parent->d_balance--;*/
                        break;
                    }

                    break;
                case 0:
                    Right_Rotation(current_node->left);
                    //current_node = current_node->parent;
                    current_node->parent->right->d_balance = -1;
                    current_node->parent->d_balance = 1;
                    break;      

                case -1:
                    Right_Rotation(current_node->left);
                    //current_node = current_node->parent;
                    /*if (current_node != root) current_node->parent->d_balance--;*/
                    current_node->parent->right->d_balance = 0;
                    current_node->parent->d_balance = 0;
                    break;
                }
                current_node = current_node->parent;
                break;
            default:
                break;
            }
        }while (current_node != root && current_node->d_balance != 0);
    }

    void Delete(Node* node) {
        Node* current_node;
        bool left_right;
        bool flag = false;
        if (node->right && node->left) {
            Node* alt = node->right;
            
            while (alt->left) {
                alt = alt->left;
            }
            
            if (alt->parent->left == alt) { left_right = true; /*alt->parent->d_balance++;*/ }
            else left_right = false;
            
            alt->d_balance = node->d_balance;

            if (alt->parent == node) {
                alt->parent = node->parent;
                alt->left = node->left;
                
                if (node->left) node->left->parent = alt;
                current_node = alt;
            }
            else{
                current_node = alt->parent;
                alt->parent->left = alt->right;
                if (alt->right) alt->right->parent = alt->parent;
                alt->left = node->left;
                alt->right = node->right;
                alt->parent = node->parent; 
                if (node->left) node->left->parent = alt;
                if (node->right) node->right->parent = alt;
            }
            
            if (node!= root)
            if (alt->parent->left == node) alt->parent->left = alt;
            else alt->parent->right = alt;      
            else {
                root = alt;
            }
            
            delete node;
        }
        else {
            if (node->left) {
                if (node != root){if (node->parent->left == node) {
                    node->parent->left = node->left; left_right = true;
                }
                else {
                    node->parent->right = node->left; left_right = false;
                }
                node->left->parent = node->parent;
                current_node = node->parent;
                delete node;
                }
                else {
                    root = node->left;
                    root->parent = nullptr;
                    root->d_balance = 0;
                    delete node;
                    return;
                }
            }
            else if (node->right) {
                if(node != root){
                if (node->parent->left == node) {
                    node->parent->left = node->right; left_right = true;
                }
                else {
                    node->parent->right = node->right; left_right = false;
                }

                node->right->parent = node->parent;
                current_node = node->parent;
                delete node;
                }
                else {
                    root = node->right;
                    root->parent = nullptr;
                    root->d_balance = 0;
                    delete node;
                    return;
                }

            }
            else {
                if(node != root){
                    if (node->parent->left == node) {
                    node->parent->left = nullptr; left_right = true;
                }
                else { node->parent->right = nullptr; left_right = false; }
                current_node = node->parent;

                delete node;
                }
                else {
                    root = nullptr;
                    return;
                }
            }
        }
         if (left_right) current_node->d_balance++;
             else current_node->d_balance--;
         
         do{
             
            if (current_node != root && flag) {
                 if (current_node->parent->left == current_node) { /*left_right = true;*/ current_node->parent->d_balance++;
                 }
                 else /*left_right = false;*/ current_node->parent->d_balance--;
                 current_node = current_node->parent;
            }
            
            switch (current_node->d_balance)
            {
            case 2:
                switch (current_node->right->d_balance)
                {
                case 1:
                    Left_Rotation(current_node->right);
                    current_node->parent->left->d_balance = 0;
                    current_node->parent->d_balance = 0;
                    break;
                case -1:
                    switch (current_node->right->left->d_balance)
                    {
                    case 1:
                        Double_Left_Rotation(current_node->right->left);
                        current_node->parent->left->d_balance = -1;
                        current_node->parent->right->d_balance = 0;
                        current_node->parent->d_balance = 0;
                        break;
                    case 0:
                        Double_Left_Rotation(current_node->right->left);
                        current_node->parent->left->d_balance = 0;
                        current_node->parent->right->d_balance = 0;
                        current_node->parent->d_balance = 0;
                        break;
                    case -1:
                        Double_Left_Rotation(current_node->right->left);
                        current_node->parent->left->d_balance = 0;
                        current_node->parent->right->d_balance = 1;
                        current_node->parent->d_balance = 0;
                        break;
                    }
                    break;
                case 0:
                    Left_Rotation(current_node->right);
                    current_node->parent->left->d_balance = 1;
                    current_node->parent->d_balance = -1;
                    break;
                }
                current_node = current_node->parent;
                break;
            case -2:
                switch (current_node->left->d_balance)
                {
                case 1:
                    switch (current_node->left->right->d_balance)
                    {
                    case -1:
                        Double_Right_Rotation(current_node->left->right);
                        current_node->parent->left->d_balance = 0;
                        current_node->parent->right->d_balance = 1;
                        current_node->parent->d_balance = 0;
                        break;
                    case 0:
                        Double_Right_Rotation(current_node->left->right);
                        current_node->parent->left->d_balance = 0;
                        current_node->parent->right->d_balance = 0;
                        current_node->parent->d_balance = 0;
                        break;
                    case 1:
                        Double_Right_Rotation(current_node->left->right);
                        current_node->parent->left->d_balance = -1;
                        current_node->parent->right->d_balance = 0;
                        current_node->parent->d_balance = 0;
                        break;
                    }

                    break;
                case 0:
                    Right_Rotation(current_node->left);
                    current_node->parent->right->d_balance = -1;
                    current_node->parent->d_balance = 1;
                    break;

                case -1:
                    Right_Rotation(current_node->left);
                    current_node->parent->right->d_balance = 0;
                    current_node->parent->d_balance = 0;
                    break;
                }
                current_node = current_node->parent;
                break;
            }
            flag = true;
        }while (current_node != root && abs(current_node->d_balance) != 1);
    }
    void Print() {
        cout << "root: " << &root->data << '; ' << root->d_balance << endl;
        cout << "root->l: " << root->left->data << '; ' << root->left->d_balance; 
       /* cout << " root->r: " << root->right->data << '; ' << root->right->d_balance << endl;*/
       /* cout << "root->l->l" << root->left->left->data << '; ' << root->left->left->d_balance;
        cout << "root->l->r" << root->left->right->data << '; ' << root->left->right->d_balance;
        cout << "root->r->l" << root->right->left->data << '; ' << root->right->left->d_balance;
        cout << "root->r->r" << root->right->right->data << '; ' << root->right->right->d_balance;*/
    }

    Node* Search(long long int elem) {
        Node* current_node = root;
        if (current_node) {
            while (current_node->data != elem) {
                if (current_node->data > elem) {
                    current_node = current_node->left;
                    if (!current_node) return nullptr;
                }
                else if (current_node->data < elem) {
                    current_node = current_node->right;
                    if (!current_node) return nullptr;
                }
            }
            return current_node;
        }
        else return nullptr;
    }
};



//int main()
//{
//    
//    AVL_Tree Test;
//    for (int i = 0; i < 10; i++) {
//        long long int x;
//        cin >> x;
//        Test.Add(x);
//    }
//    /*Test.Delete(Test.root->right);
//    Test.Delete(Test.root->right->left);
//    Test.Delete(Test.root->left);
//    Test.Delete(Test.root->left);
//    Test.Delete(Test.root->left->left);
//    Test.Delete(Test.root->left->right);
//
//    int searched_elem;
//    cin >> searched_elem;
//
//    cout <<"\nSearch: "<< Test.Search(searched_elem)->data << endl;*/
//
//    Test.Print();
//    return 0;
//}

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

    ifstream inFile("test.txt");
    ofstream outFile("result.txt");

    int n;
    AVL_Tree Tree;
    inFile >> n;
    for (int i = 0; i < n; i++) {
        string command;
        long long int x;
        inFile >> command;
        
        switch (command[0])
        {
        case '+':
            inFile >> x;
            Tree.Add(x);
            if (Tree.root)outFile << Tree.root->d_balance<< '\n';
            else outFile << 0 << endl;
            break;
        case '-':
            inFile >> x;
            
            Tree.Delete(Tree.Search(x));
            if (Tree.root) outFile << Tree.root->d_balance << '\n';
            else outFile << 0 << endl;
            break;
        case '?':
            inFile >> x;
            if (Tree.Search(x)) outFile << "true" << endl;
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
