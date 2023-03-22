// B-Tree.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include<fstream>
#include<string>

using namespace std;

class BNode {
public:
	BNode** child;					
	long long int* key;				
	int size;					
	bool leaf;						
};

class BTree {
public:
	BNode* root;
	unsigned minDegree;
	
	BTree(int t) {
		minDegree = t;
		root = new BNode;
		CreateNode(root);
		root->leaf = true;
		
	}

	~BTree() {
		deleteNode(root);
	}
		
	
	
	void CreateNode(BNode* x) {
		x->size = 0;
		x->key = new long long int[2 * minDegree - 1];
		x->child = new BNode * [2 * minDegree];
	}

	void deleteNode(BNode* x) {
		if (!x->leaf) {
			for (int i = 0; i <= x->size; i++) {
				deleteNode(x->child[i]);
			}
		}
		delete[](x->child);
		delete[](x->key);
		delete(x);
	}

	int nodeInsert(BNode* x, long long int k) {

		int i;
		
		for (i = x->size; i > 0 && k < x->key[i - 1]; i--) {
			x->key[i] = x->key[i - 1];
			x->child[i + 1] = x->child[i];
		}

		
		x->child[i + 1] = x->child[i];
		x->key[i] = k;
		x->size++;

		return i;
	}

	void splitChild(BNode* x, int i) {

		
		BNode* toSplit = x->child[i];
		BNode* newNode = new BNode;
		CreateNode(newNode);
		newNode->leaf = toSplit->leaf;
		newNode->size = minDegree - 1;

		
		for (int j = 0; j < minDegree - 1; j++) {
			newNode->key[j] = toSplit->key[j + minDegree];
		}
		if (!toSplit->leaf) {
			for (int j = 0; j < minDegree; j++) {
				newNode->child[j] = toSplit->child[j + minDegree];
			}
		}
		toSplit->size = minDegree - 1;

		nodeInsert(x, toSplit->key[minDegree - 1]);
		
		x->child[i + 1] = newNode;
	}

	
	void insert(long long int k) {
		if (root->size == 2 * minDegree - 1) {
			BNode* newRoot = new BNode;
			CreateNode(newRoot);
			newRoot->leaf = false;
			newRoot->child[0] = root;
			root = newRoot;
			splitChild(newRoot, 0);
		}

		
		BNode* current_node = root;
		while (!current_node->leaf) {

			int index = current_node->size - 1;
			while (index >= 0 && k <current_node->key[index]) {
				index--;
			}
			index++;

			if (current_node->child[index]->size == 2 * minDegree - 1) {
				splitChild(current_node, index);
				if (k > current_node->key[index]) {
					index++;
				}
			}
			current_node = current_node->child[index];
		}

		nodeInsert(current_node, k);

	}

	
	long long int searchKey(long long int k) {
		pair<BNode*, int> pair_node = search(k);
		if (pair_node.first == nullptr) {
			//
		}
		return pair_node.first->key[pair_node.second];
	}

	pair<BNode*, int> search(long long int k) {

		BNode* x = root;

		while (true) {

			int i = findIndex(x, k);
			if (i < x->size && k == x->key[i]) {
				return pair<BNode*, int>(x, i);
			}

			else if (x->leaf) {
				return pair<BNode*, int>(nullptr, 0);
			}

			else {
				x = x->child[i];
			}
		}
	}

	int findIndex(BNode* x, long long int k) {
		int i = 0;
		while (i < x->size && x->key[i] < k) {
			i++;
		}
		return i;
	}
	
	long long int remove(long long int k) {
		BNode* current_node = root;
		while (true){
			int i = findIndex(current_node, k);
			if (i < current_node->size && current_node->key[i] == k) {
				long long int toReturn = current_node->key[i];

				if (current_node->leaf) {
					nodeDelete(current_node, i);
				}

				else {
					BNode* left_child = current_node->child[i];
					BNode* right_child = current_node->child[i + 1];

					
					if (left_child->size >= minDegree) {
						while (!(left_child->leaf)) {
							fixChildSize(left_child, left_child->size);
							left_child = left_child->child[left_child->size];
						}
						current_node->key[i] = nodeDelete(left_child, left_child->size - 1);
					}
					else if (right_child->size >= minDegree) {
						while (!(right_child->leaf)) {
							fixChildSize(right_child, 0);
							right_child = right_child->child[0];
						}
						current_node->key[i] = nodeDelete(right_child, 0);
					}
					else {
						mergeChildren(current_node, i);
						current_node = left_child;
						continue;
					}
				}
				return toReturn;
			}

			else {

				char result = fixChildSize(current_node, i);
				if (result == 2) {
					current_node = root;
				}
				else {
					current_node = current_node->child[findIndex(current_node, k)];
				}
			}
		}
	}

	char mergeChildren(BNode* parent, int i) {

		BNode* left_child = parent->child[i];
		BNode* right_child = parent->child[i + 1];

		left_child->key[left_child->size] = nodeDelete(parent, i);
		int j = ++(left_child->size);

		for (int k = 0; k < right_child->size; k++) {
			left_child->key[j + k] = right_child->key[k];
			left_child->child[j + k] = right_child->child[k];
		}

		left_child->size += right_child->size;
		left_child->child[left_child->size] = right_child->child[right_child->size];

		
		delete[] right_child->child;
		delete[] right_child->key;
		delete right_child;

		if (parent->size == 0) {
			root = left_child;
			delete[] parent->child;
			delete[] parent->key;
			delete parent;
			return 2;
		}

		return 1;
	}

	char fixChildSize(BNode* parent, int index) {
		BNode* child = parent->child[index];

		if (child->size < minDegree) {

			if(index == 0 && index != parent->size && parent->child[index + 1]->size >= minDegree){
				BNode* right_child = parent->child[index + 1];
				nodeInsert(child, parent->key[index]);
				child->child[child->size] = right_child->child[0];
				
				right_child->child[0] = right_child->child[1];
				parent->key[index] = nodeDelete(right_child, 0);
			}
			
			
			else if (index != parent->size && parent->child[index + 1]->size >= minDegree && (parent->child[index + 1]->size > parent->child[index -1]->size)) {
				BNode* right_child = parent->child[index + 1];
				nodeInsert(child, parent->key[index]);
				
				child->child[child->size] = right_child->child[0];
				right_child->child[0] = right_child->child[1];
				parent->key[index] = nodeDelete(right_child, 0);
			}
			else if (index != 0 && parent->child[index - 1]->size >= minDegree) {
				BNode* left_child = parent->child[index - 1];
				for (int i = nodeInsert(child, parent->key[index - 1]); i >= 0; i--) {
					//child->child[i] = child->child[i - 1];
					child->child[i + 1] = child->child[i];
				}
				child->child[0] = left_child->child[left_child->size];
				parent->key[index - 1] = nodeDelete(left_child, left_child->size - 1);
			}

			else if (index != 0) {
				return mergeChildren(parent, index - 1);
			}
			else {
				return mergeChildren(parent, index);
			}
			return 1;
		}

		return 0;
	}

	long long int nodeDelete(BNode* x, int index) {

		long long int toReturn = x->key[index];

		x->size--;
		while (index < x->size) {
			x->key[index] = x->key[index + 1];
			x->child[index + 1] = x->child[index + 2];
			index++;
		}
		return toReturn;
	}

	

		
};


int main(/*int argc, char* argv[]*/)
{
	/*

	for (int i = 0; i < 15; i++) {
		int a;
		cin >> a;
		Test.insert(a);
	}
	*/

	

	
	/*
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
	}*/

	ifstream inFile("test.txt");
	ofstream outFile("result.txt");

	int n, t;
	inFile >> t >> n;
	
	

	BTree Tree(t);

	for (int i = 0; i < n; i++) {
		string command;
		long long int x;
		inFile >> command;

		switch (command[0])
		{
		case '+':
			inFile >> x;
			Tree.insert(x);
			outFile << Tree.root->size;
			for (int i = 0; i < Tree.root->size; i++)
				outFile << ' ' << Tree.root->key[i];
			outFile << '\n';
			break;
		case '-':
			inFile >> x;
			
			Tree.remove(x);
			if (Tree.root == nullptr) cout << '0' << '\n';
			else {
				outFile << Tree.root->size;
				for (int i = 0; i < Tree.root->size; i++)
					outFile << ' ' << Tree.root->key[i];
				outFile << '\n';
			}
			break;
		case '?':
			inFile >> x;
			if (Tree.search(x).first)  outFile << "true\n";
			else outFile << "false\n";
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
