// Fast_Substr.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//



#include <fstream>

int* pref(char* T, int subsize) {

    int n = subsize;
    int t;
    int* p = new int[n + 1];
    p[0] = 0;
    for (int i = 1; i < n; i++) {
        int j = p[i - 1];
        while (j > 0 && T[i] != T[j])
            j = p[j - 1];
        if (T[i] == T[j])
            j++;
        p[i] = j;
    }

    return p;
}

int* find(char* str, char* sub, int strsize, int subsize) {

    int* p = pref(sub, subsize);
    int* solution = new int[strsize];
    solution[0] = 0;
    int position = 0;
    
    for (int i = 0; i < strsize; i++) {
        while (position == subsize || (position > 0 && sub[position] != str[i])) {
            position = p[position - 1];
            if (subsize - position > strsize - i)
                break;
        }
        if (str[i] == sub[position]) {
            position++;
        }
        if (position == subsize) {
            solution[0]++;
            solution[solution[0]] = (i - position + 1);
        }
    }
    delete[] p;
    return solution;
}


int main(int argc, char* argv[]) {
    if (argc < 3) {
        return 1;
    }
    std::ifstream inFile(argv[1]);
    if (!inFile) {
        return 2;
    }
    std::ofstream outFile(argv[2]);
    if (!outFile) {
        return 3;
    }
    
    
    /*std::ifstream inFile("C:/Users/User/Documents/GitHub/task10-fastsubstr-tottttttt/tests/in/input9.txt");
    std::ofstream outFile("result.txt");
    */
    int size_substr = -1;
    int size_str = -1;
    char a = '0';

    bool flag = true;
    while (!inFile.eof()) {
        if (flag) {
            inFile.get(a);
            size_substr++;
            if (a == '\n') flag = false;
        }
        else {
            inFile.get(a);
            size_str++;
        }
    }

    char* substr = new char[size_substr + 1];;
    char* str = new char[size_str + 1];
    inFile.clear();
    inFile.seekg(0, std::ios::beg);
    inFile.getline(substr, size_substr + 1);
    inFile.getline(str, size_str + 1);


    int* result = find(str, substr, size_str, size_substr);
    outFile << result[0] << '\n' << result[1] + 1;
    for (int i = 2; i <= result[0]; i++) outFile << ' ' << result[i] + 1;

    delete[] substr;
    delete[] str;
    delete[] result;

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
