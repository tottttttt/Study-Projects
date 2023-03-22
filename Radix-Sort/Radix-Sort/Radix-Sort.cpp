#include <fstream>
#include<string>
#include<iostream>
#include<vector>
using namespace std;

int main(/*int argc, char* argv[]*/) {
    /*if (argc < 3) {
        return -1;
    }

    ifstream inFile(argv[1]);
    if (!inFile) {
        return -2;
    }

    ofstream outFile(argv[2]);
    if (!outFile) {
        return -3;
    }*/
    ifstream inFile("input.txt");
    ofstream outFile("result.txt");

    int N, M, K;
    vector<int> ind;
    vector<string> text;
    inFile >> N >> M >> K;
    for (int i = 0; i < N; i++) {
        ind.push_back(i);
    }
    int count = 0;
    int simb[256];
    //vector<int> ind(N);
    for (int i = 0; i < 256; i++) {
        simb[i] = 0;
    }
    for (int i = 0; i < M; i++) {
        string str;
        inFile >> str;
        text.push_back(str);
    }
    while (count < K) {
        count++;
        string str = text[text.size() - count];
        string tempstr = str;
        for (int i = 0; i < N; i++) {
            str[i] = tempstr[ind[i]];
        }
        for (int i = 0; i < N; i++) {
            int code = str[i] - 'A';
            //if ((char)tolower(str[i]) == str[i]) { code -= 6; }
            simb[code]++;
        }
        int tmp = simb[0];
        simb[0] = 0;
        for (int i = 1; i < 52; i++) {
            int tmp2 = tmp + simb[i];
            simb[i] = tmp;
            tmp = tmp2;
        }
        vector<int> tmpp;
        for (int i = 0; i < N; i++) {
            tmpp.push_back(i);
        }

        for (int i = 0; i < N; i++) {
            int code = str[i] - 'A';
            //if ((char)tolower(str[i]) == str[i]) { code -= 6; }
            tmpp[simb[code]] = ind[i];
            simb[code]++;
        }
        for (int i = 0; i < N; i++) {
            ind[i] = tmpp[i];
        }

        for (int i = 0; i < 52; i++) {
            simb[i] = 0;
        }

    }
    outFile << ind[0] + 1;
    for (int i = 1; i < N; i++) {
        outFile << ' ' << ind[i] + 1;
    }
    inFile.close();
    outFile.close();
    return 0;
}