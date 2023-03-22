// SLE_matrix_transformation.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <bits/stdc++.h>
//#include "C:\Users\User\Desktop\admin\Университет\Вычислительная математика\SLE_matrix_transformation\SLE_martix_transformation.h"
using namespace std;

int count_oper;

class Matrix
{
private:
	
	Matrix* U = nullptr;
	Matrix* P = nullptr;
	Matrix* L = nullptr;
	Matrix* R = nullptr;
	Matrix* Q = nullptr;
	Matrix* Qr = nullptr;
	Matrix* Rr = nullptr;

public: 
	vector<vector<double>> matrix_array;
	int Dim_n;
	int Dim_m;
	double Determine = 1.;
	int rank = 0;
	Matrix(int n, int m, int special_elem);
	Matrix(int n, int m, vector<vector<double>> arr);
	Matrix(vector<double> arr);
	~Matrix();
	
	
	double elem(int i, int j) {
		return matrix_array[i][j];
	}
	void show_Matrix() {
		for (int i = 0; i < Dim_n; ++i){
			for (int j = 0; j < Dim_m; ++j) {
				/*if (abs(matrix_array[i][j]) <= 1e-10) cout << 0 << " ";
				else cout << matrix_array[i][j] << " ";*/
				cout << setprecision(8) << matrix_array[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
	void show_Matrix(char Spec_Symbol) {
		
		switch (Spec_Symbol)
		{
		case 'L':
			L->show_Matrix();
			break;
		case 'P':
			P->show_Matrix();
			break;

		case 'U':
			U->show_Matrix();
			break;
		case 'R':
			R->show_Matrix();
			break;
		case 'Q':
			Q->show_Matrix();
			break;
		default:
			cout << "Incorrekt type of matrix. Correct: L, P, U, R, Q" << endl;
			break;
		}
		
	}

	Matrix operator * (Matrix M1) {
		
		if (!(this->Dim_m - M1.Dim_n)) {
			Matrix result(this->Dim_n, M1.Dim_m, 0);
			vector<double> temp_row;
			double temp_elem = 0;
			for (int i = 0; i < this->Dim_n; i++) {
				for (int j = 0; j < M1.Dim_m; j++) {
					for (int r = 0; r < this->Dim_m; r++) {

						temp_elem += this->matrix_array[i][r] * M1.matrix_array[r][j];
					}
					temp_row.push_back(temp_elem);
					temp_elem = 0;
				}
				result.matrix_array[i] = temp_row;
				temp_row.clear();
			}
			 
			return result;
		}
		/*else {
		}*/
	}

	Matrix operator - (Matrix M1) {

		if (!(Dim_n-M1.Dim_n)&&!(Dim_m-M1.Dim_m)) {
			Matrix result(Dim_n, Dim_m, 0);
			for (int i = 0; i < Dim_n; i++) {
				for (int j = 0; j < Dim_m; j++) {
					result.matrix_array[i][j] = matrix_array[i][j] - M1.matrix_array[i][j];
				}
			}

			return result;
		}
		/*else {
		}*/
	}

	Matrix operator + (Matrix M1) {

		if (!(Dim_n - M1.Dim_n) && !(Dim_m - M1.Dim_m)) {
			Matrix result(Dim_n, Dim_m, 0);
			for (int i = 0; i < Dim_n; i++) {
				for (int j = 0; j < Dim_m; j++) {
					result.matrix_array[i][j] = matrix_array[i][j] + M1.matrix_array[i][j];
				}
			}

			return result;
		}
		/*else {
		}*/
	}

	void swap_row(int f_row, int s_row) {
		vector<double> tmp;
		
		tmp = matrix_array[s_row];
		matrix_array[s_row] = matrix_array[f_row];
		matrix_array[f_row] = tmp;
	}

	void swap_col(int f_col, int s_col) {
		double tmp;
		
		for (int i = 0; i < this->Dim_n; i++) {
			tmp = matrix_array[i][s_col];
			matrix_array[i][s_col] = matrix_array[i][f_col];
			matrix_array[i][f_col] = tmp;
		}
	}

	int index_max_col_elem(int i_col, int i_row_1, int i_row_2) {
		int max_i = i_row_1;
		
		for (int i = i_row_1; i < i_row_2; i++) {
			if (matrix_array[i][i_col] > matrix_array[max_i][i_col]) {
				max_i = i;
			}
		}
		return max_i;
	}

	pair<int, int> i_max_elem( int i_row_1, int i_row_2, int i_col_1, int i_col_2) {
		int i_max = i_row_1;
		int j_max = i_col_1;
		double tmp_max = abs(matrix_array[i_row_1][i_col_1]);
		for (int i = i_row_1; i < i_row_2; i++) {
			for (int j = i_col_1; j < i_col_2; j++) {
				if (abs(matrix_array[i][j]) > tmp_max) {
					i_max = i;
					j_max = j;
					
					tmp_max = abs(matrix_array[i][j]);
				}
			}
		}
		
		pair<int, int> max_elem = make_pair(i_max, j_max); 
		return max_elem;
	}

	double max_elem() {
		pair<int, int> i_max = i_max_elem(0, Dim_n, 0, Dim_m);
		return matrix_array[i_max.first][i_max.second];
	}
	
	void LU_Decomposition() {

			P = new Matrix(Dim_n, Dim_m, 1);
			Q = new Matrix(Dim_n, Dim_m, 1);
			L = new Matrix(Dim_n, Dim_m, 0);
			U = new Matrix(Dim_n, Dim_m, 0);
		
			U->matrix_array = this->matrix_array;

			int deg_p = 1;
			double k;
			int swap_index;
		
		

			for (int n = 0; n < Dim_n - 1; n++) {


				
				pair<int, int> coord_max = U->i_max_elem(n, Dim_n, n, Dim_m);
				int i_max = coord_max.first;
				int j_max = coord_max.second;
			
				if (n != i_max) {
					deg_p *= -1;
					
					U->swap_row(n, i_max);
					P->swap_row(n, i_max);
					L->swap_row(n, i_max);

				}
				if (n != j_max) {
					deg_p *= -1;
					
					U->swap_col(n, j_max);
					Q->swap_col(n, j_max);
					L->swap_col(n, j_max);
				}


						
				for (int i = n+1; i < Dim_n; i++) {
				
					k = U->matrix_array[i][n] / U->matrix_array[n][n];
					count_oper++; // для подсчёта арифметических операций в методе Ньютона
					L->matrix_array[i][n] = k;

					for (int j = n; j < Dim_n; j++) {
						U->matrix_array[i][j] -= U->matrix_array[n][j] * k;
						count_oper++; // для подсчёта арифметических операций в методе Ньютона
					}
				}
						
			}
		
			/*cout << "\n" << "U" << endl;
			U->show_Matrix();

			cout << "\n" << "Q" << endl;
			Q->show_Matrix();*/

			for (int i = 0; i < Dim_n; i++) {
				Determine *= U->matrix_array[i][i];
				L->matrix_array[i][i] = 1;
				if (abs(U->matrix_array[i][i]) > 1e-10) this->rank++;
			}

			Determine *= deg_p;

			/*cout << " \n rank= " <<rank<< endl;

			cout <<"\n"<< "L*U" << endl;
			((*L) * (*U)).show_Matrix();

			cout << "\n Check LU:\n";

			((*L) * (*U)  - ((*P) * (*this)) * (*Q)).show_Matrix();*/

		}

	void QR_Decomposition() {
		Qr = new Matrix(Dim_n, Dim_m, 1);
		Rr = new Matrix(Dim_n, Dim_m, 0);
		Matrix TMP(Dim_n, Dim_m, 0);

		TMP.matrix_array = this->matrix_array;

		

		double c;
		double s;
		double sqr;
		for (int j = 0; j < Dim_n - 1; j++) {
			for (int i = j + 1; i < Dim_n; i++) {
				sqr = sqrt((TMP.matrix_array[j][j] * TMP.matrix_array[j][j]) + (TMP.matrix_array[i][j] * TMP.matrix_array[i][j]));
				s = TMP.matrix_array[i][j] / sqr;
				c = TMP.matrix_array[j][j] / sqr;

				Matrix Q_i_j(Dim_n, Dim_m, 1);
				
				Q_i_j.matrix_array[i][j] = -s;
				Q_i_j.matrix_array[i][i] = c;
				Q_i_j.matrix_array[j][i] = s;
				Q_i_j.matrix_array[j][j] = c;

				*Qr = (*Qr) * Q_i_j.Transp();

				/*TMP = Q_i_j * TMP;*/  //
				
				for(int k = 0; k < Dim_m; k++){
					double	elem_jj = TMP.matrix_array[j][k],  
							elem_ij = TMP.matrix_array[i][k],  
							elem_ji = TMP.matrix_array[j][k],  
							elem_ii = TMP.matrix_array[i][k];

					TMP.matrix_array[j][k] = elem_jj * c + elem_ij * s;
					TMP.matrix_array[j][k] = elem_ji * c + elem_ii * s;
					TMP.matrix_array[i][k] = -elem_jj * s + elem_ij * c;
					TMP.matrix_array[i][k] = -elem_ji * s + elem_ii * c;
				}				
					
				/*cout << "\n" << "Q_" <<i << ","<<j<< endl;
				Q_i_j.show_Matrix();

				cout << "\n" << "TMP" << endl;
				TMP.show_Matrix();*/

			}
		}
		
		*Rr = Qr->Transp() * (*this);

		/*cout << "\n" << "-------QR _ DECOMPOSITION--------" << endl;
		
		cout << "\n" << "Qr" << endl;
		Qr->show_Matrix(); 
		
		Matrix Qr_t = Qr->Transp();
		cout << "\n" << "Q*Q'" << endl;
		((*Qr) * Qr_t).show_Matrix();
		
		cout << "\n" << "Q^-1 - Q'" << endl;
		
		Qr->LU_Decomposition();
		(Qr->Reverse() - Qr_t).show_Matrix();
				
		cout << "\n" << "TMP Matrix Rr" << endl;
		TMP.show_Matrix();

		cout << "\n" << "Rr" << endl;
		*Rr = (Qr_t) * (*this);
		Rr->show_Matrix();

		cout << "\n" << "Check Q*R " << endl;

		((*this) - ((*Qr) * *Rr)).show_Matrix();

		cout << "\n" << "Check Rr " << endl;

		(*Rr - TMP).show_Matrix();

		cout << "\n" << "Check TMP " << endl;

		((*this) - ((*Qr) * TMP)).show_Matrix();*/

		
	}

	Matrix QR_SLE_Decision(Matrix b1) {
		
		Matrix b = Qr->Transp() * b1;
		Matrix decision(Dim_n, 1, 0);
		vector<double> x(Dim_n);
		x[Dim_n - 1] = b.matrix_array[Dim_n-1][0] / Rr->matrix_array[Dim_n - 1][Dim_n - 1];
		double sum_x = 0;
		
		for (int i = Dim_n - 2; i >= 0; i--) {
			for (int j = i; j < Dim_m; j++) {
				sum_x += Rr->matrix_array[i][j] * x[j];
			}
			if (abs(Rr->matrix_array[i][i]) < 1e-10) x[i] = 0;
			else  x[i] = (1 / Rr->matrix_array[i][i]) * (b.matrix_array[i][0]- sum_x);
			sum_x = 0;
		}

		for (int i = 0; i < Dim_n; i++) {
			decision.matrix_array[i][0] = x[i];
		}

		return decision;
	}

	Matrix SLE_Decision(Matrix b1) {
		
		Matrix b = *P * b1 ;
		Matrix decision_b(Dim_n, 1, 0);
		vector<double> y(Dim_n);
		vector<double> x(Dim_n);
		double b_i = 0;
		double sum_x = 0;
		
		y[0] = b.matrix_array[0][0];
		
		for (int i = 1; i < Dim_n; i++) {
			for (int j = 0; j < i; j++) {
				b_i += L->matrix_array[i][j] * y[j];
				count_oper++; // для подсчёта арифметических операций в методе Ньютона
			}
			y[i] = b.matrix_array[i][0] - b_i;
			count_oper++; // для подсчёта арифметических операций в методе Ньютона
			b_i = 0;
		}
		
		
		x[Dim_n-1] = y[Dim_n-1]/ U->matrix_array[Dim_n-1][Dim_n-1];

		for (int i = Dim_n - 2; i >= 0; i--) {
			for (int j = i; j < Dim_m; j++) {
				sum_x += U->matrix_array[i][j] * x[j];
				count_oper++; // для подсчёта арифметических операций в методе Ньютона
			}
			if (abs(U->matrix_array[i][i])< 1e-10) x[i] = 0;
			else {
				x[i] = (1 / U->matrix_array[i][i]) * (y[i] - sum_x);
				count_oper++; // для подсчёта арифметических операций в методе Ньютона
				count_oper++; // для подсчёта арифметических операций в методе Ньютона
			}
			sum_x = 0;
		}

		for (int i = 0; i < Dim_n; i++) {
			decision_b.matrix_array[i][0] = x[i];
		}
		
		return (*Q)*decision_b;

	}

	Matrix Reverse() {
		R = new Matrix(Dim_n, Dim_m, 0);
		Matrix temp_col(Dim_n, 1, 0);
		Matrix temp_sol(Dim_n, 1, 0);
		for (int j = 0; j < Dim_m; j++) {
			temp_col.matrix_array[j][0] = 1;
			temp_sol = SLE_Decision(temp_col);
			for (int i = 0; i < Dim_n; i++) {
				
				R->matrix_array[i][j] = temp_sol.matrix_array[i][0];
			}
			temp_col.matrix_array[j][0] = 0;
		}
		return *R;
	}

	Matrix Transp() {
		Matrix* T = new Matrix(Dim_m, Dim_n, 0);


		for (int i = 0; i < Dim_n; i++) {
			for (int j = 0; j < Dim_m; j++) {
				T->matrix_array[j][i] = matrix_array[i][j];
			}
		}
		return *T;
	}

	Matrix Jakobi(Matrix b, double eps){
		Matrix B(Dim_n, Dim_m, 0);
		Matrix d(Dim_n, 1, 0);
		for (int i = 0; i < Dim_n; i++) {
			for (int j = 0; j < Dim_m; j++) {
				if (i - j) B.matrix_array[i][j] = -matrix_array[i][j] / matrix_array[i][i];
			}
			d.matrix_array[i][0] = b.matrix_array[i][0] / matrix_array[i][i];
		}

		/*cout << "\n" << "B:" << endl;
		B.show_Matrix();

		cout << "\n" << "d:" << endl;
		d.show_Matrix();*/
		
		return Method_Simple_Iter(B, d, eps);
			
		
	}

	Matrix Method_Simple_Iter(Matrix B, Matrix d, double eps) {
		
		Matrix x_n(Dim_n, 1, d.matrix_array);
		Matrix x_n1(Dim_n, 1, 0);
		int count_iter = 0;
		double epsilon = 1;
		double epsilon2 = 1;
		while (abs(epsilon2) > eps)
		{
			if (count_iter < 100) {

				x_n1 = B * x_n + d;
				epsilon = abs((x_n1 - x_n).Norm());
				x_n = x_n1;
				

				if (count_iter == 0) cout << "\nApriorn Jakobi" << log(eps * (1 - B.Norm()) / epsilon) / log(B.Norm()) << endl;
				count_iter++;
				epsilon2 = (B.Norm() / (1 - B.Norm())) * epsilon;	
			}
			else break;
		}

		

		cout << "\n" << "x_n MSI:" << endl;
		x_n.show_Matrix();

		cout << "\n" << "count of iteration Jacobi= " << count_iter << endl;
		return x_n;
	}

	Matrix Seidel_Method(Matrix b, double eps) {
		
		Matrix B(Dim_n, Dim_m, 0);
		Matrix c(Dim_n, 1, 0);
		Matrix x_n(Dim_n, 1, c.matrix_array);
		Matrix x_n1(Dim_n, 1, 0);
		
		int count_iter = 0;
		double epsilon = 1;
		double epsilon2 = 1;
		double tmp1;
		double tmp2;
		double NormB = -10000000;
		double NormBR= -10000000;

		for (int i = 0; i < Dim_n; i++) {
			for (int j = 0; j < Dim_m; j++) {
				if (i - j) B.matrix_array[i][j] = -matrix_array[i][j] / matrix_array[i][i];
			}
			c.matrix_array[i][0] = b.matrix_array[i][0] / matrix_array[i][i];
		}

		for (int i = 0; i < Dim_n; i++) {
			tmp1 = 0;
			tmp2 = 0;

			for (int j = 0; j < Dim_m; j++) {
				if (i != j) tmp1 += abs(B.matrix_array[i][j]);
			}
			if (tmp1 > NormB) {
				NormB = tmp1;
			}
			
			for (int j = i; j < Dim_m; j++) {
				if (i != j) tmp2 += abs(B.matrix_array[i][j]);
			}
			if (tmp2 > NormBR) {
				NormBR = tmp2;
			}

		}
		cout << "\n Norm B  " << NormB << "\n Norm BR  " << NormBR << endl;
		while(abs(epsilon2)> eps)
		{
			if (count_iter < 10000) {
				
				


				for (int i = 0; i < Dim_n; i++) {
					x_n1.matrix_array[i][0] = c.matrix_array[i][0];
					for (int j = 0; j < i; j++) {
						 x_n1.matrix_array[i][0] += (B.matrix_array[i][j] * x_n1.matrix_array[j][0]);
					}
					for (int j = i+1; j < Dim_m; j++) {
						x_n1.matrix_array[i][0] += (B.matrix_array[i][j] * x_n.matrix_array[j][0]);
					}
				}
				
				/*for (int i = 0; i < Dim_n; i++) {
					x_n1.matrix_array[i][0] = b.matrix_array[i][0] / A.matrix_array[i][i];
					for (int j = 0; j < i; j++) {
						x_n1.matrix_array[i][0] -= A.matrix_array[i][j] * x_n1.matrix_array[j][0] / A.matrix_array[i][i];
					}
					for (int j = i + 1; j < Dim_m; j++) {
						x_n1.matrix_array[i][0] -= A.matrix_array[i][j] * x_n.matrix_array[j][0] / A.matrix_array[i][i];
					}
				}*/
								
				epsilon = (x_n1 - x_n).Norm();
				x_n = x_n1;
				
				if (!count_iter) cout << "apriorn Seidel: " << log(eps * (1 - NormB) / epsilon) / log(NormBR) << endl;
				count_iter++;
				epsilon2 = (NormBR / (1 - NormB)) * epsilon;
			}
			else break;
		}

	
		cout << "\n" << "count of iteration Seidel= " << count_iter << endl;
		return x_n;
	}

	double Norm() {
		
		double max_sum = 0;
		double tmp_sum = 0;
		for (int i = 0; i < Dim_n; i++) {
			for (int j = 0; j < Dim_m; j++) {
				tmp_sum += abs(matrix_array[i][j]);
			}

			if (max_sum < tmp_sum) max_sum = tmp_sum;
			tmp_sum = 0;
		}
		return max_sum;
	}
	
	void diag_transf() {

		for (int i = 0; i < Dim_n; i++)
		{
			double diag_k = 0;
			for (int j = 0; j < Dim_m; j++) {
				diag_k += abs(matrix_array[i][j]);

				matrix_array[i][i] = diag_k;

			}
		}
	}
};

double double_rand(double min, double max) {
	return (double)(rand()) / RAND_MAX * (max - min) + min;
}

Matrix::Matrix(int n, int m, int special_elem)
{	
	this->Dim_n =n;
	this->Dim_m = m;
	
	switch (special_elem)
	{
	case 0:
	{	for (int i = 0; i < n; ++i)
	{
		vector<double> temp;
		for (int j = 0; j < m; ++j)
			temp.push_back(0);
		this->matrix_array.push_back(temp);

	}
	break;
	}
	case 1:
	{
		for (int i = 0; i < n; ++i)
		{
			vector<double> temp;
			for (int j = 0; j < m; ++j) {
				if (j - i) temp.push_back(0);
				else temp.push_back(1);
			}
			this->matrix_array.push_back(temp);
		}
		break;
	}
	case 2:
	{for (int i = 0; i < n; ++i)
		{
			vector<double> temp;
			for (int j = 0; j < m; ++j)
				temp.push_back(double_rand(-10, 10));
			/*temp.push_back(rand() % 20+1);*/
			this->matrix_array.push_back(temp);
		}
		break;
	}
	case 3:
	{for (int i = 0; i < n; ++i)
	{
		vector<double> temp;
		for (int j = 0; j < m; ++j)
			temp.push_back(double_rand(1, 10));
		/*temp.push_back(rand() % 10);*/
		this->matrix_array.push_back(temp);
	}
	double k1 = double_rand(0, 10);
	double k2 = double_rand(0, 10);
	double k3 = double_rand(0, 10);
	double k4 = double_rand(0, 10);


	for (int i = 0; i < n; i++) {
		matrix_array[i][2] = k1 * matrix_array[i][0] + k2 * this->matrix_array[i][1];
		matrix_array[i][3] = k3 * matrix_array[i][0] + k4 * this->matrix_array[i][1];
	}
	cout << endl << "k1= " << k1 << " k2= " << k2 << " k3= " << k3 << " k4= " << k4 << endl;
	break;
	}
	case 4:
	{for (int i = 0; i < n; ++i)
	{
		vector<double> temp;
		for (int j = 0; j < m; ++j)
			temp.push_back(abs(double_rand(0, 10)));
		/*temp.push_back(rand() % 20+1);*/
		this->matrix_array.push_back(temp);
	}

	*this = (*this) * this->Transp();
	/*for (int i = 0; i < Dim_n; i++) {
		matrix_array[i][i] = matrix_array[i][i]*sqrt(max_elem());
	}*/
	break; }
	}
	/*for (int i = 0; i < n; ++i)
	{
		vector<double> temp;
		for (int j = 0; j < m; ++j)
			temp.push_back(0);
		this->L_matrix_array.push_back(temp);
		this->U_matrix_array.push_back(temp);

	}

	for (int i = 0; i < n; ++i)
	{
		vector<double> temp;
		for (int j = 0; j < m; ++j)
			if (j - i) temp.push_back(0);
			else temp.push_back(1);
		this->P_matrix_array.push_back(temp);
	}*/


}
Matrix::Matrix(int n, int m, vector<vector<double>> arr) {
	
	
	this->Dim_n = n;
	this->Dim_m = m;
	
	

	for (int i = 0; i < n; i++) {
		
		vector<double> temp;
		for (int j = 0; j < m; j++) {
			
			temp.push_back(arr[i][j]);
		}
		matrix_array.push_back(temp);
	}


}
Matrix::Matrix(vector<double> arr) {


	this->Dim_n = arr.size();
	this->Dim_m = 1;



	for (int i = 0; i < arr.size(); i++) {

		vector<double> temp;
		temp.push_back(arr[i]);
		matrix_array.push_back(temp);
	}


}
Matrix::~Matrix()
{

}

void Test_SLE(Matrix test) {
	
	cout << "\n" << "Matrix:" << endl;

	test.show_Matrix();
	
	test.LU_Decomposition();




	//cout << "\n" << "L:" << endl;
	//test.show_Matrix('L');

	//cout << "\n" << "U:" << endl;
	//test.show_Matrix('U');

	//cout << "\n" << "P" << endl;

	//test.show_Matrix('P');


	cout << "\n" << "------Decision of SLE Ax =b------" << endl;
	{
		cout << "\n" << "b:" << endl;
		Matrix b(test.Dim_n, 1, 2);
		b.show_Matrix();
		test.SLE_Decision(b);
		vector<vector<double>>  array_b;
	double tmp_elem;

	for (int i = 0; i < test.Dim_n; ++i)
	{
		vector<double> temp;
		for (int j = 0; j < 1; ++j) {

			temp.push_back(test.elem(i, 0));
		}
		array_b.push_back(temp);

	}

		test.QR_Decomposition();
		cout << "\n" << "Solition QR:" << endl;

		test.QR_SLE_Decision(b).show_Matrix();
		cout << "\n" << "Check Solution QR:" << endl;
		((test * test.QR_SLE_Decision(b)) - b).show_Matrix();

		cout << "\n" << "Solition LU:" << endl;

		test.SLE_Decision(b).show_Matrix();

		cout << "\n" << "Check Solution LU:" << endl;
		((test * test.SLE_Decision(b)) - b).show_Matrix();

	}


	/*cout << "\n" << "Reverse:" << endl;
	Matrix test_reverse = test.Reverse();
	test.show_Matrix('R');

	cout << "\n" << "Check Reverse A*A^-1" << endl;
	(test * test_reverse).show_Matrix();

	cout << "\n" << "Check Reverse A^-1*A" << endl;
	(test_reverse * test).show_Matrix();

	cout << "\n" << "Condition Number: " << test.Norm() * test_reverse.Norm() << endl;

	cout << "Determine = " << test.Determine << endl;*/

}

void test_S_J_MSI(Matrix test) {
	
	
	cout << "\n" << "diag test:" << endl;

	test.diag_transf();
	test.QR_Decomposition();
	test.show_Matrix();

	cout << "\n" << "b:" << endl;
	Matrix b(test.Dim_n, 1, 2);

	/*vector<vector<double>>  array_b;
double tmp_elem;

for (int i = 0; i < Dim_n; ++i)
{
	vector<double> temp;
	for (int j = 0; j < 1; ++j) {

		temp.push_back(test.elem(i, 0));
	}
	array_b.push_back(temp);

}

Matrix b(Dim_n, 1, array_b);*/

	b.show_Matrix();

	double eps;
	cout << "Enter eps:";
	cin >> eps;

	Matrix TestS = test.Seidel_Method(b, eps);
	Matrix TestQR = test.QR_SLE_Decision(b);
	cout << "\n" << "Solition QR:" << endl;
	TestQR.show_Matrix();

	cout << "\nSeidel solution:\n";
	TestS.show_Matrix();

	cout << "\n" << "Check Solution Seidel:" << endl;
	(TestS - TestQR).show_Matrix();

	cout << "\nJakobi solution:\n";
	Matrix testJ = test.Jakobi(b, eps);


	cout << "\n" << "Check Solution Jakobi:" << endl;
	(testJ - TestQR).show_Matrix();
}

Matrix Jakobian(Matrix x) {
	
	vector<vector<double>> J;
	for (int i = 0; i < 10; i++) {

		vector<double> temp;
		for (int j = 0; j < 10; j++) {
			
			temp.push_back(0);
		}
		cout << "3-";
		J.push_back(temp);
	}
	
	double x1 = x.matrix_array[0][0];
	double x2 = x.matrix_array[1][0];
	double x3 = x.matrix_array[2][0];
	double x4 = x.matrix_array[3][0];
	double x5 = x.matrix_array[4][0];
	double x6 = x.matrix_array[5][0];
	double x7 = x.matrix_array[6][0];
	double x8 = x.matrix_array[7][0];
	double x9 = x.matrix_array[8][0];
	double x10 = x.matrix_array[9][0];
	
	{J[0][0] = -x2 * sin(x2 * x1);
	J[0][1] = -x1 * sin(x2 * x1);
	J[0][2] = 3.0 * exp(-(3.0 * x3));
	J[0][3] = x5 * x5;
	J[0][4] = 2.0 * x4 * x5;
	J[0][5] = -1.0;
	J[0][6] = 0.0;
	J[0][7] = -2.0 * cosh(2.0 * x8) * x9;
	J[0][8] = -sinh(2.0 * x8);
	J[0][9] = 2.0;
	J[1][0] = x2 * cos(x2 * x1);
	J[1][1] = x1 * cos(x2 * x1);
	J[1][2] = x9 * x7;
	J[1][3] = 0.0;
	J[1][4] = 6.0 * x5;
	J[1][5] = -exp(-x10 + x6) - x8 - 1.0;
	J[1][6] = x3 * x9;
	J[1][7] = -x6;
	J[1][8] = x3 * x7;
	J[1][9] = exp(-x10 + x6);
	J[2][0] = 1;
	J[2][1] = -1;
	J[2][2] = 1;
	J[2][3] = -1;
	J[2][4] = 1;
	J[2][5] = -1;
	J[2][6] = 1;
	J[2][7] = -1;
	J[2][8] = 1;
	J[2][9] = -1;
	J[3][0] = -x5 * pow(x3 + x1, -2.0);
	J[3][1] = -2.0 * x2 * cos(x2 * x2);
	J[3][2] = -x5 * pow(x3 + x1, -2.0);
	J[3][3] = -2.0 * sin(-x9 + x4);
	J[3][4] = 1.0 / (x3 + x1);
	J[3][5] = 0;
	J[3][6] = -2.0 * cos(x7 * x10) * x10 * sin(x7 * x10);
	J[3][7] = -1;
	J[3][8] = 2.0 * sin(-x9 + x4);
	J[3][9] = -2.0 * cos(x7 * x10) * x7 * sin(x7 * x10);
	J[4][0] = 2 * x8;
	J[4][1] = -2.0 * sin(x2);
	J[4][2] = 2 * x8;
	J[4][3] = pow(-x9 + x4, -2.0);
	J[4][4] = cos(x5);
	J[4][5] = x7 * exp(-x7 * (-x10 + x6));
	J[4][6] = -(x10 - x6) * exp(-x7 * (-x10 + x6));
	J[4][7] = (2 * x3) + 2.0 * x1;
	J[4][8] = -pow(-x9 + x4, -2.0);
	J[4][9] = -x7 * exp(-x7 * (-x10 + x6));
	J[5][0] = exp(x1 - x4 - x9);
	J[5][1] = -3.0 / 2.0 * x10 * sin(3.0 * x10 * x2);
	J[5][2] = -x6;
	J[5][3] = -exp(x1 - x4 - x9);
	J[5][4] = 2 * x5 / x8;
	J[5][5] = -x3;
	J[5][6] = 0;
	J[5][7] = -x5 * x5 * pow(x8, (-2));
	J[5][8] = -exp(x1 - x4 - x9);
	J[5][9] = -3.0 / 2.0 * x2 * sin(3.0 * x10 * x2);
	J[6][0] = cos(x4);
	J[6][1] = 3.0 * x2 * x2 * x7;
	J[6][2] = 1;
	J[6][3] = -(x1 - x6) * sin(x4);
	J[6][4] = x10 * pow(x5, (-2)) * cos(x10 / x5 + x8);
	J[6][5] = -cos(x4);
	J[6][6] = pow(x2, 3.0);
	J[6][7] = -cos(x10 / x5 + x8);
	J[6][8] = 0;
	J[6][9] = -1.0 / x5 * cos(x10 / x5 + x8);
	J[7][0] = 2.0 * x5 * (x1 - 2.0 * x6);
	J[7][1] = -x7 * exp(x2 * x7 + x10);
	J[7][2] = -2.0 * cos(-x9 + x3);
	J[7][3] = 1.5;
	J[7][4] = pow(x1 - 2.0 * x6, 2.0);
	J[7][5] = -4.0 * x5 * (x1 - 2.0 * x6);
	J[7][6] = -x2 * exp(x2 * x7 + x10);
	J[7][7] = 0;
	J[7][8] = 2.0 * cos(-x9 + x3);
	J[7][9] = -exp(x2 * x7 + x10);
	J[8][0] = -3;
	J[8][1] = -2.0 * x8 * x10 * x7;
	J[8][2] = 0;
	J[8][3] = exp((x5 + x4));
	J[8][4] = exp((x5 + x4));
	J[8][5] = -7.0 * pow(x6, -2.0);
	J[8][6] = -2.0 * x2 * x8 * x10;
	J[8][7] = -2.0 * x2 * x10 * x7;
	J[8][8] = 3;
	J[8][9] = -2.0 * x2 * x8 * x7;
	J[9][0] = x10;
	J[9][1] = x9;
	J[9][2] = -x8;
	J[9][3] = cos(x4 + x5 + x6) * x7;
	J[9][4] = cos(x4 + x5 + x6) * x7;
	J[9][5] = cos(x4 + x5 + x6) * x7;
	J[9][6] = sin(x4 + x5 + x6);
	J[9][7] = -x3;
	J[9][8] = x2;
	J[9][9] = x1;
	}
	Matrix J_x(10, 10, J);
	return J_x;
}
Matrix F(Matrix x) {
	
	
	vector<double> F_x;

	for (int i = 0; i < 10; i++) {
		F_x.push_back(0);
	}

	double x1 = x.matrix_array[0][0];
	double x2 = x.matrix_array[1][0];
	double x3 = x.matrix_array[2][0];
	double x4 = x.matrix_array[3][0];
	double x5 = x.matrix_array[4][0];
	double x6 = x.matrix_array[5][0];
	double x7 = x.matrix_array[6][0];
	double x8 = x.matrix_array[7][0];
	double x9 = x.matrix_array[8][0];
	double x10 = x.matrix_array[9][0];
	
	F_x[0] = cos(x2 * x1) - exp(-(3.0 * x3)) + x4 * x5 * x5 - x6 - sinh((2.0 * x8)) * x9 + (2.0 * x10) + 2.000433974165385440,
	F_x[1] = sin(x2 * x1) + x3 * x9 * x7 - exp(-x10 + x6) + 3.0 * x5 * x5 - x6 * (x8 + 1.0) + 10.886272036407019994;
	F_x[2] = x1 - x2 + x3 - x4 + x5 - x6 + x7 - x8 + x9 - x10 - 3.1361904761904761904;
	F_x[3] = 2.0 * cos(-x9 + x4) + x5 / (x3 + x1) - sin(x2 * x2) + pow(cos(x7 * x10), 2.0) - x8 - 0.1707472705022304757;
	F_x[4] = sin(x5) + 2.0 * x8 * (x3 + x1) - exp(-x7 * (-x10 + x6)) + 2.0 * cos(x2) - 1.0 / (-x9 + x4) - 0.3685896273101277862;
	F_x[5] = exp(x1 - x4 - x9) + x5 * x5 / x8 + cos(3.0 * x10 * x2) / 2.0 - x6 * x3 + 2.0491086016771875115;
	F_x[6] = pow(x2, 3.0) * x7 - sin(x10 / x5 + x8) + (x1 - x6) * cos(x4) + x3 - 0.7380430076202798014;
	F_x[7] = x5 * pow(x1 - 2.0 * x6, 2.0) - 2.0 * sin(-x9 + x3) + 1.5 * x4 - exp(x2 * x7 + x10) + 3.5668321989693809040;
	F_x[8] = 7.0 / x6 + exp(x5 + x4) - 2.0 * x2 * x8 * x10 * x7 + 3.0 * x9 - 3.0 * x1 - 8.4394734508383257499;
	F_x[9] = x10 * x1 + x9 * x2 - (x8 * x3) + sin(x4 + x5 + x6) * x7 - 0.78238095238095238096;

	Matrix Fx(F_x);
	return Fx;	
}

Matrix Newton(Matrix x, double eps) {
	
	Matrix x_k = x;
	Matrix x_k1(10, 1, 0);
	Matrix J(10, 10, 0);
	Matrix b_x(10, 1, 0);
	double epsilon = 1;
	
	int count_iter = 0;
	count_oper = 0;
	double start_time = clock();
	while (abs(epsilon) > eps ) {
		cout << "1-";
		if (count_iter < 1000) {
			cout << "2-";
			J = Jakobian(x_k);
			
			J.LU_Decomposition();

			b_x = J * x_k - F(x_k);
			cout << "4-";
			x_k1 = J.SLE_Decision(b_x);
			cout << "6-";
			cout << "7-";
			epsilon = abs((x_k1 - x_k).Norm());			
			x_k = x_k1;
			count_iter++;
		}
		else break;
	}
	
	/*cout << "\n" << "count of iteration Newton= " << count_iter << endl;
	cout << "\n" << "count of oper Newton= " << count_oper<< endl;
	cout << "\n" << "Time Newton= " << (clock()-start_time)/1000 <<" s"<< endl;*/
	return x_k;
}

Matrix Newton(Matrix x, double eps, int k) {

	Matrix x_k = x;
	Matrix x_k1(10, 1, 0);
	Matrix J(10, 10, 0);
	Matrix b_x(10, 1, 0);
	double epsilon = 1;
	bool Disc = true;

	int count_iter = 0;
	count_oper = 0;
	double start_time = clock();
	
	while (abs(epsilon) > eps || Disc) {
		if (count_iter < 100000) {
			
			if (count_iter < k) {
				J = Jakobian(x_k);
				J.LU_Decomposition();
				cout << "\ncheck J on iter " << count_iter << endl;
			}
			else if (abs(F(x_k).Norm()) - eps < 1e-12) Disc = false;
				
			b_x = J * x_k - F(x_k);
			x_k1 = J.SLE_Decision(b_x);
			epsilon = abs((x_k1 - x_k).Norm());
			x_k = x_k1;
			
			count_iter++;
		}
		else break;
	}
	
	cout << "\n" << "count of iteration K_Newton= " << count_iter << endl;
	cout << "\n" << "count of oper k_Newton= " << count_oper << endl;
	cout << "\n" << "Time k_Newton= " << (clock() - start_time)/1000  << " s" << endl;
	return x_k;
}

Matrix M_Newton(Matrix x, double eps) {

	Matrix x_k = x;
	Matrix x_k1(10, 1, 0);
	Matrix J(10, 10, 0);
	Matrix b_x(10, 1, 0);
	double epsilon = 1;

	int count_iter = 0;
	count_oper = 0;
	double start_time = clock();

	J = Jakobian(x);
	J.LU_Decomposition();
	
	while (abs(epsilon) - eps > 1e-12 || abs(F(x_k).Norm()) - eps > 1e-12) {
		if (count_iter < 1000) {
						
			b_x = J * x_k - F(x_k);
			x_k1 = J.SLE_Decision(b_x);
			epsilon = abs((x_k1 - x_k).Norm());
			x_k = x_k1;
			count_iter++;

		}
		else break;
	}

	cout << "\n" << "count of iteration M_Newton= " << count_iter << endl;
	cout << "\n" << "count of oper M_Newton= " << count_oper << endl;
	cout << "\n" << "Time M_Newton= " << (clock() - start_time) / 1000 << " s" << endl;
	return x_k;
}

Matrix G_Newton(Matrix x, double eps, int k) {

	Matrix x_k = x;
	Matrix x_k1(10, 1, 0);
	Matrix J(10, 10, 0);
	Matrix b_x(10, 1, 0);
	double epsilon = 1;

	int count_iter = 0;
	count_oper = 0;
	double start_time = clock();
	while (abs(epsilon) - eps > 1e-12) {
		if (count_iter < 1000) {

			if (count_iter % k == 0) {
				J = Jakobian(x_k);
				J.LU_Decomposition();
				cout << "\n Check J on " << count_iter << " iter\n";
			}

			b_x = J * x_k - F(x_k);
			x_k1 = J.SLE_Decision(b_x);

			/*cout << "\n" << "x_k1" << endl;
			x_k1.show_Matrix();*/

			epsilon = abs((x_k1 - x_k).Norm());
			/*cout << "\n" << "epsilon= " << epsilon << endl;*/
			x_k = x_k1;
			count_iter++;
		}
		else break;
	}

	cout << "\n" << "count of iteration G_Newton= " << count_iter << endl;
	cout << "\n" << "count of oper G_Newton= " << count_oper << endl;
	cout << "\n" << "Time G_Newton= " << (clock() - start_time) / 1000 << " s" << endl;
	return x_k;
}



int main()
{
	srand(time(NULL));
	//

	
	/*Matrix test = Matrix(2, 2, 2);
	test.LU_Decomposition();

	cout << endl;
	cout << endl; cout << endl; cout << endl; cout << endl;
	test.Reverse();*/
	
	//int Dim_n, Dim_m, Special_index;
	//
	//cout << "Enter Dim_n, Dim_n and special symbol: \n 0 - Zero Matrix \n 1 - Unit Matrix \n 2 - Random Matrix \n 3 - Singular Matrix \n 4 - Symmetrical Matrix \n *Also you can enter your array, what will be converted from object of class Matrix \n";
	//
	//cin >> Dim_n >> Dim_m >> Special_index;
	//
	//Matrix test = Matrix(Dim_n, Dim_m, 2);
	//Test_SLE(test);

	//cout << count_oper;

	Matrix start_x({ 0.234656985, 0.5986543, 1.554237868, -1.012345766, -0.2235465, 1.5065423, 0.543452, -0.5464664, 1.1241, -1.5991034 });
	
	//start_x.show_Matrix();

	cout << "\n Newton: \n";
	Newton(start_x, 1e-7).show_Matrix();
	
	/*cout << "\n M_Newton: \n";
	M_Newton(start_x, 1e-12).show_Matrix();*/


	/*cout << "\n M_Newton 1: \n";
	Newton(start_x, 1e-12, 1).show_Matrix();*/

	/*cout << "\n k_Newton: \n";

	int k;
	cout << "Enter k: \n";
	cin >> k;
	Newton(start_x, 1e-12, k).show_Matrix();

	cout << "\n g_Newton: \n";

	cout << "Enter k1: \n";
	int k1;
	cin >> k1;
	G_Newton(start_x, 1e-12, k1).show_Matrix();*/
	
}
//
