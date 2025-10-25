#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <iomanip>
#include <variant>


using namespace std;


struct Max_abs{
    double value = 0;
    int i = 0;
    int j = 0;
};


void Print_Matrix(const vector<vector<double>>& matrix, string name){
    cout << "The matrix " << name << ": " << endl;
    cout << "---------------------------------------------------" << endl;
    for (int i = 0; i < matrix.size(); i++){
        for (int j = 0; j < matrix.size(); j++){
            cout << fixed << setprecision(4) << matrix[i][j] << "  ";
        }
        cout << endl;
    }
    cout << "---------------------------------------------------" << endl;
}


Max_abs Find_Max(const vector<vector<double>>& A){
    Max_abs max_abs;
    for (int i = 0; i < A.size(); i++){
        for (int j = 0; j < A.size(); j++){
            if (i != j){
                if (max_abs.value < abs(A[i][j])){
                    max_abs.value = abs(A[i][j]);
                    max_abs.i = i;
                    max_abs.j = j;
                }
            }
        }
    }
    return max_abs;
}


tuple<double, double, double> Calculate_p_c_s(const vector<vector<double>>& A, const Max_abs max_abs){
    double p = (2 * max_abs.value) / (A[max_abs.i][max_abs.i] - A[max_abs.j][max_abs.j]);
    double c = sqrt(0.5 * (1 + (1 / sqrt(1 + pow(p, 2)))));
    double s = copysign(1, p) * sqrt(0.5 * (1 - (1 / sqrt(1 + pow(p, 2)))));
    return make_tuple(p, c, s);
}


vector<vector<double>> Create_T(const vector<vector<double>>& A, const Max_abs max_abs, double p, double c, double s){
    vector<vector<double>> T(A.size(), vector<double>(A.size(), 0));
    for (int i = 0; i < T.size(); i++) T[i][i] = 1;

    T[max_abs.i][max_abs.i] = c;
    T[max_abs.j][max_abs.j] = c;
    T[max_abs.i][max_abs.j] = -s;
    T[max_abs.j][max_abs.i] = s;

    return T;
}


vector<vector<double>> Transp_Matrix(vector<vector<double>> matrix){
    vector<vector<double>> matrix_t(matrix.size(), vector<double>(matrix.size()));
    for (int i = 0; i < matrix.size(); i++){
        for (int j = 0; j < matrix.size(); j++){
            matrix_t[i][j] = matrix[j][i];
        }
    }
    return matrix_t;
}


vector<vector<double>> Mat_Mul(const vector<vector<double>>& mat1, const vector<vector<double>>& mat2){
    vector<vector<double>> res(mat1.size(), vector<double>(mat2.size()));
    
    for (int i = 0; i < mat1.size(); i++){
        for (int j = 0; j < mat1.size(); j++){
            for (int k = 0; k < mat1.size(); k++){
                res[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
    return res;
}


vector<vector<double>> To_zero_matrix(vector<vector<double>> B, double eps){
    for (int i = 0; i < B.size(); i++){
        for (int j = 0; j < B.size(); j++){
            if (i != j){
                if (B[i][j] < eps) B[i][j] = 0;
            }
        }
    }
    return B;
}


bool Check(const vector<vector<double>> B, double eps){
    vector<double> vec;

    for (int i = 0; i < B.size(); i++){
        for (int j = 0; j < B.size(); j++){
            if (i != j && i != B.size() - 1 && j != B.size() - 1){
                vec.push_back(pow(B[i][j], 2));
            }
        }
    }
    
    double check;
    for (double elem: vec){
        check += elem;
    }
    check = sqrt(check);

    if (sqrt(check) > eps) return true;
    else return false;
}


int main(){
    vector<vector<double>> A = {
        {0.0012, 0.5636, 0.1933, 0.8087},
        {0.5636, 0.4798, 0.3503, 0.8959},
        {0.1933, 0.3503, 0.7466, 0.1741},
        {0.8087, 0.8959, 0.1741, 0.7105}
    };

    double eps = 0.1;

    Print_Matrix(A, "A");

    bool check = true;
    int iterations = 1;
    vector<vector<vector<double>>> T_matrices;
    while (check == true){
        cout << "Iteration: " << iterations << endl;
        // ищем макс элем по модулю
        Max_abs max_abs = Find_Max(A);
        cout << "Max value outside the diagonal = " << max_abs.value << endl;
        cout << "i = " << max_abs.i << endl;
        cout << "j = " << max_abs.j << endl;

        //ищем 'p', 'c' и 's'
        auto [p, c, s] = Calculate_p_c_s(A, max_abs);
        cout << "---------------------------------------------------" << endl;
        cout << "p = " << p << endl;
        cout << "c = " << c << endl;
        cout << "s = " << s << endl;

        //создаём матрицу T
        vector<vector<double>> T = Create_T(A, max_abs, p, c, s);
        cout << "---------------------------------------------------" << endl;
        Print_Matrix(T, "T");
        T_matrices.push_back(T);

        //транспонируем T
        vector<vector<double>> T_t = Transp_Matrix(T);
        Print_Matrix(T_t, "T_t");

        //находим B
        vector<vector<double>> B;
        B = Mat_Mul(T_t, A);
        B = Mat_Mul(B, T);
        Print_Matrix(B, "B");

        //зануляем элемы B меньшие eps
        B = To_zero_matrix(B, eps);
        Print_Matrix(B, "B");

        //проверка условия остановки
        check = Check(B, eps);
        
        A = B;
        iterations++;
    }

    vector<vector<double>> U;
    U = Mat_Mul(T_matrices[0], T_matrices[1]);
    for (int i = 2; i < T_matrices.size(); i++){
        U = Mat_Mul(U, T_matrices[i]);
    }
    Print_Matrix(U, "U");

    vector<double> lambdas;
    for (int i = 0; i < A.size(); i++) lambdas.push_back(A[i][i]);

    vector<vector<double>> x(U[0].size(), vector<double>(U.size()));
    for (int i = 0; i < U.size(); i++){
        for (int j = 0; j < U.size(); j++){
            x[i][j] = U[j][i];
        }
    }

    for (int i = 0; i < lambdas.size(); i++){
        cout << "lambda " << i + 1 << " = " << lambdas[i] << endl;
        cout << "x" << i + 1 << " = ( ";
        for (int j = 0; j < x[i].size(); j++){
            if (j == x[i].size() - 1) cout << x[i][j] << " )" << endl;
            else{
                cout << x[i][j] << ", ";
            }
        }
    }

    return 0;
}