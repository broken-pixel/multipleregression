#include <iostream>
#include <iomanip>
#include <cmath>
double** create(int, int);
void destroy(double**, int);
bool inverse(double **matrix, double **result, int size);
void print(double**, int, int, int);
void load1(double**, int, int);
double ** mult(double **X_tr,double **X, int rows, int columns);
double ** transpose(double **matrix, int rows, int columns);
void load(double**, int, int);

int main(){
int n, m = 4,y = 1;
double **X,**X_tr,**XtrX, **Y;

do{     
	std::cout<<"num of datapoints: ";
    std::cin>>n;
} while (n<=0);
	X=create(n, m);
	Y=create(n, y);
	load1(X, n, m);
	std::cout<<"\nX1 data:\n";
	load(X, n, 1);
	std::cout<<"\nX2 data:\n";
	load(X, n, 2);
	std::cout<<"\nX3 data:\n";
	load(X, n, 3);
	std::cout<<"\nY data:\n";
	load(Y, n, 0);

/////////////////traspose X matrix/////////////////
	X_tr=transpose(X,n,m);

/////////////////Xtr * X/////////////////
	double **Xt_x;
	Xt_x = mult(X_tr, X, n, m);

/////////////////inverse/////////////////
	double **xX;
	xX=create(m, m);
	inverse(Xt_x,xX,m);

/////////////////Xtr * Y matrix/////////////////
	double **a;
	a = new double*[m];
	  for (int i = 0; i<m; i++){
	    a[i] = new double[1];
	     for(int j = 0; j<1; j++){   
	      a[i][j] = 0;
	        for(int k = 0; k<n;k++)
	        a[i][j] += X_tr[i][k] * Y[k][j];
	    }
	}

/////////////////final matrix/////////////////
	double **result_matrix;
	result_matrix = new double*[m];
	  for (int i = 0; i<m; i++){
	    result_matrix[i] = new double[1];
	     for(int j = 0; j<1; j++){   
	      result_matrix[i][j] = 0;
	        for(int k = 0; k<m;k++)
	        result_matrix[i][j] += xX[i][k] * a[k][j];
	    }
	}

	std::cout<<"Result: \n";
	std::cout<<"intercept = "<<result_matrix[0][0]
			 <<"\n x1 = "<<result_matrix[1][0]
			 <<"\n x2 = "<<result_matrix[2][0]
			 <<"\n x3 = "<<result_matrix[3][0];
return 0;
}


//---------------------------------------------------------------------------
double ** transpose(double **matrix, int rows, int columns){
    double ** trans;                
    trans=new double *[columns];        
    for(int i=0;i<columns;i++){
        trans[i]=new double[rows];
        for(int j=0;j<rows;j++)
            trans[i][j]=matrix[j][i];
    }
    return trans;
}
double ** mult(double **X_tr,double **X, int rows, int columns){
    double **c;
	c = new double*[columns];        
    for (int i = 0; i<columns; i++){
	    c[i] = new double[columns];
	     for(int j = 0; j<columns; j++){   
	      c[i][j] = 0;
	        for(int k = 0; k<rows;k++)
	        c[i][j] += X_tr[i][k] * X[k][j];
	    }
	}
	return c;
}

double** create(int rows, int columns){
        double** matrix;
        matrix = new double* [rows];
        for(int i=0; i<rows; i++)
                matrix[i] = new double [columns];
        return matrix;
}

void destroy(double** matrix, int rows){       
        for(int i=0; i<rows; i++)
                delete[] matrix[i];
        delete[] matrix;
        matrix = NULL;
}

void load(double** matrix, int rows, int columns){       
        for(int i=0; i<rows; i++)
            std::cin>>matrix[i][columns];
}

void load1(double** matrix, int rows, int columns){       
        for(int i=0; i<rows; i++)
            matrix[i][0]=1;
        std::cout<<"OK";
}

void print(double** matrix, int rows, int columns, int width){      
std::cout<<std::endl;
       for(int i=0; i<rows; i++)
       {    for(int j=0; j<columns; j++)
                       std::cout<<std::setw(width)<<matrix[i][j];
             std::cout<<std::endl;
       }
}

// ���������:
//     matrix - ������� ��� ���������
//     result - ������� ������������ ������� ��� �������� ����������
//     size   - ����������� �������
// ����������:
//     true � ������ ��������� ���������, false � ��������� ������
bool inverse(double **matrix, double **result, int size){   
    // ���������� �������������� ������� �������� ���������
    // ��������� ��������� �������
    for (int i = 0; i < size; ++i){
        for (int j = 0; j < size; ++j)
            result[i][j] = 0.0;       
        result[i][i] = 1.0;
    }
    
    // ����� �������� �������
    double **copy = new double *[size]();
    
    // ��������� ����� �������� �������
    for (int i = 0; i < size; ++i){
        copy[i] = new double [size];     
        for (int j = 0; j < size; ++j)
            copy[i][j] = matrix[i][j];
    }
    
    // �������� �� ������� ������� 
    // ������ ����. �� ������ ����� ���������� ������ ���
    // � �������� ������� ������������ � ������� �����������
    for (int k = 0; k < size; ++k){
        // ���� ������� �� ������� ��������� � ��������
        // ������ - ����, �� ���� ������, ��� �������
        // ���� �� ������� �� �������, � ������ ������
        // �������
        if (fabs(copy[k][k]) < 1e-8){
            // ����, ��������� � ���, ��� ��� ��������� ����� �����
            bool changed = false;
            // ��� �� �������, ������������� ���� ��������
            for (int i = k + 1; i < size; ++i){
                // ���� ����� ������, ��� � ��� �� �������
                // ������� ��������� �������
                if (fabs(copy[i][k]) > 1e-8){
                    // ������ ��������� � �������� ������ �������
                    // ��� � �������� �������, ��� � � ���������
                    std::swap(copy[k],   copy[i]);
                    std::swap(result[k], result[i]);
                    // ������� ���� - �������� � ������������ ������ �����
                    changed = true;
                    break;
                }
            }
            
            // ���� ����� ����� ��������� �� ��� - ������� �� ����� ����
            // ��������
            if (!changed){
                // ������ ������
                for (int i = 0; i < size; ++i)
                    delete [] copy[i];
                delete [] copy;
                // �������� � ������� ���������
                return false;
            }
        }
        
        // ���������� �������� - ������������ �������
        double div = copy[k][k];
        
        // ��� �������� �������� ������ ����� �� ������������
        // ������� ��� � �������� �������, ��� � � ���������
        for (int j = 0; j < size; ++j){
            copy[k][j]   /= div;
            result[k][j] /= div;
        }
        
        // ��� �� �������, ������� ����������� ���� ��������
        for (int i = k + 1; i < size; ++i){
            // ���������� ��������� - ������� ��������� ������,
            // ������������� ��� ������������ ��������� ��������
            // ������
            double multi = copy[i][k];
            
            // �������� �� ��������� ������ ��������, ����������
            // �� ���������� ����� ��������� ��� � ��������,
            // ��� � � ��������� �������
            for (int j = 0; j < size; ++j){
                copy[i][j]   -= multi * copy[k][j];
                result[i][j] -= multi * result[k][j];
            }
        }
    }
    
    // �������� �� �������� ����������� �������, ����������
    // �� ������ ����, ����� �����
    // �� ������ ����� ���������� �������� ���, � �� ��������
    // ������� ������������ ����������� ���������, � �� ��������� -
    // ��������
    for (int k = size - 1; k > 0; --k){
        // ��� �� �������, ������� ����������� ���� ��������
        for (int i = k - 1; i + 1 > 0; --i){
            // ���������� ��������� - ������� ��������� ������,
            // ������������� ��� ������������ ��������� ��������
            // ������
            double multi = copy[i][k];
            // �������� �� ��������� ������ ��������, ����������
            // �� ���������� ����� ��������� ��� � ��������,
            // ��� � � ��������� �������
            for (int j = 0; j < size; ++j){
                copy[i][j]   -= multi * copy[k][j];
                result[i][j] -= multi * result[k][j];
            }
        }
    }
    for (int i = 0; i < size; ++i)
        delete [] copy[i];
    delete [] copy;
    return true;
}
