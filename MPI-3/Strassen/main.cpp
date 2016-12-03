#include <mpi.h>
#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>

using namespace std;

#define MATRIX_SIZE 128 // размер матарицы по умолчанию
#define THRESHOLD 1 // размер матриц для перехода к обычному умножению по умолчанию

short log2(const short Size) // логарифм по основанию 2
{
	return short(log(Size) / log(2));
}

bool ShowPower2ornot(const short Size) // проверяем является ли размер матрицы степенью двойки
{
	return (short(pow(2,log2(Size)))) == Size;
}

bool ShowFullSquareornot(const short Size) // проверяем является ли размер матрицы степенью двойки
{
	return (short(sqrt(Size))*short(sqrt(Size)) == Size);
}

short Increase (const short Size) // увеличиваем размер матрицы до ближайшей степени двойки
{
	short res;
	if (!ShowPower2ornot(Size)) res = 1 << (log2(Size) + 1);
	else res = 1 << (log2(Size));
	return res;
}

void Input(short *Matrix1, short *Matrix2, const short N) // формирование матриц
{
	cout << "\t\t***Multiplication dense matrices. Strassen Algorithm***\n\n";

	for (short i = 0; i < Increase(N)*Increase(N); i++)
	{
		Matrix1[i] = 0;
		Matrix2[i] = 0;
	}
	srand(unsigned short(time(NULL)));
	for (short i = 0; i < N; i++)
	{
		for (short j = 0; j < N; j++)
		{
			Matrix1[i*Increase(N)+j] = rand() % 10 - 5;
			Matrix2[i*Increase(N)+j] = rand() % 10 - 5;
		}
	}
	if (N <= 10)
	{
		cout << "Matrix1: " << endl;
		for (short i = 0; i < N; i++, cout << endl)
			for (short j = 0; j < N; j++)
				cout << Matrix1[i*Increase(N) + j] << " ";
		cout << endl << "Matrix2: " << endl;
		for (short i = 0; i < N; i++, cout << endl)
			for (short j = 0; j < N; j++)
				cout << Matrix2[i*Increase(N) + j] << " ";
		cout << endl;
	}
}

void Output(const short *Result, const short N) // печать матрицы
{
	if (N <= 10)
	{
		cout << "Result: " << endl;
		for (short i = 0; i < N; i++)
		{
			for (short j = 0; j < N; j++)
				cout << Result[i*Increase(N) + j] << " ";
			cout << endl;
		}
		cout << endl;
	}
}

short *Add(short *Matrix1, short *Matrix2, short N) // сложение матриц
{
	short *result = new short[N*N];
	for (short i = 0; i < N; i++)
		for (short j = 0; j < N; j++)
			result[i*N+j] = Matrix1[i*N+j] + Matrix2[i*N+j];
	return result;
}

short *Sub(short *Matrix1, short *Matrix2, short N) // вычитание матриц
{
	short *result = new short[N*N];
	for (short i = 0; i < N; i++)
		for (short j = 0; j < N; j++)
			result[i*N+j] = Matrix1[i*N+j] - Matrix2[i*N+j];
	return result;
}

/* Обычный метод умножения матриц*/
short *StandardAlgorithm(short *Matrix1, short *Matrix2, short N) 
{
	short *result = new short[N*N];
	fill(result,result+N*N, 0);
	for (short i = 0; i < N; i++)
		for (short j = 0; j < N; j++)
			for (short k = 0; k < N; k++)
				result[i*N+j] += Matrix1[i*N+k] * Matrix2[k*N+j];
	return result;
}

//------------------------------------------------------------------------
//      Последовательная реализация
//------------------------------------------------------------------------

/* Алгоритм Штрассена умножения матриц*/
short *StrassenAlgorithm(short *Matrix1, short *Matrix2, short N, short threshold)
{
	short *result;
	/* Когда размер матриц становится достаточно малым, 
	используем обычный метод умножения матриц, 
	так как алгоритм Штрассена теряет эффективность */
	if (N <= threshold)
		result = StandardAlgorithm(Matrix1, Matrix2, N);
	else
	{
		result = new short[N*N];
		N /= 2; // Матрица разбивается на 4 блока размерами N/2 * N/2

		/* Создаём вспомогательные матрицы, выделяем память */
		short *A[4], *B[4], *C[4], *P[7];
		for (short i = 0; i < 4; i++)
		{
			A[i] = new short [N*N];
			B[i] = new short [N*N];
		}
		/* Разбиваем матрицы на 4 блока */
		for (short i =0; i < N; i++)
		{
			for (short j = 0; j < N; j++)
			{
				//short index_new = i*N+j,index_old = 2*i*N+j, N_N = 2*N*N;
				A[0][i*N+j] = Matrix1[2*i*N+j];
				A[1][i*N+j] = Matrix1[2*i*N+j+N];
				A[2][i*N+j] = Matrix1[2*i*N+j+2*N*N];
				A[3][i*N+j] = Matrix1[2*i*N+j+2*N*N+N];

				B[0][i*N+j] = Matrix2[2*i*N+j];
				B[1][i*N+j] = Matrix2[2*i*N+j+N];
				B[2][i*N+j] = Matrix2[2*i*N+j+2*N*N];
				B[3][i*N+j] = Matrix2[2*i*N+j+2*N*N+N];
			}
		}

		/* Производим рекурсивный процесс N раз, до тех пор
		пока размер матриц не станет достаточно малым*/
		P[0] = StrassenAlgorithm(Add(A[0],A[3],N),Add(B[0],B[3],N),N,threshold);
		P[1] = StrassenAlgorithm(Add(A[2],A[3],N),B[0],N,threshold);
		P[2] = StrassenAlgorithm(A[0],Sub(B[1],B[3],N),N,threshold);
		P[3] = StrassenAlgorithm(A[3],Sub(B[2],B[0],N),N,threshold);
		P[4] = StrassenAlgorithm(Add(A[0],A[1],N),B[3],N,threshold);
		P[5] = StrassenAlgorithm(Sub(A[2],A[0],N),Add(B[0],B[1],N),N,threshold);
		P[6] = StrassenAlgorithm(Sub(A[1],A[3],N),Add(B[2],B[3],N),N,threshold);

		/* Находим результирующие блоки */
		C[0] = Sub(Add(Add(P[0],P[3],N),P[6],N),P[4],N); // P[0] + P[3] + P[6] - P[4]
		C[1] = Add(P[2],P[4],N); // P[2] + P[4];
		C[2] = Add(P[1],P[3],N); // P[1] + P[3];
		C[3] = Sub(Add(Add(P[0],P[2],N),P[5],N),P[1],N); // P[0] + P[2] + P[5] - P[1]

		/*Формируем результирующую матрицу*/
		for(short i = 0; i < N; i++)
			for(short j = 0; j < N; j++)
			{
				result[i*2*N+j] = C[0][i*N+j];
				result[i*2*N+j+N] = C[1][i*N+j];
				result[i*2*N+j+2*N*N] = C[2][i*N+j];
				result[i*2*N+j+2*N*N+N] = C[3][i*N+j];
			}

			for(short i = 0; i < 4; i++)
			{
				delete[] A[i];
				delete[] B[i];
				delete[] C[i];
			}
			for(short i = 0; i < 7; i++)
				delete[] P[i];
	}
	return result;
}

void Test(short *Matrix1, short *Matrix2, short *strassen_result, short N)
{
	short* res = StandardAlgorithm(Matrix1,Matrix2, N);
	bool f = true;
	for (int i = 0; i < N; i++)
	{
		if (res[i] != strassen_result[i])
		{
			f = false;
			break;
		}
	}
	if (f == true)
		cout << "\nThe results of the standard algorithm and Strassen algorithm are equal\n" << endl;
	else
		cout << "\nThe results of the standard algorithm and Strassen algorithm are NOT equal\n" << endl;
}

void MegaTest(short Size)
{
	for (int i = 0; i < Size; i++)
	{
		short *res1 = new short [i];
		short *res2 = new short [i];
	}
}

short main(int argc, char *argv[])
{
	int rank, procnum;
	short sqr, new_N, Size = MATRIX_SIZE, threshold = THRESHOLD, 
		**A, **B, **tmp, *Matrix1, *Matrix2, *result_s, *result_p;
	if (argc >= 2) Size = atoi(argv[1]);
	if (argc >= 3) threshold = atoi(argv[2]);
	short new_Size = Increase(Size);
	double start_time, end_time, serial_time, parallel_time;
	MPI_Status Status;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &procnum);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	try
	{

		if (rank == 0)
		{
			if (!ShowPower2ornot(procnum) || !ShowFullSquareornot(procnum))
				throw "The number of processes must be a power of 2 and a full square";
			Matrix1 = new short [new_Size*new_Size]; 
			Matrix2 = new short [new_Size*new_Size];
			result_p = new short [new_Size*new_Size]; 
			Input(Matrix1, Matrix2, Size);
			start_time = MPI_Wtime();
			result_s = StrassenAlgorithm(Matrix1,Matrix2,Increase(Size),threshold);
			end_time = MPI_Wtime();
			cout << "Serial realisation: " << endl;
			Output(result_s, Size); // результат умножения матриц обычным методом
			serial_time = end_time - start_time;
			cout << "Serial time: " << fixed << serial_time << endl;
			//Test(Matrix1,Matrix2,result_s,Increase(Size));

			//------------------------------------------------------------------------
			//      Параллельная реализация
			//------------------------------------------------------------------------

			start_time = MPI_Wtime();
			MPI_Bcast(&new_Size, 1, MPI_SHORT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&threshold, 1, MPI_SHORT, 0, MPI_COMM_WORLD);
			sqr = (short)sqrt((double)procnum), new_N = new_Size/sqr;
			A = new short*[procnum], B = new short*[procnum];
			for(short i = 0; i < procnum; i++)
			{
				A[i] = new short [new_N*new_N];
				B[i] = new short [new_N*new_N];
			}
			for(short i = 0; i < new_Size; i++)
			{
				for(short j = 0; j < new_Size; j++)
				{
					A[sqr*(i/new_N)+j/new_N][(i%new_N)*new_N+(j%new_N)] = Matrix1[i*new_Size+j];
					B[sqr*(i/new_N)+j/new_N][(i%new_N)*new_N+(j%new_N)] = Matrix2[i*new_Size+j];
				}
			}
			/*Рассылка данных другим процессам*/
			MPI_Bcast(&new_Size, 1, MPI_SHORT, 0, MPI_COMM_WORLD); 
			MPI_Bcast(&threshold, 1, MPI_SHORT, 0, MPI_COMM_WORLD);
			for(short i = 1; i < procnum; i++)
			{
				short coef_A = sqr*(i / sqr), coef_B = i % sqr;
				for(short j = 0; j < sqr; j++)
				{
					MPI_Send(A[coef_A], new_N*new_N, MPI_SHORT, i , 0, MPI_COMM_WORLD); 
					MPI_Send(B[coef_B], new_N*new_N, MPI_SHORT, i , 0, MPI_COMM_WORLD); 
					coef_A++;
					coef_B += sqr;
				}
			}

			/*Swap указателей для однообразных вычислений*/
			for(short i = 0; i < sqr; i++)
			{
				short* С = B[i];
				B[i] = B[i*sqr];
				B[i*sqr] = С;
			}
		}
		/*Прием данных от процесса-root и формировка нужных данных*/
		if(rank != 0)
		{
			MPI_Bcast(&new_Size, 1, MPI_SHORT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&threshold, 1, MPI_SHORT, 0, MPI_COMM_WORLD);
			sqr = (short)sqrt((double)procnum), new_N = new_Size/sqr;
			A = new short*[sqr], B = new short*[sqr];
			for(short i = 0; i < sqr; i++)
			{
				A[i] = new short [new_N*new_N];
				B[i] = new short [new_N*new_N];
			}
			for(short i = 0; i < sqr; i++)
			{
				MPI_Recv(A[i], new_N*new_N, MPI_SHORT, 0, 0, MPI_COMM_WORLD, &Status);
				MPI_Recv(B[i], new_N*new_N, MPI_SHORT, 0, 0, MPI_COMM_WORLD, &Status);
			}
		}
		/*вычисление каждым процессом своего куска матрицы*/
		tmp = new short*[sqr+1];
		for(short i = 0; i < sqr; i++)
			tmp[i+1] = StandardAlgorithm(A[i], B[i], new_N);
		if (procnum == 1)
			tmp[0] = tmp[1];
		if(procnum == 4)
			tmp[0] = Add(tmp[1], tmp[2], new_N);
		if(procnum == 16)
			tmp[0] = Add(Add(Add(tmp[1], tmp[2], new_N),tmp[3],new_N),tmp[4],new_N);
		/*Освобождение вспомогательной памяти*/
		if(procnum == 0)
		{
			for(short i = 0; i < procnum; i++)
			{
				delete[] A[i];
				delete[] B[i];
			}
		}
		else
		{
			for(short i = 0; i < sqr; i++)
			{
				delete[] A[i];
				delete[] B[i];
			}
		}
		for(short i = 1; i < sqr+1; i++)
			delete[] tmp[i];
		delete[] A;
		delete[] B;

		/*Отправка результата на 0 процесс*/
		if (rank != 0)	
		{
			MPI_Send(tmp[0], new_N*new_N, MPI_SHORT, 0, 0, MPI_COMM_WORLD);
			delete[] tmp[0];
		}
		if (rank == 0)
		{
			short coef = (short)sqrt((double)procnum);
			/*записываем результат совей работы*/
			for(short i = 0; i < new_N; i++)
				for(short j = 0; j < new_N; j++)
					result_p[coef*i*new_N+j] = tmp[0][i*new_N+j];
			for(short k = 1; k < procnum; k++)
			{
				/*принимаем и записываем результаты работы других процессов*/
				MPI_Recv(tmp[0], new_N*new_N, MPI_SHORT, k, 0, MPI_COMM_WORLD, &Status);
				for(short i = 0; i < new_N; i++)
					for(short j = 0; j < new_N; j++)
						result_p[(k/coef)*new_N*new_Size+(k%coef)*new_N+coef*i*new_N+j] = tmp[0][i*new_N+j];
			}
		}
		if (rank == 0)
		{
			end_time = MPI_Wtime();
			cout << "Parallel realisation: " << endl;
			Output(result_p, Size); // результат умножения матриц обычным методом
			parallel_time = end_time - start_time;
			cout << "Parallel time: " << fixed << parallel_time << endl;
			bool flag = true;
			for (short i = 0; (i < new_Size * new_Size) && (flag == true); i++)
				if (result_s[i] != result_p[i])
					flag = false;
			if(flag) cout << "The results are equal" << endl;
			else cout << "The results are not equal" << endl;
			delete [] result_p;
			delete [] result_s;
			delete [] Matrix1;
			delete [] Matrix2;
			cout << fixed << "Boost: "<< serial_time / parallel_time << endl;
		}

		MPI_Finalize();
		return 0;
	}
	catch(const char* error)
	{
		cout << error;
	}
}