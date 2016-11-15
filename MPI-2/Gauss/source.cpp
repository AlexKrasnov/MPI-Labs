#include "mpi.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

using namespace std;

#define MAX_VALUE 10
#define SIZE_MATRIX 2
//#define FILE_INPUT
//#define FILE_OUTPUT
//#define KEYBOARD_INPUT

int* LeadingRowPos;  // Номера ведущих строк прямого хода Гаусса
int* LeadingRowIter; // Номера итераций прямого хода, на которых строки выбирались в качестве ведущих

void Input(double*& Matrix, double*& Result, int& Size) // ввод системы уравнений
{
	cout << "\t\t***The Gauss method - horizontal band scheme***\n\n";
	Matrix = new double[Size * (Size+1)];
	Result = new double[Size];
#ifdef FILE_INPUT
	ifstream fin("input.txt");
	for (int i = 0; i < Size; i++) 
		for (int j = 0; j < Size + 1; j++)
			fin >> Matrix [i * (Size + 1) + j];
	fin.close();
	/*#ifdef KEYBOARD_INPUT
	for (int i = 0; i < Size; i++) 
	for (int j = 0; j < Size + 1; j++)
	scanf("%lf", Matrix + i * (Size + 1) + j);*/
#else
	srand(unsigned int(time(NULL)));
	cout << "\nThe system of equation:\n";
	for(int i = 0; i < Size; i++, cout << endl)
		for(int j = 0; j < Size+1; j++)
		{
			Matrix[(Size+1) * i + j] = rand() % MAX_VALUE - MAX_VALUE / 2;
			cout << Matrix[(Size+1) * i + j] << " ";
			if (j == Size - 1) cout << "| ";
		}
#endif
}

void Output(double* Result, int Size) // вывод результата
{
#ifdef FILE_OUTPUT
	ofstream fout("output.txt");
	for (int i = 0; i < Size; i++) 
		fout << Result[i] << " ";
	fout.close();
#else
	cout << "\nResult:";
	for(int i = 0; i < Size; i++)
		cout << Result[i] << " ";
	cout << endl;
#endif
}

//------------------------------------------------------------------------
//      Последовательная реализация
//------------------------------------------------------------------------

int SerialFindPosLeadingRow(double* pMatrix, int Size, int Iter) // Определение позиции ведущей строки
{
	int LeadingRow = -1;
	double MaxValue = 0;
	for(int i = 0; i < Size; i++)
	{
		if( (LeadingRowIter[i] == -1) && (fabs(pMatrix[i*(Size+1) + Iter]) > MaxValue) )
		{
			LeadingRow = i;
			MaxValue =  pMatrix[i*(Size+1) + Iter];
		}
	}
	return LeadingRow;
}
void SerialGaussTransformation(double* pMatrix, int Size, int Iter, int PivotRow) // Преобразования метода Гаусса
{
	double LeadingRow, Koef;
	LeadingRow = pMatrix[PivotRow*(Size+1) + Iter];
	for (int i = 0; i < Size; i++)
		if (LeadingRowIter[i] == -1)
		{
			Koef = pMatrix[i*(Size+1) + Iter] / LeadingRow;
			for(int j = Iter; j < Size + 1; j++)
				pMatrix[i*(Size+1) + j] -= Koef * pMatrix[PivotRow*(Size+1) + j];
		}
}
void SerialGaussPart1(double* pMatrix, int Size) // Прямой ход
{
	int LeadingRow; // Номер ведущей строки на данной итерации
	for(int i = 0; i < Size; i++)
	{
		LeadingRow = SerialFindPosLeadingRow(pMatrix, Size, i); // Выбираем ведущую строку
		LeadingRowPos[i] = LeadingRow;
		LeadingRowIter[LeadingRow] = i;
		SerialGaussTransformation(pMatrix, Size, i, LeadingRow); // Делаем преобразования
	}
}
void SerialGaussPart2(double* pMatrix, double* pResult, int Size) // Обратный ход
{
	int RowIndex, Row;
	for(int i = Size - 1; i >= 0; --i)
	{
		RowIndex = LeadingRowPos[i];
		pResult[i] = pMatrix[RowIndex*(Size+1) + Size] / pMatrix[(Size+1) * RowIndex + i];
		pMatrix[RowIndex*(Size+1) + i] = 1;
		for(int j = 0; j < i; ++j)
		{
			Row = LeadingRowPos[j];
			pMatrix[Row*(Size+1) + Size] -= pMatrix[Row*(Size+1) + i] * pResult[i];
			pMatrix[Row*(Size+1) + i] = 0;
		}
	}
}

//------------------------------------------------------------------------
//      Параллельная реализация
//------------------------------------------------------------------------

typedef struct {int Row; double Value;} sPivotRow; 

void ParallelFindPosLeadingRow(double* pMatrixp, double* pSubMatr, int Size, int Iter, int ProcRank, int ProcNum, MPI_Datatype TPivot, int &PivotRow, int* pPivotIterL, int* pPivotIterp)
{   // Определение позиции ведущей строки
	sPivotRow* recvtmp = new sPivotRow[ProcNum];
	sPivotRow A; A.Row = -1; A.Value = 0;
	int k = Size/ProcNum;
	MPI_Scatter(pMatrixp, (Size+1)*k, MPI_DOUBLE, pSubMatr, (Size+1)*k, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	for(int i = 0; i < k; i++)
	{
		if( (pPivotIterL[i] == -1) && fabs(pSubMatr[i*(Size+1)+Iter]) > A.Value )
		{
			A.Value = pSubMatr[i*(Size+1)+Iter];
			A.Row = ProcRank*k + i;
		}
	} 
	MPI_Gather(&A, 1, TPivot, recvtmp, 1, TPivot,0, MPI_COMM_WORLD); 
	if( ProcRank == 0 )
	{
		int MaxValue = 0;
		for(int i = 0; i < ProcNum; i++)
			if( fabs(recvtmp[i].Value) > MaxValue && recvtmp[i].Row != -1 )
			{
				PivotRow = recvtmp[i].Row;
				MaxValue = int(recvtmp[i].Value);
			}	
	}
}

int main(int argc, char* argv[])
{
	int ProcNum, ProcRank; // число процессов, ранг процесса
	int Size = SIZE_MATRIX, PivotRow; // размер матрицы, ведущая строка
	double start_time, end_time, serial_time, parallel_time; // время работы
	double* pMatrixs(NULL), *pMatrixp(NULL), *pResults(NULL), *pResultp(NULL), *pSubMatr(NULL);
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Datatype TPivot;
	int blocklens[2] = {1, 1};
	MPI_Aint indices[2] = {offsetof(sPivotRow,Row), offsetof(sPivotRow, Value)};
	MPI_Datatype oldtypes[2] = {MPI_INT, MPI_DOUBLE};
	MPI_Type_create_struct(2, blocklens, indices, oldtypes, &TPivot);
	MPI_Type_commit(&TPivot);
	LeadingRowPos = new int[Size];
	LeadingRowIter = new int[Size];
	if (ProcRank == 0)
	{
	    fill(LeadingRowPos, LeadingRowPos + Size, 0);
	    fill(LeadingRowIter, LeadingRowIter + Size, -1);
		Input(pMatrixs, pResults, Size);
		pMatrixp = new double[(Size+1)*Size];
		for(int i = 0; i < Size; i++)
			for(int j = 0; j < Size+1; j++)
				pMatrixp[i*(Size+1)+j] = pMatrixs[i*(Size+1)+j]; 
		cout << "Serial realisation\n";
		start_time = MPI_Wtime();
		SerialGaussPart1(pMatrixs, Size);
		SerialGaussPart2(pMatrixs, pResults, Size);
		Output(pResults, Size);	
		end_time = MPI_Wtime();
		serial_time = end_time - start_time;
		cout << "Time: " << serial_time << "\n";
	} 
	MPI_Barrier(MPI_COMM_WORLD);
	fill(LeadingRowPos, LeadingRowPos + Size, 0);
	fill(LeadingRowIter, LeadingRowIter + Size, -1);
	int k = Size / ProcNum; // Число лент при распараллеливании, желательно, чтобы Size было кратно ProcNum
	int Iter; //Номер текущей итерации прямого хода
	int* LeadingRowIterL = new int[k];
	fill(LeadingRowIterL, LeadingRowIterL + Size, -1);
	pSubMatr = new double[(Size+1)*k];
	double* pPivotStr = new double[Size+1];
	start_time = MPI_Wtime();
	for(Iter = 0; Iter < Size; Iter++) // прямой ход
	{
		ParallelFindPosLeadingRow(pMatrixp, pSubMatr, Size, Iter, ProcRank, ProcNum, TPivot, PivotRow, LeadingRowIterL, LeadingRowIter);
		if(ProcRank == 0)
		{
			for(int i = 0; i < Size+1; i++)
				pPivotStr[i] = pMatrixp[PivotRow*(Size+1)+i];
			LeadingRowPos[Iter] = PivotRow;
			LeadingRowIter[PivotRow] = Iter;
		} 
		MPI_Scatter(LeadingRowIter, k, MPI_INT, LeadingRowIterL, k, MPI_INT, 0, MPI_COMM_WORLD); 
		MPI_Bcast(pPivotStr, Size+1, MPI_DOUBLE, 0, MPI_COMM_WORLD); //Собираем ведущую строку на всех процессах
		double LeadingRow = pPivotStr[Iter], koef = 1;
		for(int i = 0; i < k; i++) // преобразования метода Гауса для каждой ленты
		{
			if( LeadingRowIterL[i] == -1 )
			{
				koef = pSubMatr[i*(Size+1) + Iter] / LeadingRow;
				for(int j = Iter; j < Size; j++)
					pSubMatr[i*(Size+1) + j] -= koef * pPivotStr[j];
				pSubMatr[i*(Size+1) + (Size)] -= koef * pPivotStr[Size];
			}
		} 
		MPI_Gather(pSubMatr, k*(Size+1), MPI_DOUBLE, pMatrixp, k*(Size+1), MPI_DOUBLE,0, MPI_COMM_WORLD); 
		//Внешний цикл выполняется последовательно, распараллеливаются лишь части внутри него
	}  
	pResultp = new double[Size];
	MPI_Bcast(LeadingRowPos, Size, MPI_INT, 0, MPI_COMM_WORLD);
	for(int i = Size-1; i >= 0; --i) // обратный ход
	{
		int StrPos = LeadingRowPos[i]; //Номер строки, соответствующей i-ой итерации
		int rank = StrPos / k; //Номер процесса, который исключается из рассылки (такая строка у него уже есть)
		if( ProcRank == rank )
		{
			int StrPosL = StrPos % k;
			pResultp[i] = pSubMatr[StrPosL*(Size+1) + Size] / pSubMatr[(Size+1) * StrPosL + i];
			pSubMatr[StrPosL*(Size+1) + i] = 1;
			for(int m = 0; m < k; m++)
				if( m != StrPosL && i > LeadingRowIterL[m]) //Итерация, на которой данная строка была ведущей
				{
					pSubMatr[m*(Size+1) + Size] -= pSubMatr[m*(Size+1) + i] * pResultp[i];
					pSubMatr[m*(Size+1) + i] = 0;
				} 
		}
		MPI_Bcast(&pResultp[i], 1, MPI_DOUBLE, rank, MPI_COMM_WORLD); 
		if(ProcRank != rank)
			for(int m = 0; m < k; m++)
				if(i > LeadingRowIterL[m]) //Итерация, на которой эта строка была ведущей(m-я строка)
				{
					pSubMatr[m*(Size+1) + Size] -= pSubMatr[m*(Size+1) + i] * pResultp[i];
					pSubMatr[m*(Size+1) + i] = 0;
				}
	}
	end_time = MPI_Wtime();
	if( ProcRank == 0 )
	{
		cout << "\nParallel realisation";
		Output(pResultp, Size);
		parallel_time = end_time - start_time;
		cout << "\nTime: " << parallel_time << endl;
		cout << "\nBoost: " << serial_time / parallel_time << endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}