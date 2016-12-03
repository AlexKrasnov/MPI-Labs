#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

using namespace std;

#define MAX_VALUE 9
#define SIZE_MATRIX 200
//#define FILE_INPUT
//#define FILE_OUTPUT
//#define KEYBOARD_INPUT

int* LeadingRowPos;  // Номера ведущих строк прямого хода Гаусса
int* LeadingRowIter; // Номера итераций прямого хода, на которых строки выбирались в качестве ведущих

void Input(double*& Matrix, double*& Result, int& Size) // ввод системы уравнений
{
	cout << "\n\t\t***The Gauss method - horizontal band scheme***\n\n";
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
	cout << "The system of equation:\n";
	for (int i = 0; i < Size; i++)
		for (int j = 0; j < Size + 1; j++)
			Matrix[(Size+1) * i + j] = rand() % MAX_VALUE + 1;
	if (Size <= 10)
		for(int i = 0; i < Size; i++, cout << endl)
			for(int j = 0; j < Size + 1; j++)
			{
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
	if(Size <= 10)
	{
		cout << "\nResult :";
		for(int i = 0; i < Size; i++)
			cout << floor(Result[i]*1000 + 0.5)/1000. << " ";
		cout << endl;
	}
#endif
}

void Validation(double* sResult, double* pResult, int N)
{
	bool f = true;
	if (sResult != NULL && pResult != NULL)
		for (int i = 0; i < N && f == true; i++)
		{
			if (floor(sResult[i]*1000 + 0.5)/1000. != floor(pResult[i]*1000 + 0.5)/1000.)
				f = false;
		}
		if (f == true)
			cout << "\nResults are equal\n" << endl;
		else 
			cout << "\nResults are not equal\n" << endl;
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

void ParallelCalculation(double* Matrix, double* Result , int Size , int ProcNum , int ProcRank)
{
	int PivotRow; // ведущая строка
	int *pPivotPosp = new int[Size]; // Номера строк, которые были выбраны ведущими
	int *pPivotIterp = new int[Size]; // Номера итераций, на которых строки процессора
	// использовались в качестве ведущих
	fill(pPivotIterp, pPivotIterp+Size,-1);
	double* pSubMatr = NULL;
	typedef struct {int Row; double Value;} sPivotRow; // Структура для выбора ведущего элемента
	// sPivotRow используется для хранения текущей строки матрицы и соответствующего жлемента вектора b
	MPI_Datatype TPivot;
	int blocklens[2] = {1, 1};
	MPI_Aint indices[2] = {offsetof(sPivotRow,Row ), offsetof(sPivotRow, Value)};
	MPI_Datatype oldtypes[2] = {MPI_INT, MPI_DOUBLE};
	MPI_Type_create_struct(2, blocklens, indices, oldtypes, &TPivot);
	MPI_Type_commit( &TPivot );
	int k = Size / ProcNum, p = k + Size % ProcNum, Iter, *pPivotIterL;
	if(ProcRank == 0)
	{
		pPivotIterL = new int[p];
		fill(pPivotIterL,pPivotIterL+p,-1);
	} 
	else 
	{
		pPivotIterL = new int[k];
		fill(pPivotIterL,pPivotIterL+k,-1);
	}
	int *send_counts, *displs, *send_countsg, *displsg;
	send_counts = new int[ProcNum];
	displs = new int[ProcNum];
	send_countsg = new int[ProcNum];
	displsg = new int[ProcNum];
	send_counts[0] = p*(Size+1);
	for(int i = 1; i < ProcNum;i++)
		send_counts[i] = k*(Size+1);
	displs[0] = 0;
	for(int i = 1; i < ProcNum;i++)
		displs[i] = p*(Size+1) + (i-1)*k*(Size+1);
	send_countsg[0] = p*(Size+1);
	for(int i = 1; i < ProcNum;i++)
		send_countsg[i] = k*(Size+1);
	displsg[0] = 0;
	for(int i = 1; i < ProcNum;i++)
		displsg[i] = p*(Size+1) + (i-1)*k*(Size+1);
	double* pPivotString = new double[Size+1];
	if(ProcRank == 0) pSubMatr = new double[(Size+1)*p];
	else pSubMatr = new double[(Size+1)*k];
	double* pPivotStr = new double[Size+1];
	MPI_Scatterv(Matrix, send_counts, displs, MPI_DOUBLE, pSubMatr, send_counts[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
	for(Iter = 0; Iter < Size; Iter++) // прямой ход метода Гауса
	{
		// Определение позиции ведущей строки
		sPivotRow* recvtmp = new sPivotRow[ProcNum];
		sPivotRow A; A.Row = -1; A.Value = 0;
		if(ProcRank == 0)
		{
			int s;
			if (ProcRank == 0) s = p;
			else s = k;
			for(int i = 0; i < p; i++)
				if( (pPivotIterL[i] == -1) && fabs(pSubMatr[i*(Size+1)+Iter]) > A.Value )
				{
					A.Value = pSubMatr[i*(Size+1)+Iter];
					A.Row = i;
				}
		} 
		else 
		{
			for(int i = 0; i < k; i++)
				if( (pPivotIterL[i] == -1) && fabs(pSubMatr[i*(Size+1)+Iter]) > A.Value )
				{
					A.Value = pSubMatr[i*(Size+1)+Iter];
					A.Row = p + (ProcRank-1)*k + i;
				}
		}
		MPI_Gather(&A, 1, TPivot, recvtmp, 1, TPivot,0, MPI_COMM_WORLD); 
		if( ProcRank == 0 )
		{
			double MaxValue = 0;
			for(int i = 0; i < ProcNum; i++)
				if( fabs(recvtmp[i].Value) > MaxValue && recvtmp[i].Row != -1 )
				{
					PivotRow = recvtmp[i].Row;
					MaxValue = recvtmp[i].Value;
				}	
		}
		MPI_Bcast(&PivotRow, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if(ProcRank == 0)
			pPivotPosp[Iter] = PivotRow;
		int r;
		if(PivotRow < p) r = 0;
		else  r = (PivotRow - p) / k + 1; 
		int Str;
		if( r == 0 ) Str = PivotRow;
		else Str = (PivotRow - p - k*(ProcRank - 1))% k;
		if( ProcRank == r )	
		{
			for(int i = 0; i < Size+1; i++)
				pPivotStr[i] = pSubMatr[Str*(Size+1)+i];
			pPivotIterL[Str] = Iter;
		}
		MPI_Bcast(pPivotStr, Size+1, MPI_DOUBLE, r, MPI_COMM_WORLD);
		double PivotValue, Koef;
		PivotValue = pPivotStr[Iter];
		int s = 0;
		if (ProcRank == 0) s = p;
		else s = k;
		for(int i = 0; i < s; i++) // преобразования прямого хода
			if( pPivotIterL[i] == -1 )
			{
				Koef = pSubMatr[i*(Size+1) + Iter] / PivotValue;
				for(int j = Iter; j < Size; j++)
					pSubMatr[i*(Size+1) + j] -= Koef * pPivotStr[j];
				pSubMatr[i*(Size+1) + (Size)] -= Koef * pPivotStr[Size];
			}
	}  
	MPI_Gatherv(pSubMatr, send_countsg[ProcRank] , MPI_DOUBLE , Matrix , send_countsg, displsg, MPI_DOUBLE,0, MPI_COMM_WORLD);
	MPI_Bcast(pPivotPosp, Size, MPI_INT, 0, MPI_COMM_WORLD);
	for(int i = Size-1; i >= 0; --i) // обратный ход метода Гауса
	{
		int StrPos = pPivotPosp[i];  // Номер строки, соответствующей i-ой итерации
		int rank = StrPos / k; // Номер процесса, который исключается из рассылки, у него есть такая строка
		if( ProcRank == rank )
		{
			int StrPosL = StrPos % k;
			Result[i] = pSubMatr[StrPosL*(Size+1) + Size] / pSubMatr[(Size+1) * StrPosL + i];
			pSubMatr[StrPosL*(Size+1) + i] = 1;
			for(int m = 0; m < k; m++)
				if( m != StrPosL && i > pPivotIterL[m]) // Итерация, на которой данная строка была ведущей
				{
					pSubMatr[m*(Size+1) + Size] -= pSubMatr[m*(Size+1) + i] * Result[i];
					pSubMatr[m*(Size+1) + i] = 0;
				} 
		}
		MPI_Bcast(&Result[i], 1, MPI_DOUBLE, rank, MPI_COMM_WORLD); 
		if(ProcRank != rank)
			for(int m = 0; m < k; m++)
				if(i > pPivotIterL[m]) //Итерация, на которой эта строка была ведущей(m-я строка)
				{
					pSubMatr[m*(Size+1) + Size] -= pSubMatr[m*(Size+1) + i] * Result[i];
					pSubMatr[m*(Size+1) + i] = 0;
				}
	}
} 

int main(int argc, char* argv[])
{
	int ProcNum, ProcRank, Size; // число процессов, ранг процесса, ведущая строка, размер матрицы
	if (argc < 2) Size = SIZE_MATRIX ;
	else Size = atoi(argv[1]);
	double start_time, end_time, serial_time, parallel_time; // время работы
	double* pMatrixs(NULL), *pMatrixp(NULL), *pResults(NULL), *pResultp(NULL), *pSubMatr(NULL);
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
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
		cout << "\nSerial realisation:\n";
		start_time = MPI_Wtime();
		SerialGaussPart1(pMatrixs, Size);
		SerialGaussPart2(pMatrixs, pResults, Size);
		Output(pResults, Size);	
		end_time = MPI_Wtime();
		serial_time = end_time - start_time;
		cout << "\nTime: " << serial_time << "\n";

		start_time = MPI_Wtime();
	} 
	MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	pResultp = new double [Size];
	ParallelCalculation(pMatrixp,pResultp,Size,ProcNum,ProcRank);
	if(ProcRank == 0)
	{
		end_time = MPI_Wtime();
		cout << "\nParallel realisation:" << endl;
		Output(pResultp, Size);
		parallel_time = end_time - start_time;
		cout << "\nTime: " << parallel_time << endl;
		cout << "\nBoost: " << serial_time / parallel_time << endl;
		Validation(pResults, pResultp, Size);
	}
	MPI_Finalize();
	return 0;
}