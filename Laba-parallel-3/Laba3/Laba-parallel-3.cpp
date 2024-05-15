#include <iostream>
#include <fstream>
#include <stdio.h>
#include <conio.h>
#include <random>
#include <chrono>
#include <string>
#include<windows.h>

#include "mpi.h"

using namespace std;

#define MASTER_TO_SLAVE_TAG 1 
#define SLAVE_TO_MASTER_TAG 4


void FillMat(int n, int*** Mat) {
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			(*Mat)[i][j] = rand() % (100 + 100 + 1) - 100;
		}
	}
}

int** Reader(const std::string filename) {
	std::fstream fin;
	fin.open(filename);

	int n;
	fin >> n;

	int** Mat = new int* [n];
	for (int i = 0; i < n; i++)
	{
		Mat[i] = new int[n];
	}

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++) {
			fin >> Mat[i][j];
		}
	}

	return Mat;
}

void Writer(int n, int** Mat, const std::string filename)
{
	std::fstream fout;
	fout.open(filename, std::ofstream::out | std::ofstream::trunc);
	fout << n << std::endl;

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			fout << Mat[i][j] << ' ';
		}
		fout << std::endl;
	}
}

int** resize_matrix(int size)
{
	int** matrix = new int* [size];
	for (int i = 0; i < size; ++i)
		matrix[i] = new int[size];
	return matrix;
}

int main(int argc, char* argv[]) {

	srand(time(0));

	CreateDirectory(L"res_matrix", NULL);
	CreateDirectory(L"1_matrix", NULL);
	CreateDirectory(L"2_matrix", NULL);
	CreateDirectory(L"times", NULL);

	int** Mat1 = NULL;
	int** Mat2 = NULL;
	int** Mat_res = NULL;

	int rank, size;
	double ti[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	std::fstream fout;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (rank == 0)
	{
		fout.open(std::string("times/time_thread_") + std::to_string(size) + std::string(".txt"), std::ofstream::out | std::ofstream::trunc);
		std::cout << "size = " << size << endl;
	}

	for (int n = 200; n <= 2000; n += 200)
	{
		if (rank == 0)
		{
			
			std::cout << n << std::endl;
		}
		Mat1 = resize_matrix(n);
		Mat2 = resize_matrix(n);
		Mat_res = resize_matrix(n);

		MPI_Barrier(MPI_COMM_WORLD);
		double wtime = MPI_Wtime();

		for (int t = 1; t <= 10; t++)
		{
			if (rank == 0)
			{
				FillMat(n, &Mat1);
				FillMat(n, &Mat2);
				//std::cout << '\t' << i;
			}
			MPI_Bcast(&(Mat2[0][0]), n * n, MPI_INT, 0, MPI_COMM_WORLD);

			int rows_per_process = n / size;
			int rows_remaining = n % size;
			int start_row = rank * rows_per_process + min(rank, rows_remaining);
			int end_row = start_row + rows_per_process + (rank < rows_remaining ? 1 : 0);

			for (int i = start_row; i < end_row; i++)
			{
				for (int j = 0; j < n; j++)
				{
					Mat_res[i][j] = 0;
					for (int k = 0; k < n; k++)
						Mat_res[i][j] += Mat2[i][k] * Mat1[k][j];
				}
			}

			MPI_Barrier(MPI_COMM_WORLD);
			if (rank == 0)
				ti[t - 1] = MPI_Wtime() - wtime;
		}

		MPI_Barrier(MPI_COMM_WORLD);

		if (rank == 0)
		{
			Writer(n, Mat_res, std::string("res_matrix/res_matrix") + std::to_string(n) + std::string(".txt"));
			Writer(n, Mat1, std::string("1_matrix/1_matrix") + std::to_string(n) + std::string(".txt"));
			Writer(n, Mat2, std::string("2_matrix/2_matrix") + std::to_string(n) + std::string(".txt"));
			fout << n << ';';
			for (int i = 0; i < 10; i++)
			{
				fout << ti[i] << ';';
				std::cout << ti[i] << ';';
			}
			fout << std::endl;
			std::cout << std::endl;

		}
		for (int i = 0; i < n; i++)
		{
			delete[] Mat_res[i];
			delete[] Mat1[i];
			delete[] Mat2[i];
		}

		delete[] Mat_res;
		delete[] Mat1;
		delete[] Mat2;

		Mat_res = NULL;
		Mat1 = NULL;
		Mat2 = NULL;

	}

	MPI_Finalize();
	fout.close();


	return 0;
}