#pragma once
#define CUDA_API __declspec(dllexport)

#ifdef DLL_EXPORT
	__declspec(dllexport) void cuda7VertexWorm(int N0, int N1, const double M_start, const int n_data, const int n_mass);
#else
	__declspec(dllimport) void cuda7VertexWorm(int N0, int N1, const double M_start, const int n_data, const int n_mass);
#endif