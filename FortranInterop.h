#pragma once

/*	Create a new Fortran DLL Project within the same solution.
*	Go to the main C++ Project Settings.
*	1) Add the .dll folder to Linker->General->Additional Library Directories
*	2) Add the .lib filename to Linker -> Input -> Additional Dependencies
*	3) Go to Build Events -> Post Build Events -> Command Line and add xcopy /y /d "C:/MyPath/MyLibrary.dll" "$(OutDir)"
*	You might also have to right-click the solution and manually add it as a dependancy.
*	Warning! Will only update if you right click on the Fortran project and click rebuild!
*/

#ifdef USING_FORTRAN_DLL
#define FORTRAN_API __declspec(dllimport)
#else
#define FORTRAN_API
#endif

extern "C" FORTRAN_API void generate7VertexWormDataFortran_v1(const int* N0, const int* N1, const float* M, const int* n_data);
extern "C" FORTRAN_API void generate7VertexWormDataFortran_v2(const int *N0, const int *N1, const float *M, const int *n_data);
extern "C" FORTRAN_API void generate7VertexWormDataFortran_v3(const int *N0, const int *N1, const float *M, const int *n_data);
extern "C" FORTRAN_API void generate7VertexWormDataFortran_v4(const int *N0, const int *N1, const float *M, const int *n_data);
extern "C" FORTRAN_API void generate7VertexWormDataFortran_v5(const int *N0, const int *N1, const float *M, const int *n_data);
extern "C" FORTRAN_API void generate7VertexWormDataFortran_v6(const int *N0, const int *N1, const float *M, const int *n_data);
extern "C" FORTRAN_API void generate7VertexWormDataFortran_v7(const int *N0, const int *N1, const float *M, const int *n_data);
extern "C" FORTRAN_API void generate7VertexWormDataFortran_v8(const int *N0, const int *N1, const float *M, const int *n_data);
extern "C" FORTRAN_API void generate7VertexWormDataFortran_v9(const int *N0, const int *N1, const float *M, const int *n_data);
extern "C" FORTRAN_API void generate7VertexWormDataFortran_v10(const int *N0, const int *N1, const float *M, const int *n_data);