// SGM_Header.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include "SGM_header.h"
using namespace std;
using namespace SGM;



int main()
{
	//TODO here code
	cout << GeneratePassword(8, 2, 5, 1) << endl;



	system("pause");
	return 0;
}

