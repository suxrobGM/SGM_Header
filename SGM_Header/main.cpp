// SGM_Header.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include "SGM_header.h"
using namespace std;
using namespace SGM;



int main()
{
	//TODO here code
	
	Graph graph(7);
	graph.AddEdge(0, 1);
	graph.AddEdge(0, 2);
	graph.AddEdge(0, 3);
	graph.AddEdge(1, 4);
	graph.AddEdge(1, 5);
	graph.AddEdge(2, 6);

	graph.DepthFirstSearch(0);

	system("pause");
	return 0;
}

