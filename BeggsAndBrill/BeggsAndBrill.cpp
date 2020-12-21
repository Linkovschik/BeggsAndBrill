// BeggsAndBrill.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include "TestBeggsAndBrillPVT.h"
#include "TestBeggsAndBrillEnvironment.h"
#include "TestBeggsAndBrillReadableTube.h"
#include "BeggsAndBrillAlgorithm.h"
int main()
{
	setlocale(LC_ALL, "Russian");
    std::cout << "Hello World!\n";

	TestBeggsAndBrillEnvironment* environment = new TestBeggsAndBrillEnvironment();
	TestBeggsAndBrillPVT* pvt_module = new TestBeggsAndBrillPVT();
	TestBeggsAndBrillReadableTube* tube = new TestBeggsAndBrillReadableTube();

	double dot_count = 10;
	TubeStreamParaeters::NUnknownParameter n_unknown_parameters = TubeStreamParaeters::NUnknownParameter::First;
	
	BeggsAndBrillAlgorithm* algorithmInstance = BeggsAndBrillAlgorithm::GetInstance(environment,tube,pvt_module,dot_count,n_unknown_parameters);
	algorithmInstance->Execute();

	delete environment;
	delete pvt_module;
	delete tube;
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
