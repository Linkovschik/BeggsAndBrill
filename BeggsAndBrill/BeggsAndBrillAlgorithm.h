#pragma once
#define _USE_MATH_DEFINES
#include <string>
#include <iostream>
#include <cmath>
#include "IReadableBeggsAndBrill.h"
#include "IPVTBeggsAndBrill.h"


class BeggsAndBrillAlgorithm
{
public:
	double Execute();
protected:
	std::string file_name;

	BeggsAndBrillAlgorithm();
	virtual ~BeggsAndBrillAlgorithm();
private:
	//параметры трубы
	struct PipeParameters {

		double d;			//диаметр
		double e;			//шероховатость

	} pipe_parameters;

	//информация о потоке
	struct FlowInfo {
		std::string name;	//название (нефть, вода, газ)
		double q;			//объёмный дебит
		double q_normal;	//объёмный дебит при н.у.
		double v_reduced;	//приведённая скорость (отсутствует у нефти = 0.0) 
							//(v_sl и v_sg у жикости и газа соответственно)

	};
	FlowInfo oilFlow{ "oil" };
	FlowInfo waterFlow{ "water" };
	FlowInfo gasFlow{ "gas" };


	double t_initial;		// температура в заданной точке
	double p_initial;		//давление в заданной точке

	double volume_coefficient_oil;				//объёмный коэффициент нефти (B_o)
	double volume_coefficient_water;			//объёмный коэффициент воды (B_w)
	double volume_coefficient_gas;				//объёмный коэффициент воды (B_w)

	double solubility_oil;			//растворимость газа в нефти (R_s)
	double solubility_water;		//растворимость газа в воде (R_sw)

	double mixture_speed;			//скорость смеси (v_m)

	//считывание начальных значений (параметров трубы и объёмных дебитов потока)
	//(не переопределяемое)
	virtual void ReadInitialParameters(const IReadableBeggsAndBrill& tube) final;

	//расчет объемного коэффициента нефти и воды
	virtual void VolumeCoefficientCount(const IPVTBeggsAndBrill& PVT);

	//расчет растворимости газа в нефти и в воде
	virtual void SolubilityCount(const IPVTBeggsAndBrill& PVT);

	//расчет объемных дебитов нефти, воды и газа в заданной точке
	virtual void VolumeCount();

	//вычисление приведенной скорости жидкости и газа
	//вычисление скорости смеси
	//Вычисление объемного содержания жидкости
	void ReducedSpeedCount();

	
};

