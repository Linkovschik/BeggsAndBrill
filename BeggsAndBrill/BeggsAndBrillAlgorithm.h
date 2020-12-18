#pragma once
#define _USE_MATH_DEFINES
#include <string>
#include <iostream>
#include <cmath>
#include "Constants.h"
#include "IReadableBeggsAndBrill.h"
#include "IPVTBeggsAndBrill.h"


class BeggsAndBrillAlgorithm
{
public:
	//Execute - шаблонный метод, не должен переопределяться
	double Execute(const IReadableBeggsAndBrill& tube);
protected:
	std::string file_name;

	//конструктор
	BeggsAndBrillAlgorithm();

	//деструктор
	virtual ~BeggsAndBrillAlgorithm();

	//параметры трубы (полученные извне)
	struct PipeParameters {

		double diameter;						//диаметр (d)
		double roughness;						//шероховатость (e)
		double temperature;						//температура в указанной точке (T_initial)
		double pressure;						//давление в указанной точке (P_initial)
	};

	//информация о потоке (полученные извне и расчитанная в алгоритме)
	struct FlowInfo {

		std::string name;						//название (нефть, вода, газ)
		double volume_debit;					//объёмный дебит (q)
		double volume_debit_normal;				//объёмный дебит при н.у. (q_normal)
		double volume_coefficient;				//объёмный коэффициент (B_o и B_w у нефти и воды соотв.)
		double gas_solubility;					//растворимость газа (R_s, R_sw в нефти и воде соотв.)
		double speed_reduced;					//приведённая скорость (v_sl и v_sg у жикости и газа соотв.)
		double viscosity;						//вязкость (mu)
		double toughness;						//плотность (ro)
		double surface_tension;					//поверхностное натяжение (delt)
	};

	//информация о других данных, расчитанных в алгоритме
	struct MixtureInfo {

		double speed;								//скорость смеси (v_m)
		double lambda;								//объёмное содержание жидкости (lambda)
		double viscosity_liquid;					//вязкость смеси (mu)
		double toughness_liquid;					//плотность смеси (ro)
		double surface_tension_liquid;				//поверхностное натяжение смеси (delt)
		double viscosity_secondphase;				//вязкость смеси  (двухфазн.) (mu)
		double toughness_secondphase;				//плотность смеси (двухфазн) (ro) 
		double froud_number;						//число Фруда (N_fr)
		double liquid_speed_index;					//показатель скорости жидкости (N_lv)
	};
private:

	//считывание начальных значений (параметров трубы)
	virtual void ReadPipeParameters(const IReadableBeggsAndBrill& tube, 
		double& out_diameter, double& out_roughness, double& out_temperature, double& out_pressure);

	//считывание параметров потока
	virtual void ReadFlowNormalDebits(const IReadableBeggsAndBrill& tube,
		double& out_oil_normal_debit, double& out_water_normal_debit, double& out_gas_normal_debit);

	//расчет объемного коэффициента нефти и воды
	virtual void VolumeCoefficientCount(const IPVTBeggsAndBrill& PVT,
		const PipeParameters& pipe_parameters,
		double& out_volume_coefficient_oil, double& out_volume_coefficient_water, double& out_volume_coefficient_gas);

	//расчет растворимости газа в нефти и в воде
	virtual void SolubilityCount(const IPVTBeggsAndBrill& PVT,
		const PipeParameters& pipe_parameters,
		double& out_gas_sollubility_oil, double& out_gas_sollubility_water);

	//расчет объемных дебитов нефти, воды и газа в заданной точке
	virtual void VolumeDebitCount(const FlowInfo& oilFlow, const FlowInfo& waterFlow, const FlowInfo& gasFlow,
		double& out_volume_debit_oil, double& out_volume_debit_water, double& out_volume_debit_gas);

	double ApCount(const double& diameter);

	//вычисление приведенной скорости жидкости и газа
	virtual void ReducedSpeedCount(const FlowInfo& oilFlow, const FlowInfo& waterFlow, const FlowInfo& gasFlow,
		const PipeParameters& pipe_parameters,
		double& out_reduced_speed_water, double& out_reduced_speed_gas);

	//вычисление скорости смеси
	virtual void MixtureSpeedCount(const FlowInfo& waterFlow, const FlowInfo& gasFlow, 
		double& out_mixture_speed);

	//вычисление объемного содержания жидкости
	virtual void VolumeContainityCount(const FlowInfo& oilFlow, const FlowInfo& waterFlow, const FlowInfo& gasFlow,
		double& out_lambda);

	//расчёт вязкости нефти, воды, газа
	virtual void ViscosityCount(const IPVTBeggsAndBrill& PVT,
		const FlowInfo & oilFlow, const FlowInfo & waterFlow, const FlowInfo & gasFlow,
		double& out_viscosity_oil, double& out_viscosity_water, double& out_viscosity_gas);

	//расчёт плотности нефти, воды, газа
	virtual void ToughnessCount(const IPVTBeggsAndBrill& PVT,
		const FlowInfo & oilFlow, const FlowInfo & waterFlow, const FlowInfo & gasFlow,
		double& out_toughness_oil, double& out_toughness_water, double& out_toughness_gas);

	//расчёт поверхностного натяжения нефти и воды
	virtual void SurfaceTensionCount(const IPVTBeggsAndBrill& PVT,
		const FlowInfo & oilFlow, const FlowInfo & waterFlow, const FlowInfo & gasFlow,
		double& out_surface_tension_oil, double& out_surface_tension_water);

	//Вычисляем вязкость, плотность и поверхностное натяжение смеси
	virtual void MixtureViscosity_Toughness_SurfaceTension_Count(const FlowInfo & oilFlow, const FlowInfo & waterFlow, const FlowInfo & gasFlow,
		const MixtureInfo& mixture,
		double& out_mixture_viscosity_secondphase, double& out_mixture_toughness_secondphase,
		double& out_mixture_viscosity_liquid, double& out_mixture_toughness_liquid, double& out_mixture_surfacetansion_liquid);

	//расчет числа Фруда на основе скорости смеси и диаметра трубы
	virtual void FroudNumberCount(const MixtureInfo& mixture, const PipeParameters& pipe_parameters,
		double& out_froud_number);

	//расчет показателя скорости жидкости
	virtual
};

