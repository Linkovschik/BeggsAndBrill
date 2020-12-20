#pragma once
#define _USE_MATH_DEFINES
#include <string>
#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
#include "IHydraulicСomputationStrategy.h"
#include "Constants.h"
#include "IReadableBeggsAndBrill.h"
#include "IPVTBeggsAndBrill.h"


class BeggsAndBrillAlgorithm : public IHydraulicСomputationStrategy
{
public:
	//Execute - шаблонный метод, не должен переопределяться
	virtual void Execute(const IReadableBeggsAndBrill& tube) override;
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
		double reduced_speed;					//приведённая скорость (v_sg у газа)
		double viscosity;						//вязкость (mu)
		double toughness;						//плотность (ro)
		double surface_tension;					//поверхностное натяжение (delt)
	};

	//информация о других данных, расчитанных в алгоритме
	struct MixtureInfo {

		double reduced_speed;					//приведённая скорость (v_sl у жидкости)
		double speed;							//скорость смеси (v_m)
		double lambda;							//объёмное содержание жидкости (lambda)
		double viscosity_liquid;				//вязкость смеси (mu_l)
		double toughness_liquid;				//плотность смеси (ro_l)
		double surface_tension_liquid;			//поверхностное натяжение смеси (delt)
		double viscosity_secondphase;			//вязкость смеси  (двухфазн.) (mu_n)
		double toughness_secondphase;			//плотность смеси (двухфазн) (ro_n) 
		double froud_number;					//число Фруда (N_fr)
		double liquid_speed_index;				//показатель скорости жидкости (N_lv)
	};

	//гранцы потока
	struct FlowBorders {
		double L1, L2, L3, L4;
	};

	//состояние потока
	class FlowState {
	public:
		//конструктор
		FlowState(const MixtureInfo* ptr_mixture_info);
		//деструктор
		virtual ~FlowState();

		//возвращает объемное содержание жидкости с поправкой на угол наклона
		double GetVolumetric_content_liquid_angled(const IReadableBeggsAndBrill& tube);

		//возвращает true, если исходные данные подохдят для этого режима
		virtual bool ModeChecking(const FlowBorders& borders) = 0;
	protected:
		const MixtureInfo* ptr_mixture;							//указатель на структуру "смесь"
		double a, b, c;											//параметры для расчёта объёмного содержания жидкости
		double e_down_up, f_down_up, g_down_up, h_down_up;		//параметры для расчёта некоего коэффициента C (снизу_вверх)
		double e_up_down, f_up_down, g_up_down, h_up_down;		//параметры для расчёта некоего коэффициента C (сверху_вниз)
	private:
		//расчёт объёмного содержания жикости без поправки на угол
		virtual double GetVolumetric_content_liquid_not_angled();
		//возвращает значение коэффициента C
		virtual double GetC_Coefficient(const Flow::FlowDirection& flow_direction);
		//расчёт поправочного коэффициента на угол наклона трубы (пси)
		virtual double GetAngleCorrectionCoefficient(const IReadableBeggsAndBrill & tube);
	};

	//разделённый режим
	class DividedFlowState : public FlowState {
	public:
		//конструктор
		DividedFlowState(const MixtureInfo* ptr_mixture_info);
		//деструктор
		virtual ~DividedFlowState();

		virtual bool ModeChecking(const FlowBorders& borders) override;
	};

	//прерывистый режим
	class IntermittentFlowState : public FlowState {
	public:
		//конструктор
		IntermittentFlowState(const MixtureInfo* ptr_mixture_info);
		//деструктор
		virtual ~IntermittentFlowState();

		virtual bool ModeChecking(const FlowBorders& borders) override;
	};

	//распределённый режим
	class DistributedFlowState : public FlowState {
	public:
		//конструктор
		DistributedFlowState(const MixtureInfo* ptr_mixture_info);
		//деструктор
		virtual ~DistributedFlowState();

		virtual bool ModeChecking(const FlowBorders& borders) override;

		virtual double GetC_Coefficient(const Flow::FlowDirection& flow_direction) override;
		virtual double GetAngleCorrectionCoefficient(const IReadableBeggsAndBrill & tube) override;
	};
private:
	virtual void Update() override;

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
		double& out_reduced_speed_liquid, double& out_reduced_speed_gas);

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
	virtual void LiquidIndexCount(const MixtureInfo& mixture,
		double& out_liquid_speed_index);

	//расчёт границ режимом потока
	virtual void FlowBordersCount(const MixtureInfo& mixture,
		double& out_L1, double& out_L2, double& out_L3, double& out_L4);

	//объемное содержание жидкости с поправкой на угол наклона
	virtual void VolumeCorrectionLiquidCount(const IReadableBeggsAndBrill& tube, const MixtureInfo& mixture, const FlowBorders& borders,
		double& out_volume_correction_liquid);

	//использовать поправку для объемного содержания жидкости (Payne)
	virtual void PayneVolumeCorrectionCount(const IReadableBeggsAndBrill& tube, const double& volume_correction_liquid,
		double &out_payne_corrected_vlume_correction_liquid);

	//Вычисление числа Рейнольдса
	virtual void RaynoldsNumberCount(const MixtureInfo& mixture, const PipeParameters& pipe_parameters,
		double& out_raynolds_number);

	//на основе шероховатости и числа Рейнольдса определяем значение нормирующего коэффициента трения
	virtual void NormalizingFrictionCoefficientCount(const IReadableBeggsAndBrill& tube, const PipeParameters& pipe_parameters, 
		const double& raynolds_number,
		double& out_normalizing_friction_coefficient);

	//вычисляем параметр y
	double GetY(const double& lambda, const double& volume_correction_liquid);

	//вычисляем параметр s
	double GetS(const double& y);

	//вычисляем коэффициент трения двухфазного потока
	virtual void FrictionCoefficientCount(const MixtureInfo& mixture, const double& volume_correction_liquid,
		const double& normilizing_friction_coefficient,
		double& out_friction_coefficient);

	//вычисляем градиент давления
	virtual void PressureGradientCount(const IReadableBeggsAndBrill & tube, 
		const PipeParameters& pip_parameters, const MixtureInfo& mixture, const FlowInfo& gasFlow,
		const double& friction_coefficient, const double& volume_correction_liquid,
		double& out_pressure_gradient);
}; 
