#include "BeggsAndBrillAlgorithm.h"


double BeggsAndBrillAlgorithm::Execute(const IReadableBeggsAndBrill& tube)
{
	PipeParameters pipe_parameters;
	ReadPipeParameters(tube,
		pipe_parameters.diameter, pipe_parameters.roughness,
		pipe_parameters.temperature, pipe_parameters.pressure);

	FlowInfo oilFlow{ "oil" };
	FlowInfo waterFlow{ "water" };
	FlowInfo gasFlow{ "gas" };

	ReadFlowNormalDebits(tube,
		oilFlow.volume_debit_normal,
		waterFlow.volume_debit_normal,
		gasFlow.volume_debit_normal);

	return 0.0;
}

void BeggsAndBrillAlgorithm::ReadPipeParameters(const IReadableBeggsAndBrill & tube,
	double & out_diameter, double & out_roughness, double & out_temperature, double & out_pressure)
{
	out_diameter = tube.GetDiameter();
	out_roughness = tube.GetRoughness();
	out_temperature = tube.GetTemperature();
	out_pressure = tube.GetPressure();
}

void BeggsAndBrillAlgorithm::ReadFlowNormalDebits(const IReadableBeggsAndBrill & tube, 
	double & out_oil_normal_debit, double & out_water_normal_debit, double & out_gas_normal_debit)
{
	out_oil_normal_debit = tube.GetNormalOilDebit();
	out_water_normal_debit = tube.GetNormalWaterDebit();
	out_gas_normal_debit = tube.GetNormalGasDebit();
}

void BeggsAndBrillAlgorithm::VolumeCoefficientCount(const IPVTBeggsAndBrill & PVT,
	const PipeParameters& pipe_parameters, 
	double & out_volume_coefficient_oil, double & out_volume_coefficient_water, double & out_volume_coefficient_gas)
{
	const double& temperature = pipe_parameters.temperature;
	const double& pressure = pipe_parameters.pressure;
	out_volume_coefficient_oil = PVT.OilVolumeCoefficientCount(temperature,pressure);
	out_volume_coefficient_water = PVT.WaterVolumeCoefficientCount(temperature, pressure);
	out_volume_coefficient_gas = PVT.GasVolumeCoefficientCount(temperature, pressure);
}

void BeggsAndBrillAlgorithm::SolubilityCount(const IPVTBeggsAndBrill & PVT, 
	const PipeParameters& pipe_parameters,
	double & out_gas_sollubility_oil, double & out_gas_sollubility_water)
{
	const double& temperature = pipe_parameters.temperature;
	const double& pressure = pipe_parameters.pressure;
	out_gas_sollubility_oil = PVT.OilSolubilityCount(temperature, pressure);
	out_gas_sollubility_water = PVT.WaterSolubilityCount(temperature, pressure);
}

void BeggsAndBrillAlgorithm::VolumeDebitCount(const FlowInfo& oilFlow, const FlowInfo& waterFlow, const FlowInfo& gasFlow,
	double& out_volume_debit_oil, double& out_volume_debit_water, double& out_volume_debit_gas)
{
	out_volume_debit_oil = oilFlow.volume_debit_normal*oilFlow.volume_coefficient / 86400;
	out_volume_debit_water = waterFlow.volume_debit_normal*waterFlow.volume_coefficient / 86400;
	out_volume_debit_gas = (gasFlow.volume_debit_normal - oilFlow.volume_debit_normal*oilFlow.gas_solubility - waterFlow.volume_debit_normal*waterFlow.gas_solubility)*gasFlow.volume_coefficient / 86400;
}

double BeggsAndBrillAlgorithm::ApCount(const double& diameter)
{
	return M_PI * pow(diameter, 2)*1.0e-5 / 4;
}

void BeggsAndBrillAlgorithm::ReducedSpeedCount(const FlowInfo & oilFlow, const FlowInfo & waterFlow, const FlowInfo & gasFlow,
	const PipeParameters& pipe_parameters,
	double& out_reduced_speed_water, double& out_reduced_speed_gas)
{
	double q_l = (oilFlow.volume_debit + waterFlow.volume_debit);
	double A_p = ApCount(pipe_parameters.diameter);
	out_reduced_speed_water = q_l / A_p;
	out_reduced_speed_gas = (gasFlow.volume_debit / A_p);
}

void BeggsAndBrillAlgorithm::MixtureSpeedCount(const FlowInfo & waterFlow, const FlowInfo & gasFlow, double & out_mixture_speed)
{
	out_mixture_speed = waterFlow.speed_reduced + gasFlow.speed_reduced;
}

void BeggsAndBrillAlgorithm::VolumeContainityCount(const FlowInfo & oilFlow, const FlowInfo & waterFlow, const FlowInfo & gasFlow, double & out_lambda)
{
	double q_l = oilFlow.volume_debit + waterFlow.volume_debit;
	out_lambda = q_l / (q_l + gasFlow.volume_debit);
}

void BeggsAndBrillAlgorithm::ViscosityCount(const IPVTBeggsAndBrill& PVT,
	const FlowInfo & oilFlow, const FlowInfo & waterFlow, const FlowInfo & gasFlow,
	double& out_viscosity_oil, double& out_viscosity_water, double& out_viscosity_gas)
{
	out_viscosity_oil = PVT.OilViscosityCount(42);
	out_viscosity_water = PVT.WaterViscosityCount(42);
	out_viscosity_gas = PVT.GasViscosityCount(42);
}

void BeggsAndBrillAlgorithm::ToughnessCount(const IPVTBeggsAndBrill & PVT, 
	const FlowInfo & oilFlow, const FlowInfo & waterFlow, const FlowInfo & gasFlow, 
	double & out_toughness_oil, double & out_toughness_water, double & out_toughness_gas)
{
	out_toughness_oil = PVT.OilToughnessCount(42);
	out_toughness_water = PVT.WaterToughnessCount(42);
	out_toughness_gas = PVT.GasToughnessCount(42);
}

void BeggsAndBrillAlgorithm::SurfaceTensionCount(const IPVTBeggsAndBrill & PVT,
	const FlowInfo & oilFlow, const FlowInfo & waterFlow, const FlowInfo & gasFlow, 
	double & out_surface_tension_oil, double & out_surface_tension_water)
{
	out_surface_tension_oil = PVT.OilSurfaceTensionCount(42);
	out_surface_tension_water = PVT.WaterSurfaceTensionCount(42);
}

void BeggsAndBrillAlgorithm::MixtureViscosity_Toughness_SurfaceTension_Count(const FlowInfo & oilFlow, const FlowInfo & waterFlow, const FlowInfo & gasFlow, 
	const MixtureInfo& mixture,
	double & out_mixture_viscosity_secondphase, double & out_mixture_toughness_secondphase,
	double & out_mixture_viscosity_liquid, double & out_mixture_toughness_liquid, double & out_mixture_surfacetansion_liquid)
{
	double oil_proportion = oilFlow.volume_debit / (oilFlow.volume_debit + waterFlow.volume_debit);	//f_o - доля нефти в жидкой фазе
	out_mixture_viscosity_liquid = 1.0e-2 * (oilFlow.viscosity*oil_proportion + waterFlow.viscosity*(1 - oil_proportion));
	out_mixture_viscosity_secondphase = out_mixture_viscosity_liquid * mixture.lambda + 1.0e-2 * gasFlow.viscosity*(1 - mixture.lambda);
	out_mixture_toughness_liquid = oilFlow.toughness*oil_proportion + waterFlow.toughness*(1 - oil_proportion);
	out_mixture_toughness_secondphase = out_mixture_toughness_liquid * mixture.lambda + gasFlow.toughness*(1 - mixture.lambda);
	out_mixture_surfacetansion_liquid = oilFlow.surface_tension*oil_proportion + waterFlow.surface_tension*(1 - oil_proportion);
}

void BeggsAndBrillAlgorithm::FroudNumberCount(const MixtureInfo & mixture, const PipeParameters & pipe_parameters, double & out_froud_number)
{
	out_froud_number = pow(mixture.speed, 2) / PhysicConstants::g*pipe_parameters.diameter * 1.0e-2;
}






