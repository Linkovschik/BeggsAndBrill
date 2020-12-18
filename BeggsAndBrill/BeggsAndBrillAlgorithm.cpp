#include "BeggsAndBrillAlgorithm.h"


double BeggsAndBrillAlgorithm::Execute()
{
	return 0.0;
}

void BeggsAndBrillAlgorithm::ReadInitialParameters(const IReadableBeggsAndBrill& tube)
{
	pipe_parameters.d = tube.GetDiameter();
	pipe_parameters.e = tube.GetRoughness();
	oilFlow.q_normal = tube.GetNormalOilDebit();
	waterFlow.q_normal = tube.GetNormalWaterDebit();
	gasFlow.q_normal = tube.GetNormalGasDebit();
	t_initial = tube.GetTemperature();
	p_initial = tube.GetPressure();
}

void BeggsAndBrillAlgorithm::VolumeCoefficientCount(const IPVTBeggsAndBrill& PVT)
{
	volume_coefficient_oil = PVT.OilVolumeCoefficientCount(this->t_initial, this->p_initial);
	volume_coefficient_water = PVT.WaterVolumeCoefficientCount(this->t_initial, this->p_initial);
	volume_coefficient_gas = PVT.GasVolumeCoefficientCount(this->t_initial, this->p_initial);
}

void BeggsAndBrillAlgorithm::SolubilityCount(const IPVTBeggsAndBrill& PVT)
{
	solubility_oil = PVT.OilSolubilityCount(this->t_initial, this->p_initial);
	solubility_water = PVT.WaterSolubilityCount(this->t_initial, this->p_initial);
}

void BeggsAndBrillAlgorithm::VolumeCount()
{
	oilFlow.q = oilFlow.q_normal*volume_coefficient_oil / 86400;
	waterFlow.q = waterFlow.q_normal*volume_coefficient_water / 86400;
	gasFlow.q = (gasFlow.q_normal - oilFlow.q_normal*solubility_oil - waterFlow.q_normal*solubility_water)*volume_coefficient_gas / 86400;
}

void BeggsAndBrillAlgorithm::ReducedSpeedCount()
{
	double q_l = (oilFlow.q + waterFlow.q);
	double A_p = M_PI * pow(pipe_parameters.d, 2)*1.0e-5 / 4;
	waterFlow.v_reduced = q_l / A_p;
	gasFlow.v_reduced = (gasFlow.q / A_p);
}



