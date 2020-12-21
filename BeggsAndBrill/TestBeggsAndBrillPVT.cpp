#include "TestBeggsAndBrillPVT.h"

double TestBeggsAndBrillPVT::OilVolumeCoefficientCount(double temperature, double pressure) const
{
	return 1.0;
}

double TestBeggsAndBrillPVT::WaterVolumeCoefficientCount(double temperature, double pressure) const
{
	return 2.0;
}

double TestBeggsAndBrillPVT::GasVolumeCoefficientCount(double temperature, double pressure) const
{
	return 3.0;
}

double TestBeggsAndBrillPVT::OilSolubilityCount(double temperaturt, double pressure) const
{
	return 13.0;
}

double TestBeggsAndBrillPVT::WaterSolubilityCount(double temperaturt, double pressure) const
{
	return 15.0;
}

double TestBeggsAndBrillPVT::OilViscosityCount(double ANY_PARAMETER) const
{
	return 21.0;
}

double TestBeggsAndBrillPVT::WaterViscosityCount(double ANY_PARAMETER) const
{
	return 31.0;
}

double TestBeggsAndBrillPVT::GasViscosityCount(double ANY_PARAMETER) const
{
	return 41.0;
}

double TestBeggsAndBrillPVT::OilToughnessCount(double ANY_PARAMETER) const
{
	return 4.0;
}

double TestBeggsAndBrillPVT::WaterToughnessCount(double ANY_PARAMETER) const
{
	return 2.0;
}

double TestBeggsAndBrillPVT::GasToughnessCount(double ANY_PARAMETER) const
{
	return 0.3;
}

double TestBeggsAndBrillPVT::OilSurfaceTensionCount(double ANY_PARAMETER) const
{
	return 1.1;
}

double TestBeggsAndBrillPVT::WaterSurfaceTensionCount(double ANY_PARAMETER) const
{
	return 0.5;
}
