#pragma once
#include "IPVTBeggsAndBrill.h"
class TestBeggsAndBrillPVT :
	public IPVTBeggsAndBrill
{
	// Унаследовано через IPVTBeggsAndBrill
	virtual double OilVolumeCoefficientCount(double temperature, double pressure) const override;
	virtual double WaterVolumeCoefficientCount(double temperature, double pressure) const override;
	virtual double GasVolumeCoefficientCount(double temperature, double pressure) const override;
	virtual double OilSolubilityCount(double temperaturt, double pressure) const override;
	virtual double WaterSolubilityCount(double temperaturt, double pressure) const override;
	virtual double OilViscosityCount(double ANY_PARAMETER) const override;
	virtual double WaterViscosityCount(double ANY_PARAMETER) const override;
	virtual double GasViscosityCount(double ANY_PARAMETER) const override;
	virtual double OilToughnessCount(double ANY_PARAMETER) const override;
	virtual double WaterToughnessCount(double ANY_PARAMETER) const override;
	virtual double GasToughnessCount(double ANY_PARAMETER) const override;
	virtual double OilSurfaceTensionCount(double ANY_PARAMETER) const override;
	virtual double WaterSurfaceTensionCount(double ANY_PARAMETER) const override;
};

