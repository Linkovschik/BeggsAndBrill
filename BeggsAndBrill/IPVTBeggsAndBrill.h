#pragma once
class IPVTBeggsAndBrill
{
public:
	virtual double OilVolumeCoefficientCount(double temperature, double pressure) const = 0;
	virtual double WaterVolumeCoefficientCount(double temperature, double pressure) const = 0;
	virtual double GasVolumeCoefficientCount(double temperature, double pressure) const = 0;

	virtual double OilSolubilityCount(double temperaturt, double pressure) const = 0;
	virtual double WaterSolubilityCount(double temperaturt, double pressure) const = 0;

	virtual double OilViscosityCount(double ANY_PARAMETER) const = 0;
	virtual double WaterViscosityCount(double ANY_PARAMETER) const = 0;
	virtual double GasViscosityCount(double ANY_PARAMETER) const = 0;

	virtual double OilToughnessCount(double ANY_PARAMETER) const = 0;
	virtual double WaterToughnessCount(double ANY_PARAMETER) const = 0;
	virtual double GasToughnessCount(double ANY_PARAMETER) const = 0;

	virtual double OilSurfaceTensionCount(double ANY_PARAMETER) const = 0;
	virtual double WaterSurfaceTensionCount(double ANY_PARAMETER) const = 0;
};

