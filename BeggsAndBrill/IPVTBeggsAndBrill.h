#pragma once
class IPVTBeggsAndBrill
{
public:
	virtual double OilVolumeCoefficientCount(double temperature, double pressure) const = 0;
	virtual double WaterVolumeCoefficientCount(double temperature, double pressure) const = 0;
	virtual double GasVolumeCoefficientCount(double temperature, double pressure) const = 0;

	virtual double OilSolubilityCount(double temperaturt, double pressure) const = 0;
	virtual double WaterSolubilityCount(double temperaturt, double pressure) const = 0;
};

