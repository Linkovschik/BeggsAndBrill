#pragma once
#include "IEnvironmentBeggsAndBrill.h"
class TestBeggsAndBrillEnvironment :
	public IEnvironmentBeggsAndBrill
{
	// Унаследовано через IEnvironmentBeggsAndBrill
	virtual double GetNormalToughnessAround() const override;
	virtual double GetNormalTemperatureAround() const override;
	virtual double GetNormalPressureAround() const override;
	virtual double GetCurrentSupercompressibilityCoefficient() const override;
	virtual double GetCurrentTemperatureAround() const override;
	virtual double GetCurrentPressureAround() const override;
};

