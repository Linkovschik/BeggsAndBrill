#pragma once
class IEnvironmentBeggsAndBrill
{
public:
	virtual double GetNormalToughnessAround() const= 0;				//плотность в затрубном прос-ве при н.у.
	virtual double GetNormalTemperatureAround() const = 0;			//температура в затрубном прос-ве при н.у.
	virtual double GetNormalPressureAround() const = 0;				//давление в затрубном прос-ве при н.у.

	virtual double GetCurrentSupercompressibilityCoefficient() const = 0;	//коэффициент сверхсжимаемости при текущих условиях
	virtual double GetCurrentTemperatureAround() const = 0;					//температура газа в затрубном пост-ве при текущих условиях
	virtual double GetCurrentPressureAround() const = 0;					//давление газа в затрубном прост-ве при текущих условиях
};

 