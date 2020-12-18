#pragma once
class IReadableBeggsAndBrill
{
public:
	virtual double GetDiameter() const = 0;					// d - диаметр
	virtual double GetRoughness() const = 0;				// e - шероховатость
	virtual double GetNormalOilDebit() const = 0;			//q_o - объёмный дебит нефти
	virtual double GetNormalWaterDebit() const = 0;			//q_w - объёмный дебит воды
	virtual double GetNormalGasDebit() const = 0;			//q_g - объёмный дебит газа
	virtual double GetTemperature() const = 0;				//t_initial - температура
	virtual double GetPressure() const = 0;					//p_initial - давление
};

