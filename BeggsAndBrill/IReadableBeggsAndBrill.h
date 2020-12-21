#pragma once
#include "Constants.h"
class IReadableBeggsAndBrill
{
public:
	virtual double GetDiameter() const = 0;											// d - диаметр
	virtual double GetRoughness() const = 0;										// e - шероховатость
	virtual double GetNormalOilDebit() const = 0;									//q_o - объёмный дебит нефти
	virtual double GetNormalWaterDebit() const = 0;									//q_w - объёмный дебит воды
	virtual double GetNormalGasDebit() const = 0;									//q_g - объёмный дебит газа
	virtual double GetTemperature(const int& dot_index, const int& dot_count) const = 0;	//t_initial - температура трубы в данной точке
	virtual double GetPressure(const int& dot_index, const int& dot_count) const = 0;		//p_initial - давление трубы в данной точке
	virtual double GetRealInclineAngle() const = 0;									//teta - фактический угол наклона трубы
	virtual Flow::FlowDirection GetFlowDirection() const = 0;						//напраление потока
	virtual TubeStreamParaeters::StreamType GetStramType() const = 0;				//тип потока
	virtual TubeStreamParaeters::TubeRoughType GetTubeRoughType() const = 0;		//тип шероховатости требы
	virtual double GetLength() const = 0;											//длина трубы
	virtual double GetPwh() const = 0;												//неподписанный неизвестный параметр (скорее всего какое-то давление)
	virtual double GetHdyn() const = 0;												//ещё один неописанный параметр
};

