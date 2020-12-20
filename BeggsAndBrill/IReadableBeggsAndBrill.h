#pragma once
#include "Constants.h"
class IHydraulicСomputationStrategy;
class IReadableBeggsAndBrill
{
public:
	virtual double GetDiameter() const = 0;									// d - диаметр
	virtual double GetRoughness() const = 0;								// e - шероховатость
	virtual double GetNormalOilDebit() const = 0;							//q_o - объёмный дебит нефти
	virtual double GetNormalWaterDebit() const = 0;							//q_w - объёмный дебит воды
	virtual double GetNormalGasDebit() const = 0;							//q_g - объёмный дебит газа
	virtual double GetTemperature() const = 0;								//t_initial - температура
	virtual double GetPressure() const = 0;									//p_initial - давление
	virtual double GetRealInclineAngle() const = 0;							//teta - фактический угол наклона трубы
	virtual Flow::FlowDirection GetFlowDirection() const = 0;				//напраление потока
	virtual TubeStreamParaeters::StreamType GetStramType() const = 0;		//тип потока
	virtual TubeStreamParaeters::TubeRoughType GetTubeRoughType() const = 0;//тип шероховатости требы
	virtual double GetLength() const = 0;									//длина трубы

	virtual void SetHydraulicComputationStrategy(IHydraulicСomputationStrategy* straegy) = 0; //установка алгоритма-обсервера гидравлического расчёта потока в трубе
private:
	virtual void NotifyObservingHydraulicComputationStrategy() = 0;		//уведомление алгоритма об обновлении состояния трубы
};

