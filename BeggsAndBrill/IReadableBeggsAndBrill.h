#pragma once
class IReadableBeggsAndBrill
{
public:
	virtual double GetDiameter() const = 0;					// d - �������
	virtual double GetRoughness() const = 0;				// e - �������������
	virtual double GetNormalOilDebit() const = 0;			//q_o - �������� ����� �����
	virtual double GetNormalWaterDebit() const = 0;			//q_w - �������� ����� ����
	virtual double GetNormalGasDebit() const = 0;			//q_g - �������� ����� ����
	virtual double GetTemperature() const = 0;				//t_initial - �����������
	virtual double GetPressure() const = 0;					//p_initial - ��������
};

