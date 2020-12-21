#pragma once
#include "Constants.h"
class IReadableBeggsAndBrill
{
public:
	virtual double GetDiameter() const = 0;											// d - �������
	virtual double GetRoughness() const = 0;										// e - �������������
	virtual double GetNormalOilDebit() const = 0;									//q_o - �������� ����� �����
	virtual double GetNormalWaterDebit() const = 0;									//q_w - �������� ����� ����
	virtual double GetNormalGasDebit() const = 0;									//q_g - �������� ����� ����
	virtual double GetTemperature(const int& dot_index, const int& dot_count) const = 0;	//t_initial - ����������� ����� � ������ �����
	virtual double GetPressure(const int& dot_index, const int& dot_count) const = 0;		//p_initial - �������� ����� � ������ �����
	virtual double GetRealInclineAngle() const = 0;									//teta - ����������� ���� ������� �����
	virtual Flow::FlowDirection GetFlowDirection() const = 0;						//���������� ������
	virtual TubeStreamParaeters::StreamType GetStramType() const = 0;				//��� ������
	virtual TubeStreamParaeters::TubeRoughType GetTubeRoughType() const = 0;		//��� ������������� �����
	virtual double GetLength() const = 0;											//����� �����
	virtual double GetPwh() const = 0;												//������������� ����������� �������� (������ ����� �����-�� ��������)
	virtual double GetHdyn() const = 0;												//��� ���� ����������� ��������
};

