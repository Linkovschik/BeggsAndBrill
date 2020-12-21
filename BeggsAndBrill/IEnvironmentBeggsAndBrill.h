#pragma once
class IEnvironmentBeggsAndBrill
{
public:
	virtual double GetNormalToughnessAround() const= 0;				//��������� � ��������� ����-�� ��� �.�.
	virtual double GetNormalTemperatureAround() const = 0;			//����������� � ��������� ����-�� ��� �.�.
	virtual double GetNormalPressureAround() const = 0;				//�������� � ��������� ����-�� ��� �.�.

	virtual double GetCurrentSupercompressibilityCoefficient() const = 0;	//����������� ���������������� ��� ������� ��������
	virtual double GetCurrentTemperatureAround() const = 0;					//����������� ���� � ��������� ����-�� ��� ������� ��������
	virtual double GetCurrentPressureAround() const = 0;					//�������� ���� � ��������� �����-�� ��� ������� ��������
};

 