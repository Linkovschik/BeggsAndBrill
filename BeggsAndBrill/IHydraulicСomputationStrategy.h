#pragma once
class IReadableBeggsAndBrill;
class IHydraulic—omputationStrategy {
public:
	virtual void Execute(const IReadableBeggsAndBrill & tube) = 0;
private:
	virtual void Update() = 0;
};

