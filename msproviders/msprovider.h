#ifndef MSPROVIDER_H
#define MSPROVIDER_H

#include "../polarization.h"

#include <casacore/casa/Arrays/Array.h>

#include <casacore/tables/Tables/ArrayColumn.h>

#include <complex>
#include <set>

namespace casacore {
	class MeasurementSet;
}
class MSSelection;

/**
 * The abstract MSProvider class is the base class for classes that read and write the visibilities.
 * An MSProvider knows which rows are selected and doesn't read or write to unselected rows. 
 * It provides the visibilities weighted with the visibility weight and converts the visibilities
 * to a requested polarization.
 * Currently, the @ref ContiguousMS and @ref PartitionedMS classes implement the MSProvider interface.
 */
class MSProvider
{
public:
	friend class MSRowProvider;
	friend class DirectMSRowProvider;

	virtual ~MSProvider() { }
	
	virtual casacore::MeasurementSet &MS() = 0;
	
	virtual size_t RowId() const = 0;
	
	virtual bool CurrentRowAvailable() = 0;
	
	virtual void NextRow() = 0;
	
	virtual void Reset() = 0;
	
	virtual void ReadMeta(double& u, double& v, double& w, size_t& dataDescId) = 0;
	
	virtual void ReadMeta(double& u, double& v, double& w, size_t& dataDescId, size_t& antenna1, size_t& antenna2) = 0;
	
	virtual void ReadData(std::complex<float>* buffer) = 0;
	
	virtual void ReadModel(std::complex<float>* buffer) = 0;
	
	virtual void WriteModel(size_t rowId, std::complex<float>* buffer) = 0;
	
	virtual void ReadWeights(float* buffer) = 0;
	
	virtual void ReadWeights(std::complex<float>* buffer) = 0;
	
	virtual void ReopenRW() = 0;
	
	virtual double StartTime() = 0;
	
	virtual void MakeIdToMSRowMapping(std::vector<size_t>& idToMSRow) = 0;
	
	virtual PolarizationEnum Polarization() = 0;
	
	static std::vector<PolarizationEnum> GetMSPolarizations(casacore::MeasurementSet& ms);
protected:
	static void copyWeightedData(std::complex<float>* dest, size_t startChannel, size_t endChannel, const std::vector<PolarizationEnum>& polsIn, const casacore::Array<std::complex<float>>& data, const casacore::Array<float>& weights, const casacore::Array<bool>& flags, PolarizationEnum polOut);
	
	template<typename NumType>
	static void copyWeights(NumType* dest, size_t startChannel, size_t endChannel, const std::vector<PolarizationEnum>& polsIn, const casacore::Array<std::complex<float>>& data, const casacore::Array<float>& weights, const casacore::Array<bool>& flags, PolarizationEnum polOut);
	
	static void reverseCopyData(casacore::Array<std::complex<float>>& dest, size_t startChannel, size_t endChannel, const std::vector<PolarizationEnum>& polsDest, const std::complex<float>* source, PolarizationEnum polSource);
	
	static void getRowRange(casacore::MeasurementSet& ms, const MSSelection& selection, size_t& startRow, size_t& endRow);
	
	static void getRowRangeAndIDMap(casacore::MeasurementSet& ms, const MSSelection& selection, size_t& startRow, size_t& endRow, const std::set<size_t>& dataDescIdMap, vector<size_t>& idToMSRow);
	
	static void copyRealToComplex(std::complex<float>* dest, const float* source, size_t n)
	{
		const float* end = source + n;
		while(source != end)
		{
			*dest = *source;
			++dest;
			++source;
		}
	}
	
	static void initializeModelColumn(casacore::MeasurementSet& ms);
	
	/**
	 * Make an arraycolumn object for the weight spectrum column if it exists and is valid.
	 * The weight spectrum column is an optional column, the weight column should be used if it doesn't exist.
	 * Moreover, some measurement sets have an empty or invalid sized weight spectrum column; this method
	 * only returns true if the column can be used.
	 */
	static bool openWeightSpectrumColumn(casacore::MeasurementSet& ms, std::unique_ptr<casacore::ROArrayColumn<float>>& weightColumn, const casacore::IPosition& dataColumnShape);
	
	static void expandScalarWeights(const casa::Array<float>& weightScalarArray, casa::Array<float>& weightSpectrumArray)
	{
		casacore::Array<float>::const_contiter src = weightScalarArray.cbegin();
		for(casacore::Array<float>::contiter i=weightSpectrumArray.cbegin(); i!=weightSpectrumArray.cend(); ++i)
		{
			*i = *src;
			++src;
			if(src == weightScalarArray.cend()) src = weightScalarArray.cbegin();
		}
	}
	
	MSProvider() { }
private:
	MSProvider(const MSProvider&) { }
	void operator=(const MSProvider&) { }
};

#endif
