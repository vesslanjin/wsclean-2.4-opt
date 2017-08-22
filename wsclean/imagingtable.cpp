#include "imagingtable.h"
#include "logger.h"
#include <map>

ImagingTableEntry::ImagingTableEntry() :
	index(0),
	lowestFrequency(0.0),
	highestFrequency(0.0),
	bandStartFrequency(0.0),
	bandEndFrequency(0.0),
	polarization(Polarization::StokesI),
	outputChannelIndex(0),
	outputIntervalIndex(0),
	msData(),
	squaredDeconvolutionIndex(0),
	joinedGroupIndex(0),
	imageCount(0),
	tmpFilePrefix(),
	imageWeight(0.0)
{
}

ImagingTable ImagingTable::GetIndependentGroup(size_t index) const
{
	ImagingTable table;
	const std::vector<ImagingTableEntryPtr>& entries = _independentGroupLookup[index];
	table._entries = entries;
	table.Update();
	return table;
}

ImagingTable ImagingTable::GetSquaredGroup(size_t index) const
{
	ImagingTable table;
	const std::vector<ImagingTableEntryPtr>& entries = _squaredGroupLookup[index];
	table._entries = entries;
	table.Update();
	return table;
}

void ImagingTable::Print()
{
	Logger::Info << "=== IMAGING TABLE ===\n"
		"       # Pol Ch JG ²G In Freq(MHz)\n";
	for(size_t i=0; i!=IndependentGroupCount(); ++i)
	{
		bool isLastGroup = ((i+1) == IndependentGroupCount());
		Logger::Info << "| Independent group:\n";
		GetIndependentGroup(i).printIndependentGroup(isLastGroup);
		if(!isLastGroup)
			Logger::Info << "|\n";
	}
}

void ImagingTable::printIndependentGroup(bool isFinal)
{
	for(size_t i=0; i!=SquaredGroupCount(); ++i)
	{
		bool isSecondFinal = (i+1)==SquaredGroupCount();
		ImagingTable table = GetSquaredGroup(i);
		for(size_t j=0; j!=table._entries.size(); ++j)
		{
			if(i == 0 && j==0)
				Logger::Info << "+-";
			else if(isFinal)
				Logger::Info << "  ";
			else
				Logger::Info << "| ";
			
			if(j == 0)
				Logger::Info << "+-";
			else if(isSecondFinal)
				Logger::Info << "  ";
			else
				Logger::Info << "| ";
			Logger::Info << "J-";
			Logger::Info << table._entries[j]->ToString() << '\n';
		}
		if(isFinal && isSecondFinal)
			Logger::Info << '\n';
		else if(isFinal)
			Logger::Info << "  |\n";
		else if(isSecondFinal)
			Logger::Info << "|\n";
		else
			Logger::Info << "| |\n";
	}
}

string ImagingTableEntry::ToString()
{
	std::ostringstream str;
	if(index < 10) str << ' ';
	str << index << ' ';
	std::string polStr = Polarization::TypeToShortString(polarization);
	if(polStr.size() < 2) str << ' ';
	str << polStr << "  ";
	if(outputChannelIndex < 10) str << ' ';
	str << outputChannelIndex
		<< "  " << joinedGroupIndex
		<< "  " << squaredDeconvolutionIndex
		<< "  " << outputIntervalIndex
		<< "  " << round(bandStartFrequency*1e-6)
		<< "-" << round(bandEndFrequency*1e-6);
	
	return str.str();
}

void ImagingTable::updateIndependentGroupLookup()
{
	_independentGroupLookup.clear();
	std::map<size_t, size_t> indexToLookupSet;
	for(ImagingTableEntryPtr e : _entries)
	{
		size_t groupIndex = e->joinedGroupIndex;
		if(indexToLookupSet.count(groupIndex) == 0)
		{
			indexToLookupSet.insert(std::make_pair(groupIndex, _independentGroupLookup.size()));
			_independentGroupLookup.push_back(std::vector<ImagingTableEntryPtr>(1, e));
		}
		else {
			_independentGroupLookup[indexToLookupSet[groupIndex]].push_back(e);
		}
	}
}

void ImagingTable::updateSquaredGroupLookup()
{
	_squaredGroupLookup.clear();
	std::map<size_t, size_t> indexToLookupSet;
	for(ImagingTableEntryPtr e : _entries)
	{
		size_t groupIndex = e->squaredDeconvolutionIndex;
		if(indexToLookupSet.count(groupIndex) == 0)
		{
			indexToLookupSet.insert(std::make_pair(groupIndex, _squaredGroupLookup.size()));
			_squaredGroupLookup.push_back(std::vector<ImagingTableEntryPtr>(1, e));
		}
		else {
			_squaredGroupLookup[indexToLookupSet[groupIndex]].push_back(e);
		}
	}
}

void ImagingTable::updateImageLookup()
{
	_imageLookup.clear();
	for(ImagingTableEntryPtr e : _entries)
	{
		bool isImaginary = false; 
		for(size_t i=0; i!=e->imageCount; ++i)
		{
			_imageLookup.push_back(std::make_pair(e, isImaginary));
			isImaginary = true;
		}
	}
}
