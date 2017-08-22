#include "partitionedms.h"

#include "averagingmsrowprovider.h"
#include "directmsrowprovider.h"
#include "msrowprovider.h"
#include "noisemsrowprovider.h"

#include "../progressbar.h"

#include "../wsclean/logger.h"
#include "../wsclean/wscleansettings.h"

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <errno.h>
#include <fcntl.h>
#include <string.h>

#include <cstdio>
#include <fstream>
#include <sstream>
#include <memory>
#include <vector>

#include <boost/filesystem/path.hpp>

#include <casacore/measures/Measures/MEpoch.h>

#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>

// #define REDUNDANT_VALIDATION 1

PartitionedMS::PartitionedMS(const Handle& handle, size_t partIndex, PolarizationEnum polarization, size_t dataDescId) :
	_handle(handle),
	_metaFile(getMetaFilename(handle._data->_msPath, handle._data->_temporaryDirectory, dataDescId)),
	_modelFileMap(0),
	_currentRow(0),
	_readPtrIsOk(true),
	_metaPtrIsOk(true),
	_weightPtrIsOk(true),
	_polarization(polarization)
{
	_metaFile.read(reinterpret_cast<char*>(&_metaHeader), sizeof(MetaHeader));
	std::vector<char> msPath(_metaHeader.filenameLength+1, char(0));
	_metaFile.read(msPath.data(), _metaHeader.filenameLength);
	Logger::Info << "Opening reordered part " << partIndex << " spw " << dataDescId << " for " << msPath.data() << '\n';
	_msPath = msPath.data();
	std::string partPrefix = getPartPrefix(msPath.data(), partIndex, polarization, dataDescId, handle._data->_temporaryDirectory);
	
	_dataFile.open(partPrefix+".tmp", std::ios::in);
	if(!_dataFile.good())
		throw std::runtime_error("Error opening temporary data file");
	_dataFile.read(reinterpret_cast<char*>(&_partHeader), sizeof(PartHeader));
	if(!_dataFile.good())
		throw std::runtime_error("Error reading header from file");
	
	if(_partHeader.hasModel)
	{
		_fd = open((partPrefix+"-m.tmp").c_str(), O_RDWR);
		if(_fd == -1)
			throw std::runtime_error("Error opening temporary model data file");
		size_t length = _partHeader.channelCount * _metaHeader.selectedRowCount * sizeof(std::complex<float>);
		if(length == 0)
			_modelFileMap = 0;
		else {
			_modelFileMap = reinterpret_cast<char*>( mmap(NULL, length, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_NORESERVE, _fd, 0) );
			if(_modelFileMap == MAP_FAILED)
			{
				int errsv = errno;
				char buffer[1024];
				const char* msg = strerror_r(errsv, buffer, 1024); 
				
				_modelFileMap = 0;
				throw std::runtime_error(std::string("Error creating memory map to temporary model file: mmap() returned MAP_FAILED with error message: ") + msg);
			}
		}
	}
	
	_weightFile.open(partPrefix+"-w.tmp", std::ios::in);
	if(!_weightFile.good())
		throw std::runtime_error("Error opening temporary data file");
	_weightBuffer.resize(_partHeader.channelCount);
	_modelBuffer.resize(_partHeader.channelCount);
}

PartitionedMS::~PartitionedMS()
{
	if(_modelFileMap!=0)
	{
		size_t length = _partHeader.channelCount * _metaHeader.selectedRowCount * sizeof(std::complex<float>);
		if(length != 0)
			munmap(_modelFileMap, length);
	}
	if(_partHeader.hasModel)
		close(_fd);
}

void PartitionedMS::Reset()
{
	_currentRow = 0;
	_metaFile.seekg(sizeof(MetaHeader) + _metaHeader.filenameLength, std::ios::beg);
	_dataFile.seekg(sizeof(PartHeader), std::ios::beg);
	_weightFile.seekg(0, std::ios::beg);
	_readPtrIsOk = true;
	_metaPtrIsOk = true;
	_weightPtrIsOk = true;
}

bool PartitionedMS::CurrentRowAvailable()
{
	return _currentRow < _metaHeader.selectedRowCount;
}

void PartitionedMS::NextRow()
{
	++_currentRow;
	if(_currentRow < _metaHeader.selectedRowCount)
	{
		if(_readPtrIsOk)
			_dataFile.seekg(_partHeader.channelCount * sizeof(std::complex<float>), std::ios::cur);
		else
			_readPtrIsOk = true;
		
		if(_metaPtrIsOk)
			_metaFile.seekg(sizeof(MetaRecord), std::ios::cur);
		else
			_metaPtrIsOk = true;
		
		if(_weightPtrIsOk && _partHeader.hasWeights)
			_weightFile.seekg(_partHeader.channelCount * sizeof(float), std::ios::cur);
		_weightPtrIsOk = true;
	}
}

void PartitionedMS::ReadMeta(double& u, double& v, double& w, size_t& dataDescId)
{
	if(!_metaPtrIsOk)
		_metaFile.seekg(-sizeof(MetaRecord), std::ios::cur);
	_metaPtrIsOk = false;
	
	MetaRecord record;
	_metaFile.read(reinterpret_cast<char*>(&record), sizeof(MetaRecord));
	u = record.u;
	v = record.v;
	w = record.w;
	dataDescId = record.dataDescId;
}

void PartitionedMS::ReadMeta(double& u, double& v, double& w, size_t& dataDescId, size_t& antenna1, size_t& antenna2)
{
	if(!_metaPtrIsOk)
		_metaFile.seekg(-sizeof(MetaRecord), std::ios::cur);
	_metaPtrIsOk = false;
	
	MetaRecord record;
	_metaFile.read(reinterpret_cast<char*>(&record), sizeof(MetaRecord));
	u = record.u;
	v = record.v;
	w = record.w;
	dataDescId = record.dataDescId;
	antenna1 = record.antenna1;
	antenna2 = record.antenna2;
}

void PartitionedMS::ReadData(std::complex<float>* buffer)
{
	if(!_readPtrIsOk)
	{
		_dataFile.seekg(-_partHeader.channelCount * sizeof(std::complex<float>), std::ios::cur);
	}
#ifdef REDUNDANT_VALIDATION
	size_t pos = size_t(_dataFile.tellg()) - sizeof(PartHeader);
	size_t fact = _partHeader.hasModel ? 2 : 1;
	if(pos != fact * _currentRow * _partHeader.channelCount * sizeof(std::complex<float>))
	{
		std::ostringstream s;
		s << "Not on right pos: " << pos << " instead of " << fact * _currentRow * _partHeader.channelCount * sizeof(std::complex<float>) <<
			" (row " << (pos / (fact * _partHeader.channelCount * sizeof(std::complex<float>))) << " instead of " << _currentRow << ")";
		throw std::runtime_error(s.str());
	}
#endif
	_dataFile.read(reinterpret_cast<char*>(buffer), _partHeader.channelCount * sizeof(std::complex<float>));
	_readPtrIsOk = false;
}

void PartitionedMS::ReadModel(std::complex<float>* buffer)
{
#ifdef REDUNDANT_VALIDATION
	if(!_partHeader.hasModel)
		throw std::runtime_error("Partitioned MS initialized without model");
#endif
	size_t rowLength = _partHeader.channelCount * sizeof(std::complex<float>);
	memcpy(reinterpret_cast<char*>(buffer), _modelFileMap + rowLength*_currentRow, rowLength);
}

void PartitionedMS::WriteModel(size_t rowId, std::complex<float>* buffer)
{
#ifdef REDUNDANT_VALIDATION
	if(!_partHeader.hasModel)
		throw std::runtime_error("Partitioned MS initialized without model");
#endif
	_weightFile.seekg(_partHeader.channelCount * sizeof(float) * rowId, std::ios::beg);
	_weightFile.read(reinterpret_cast<char*>(_weightBuffer.data()), _partHeader.channelCount * sizeof(float));
	for(size_t i=0; i!=_partHeader.channelCount; ++i)
		buffer[i] *= _weightBuffer[i];
	
	size_t rowLength = _partHeader.channelCount * sizeof(std::complex<float>);
	std::complex<float>* modelWritePtr = reinterpret_cast<std::complex<float>*>(_modelFileMap + rowLength*rowId);
	
	// In case the value was not sampled in this pass, it will be set to infinite and should not overwrite the current
	// value in the set.
	for(size_t i=0; i!=_partHeader.channelCount; ++i)
	{
		if(std::isfinite(buffer[i].real()))
			modelWritePtr[i] = buffer[i];
	}
}

void PartitionedMS::ReadWeights(std::complex<float>* buffer)
{
	if(!_weightPtrIsOk)
		_weightFile.seekg(-_partHeader.channelCount * sizeof(float), std::ios::cur);
	float* displacedBuffer = reinterpret_cast<float*>(buffer)+_partHeader.channelCount;
	_weightFile.read(reinterpret_cast<char*>(displacedBuffer), _partHeader.channelCount * sizeof(float));
	_weightPtrIsOk = false;
	copyRealToComplex(buffer, displacedBuffer, _partHeader.channelCount);
}

void PartitionedMS::ReadWeights(float* buffer)
{
	if(!_weightPtrIsOk)
		_weightFile.seekg(-_partHeader.channelCount * sizeof(float), std::ios::cur);
	_weightFile.read(reinterpret_cast<char*>(buffer), _partHeader.channelCount * sizeof(float));
	_weightPtrIsOk = false;
}

std::string PartitionedMS::getPartPrefix(const std::string& msPathStr, size_t partIndex, PolarizationEnum pol, size_t dataDescId, const std::string& tempDir)
{
	boost::filesystem::path
		msPath(msPathStr),
		prefixPath;
	if(tempDir.empty())
		prefixPath = msPath;
	else
		prefixPath = boost::filesystem::path(tempDir) / msPath.filename();
		
	std::string prefix(prefixPath.string());
	while(!prefix.empty() && *prefix.rbegin() == '/')
		prefix.resize(prefix.size()-1);
	
	std::ostringstream partPrefix;
	partPrefix << prefix << "-part";
	if(partIndex < 1000) partPrefix << '0';
	if(partIndex < 100) partPrefix << '0';
	if(partIndex < 10) partPrefix << '0';
	partPrefix << partIndex;
	partPrefix << "-";
	partPrefix << Polarization::TypeToShortString(pol);
	partPrefix << "-b" << dataDescId;
	return partPrefix.str();
}

string PartitionedMS::getMetaFilename(const string& msPathStr, const std::string& tempDir, size_t dataDescId)
{
	boost::filesystem::path
		msPath(msPathStr),
		prefixPath;
	if(tempDir.empty())
		prefixPath = msPath;
	else
		prefixPath = boost::filesystem::path(tempDir) / msPath.filename();
	std::string prefix(prefixPath.string());
	while(!prefix.empty() && *prefix.rbegin() == '/')
		prefix.resize(prefix.size()-1);
	std::ostringstream s;
	s << prefix << "-spw" << dataDescId << "-parted-meta.tmp";
	return s.str();
}

// should be private but is not allowed on older compilers
struct PartitionFiles
{
	std::ofstream
		*data,
		*weight,
		*model;
};

/*
 * When partitioned:
 * One global file stores:
 * - Metadata:
 *   * Number of selected rows
 *   * Filename length + string
 *   * [ UVW, dataDescId ]
 * The binary parts store the following information:
 * - Number of channels
 * - Start channel in MS
 * - Total weight in part
 * - Data    (single polarization, as requested)
 * - Weights (single, only needed when imaging PSF)
 * - Model, optionally
 */
PartitionedMS::Handle PartitionedMS::Partition(const string& msPath, const std::vector<ChannelRange>& channels, MSSelection& selection, const string& dataColumnName, bool includeModel, bool initialModelRequired, const WSCleanSettings& settings)
{
	const bool modelUpdateRequired = settings.modelUpdateRequired;
	std::set<PolarizationEnum> polsOut;
	if(settings.useIDG)
		polsOut.insert(Polarization::Instrumental);
	else
		polsOut = settings.polarizations;
	const std::string& temporaryDirectory = settings.temporaryDirectory;
	
	size_t channelParts = channels.size();
	
	if(channelParts != 1)
	{
		Logger::Debug << "Partitioning in " << channels.size() << " channels:";
		for(size_t i=0; i!=channels.size(); ++i)
			Logger::Debug << ' ' << channels[i].dataDescId << ':' << channels[i].start << '-' << channels[i].end;
	}
	Logger::Debug << '\n';
	
	// We need to enumerate the data desc ids, because each one needs a separate
	// meta file because they can have different uvws and other info
	std::map<size_t,size_t> selectedDataDescIds;
	getDataDescIdMap(selectedDataDescIds, channels);
	
	// Ordered as files[pol x channelpart]
	std::vector<PartitionFiles>
		files(channelParts*polsOut.size());
	
	size_t fileIndex = 0;
	for(size_t part=0; part!=channelParts; ++part)
	{
		for(std::set<PolarizationEnum>::const_iterator p=polsOut.begin(); p!=polsOut.end(); ++p)
		{
			PartitionFiles& f = files[fileIndex];
			std::string partPrefix = getPartPrefix(msPath, part, *p, channels[part].dataDescId, temporaryDirectory);
			f.data = new std::ofstream(partPrefix + ".tmp");
			f.weight = new std::ofstream(partPrefix + "-w.tmp");
			if(initialModelRequired)
				f.model = new std::ofstream(partPrefix + "-m.tmp");
			f.data->seekp(sizeof(PartHeader), std::ios::beg);
			
			++fileIndex;
		}
	}
	
	std::unique_ptr<MSRowProvider> rowProvider;
	if(settings.baselineDependentAveragingInWavelengths == 0.0)
	{
		if(settings.simulateNoise)
			rowProvider.reset(new NoiseMSRowProvider(settings.simulatedNoiseStdDev, msPath, selection, selectedDataDescIds, dataColumnName, initialModelRequired));
		else
			rowProvider.reset(new DirectMSRowProvider(msPath, selection, selectedDataDescIds, dataColumnName, initialModelRequired));
	}
	else
		rowProvider.reset(new AveragingMSRowProvider(settings.baselineDependentAveragingInWavelengths, msPath, selection, selectedDataDescIds, dataColumnName, initialModelRequired));
	
	std::vector<PolarizationEnum> msPolarizations = GetMSPolarizations(rowProvider->MS());
	
	const casacore::IPosition shape(rowProvider->DataShape());
	const size_t polarizationCount = shape[0];
	size_t channelCount = shape[1];
	
	Logger::Info << "Reordering " << msPath << " into " << channelParts << " x " << polsOut.size() << " parts.\n";

	// Write header of meta file, one meta file for each data desc id
	// TODO rather than writing we can just skip and write later
	std::vector<std::unique_ptr<std::ofstream>> metaFiles(selectedDataDescIds.size());
	for(std::map<size_t,size_t>::const_iterator i=selectedDataDescIds.begin();
			i!=selectedDataDescIds.end(); ++i)
	{
		size_t dataDescId = i->first;
		size_t spwIndex = i->second;
		std::string metaFilename = getMetaFilename(msPath, temporaryDirectory, dataDescId);
		metaFiles[spwIndex].reset(new std::ofstream(metaFilename));
		MetaHeader metaHeader;
		memset(&metaHeader, 0, sizeof(MetaHeader));
		metaHeader.selectedRowCount = 0; // not yet known
		metaHeader.filenameLength = msPath.size();
		metaHeader.startTime = rowProvider->StartTime();
		metaFiles[spwIndex]->write(reinterpret_cast<char*>(&metaHeader), sizeof(metaHeader));
		metaFiles[spwIndex]->write(msPath.c_str(), msPath.size());
	}
	
	// Write actual data
	size_t polarizationsPerFile = settings.useIDG ? 4 : 1;
	std::vector<std::complex<float>> dataBuffer(polarizationCount * channelCount);
	std::vector<float> weightBuffer(polarizationCount * channelCount);
	
	casacore::Array<std::complex<float>> dataArray(shape), modelArray(shape);
	casacore::Array<float> weightSpectrumArray(shape);
	casacore::Array<bool> flagArray(shape);
	ProgressBar progress1("Reordering");
	
	size_t selectedRowsTotal = 0;
	ao::uvector<size_t> selectedRowCountPerSpwIndex(selectedDataDescIds.size(), 0);
	while(!rowProvider->AtEnd())
	{
		progress1.SetProgress(rowProvider->CurrentProgress(), rowProvider->TotalProgress());
		
		MetaRecord meta;
		memset(&meta, 0, sizeof(MetaRecord));

		uint32_t dataDescId, antenna1, antenna2;
		rowProvider->ReadData(dataArray, flagArray, weightSpectrumArray, meta.u, meta.v, meta.w, dataDescId, antenna1, antenna2);
		meta.dataDescId = dataDescId;
		meta.antenna1 = antenna1;
		meta.antenna2 = antenna2;
		size_t spwIndex = selectedDataDescIds[meta.dataDescId];
		++selectedRowCountPerSpwIndex[spwIndex];
		++selectedRowsTotal;
		std::ofstream& metaFile = *metaFiles[spwIndex];
		metaFile.write(reinterpret_cast<char*>(&meta), sizeof(MetaRecord));
		if(!metaFile.good())
			throw std::runtime_error("Error writing to temporary file");
		
		if(initialModelRequired)
			rowProvider->ReadModel(modelArray);
			
		fileIndex = 0;
		for(size_t part=0; part!=channelParts; ++part)
		{
			if(channels[part].dataDescId == int(meta.dataDescId))
			{
				size_t
					partStartCh = channels[part].start,
					partEndCh = channels[part].end;
				
				for(std::set<PolarizationEnum>::const_iterator p=polsOut.begin(); p!=polsOut.end(); ++p)
				{
					PartitionFiles& f = files[fileIndex];
					copyWeightedData(dataBuffer.data(), partStartCh, partEndCh, msPolarizations, dataArray, weightSpectrumArray, flagArray, *p);
					f.data->write(reinterpret_cast<char*>(dataBuffer.data()), (partEndCh - partStartCh) * sizeof(std::complex<float>) * polarizationsPerFile);
					if(!f.data->good())
						throw std::runtime_error("Error writing to temporary data file");
					
					if(initialModelRequired)
					{
						copyWeightedData(dataBuffer.data(), partStartCh, partEndCh, msPolarizations, modelArray, weightSpectrumArray, flagArray, *p);
						f.model->write(reinterpret_cast<char*>(dataBuffer.data()), (partEndCh - partStartCh) * sizeof(std::complex<float>) * polarizationsPerFile);
						if(!f.model->good())
							throw std::runtime_error("Error writing to temporary data file");
					}
					
					copyWeights(weightBuffer.data(), partStartCh, partEndCh, msPolarizations, dataArray, weightSpectrumArray, flagArray, *p);
					f.weight->write(reinterpret_cast<char*>(weightBuffer.data()), (partEndCh - partStartCh) * sizeof(float) * polarizationsPerFile);
					if(!f.weight->good())
						throw std::runtime_error("Error writing to temporary weights file");
					++fileIndex;
				}
			} else {
				fileIndex += polsOut.size();
			}
		}
		
		rowProvider->NextRow();
	}
	progress1.SetProgress(rowProvider->TotalProgress(), rowProvider->TotalProgress());
	Logger::Debug << "Total selected rows: " << selectedRowsTotal << '\n';
	rowProvider->OutputStatistics();
	
	// Rewrite meta headers to include selected row count
	for(std::map<size_t,size_t>::const_iterator i=selectedDataDescIds.begin();
			i!=selectedDataDescIds.end(); ++i)
	{
		size_t spwIndex = i->second;
		MetaHeader metaHeader;
		memset(&metaHeader, 0, sizeof(MetaHeader));
		metaHeader.selectedRowCount = selectedRowCountPerSpwIndex[spwIndex];
		metaHeader.filenameLength = msPath.size();
		metaHeader.startTime = rowProvider->StartTime();
		metaFiles[spwIndex]->seekp(0);
		metaFiles[spwIndex]->write(reinterpret_cast<char*>(&metaHeader), sizeof(metaHeader));
		metaFiles[spwIndex]->write(msPath.c_str(), msPath.size());
	}
	
	// Write header to parts and write empty model files (if requested)
	PartHeader header;
	memset(&header, 0, sizeof(PartHeader));
	header.hasModel = includeModel;
	header.hasWeights = true;
	fileIndex = 0;
	dataBuffer.assign(channelCount, 0.0);
	std::unique_ptr<ProgressBar> progress2;
	if(includeModel && !initialModelRequired)
		progress2.reset(new ProgressBar("Initializing model visibilities"));
	for(size_t part=0; part!=channelParts; ++part)
	{
		header.channelStart = channels[part].start,
		header.channelCount = channels[part].end - header.channelStart;
		header.dataDescId = channels[part].dataDescId;
		for(std::set<PolarizationEnum>::const_iterator p=polsOut.begin(); p!=polsOut.end(); ++p)
		{
			PartitionFiles& f = files[fileIndex];
			f.data->seekp(0, std::ios::beg);
			f.data->write(reinterpret_cast<char*>(&header), sizeof(PartHeader));
			if(!f.data->good())
				throw std::runtime_error("Error writing to temporary data file");
			
			delete f.data;
			delete f.weight;
			if(initialModelRequired)
				delete f.model;
			++fileIndex;
			
			// If model is requested, fill model file with zeros
			if(includeModel && !initialModelRequired)
			{
				std::string partPrefix = getPartPrefix(msPath, part, *p, header.dataDescId, temporaryDirectory);
				std::ofstream modelFile(partPrefix + "-m.tmp");
				const size_t selectedRowCount = selectedRowCountPerSpwIndex[selectedDataDescIds[channels[part].dataDescId]];
				for(size_t i=0; i!=selectedRowCount; ++i)
				{
					modelFile.write(reinterpret_cast<char*>(dataBuffer.data()), header.channelCount * sizeof(std::complex<float>) * polarizationsPerFile);
					progress2->SetProgress(part*selectedRowCount + i, channelParts*selectedRowCount);
				}
			}
		}
	}
	progress2.reset();
	
	return Handle(msPath, dataColumnName, temporaryDirectory, channels, initialModelRequired, modelUpdateRequired, polsOut, selection);
}

void PartitionedMS::unpartition(const PartitionedMS::Handle& handle)
{
	const std::set<PolarizationEnum> pols = handle._data->_polarizations;
	
	std::map<size_t,size_t> dataDescIds;
	getDataDescIdMap(dataDescIds, handle._data->_channels);
	
	std::vector<MetaHeader> metaHeaders(dataDescIds.size());
	for(const std::pair<size_t,size_t>& dataDescId : dataDescIds)
	{
		std::ifstream metaFile(getMetaFilename(handle._data->_msPath, handle._data->_temporaryDirectory, dataDescId.first));
		MetaHeader& metaHeader = metaHeaders[dataDescId.second];
		metaFile.read(reinterpret_cast<char*>(&metaHeader), sizeof(MetaHeader));
		std::vector<char> msPath(metaHeader.filenameLength+1, char(0));
		metaFile.read(msPath.data(), metaHeader.filenameLength);
	}
	
	ChannelRange firstRange = handle._data->_channels[0];
	std::ifstream firstDataFile(getPartPrefix(handle._data->_msPath, 0, *pols.begin(), firstRange.dataDescId, handle._data->_temporaryDirectory)+".tmp", std::ios::in);
	if(!firstDataFile.good())
		throw std::runtime_error("Error opening temporary data file");
	PartHeader firstPartHeader;
	firstDataFile.read(reinterpret_cast<char*>(&firstPartHeader), sizeof(PartHeader));
	if(!firstDataFile.good())
		throw std::runtime_error("Error reading from temporary data file");
	
	if(firstPartHeader.hasModel)
	{
		const size_t channelParts = handle._data->_channels.size();
		
		// Open the temporary files
		std::vector<std::ifstream*> modelFiles(channelParts*pols.size()), weightFiles(channelParts*pols.size());
		size_t fileIndex = 0;
		for(size_t part=0; part!=channelParts; ++part)
		{
			size_t dataDescId = handle._data->_channels[part].dataDescId;
			for(std::set<PolarizationEnum>::const_iterator p=pols.begin(); p!=pols.end(); ++p)
			{
				std::string partPrefix = getPartPrefix(handle._data->_msPath, part, *p, dataDescId, handle._data->_temporaryDirectory);
				modelFiles[fileIndex] = new std::ifstream(partPrefix + "-m.tmp");
				if(firstPartHeader.hasWeights)
					weightFiles[fileIndex] = new std::ifstream(partPrefix + "-w.tmp");
				++fileIndex;
			}
		}
		
		casacore::MeasurementSet ms(handle._data->_msPath, casacore::Table::Update);
		const std::vector<PolarizationEnum> msPolarizations = GetMSPolarizations(ms);
		initializeModelColumn(ms);
		casacore::ROScalarColumn<int> antenna1Column(ms, ms.columnName(casacore::MSMainEnums::ANTENNA1));
		casacore::ROScalarColumn<int> antenna2Column(ms, ms.columnName(casacore::MSMainEnums::ANTENNA2));
		casacore::ROScalarColumn<int> fieldIdColumn(ms, ms.columnName(casacore::MSMainEnums::FIELD_ID));
		casacore::ROScalarColumn<double> timeColumn(ms, ms.columnName(casacore::MSMainEnums::TIME));
		casacore::ROScalarColumn<int> dataDescIdColumn(ms, ms.columnName(casacore::MSMainEnums::DATA_DESC_ID));
		casacore::ROArrayColumn<casacore::Complex> dataColumn(ms, handle._data->_dataColumnName);
		casacore::ArrayColumn<casacore::Complex> modelColumn(ms, ms.columnName(casacore::MSMainEnums::MODEL_DATA));
		casacore::ROArrayColumn<double> uvwColumn(ms, ms.columnName(casacore::MSMainEnums::UVW));
		
		const casacore::IPosition shape(dataColumn.shape(0));
		size_t channelCount = shape[1];
		
		std::vector<std::complex<float>> modelDataBuffer(channelCount);
		std::vector<float> weightBuffer(channelCount);
		casacore::Array<std::complex<float>> modelDataArray(shape);
	
		ProgressBar progress(std::string("Writing changed model back to ") + handle._data->_msPath);
		size_t startRow, endRow;
		getRowRange(ms, handle._data->_selection, startRow, endRow);
		size_t timestep = handle._data->_selection.HasInterval() ? handle._data->_selection.IntervalStart() : 0;
		double time = timeColumn(startRow);
		size_t selectedRowCountForDebug = 0;
		for(size_t row=startRow; row!=endRow; ++row)
		{
			progress.SetProgress(row - startRow, startRow - endRow);
			const int
				a1 = antenna1Column(row), a2 = antenna2Column(row),
				fieldId = fieldIdColumn(row), dataDescId = dataDescIdColumn(row);
			casacore::Vector<double> uvw = uvwColumn(row);
				
			if(time != timeColumn(row))
			{
				++timestep;
				time = timeColumn(row);
			}
			if(handle._data->_selection.IsSelected(fieldId, timestep, a1, a2, uvw))
			{
				std::map<size_t,size_t>::const_iterator dataDescIdIter = dataDescIds.find(dataDescId);
				if(dataDescIdIter != dataDescIds.end())
				{
					modelColumn.get(row, modelDataArray);
					size_t fileIndex = 0;
					for(size_t part=0; part!=channelParts; ++part)
					{
						size_t
							dataDescId = handle._data->_channels[part].dataDescId,
							partStartCh = handle._data->_channels[part].start,
							partEndCh = handle._data->_channels[part].end;
						if(dataDescId == dataDescIdIter->second)
						{
							for(std::set<PolarizationEnum>::const_iterator p=pols.begin(); p!=pols.end(); ++p)
							{
								modelFiles[fileIndex]->read(reinterpret_cast<char*>(modelDataBuffer.data()), (partEndCh - partStartCh) * sizeof(std::complex<float>));
								if(!modelFiles[fileIndex]->good())
									throw std::runtime_error("Error reading from temporary model data file");
								if(firstPartHeader.hasWeights)
								{
									weightFiles[fileIndex]->read(reinterpret_cast<char*>(weightBuffer.data()), (partEndCh - partStartCh) * sizeof(float));
									if(!weightFiles[fileIndex]->good())
										throw std::runtime_error("Error reading from temporary weight data file");
									for(size_t i=0; i!=partEndCh - partStartCh; ++i)
										modelDataBuffer[i] /= weightBuffer[i];
								}
								reverseCopyData(modelDataArray, partStartCh, partEndCh, msPolarizations, modelDataBuffer.data(), *p);
								
								++fileIndex;
							}
						}
						else {
							fileIndex += pols.size();
						}
					}
					modelColumn.put(row, modelDataArray);
					selectedRowCountForDebug++;
				}
			}
		}
		progress.SetProgress(ms.nrow(),ms.nrow());
		
		Logger::Debug << "Row count during unpartitioning: " << selectedRowCountForDebug << '\n';
		
		fileIndex = 0;
		for(size_t part=0; part!=channelParts; ++part)
		{
			for(std::set<PolarizationEnum>::const_iterator p=pols.begin(); p!=pols.end(); ++p)
			{
				delete modelFiles[fileIndex];
				if(firstPartHeader.hasWeights)
					delete weightFiles[fileIndex];
				++fileIndex;
			}
		}
	}
}

void PartitionedMS::Handle::decrease()
{
	--(_data->_referenceCount);
	if(_data->_referenceCount == 0)
	{
		if(_data->_modelUpdateRequired && !_data->_initialModelRequired)
			PartitionedMS::unpartition(*this);
		
		Logger::Info << "Cleaning up temporary files...\n";
		
		std::set<size_t> removedMetaFiles;
		for(size_t part=0; part!=_data->_channels.size(); ++part)
		{
			for(std::set<PolarizationEnum>::const_iterator p=_data->_polarizations.begin(); p!=_data->_polarizations.end(); ++p)
			{
				std::string prefix = getPartPrefix(_data->_msPath, part, *p, _data->_channels[part].dataDescId, _data->_temporaryDirectory);
				std::remove((prefix + ".tmp").c_str());
				std::remove((prefix + "-w.tmp").c_str());
				std::remove((prefix + "-m.tmp").c_str());
			}
			size_t dataDescId = _data->_channels[part].dataDescId;
			if(removedMetaFiles.count(dataDescId) == 0)
			{
				removedMetaFiles.insert(dataDescId);
				std::string metaFile = getMetaFilename(_data->_msPath, _data->_temporaryDirectory, dataDescId);
				std::remove(metaFile.c_str());
			}
		}
		delete _data;
	}
}

void PartitionedMS::openMS()
{
	if(_ms == 0)
		_ms.reset(new casacore::MeasurementSet(_msPath.data()));
}

void PartitionedMS::MakeIdToMSRowMapping(vector<size_t>& idToMSRow)
{
	openMS();
	const MSSelection& selection = _handle._data->_selection;
	std::map<size_t, size_t> dataDescIds;
	getDataDescIdMap(dataDescIds, _handle._data->_channels);
	std::set<size_t> dataDescIdSet;
	for(std::map<size_t, size_t>::const_iterator i=dataDescIds.begin(); i!=dataDescIds.end(); ++i)
		dataDescIdSet.insert(i->first);
	size_t startRow, endRow;
	getRowRangeAndIDMap(*_ms, selection, startRow, endRow, dataDescIdSet, idToMSRow);
}

void PartitionedMS::getDataDescIdMap(std::map<size_t, size_t>& dataDescIds, const vector<PartitionedMS::ChannelRange>& channels)
{
	size_t spwIndex = 0;
	for(size_t i=0; i!=channels.size(); ++i)
	{
		if(dataDescIds.count(channels[i].dataDescId) == 0)
		{
			dataDescIds.insert(std::make_pair(channels[i].dataDescId, spwIndex));
			++spwIndex;
		}
	}
}
