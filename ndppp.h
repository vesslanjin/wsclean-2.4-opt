#ifndef NDPPP_H
#define NDPPP_H

#include <fstream>
#include <boost/filesystem/operations.hpp>

#include "model/model.h"
#include "model/powerlawsed.h"

#include "application.h"

class NDPPP
{
public:
	static void ConvertSkyModelToSourceDB(const std::string& destination, const std::string& input)
	{
		boost::filesystem::remove_all(destination);
		Application::Run("makesourcedb in=" + input + " out=" + destination + " format='<'");
	}
	
	static void WriteStandardHeader(std::ostream& stream, double refFrequency)
	{
		stream.precision(15);
		stream << "# (Name, Type, Ra, Dec, I, Q, U, V, ReferenceFrequency='" << refFrequency << "', SpectralIndex, MajorAxis, MinorAxis, Orientation) = format\n\n";
	}
	
	static void WriteHeaderForSpectralTerms(std::ostream& stream, double refFrequency)
	{
		stream.precision(15);
		stream
			<< "Format = Name, Type, Ra, Dec, I, SpectralIndex, LogarithmicSI, ReferenceFrequency='" << refFrequency << "', MajorAxis, MinorAxis, Orientation\n";
	}
	
	static void WriteOldHeaderForSpectralTerms(std::ostream& stream, double refFrequency, const std::string& spectralFunction)
	{
		stream.precision(15);
		stream
			<< "# (Name, Type, Ra, Dec, SpectralTerms, MajorAxis, MinorAxis, Orientation) = format\n"
			<< "# ReferenceFrequency = " << refFrequency << '\n'
			<< "# SpectralFunction = " << spectralFunction << '\n';
	}
	
	static void addSITerms(std::ostream& stream, const ao::uvector<double>& siTerms)
	{
		stream << '[';
		if(!siTerms.empty())
		{
			stream << siTerms[0];
			for(size_t i=1; i!=siTerms.size(); ++i)
			{
				stream << ',' << siTerms[i];
			}
		}
		stream << ']';
	}
	
	static void WritePointComponent(std::ostream& stream, const std::string& name, long double ra, long double dec, double i, double q, double u, double v, double freq, const ao::uvector<double>& siTerms)
	{
		stream << name << ",POINT,"
			<< RaDecCoord::RAToString(ra, ':') << ','
			<< RaDecCoord::DecToString(dec, '.') << ','
			<< i << ',' << q << ',' << u << ',' << v << ','
			<< freq << ",";
		addSITerms(stream, siTerms);
		stream << ",,,\n";
	}
	
	static void WriteGaussianComponent(std::ostream& stream, const std::string& name, long double ra, long double dec, double i, double q, double u, double v, double freq, const ao::uvector<double>& siTerms, double maj, double min, double posangle)
	{
		stream << name << ",GAUSSIAN,"
			<< RaDecCoord::RAToString(ra, ':') << ','
			<< RaDecCoord::DecToString(dec, '.') << ','
			<< i << ',' << q << ',' << u << ',' << v << ','
			<< freq << ",";
		addSITerms(stream, siTerms);
		stream << "," << maj << ',' << min << ',' << posangle << "\n";
	}
	
	static void WritePolynomialPointComponent(std::ostream& stream, const std::string& name, long double ra, long double dec, double i, bool useLogSI, const ao::uvector<double>& polTerms, double referenceFrequencyHz)
	{
		stream << name << ",POINT,"
			<< RaDecCoord::RAToString(ra, ':') << ','
			<< RaDecCoord::DecToString(dec, '.') << ','
      << i << ',';
		addSITerms(stream, polTerms);
    stream
      << "," << (useLogSI ? "true" : "false")
      << "," << referenceFrequencyHz
      << ",,,\n";
	}
	
	static void WritePolynomialGaussianComponent(std::ostream& stream, const std::string& name, long double ra, long double dec, double i, bool useLogSI, const ao::uvector<double>& polTerms, double referenceFrequencyHz, double maj, double min, double posangle)
	{
		stream << name << ",GAUSSIAN,"
			<< RaDecCoord::RAToString(ra, ':') << ','
			<< RaDecCoord::DecToString(dec, '.') << ','
			<< i << ',';
		addSITerms(stream, polTerms);
    stream
      << "," << (useLogSI ? "true" : "false")
      << "," << referenceFrequencyHz
      << "," << maj << ',' << min << ',' << posangle << "\n";
	}
	
	static void WriteOldPolynomialPointComponent(std::ostream& stream, const std::string& name, long double ra, long double dec, const ao::uvector<double>& polTerms)
	{
		stream << name << ",POINT,"
			<< RaDecCoord::RAToString(ra, ':') << ','
			<< RaDecCoord::DecToString(dec, '.') << ',';
		addSITerms(stream, polTerms);
		stream << ",,,\n";
	}
	
	static void WriteOldPolynomialGaussianComponent(std::ostream& stream, const std::string& name, long double ra, long double dec, const ao::uvector<double>& polTerms, double maj, double min, double posangle)
	{
		stream << name << ",GAUSSIAN,"
			<< RaDecCoord::RAToString(ra, ':') << ','
			<< RaDecCoord::DecToString(dec, '.') << ',';
		addSITerms(stream, polTerms);
		stream << "," << maj << ',' << min << ',' << posangle << "\n";
	}
	
	static void SaveSkyModel(const std::string& destination, const Model& model, bool convertClustersToPatches)
	{
		std::ofstream file(destination);
		file << "# (Name, Patch, Type, Ra, Dec, I, Q, U, V, ReferenceFrequency='150.e6', SpectralIndex) = format\n\n";
		std::string patchName = ",";
		for(Model::const_iterator s=model.begin(); s!=model.end(); ++s)
		{
			// Define a patch for this source
			// A patch is created by not giving a source name
			if(convertClustersToPatches)
			{
				if(patchName != s->ClusterName())
				{
					patchName = s->ClusterName();
					file << ", " << s->ClusterName() << ", POINT, , , , , , , ,\n";
				}
			}
			else {
				file << ", " << s->Name() << ", POINT, , , , , , , ,\n";
				patchName = s->Name();
			}
					
			for(size_t ci=0; ci!=s->ComponentCount(); ++ci)
			{
				const ModelComponent& c = s->Component(ci);
				file << s->Name() << '_' << ci << ", " << patchName << ", POINT, "
					<< RaDecCoord::RAToString(c.PosRA(), ':') << ", "
					<< RaDecCoord::DecToString(c.PosDec(), '.') << ", ";
				if(c.HasMeasuredSED())
				{
					const MeasuredSED& sed = c.MSED();
					if(sed.MeasurementCount() != 1)
						throw std::runtime_error("Can only predict single-measurement sky models");
					double
						refFreq = sed.ReferenceFrequencyHz(),
						i = sed.FluxAtFrequency(refFreq, Polarization::StokesI),
						q = sed.FluxAtFrequency(refFreq, Polarization::StokesQ),
						u = sed.FluxAtFrequency(refFreq, Polarization::StokesU),
						v = sed.FluxAtFrequency(refFreq, Polarization::StokesV);
					file << i << ", " << q << ", " << u << ", " << v << ", " << refFreq << ", []\n";
				}
				else {
					const PowerLawSED& sed = static_cast<const PowerLawSED&>(c.SED());
					double refFreq, flux[4];
					std::vector<double> siterms;
					sed.GetData(refFreq, flux, siterms);
					file << flux[0] << ", " << flux[1] << ", " << flux[2] << ", " << flux[3] << ", " << refFreq << ", [" ;
					if(siterms.empty())
						file << "0";
					else {
						file << siterms[0];
						for(size_t i=1; i!=siterms.size(); ++i)
							file << ", " << siterms[i];
					}
					file << "]\n";
				}
			}
		}
	}
	
	static void Predict(const std::string& msName, bool predictWithBeam)
	{
		std::ofstream ndpppParset("wsclean-prediction.parset");
		ndpppParset <<
			"msin = " << msName << "\n"
			"msin.datacolumn = DATA\n"
			"msout = .\n"
			"msout.datacolumn = MODEL_DATA\n"
			"steps = [predict]\n"
			"predict.sourcedb=wsclean-prediction.sourcedb\n"
			"predict.usebeammodel=";
		if(predictWithBeam)
			ndpppParset << "true\n";
		else
			ndpppParset << "false\n";
		ndpppParset.close();
		Application::Run("NDPPP wsclean-prediction.parset");
	}
};

#endif
