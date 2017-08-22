#include "../imagebufferallocator.h"
#include "../wstackinggridder.h"

#include "../../model/model.h"
#include "../../dftpredictionalgorithm.h"
#include "../../fitswriter.h"
#include "../../stopwatch.h"

#include <iostream>
#include <random>

int main(int argc, char* argv[])
{
	ImageBufferAllocator allocator;
	double
		pixelSizeX = 1.0 * (M_PI / 180.0 / 60.0), // 1'
		pixelSizeY = 1.0 * (M_PI / 180.0 / 60.0), // 1'
		phaseRA = RaDecCoord::ParseRA("10h12m47.11s"),
		phaseDec = RaDecCoord::ParseDec("-26d37m20.47s");
	size_t width = 2048, height = 2048;
	const double frequencies[1] = { 150.0e6 } ;
	BandData band(1, frequencies);
	
	Model model;
	for(size_t i=0; i!=8; ++i)
	{
		ModelComponent comp;
		double ra = RaDecCoord::ParseRA("10h12m47.11s");
		double dec = RaDecCoord::ParseDec("-26d37m20.47s");
		ra += double(i)*(M_PI/12.0); // add i hours to RA
		dec -= double(i)*10.0*(M_PI/180.0/60.0); // add i*10 aminutes
		comp.SetPosRA(ra);
		comp.SetPosDec(dec);
		comp.SetSED(MeasuredSED((i>=2 ? 100.0 : 1.0), 150.0e6));
		comp.SetType(ModelComponent::PointSource);
		ModelSource source;
		source.AddComponent(comp);
		model.AddSource(source);
	}
	
	DFTPredictionInput input;
	input.InitializeFromModel(model, phaseRA, phaseDec, band);
	DFTPredictionAlgorithm predicter(input, band);
	
	for(size_t methodIter=0; methodIter!=8; ++methodIter)
	{
		std::unique_ptr<WStackingGridder> gridder;
		
		switch(methodIter)
		{
			case 0:
				gridder.reset(new WStackingGridder(width, height, pixelSizeX, pixelSizeY, 1, &allocator));
				gridder->SetGridMode(WStackingGridder::NearestNeighbour);
				break;
			case 1:
				gridder.reset(new WStackingGridder(width, height, pixelSizeX, pixelSizeY, 1, &allocator, 7, 63));
				break;
			case 2:
				gridder.reset(new WStackingGridder(width, height, pixelSizeX, pixelSizeY, 1, &allocator, 15, 63));
				break;
			case 3:
				gridder.reset(new WStackingGridder(width, height, pixelSizeX, pixelSizeY, 1, &allocator, 7, 127));
				break;
			case 4:
				gridder.reset(new WStackingGridder(width, height, pixelSizeX, pixelSizeY, 1, &allocator, 3, 31));
				break;
			case 5:
				gridder.reset(new WStackingGridder(width, height, pixelSizeX, pixelSizeY, 1, &allocator, 7, 63));
				break;
			case 6:
				gridder.reset(new WStackingGridder(width, height, pixelSizeX, pixelSizeY, 1, &allocator, 11, 63));
				break;
			case 7:
				gridder.reset(new WStackingGridder(width, height, pixelSizeX, pixelSizeY, 1, &allocator, 15, 15));
				break;
		}
		
		gridder->PrepareWLayers(1, 1e12, -1.0, 1.0);
		gridder->StartInversionPass(0);
		
		double maxUV = 0.1 * (1.0/pixelSizeX);
		std::mt19937 rng;
		std::uniform_real_distribution<double> dist(-maxUV, maxUV);
		if(methodIter == 5)
			dist = std::uniform_real_distribution<double>(0, maxUV);
		
		size_t sampleCount = 10000000;
		Stopwatch watch(true);
		for(size_t i=0; i!=sampleCount; ++i)
		{
			double
				u = dist(rng),
				v = dist(rng),
				w = 0.0;
			
			MC2x2 m;
			predicter.Predict(m, u, v, w, 0, 0, 0);
			
			std::complex<float> stokesISample = 0.5f * std::complex<float>(m[0].real() + m[3].real(), m[0].imag() + m[3].imag());
			gridder->AddDataSample(stokesISample, u, v, w);
		}
		
		FitsWriter writer;
		writer.SetImageDimensions(width, height, phaseRA, phaseDec, pixelSizeX, pixelSizeY);
		
		std::string filename;
		switch(methodIter)
		{
			case 0: filename = "aliastest-nn"; break;
			case 1: filename = "aliastest-kb-7-63"; break;
			case 2: filename = "aliastest-kb-15-63"; break;
			case 3: filename = "aliastest-kb-7-127"; break;
			case 4: filename = "aliastest-kb-3-31"; break;
			case 5: filename = "aliastest-kb-7-63-positive"; break;
			case 6: filename = "aliastest-kb-11-63"; break;
			case 7: filename = "aliastest-kb-15-15"; break;
		}
		
		std::cout << "Gridding time " << filename << ": " << watch.ToString() << '\n';
		
		gridder->FinishInversionPass();
		const std::complex<double>* uvlayer = gridder->GetGriddedUVLayer(0);
		ImageBufferAllocator::Ptr uvImage;
		allocator.Allocate(width*height, uvImage);
		for(size_t y=0; y!=height; ++y) {
			for(size_t x=0; x!=width; ++x) {
				uvImage[y*width + x] = std::abs(uvlayer[y*width + x]);
			}
		}
		writer.Write(filename + "-uv.fits", uvImage.data());

		gridder->FinalizeImage(1.0/sampleCount, false);
		
		const double* image = gridder->RealImage();
		writer.Write(filename + ".fits", image);
	}
}
