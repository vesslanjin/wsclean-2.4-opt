#ifndef IUWT_DECONVOLUTION_H
#define IUWT_DECONVOLUTION_H

#include <memory>
#include <string>

#include "../uvector.h"

#include "deconvolutionalgorithm.h"
#include "imageset.h"

#include "../iuwt/iuwtdeconvolutionalgorithm.h"

#include "../wsclean/imagingtable.h"

class IUWTDeconvolution : public DeconvolutionAlgorithm
{
public:
	IUWTDeconvolution() : _useSNRTest(false) { }
	
	virtual void ExecuteMajorIteration(ImageSet& dataImage, ImageSet& modelImage, const ao::uvector<const double*>& psfImages, size_t width, size_t height, bool& reachedMajorThreshold) final override
	{
		IUWTDeconvolutionAlgorithm algorithm(width, height, _gain, _mGain, _cleanBorderRatio, _allowNegativeComponents, _cleanMask, _threshold, _useSNRTest);
		algorithm.PerformMajorIteration(_iterationNumber, MaxNIter(), modelImage, dataImage, psfImages, reachedMajorThreshold);
		if(_iterationNumber >= MaxNIter())
			reachedMajorThreshold = false;
	}
	
	void SetUseSNRTest(bool useSNRTest) { _useSNRTest = useSNRTest; }
	
private:
	bool _useSNRTest;
};

#endif
