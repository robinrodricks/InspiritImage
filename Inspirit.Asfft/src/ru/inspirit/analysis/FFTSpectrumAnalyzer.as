package ru.inspirit.analysis 
{
	import flash.utils.ByteArray;
	import flash.utils.Endian;

	/**
	 * FFT Spectrum Analyzer is a helper class to simplify
	 * workflow with sound spectrum
	 * 
	 * I have used Minim Audio Lib as reference for some of the presented methods
	 * @see http://code.compartmental.net/tools/minim/
	 *  
	 * @author Eugene Zatepyakin
	 */
	public final class FFTSpectrumAnalyzer 
	{
		public static const AVERAGE_NO:int = 0;
		public static const AVERAGE_LINEAR:int = 1;
		public static const AVERAGE_LOGARITHMIC:int = 2;
		
		public var fft:FFT;
		
		private var _sampleRate:int;
		private var _bandWidth:Number;
		private var _numberOfBands:int = 256;
		private var _averageMode:int = 0;
		private var _memory:ByteArray;
		
		public var spectrumData:ByteArray;
		public var averagesData:ByteArray;
		
		public function FFTSpectrumAnalyzer(fft:FFT, sampleRate:int = 44100)
		{
			this.fft = fft;
			
			_sampleRate = sampleRate;
			_bandWidth = (2.0 / fft.signalLength) * (_sampleRate / 2.0);
			
			averagesData = new ByteArray();
			averagesData.endian = Endian.LITTLE_ENDIAN;
			
			_memory = fft.coreMemory;
		}

		/**
		 * Logarithmic averages from Spectrum data
		 * @param minBandwidth 		min bandwidth to use
		 * @param bandsPerOctave	number of bands to divide each octave to
		 */
		public function initLogarithmicAverages(minBandwidth:int = 11, bandsPerOctave:int = 1):int
		{
			_numberOfBands = fft.core.initLogarithmicAverages(minBandwidth, bandsPerOctave, _sampleRate);
			_averageMode = AVERAGE_LOGARITHMIC;
			return _numberOfBands;
		}
		
		/**
		 * Simple Linear everages
		 * @param numberOfBands 	number ob bands to divide spectrum to
		 */
		public function initLinearAverages(numberOfBands:int):void
		{
			fft.core.initLinearAverages(numberOfBands);
			_numberOfBands = numberOfBands;
			_averageMode = AVERAGE_LINEAR;
		}
		
		/**
		 * Analyze Spectrum
		 * if you specify any of average modes it will return average data
		 * if no everages selected then original sepctrum returned
		 * you can also choose if you want normalized (using max value) data
		 */
		public function analyzeSpectrum(normalizeSpectrum:Boolean = false):ByteArray
		{
			fft.core.analyzeSpectrum(normalizeSpectrum ? 1 : 0);
			
			spectrumData = fft.getSpectrumData();
			
			if(_averageMode > 0)
			{
				var pos:int = fft.getDataMemoryOffset(FFT.TEMP_DATA);
				
				averagesData.clear();
				
				averagesData.writeBytes(_memory, pos, (_numberOfBands * fft.numberOfChannels) << 2);
				averagesData.position = 0;
				
				return averagesData;
			}
			else
			{
				return spectrumData;
			}			
		}
		
		public function set averageMode(mode:int):void
		{
			if(mode == AVERAGE_LINEAR) 
			{
				initLinearAverages(_numberOfBands);
			}
			else if(mode == AVERAGE_LOGARITHMIC)
			{
				initLogarithmicAverages();
			}
			else 
			{
				fft.core.initNoAverages();
			}
			_averageMode = mode;
		}
		public function get averageMode():int
		{
			return _averageMode;
		}

		public function frequencyToIndex(frequency:Number):int
		{
			if (frequency < _bandWidth * 0.5) return 0;
			
			if (frequency > (_sampleRate >> 1) - _bandWidth * 0.5) return fft.spectrumLength - 1;
			
			var fraction:Number = frequency / _sampleRate;
			
			return int(fft.signalLength * fraction + 0.5);
		}
		
		public function indexToFrequency(ind:int):Number
		{			
			if ( ind == 0 ) return _bandWidth * 0.25;
			
			if ( ind == fft.spectrumLength - 1 ) 
			{
				return ((_sampleRate >> 1) - (_bandWidth * 0.5)) + _bandWidth * 0.25;
			}
			
			return ind * _bandWidth;
		}
		
		public function getBand(i:int, channel:String = 'left'):Number
		{
			var maxI:int = fft.spectrumLength - 1;
			
			if (i < 0) i = 0;
			if (i > maxI) i = maxI;
			
			if(channel == 'right') i += fft.spectrumLength;
			
			spectrumData.position = i << 2;
			
			return spectrumData.readFloat();
		}
		
		public function getFrequency(frequency:Number, channel:String = 'left'):Number
		{
			return getBand(frequencyToIndex(frequency), channel);
		}
		
		public function getAverage(ind:int, channel:String = 'left'):Number
		{
			var bind:int = (ind < 0 ? 0 : (ind > _numberOfBands ? _numberOfBands : ind));
			if(channel == 'right') bind += _numberOfBands;
			
			averagesData.position = bind << 2;
			
			return averagesData.readFloat();
		}

		public function get numberOfAverageBands():int
		{
			return _numberOfBands;
		}
		
		public function get sampleRate():int
		{
			return _sampleRate;
		}
		public function set sampleRate(sr:int):void
		{
			_sampleRate = sr;
			_bandWidth = (2.0 / fft.signalLength) * (_sampleRate >> 1);
		}
		
		public function get bandWidth():Number
		{
			return _bandWidth;
		}
	}
}
